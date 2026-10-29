# FGMFoam agent — a code-aware CLI assistant for the FGMFoam solver.
#
# Copyright (C) 2025-2026 Stijn Schepers, Eindhoven University of Technology
#
# This file is part of FGMFoam. FGMFoam is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# FGMFoam is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# FGMFoam. If not, see <https://www.gnu.org/licenses/>.

"""Run the agent evaluation campaign against one or more chat models.

    agent/.venv/bin/python agent/evals/run_eval.py                       # mistral + local-3b, all
    agent/.venv/bin/python agent/evals/run_eval.py --models local-3b --suite qa
    agent/.venv/bin/python agent/evals/run_eval.py --only S05 Q03 --repeat 3
    agent/.venv/bin/python agent/evals/run_eval.py --list
    agent/.venv/bin/python agent/evals/run_eval.py --guards shadow       # the model on its own

Every scenario runs in a fresh ChatAgent through the real chat loop (retrieval,
forced reads, guards, tool dispatch), so what is measured is the agent a user
gets. With --guards shadow the chat-loop guards only report what they would
have done (FGM_GUARDS=shadow), so the score is the model's own; every guard
detection is recorded either way. Two things are sandboxed:

  * case and campaign output goes to agent/evals/.sandbox/, never testcases/ or
    agent/campaigns/ — an eval can not overwrite a real case or campaign;
  * the solver is never started: run_case / run_cases / restart_case /
    start_campaign / kill_jobs are replaced by stubs that report success.

Results: agent/evals/results/<timestamp>/results.json (full transcripts, tool
calls, per-check outcomes) and report.md (the tables).
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

EVAL_DIR = Path(__file__).resolve().parent
AGENT_DIR = EVAL_DIR.parent
sys.path.insert(0, str(AGENT_DIR))
sys.path.insert(0, str(EVAL_DIR))

import campaign as campaign_module  # noqa: E402
import config as config_module  # noqa: E402
import of_writer  # noqa: E402
import tools as tools_module  # noqa: E402
import validator  # noqa: E402
from chat import ChatAgent  # noqa: E402
from rag.index import VectorIndex  # noqa: E402
from scenarios import SUITES  # noqa: E402

SANDBOX = EVAL_DIR / ".sandbox"
RESULTS_DIR = EVAL_DIR / "results"
REPO_ROOT_TESTCASES = AGENT_DIR.parent / "testcases"
# The release default and the measured local model; the 7B is slow on a small GPU.
DEFAULT_MODELS = ("mistral", "local-3b")
REL_TOL = 0.01

_STUBBED_RUN_TOOLS = ("run_case", "run_cases", "restart_case", "start_campaign", "kill_jobs")


# --------------------------------------------------------------------------- #
# Sandbox                                                                     #
# --------------------------------------------------------------------------- #

def _stub(name: str):
    def run(args: dict[str, Any]) -> dict[str, Any]:
        return {"status": "started",
                "note": "evaluation sandbox: the solver was not actually launched",
                "tool": name}
    return run


def install_sandbox() -> None:
    """Point every write at the sandbox and replace the solver launchers."""
    if SANDBOX.exists():
        shutil.rmtree(SANDBOX)
    (SANDBOX / "testcases").mkdir(parents=True)
    (SANDBOX / "campaigns").mkdir()
    of_writer.TESTCASES_DIR = SANDBOX / "testcases"
    campaign_module.CAMPAIGN_DIR = SANDBOX / "campaigns"
    for name in _STUBBED_RUN_TOOLS:
        stub = _stub(name)
        tools_module._DISPATCH[name] = stub
        # run_campaign / start_campaign call the launcher directly, not via _DISPATCH.
        setattr(tools_module, f"_run_{name}", stub)


def reset_sandbox() -> None:
    for sub in ("testcases", "campaigns"):
        path = SANDBOX / sub
        shutil.rmtree(path, ignore_errors=True)
        path.mkdir(parents=True)


# --------------------------------------------------------------------------- #
# Recording                                                                   #
# --------------------------------------------------------------------------- #

@dataclass
class Call:
    name: str
    args: dict[str, Any]
    result: dict[str, Any] = field(default_factory=dict)
    forced: bool = False      # the chat loop's deterministic read, not the model's choice


@dataclass
class TurnLog:
    calls: list[Call] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)   # guard-blocked tool names
    events: list[str] = field(default_factory=list)
    guards: list[dict[str, Any]] = field(default_factory=list)  # every detection


class Recorder:
    """Hooks into ChatAgent and the tool layer to see what really happened."""

    def __init__(self) -> None:
        self.turn = TurnLog()
        self._forced_pending = False
        self._real_execute = tools_module.execute_tool
        tools_module.execute_tool = self._execute

    def new_turn(self) -> TurnLog:
        self.turn = TurnLog()
        self._forced_pending = False
        return self.turn

    def _execute(self, name: str, args: dict[str, Any]) -> dict[str, Any]:
        result = self._real_execute(name, args)
        self.turn.calls.append(Call(name, dict(args), result))
        return result

    def on_event(self, message: str) -> None:
        self.turn.events.append(message)
        if message.startswith("implementation question"):
            self._forced_pending = True
        match = re.match(r"skipped (\w+)", message)
        if match:
            self.turn.skipped.append(match.group(1))

    def on_guard(self, kind: str, name: str, message: str, acted: bool) -> None:
        self.turn.guards.append({"kind": kind, "name": name, "message": message,
                                 "acted": acted})

    def on_tool_call(self, name: str, args: dict[str, Any]) -> None:
        if self._forced_pending and name == "read_case_file":
            # The forced read executed (and was recorded) just before this hook;
            # it may also have probed a sibling .C first — mark all of them.
            for call in reversed(self.turn.calls):
                if call.name != "read_case_file" or call.forced:
                    break
                call.forced = True
            self._forced_pending = False

    def restore(self) -> None:
        tools_module.execute_tool = self._real_execute


# --------------------------------------------------------------------------- #
# Checks                                                                      #
# --------------------------------------------------------------------------- #

def _num(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _value_matches(expected: Any, actual: Any) -> bool:
    if isinstance(expected, bool) or isinstance(actual, bool):
        if isinstance(actual, str):
            actual = actual.strip().lower() in ("true", "1", "yes")
        return bool(expected) == bool(actual)
    if isinstance(expected, (int, float)):
        got = _num(actual)
        return got is not None and abs(got - expected) <= REL_TOL * max(abs(expected), 1e-30)
    if isinstance(expected, list):
        got = [_num(v) for v in actual] if isinstance(actual, list) else []
        if len(got) != len(expected) or None in got:
            return False
        return all(_value_matches(e, a) for e, a in zip(sorted(expected), sorted(got)))
    return str(expected).strip().lower() == str(actual).strip().lower()


def _args_match(expected: dict[str, Any], args: dict[str, Any]) -> bool:
    canon = {validator.canonical_name(k): v for k, v in args.items()}
    for key, value in expected.items():
        if key == "case_type":
            got = canon.get("case_type")
            if got is None or value.lower() not in str(got).lower().replace(" ", ""):
                return False
        elif key not in canon or not _value_matches(value, canon[key]):
            return False
    return True


def _status(call: Call) -> str:
    return str(call.result.get("status") or ("error" if "error" in call.result else ""))


def _matches(pattern: str | list[str], text: str) -> bool:
    patterns = pattern if isinstance(pattern, list) else [pattern]
    return any(re.search(p, text, re.MULTILINE) for p in patterns)


def evaluate_turn(spec: dict[str, Any], log: TurnLog, reply: str,
                  last_built: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Return one {check, ok, detail} per check in the turn spec."""
    out: list[dict[str, Any]] = []
    model_calls = [c for c in log.calls if not c.forced]

    if "tool" in spec:
        want = spec["tool"]
        # A name or a list of acceptable names (the outcome matters, e.g. a
        # rebuild may come from build_case or from change_settings).
        names = want["name"] if isinstance(want["name"], list) else [want["name"]]
        candidates = [c for c in model_calls if c.name in names]
        hit = next((c for c in candidates
                    if _args_match(want.get("args", {}), c.args)
                    and ("status" not in want or _status(c) == want["status"])), None)
        detail = "; ".join(f"{c.name}({json.dumps(c.args)}) -> {_status(c)}"
                           for c in candidates) or "not called"
        out.append({"check": f"tool {'/'.join(names)}", "ok": hit is not None,
                    "detail": detail})

    for name, status in spec.get("no_status", []):
        bad = [c for c in model_calls if c.name == name and _status(c) == status]
        out.append({"check": f"no {name}->{status}", "ok": not bad,
                    "detail": f"{len(bad)} call(s)" if bad else ""})

    for name in spec.get("no_tool", []):
        bad = [c for c in model_calls if c.name == name]
        blocked = log.skipped.count(name)
        out.append({"check": f"no {name}", "ok": not bad,
                    "detail": (f"executed {len(bad)}x" if bad else
                               f"attempted, blocked by guard {blocked}x" if blocked else "")})

    if "case" in spec:
        params = (last_built or {}).get("parameters", {})
        wrong = {k: params.get(k) for k, v in spec["case"].items()
                 if not _value_matches(v, params.get(k))}
        out.append({"check": "case parameters", "ok": last_built is not None and not wrong,
                    "detail": "no case built" if last_built is None else
                              (f"got {wrong}" if wrong else "")})

    for pattern in spec.get("reply_all", []):
        out.append({"check": f"reply ~ {pattern}", "ok": _matches(pattern, reply), "detail": ""})
    for pattern in spec.get("reply_none", []):
        found = re.search(pattern, reply, re.MULTILINE)
        out.append({"check": f"reply !~ {pattern}", "ok": found is None,
                    "detail": found.group(0) if found else ""})

    reads = [c for c in log.calls if c.name == "read_case_file" and _status(c) == "ok"]
    if "grounded" in spec:
        paths = [str(c.args.get("path")) for c in reads]
        ok = any(re.search(spec["grounded"], p) for p in paths)
        out.append({"check": "grounded (right file read)", "ok": ok,
                    "detail": ", ".join(paths) or "nothing read"})
        model_paths = [str(c.args.get("path")) for c in reads if not c.forced]
        out.append({"check": "model chose to read", "ok": bool(model_paths),
                    "detail": ", ".join(model_paths), "soft": True})
    return out


# --------------------------------------------------------------------------- #
# Running                                                                     #
# --------------------------------------------------------------------------- #

def run_scenario(model: str, scenario: dict[str, Any], index: VectorIndex,
                 recorder: Recorder, guards: str = "on",
                 on_agent: Callable[[ChatAgent], None] | None = None) -> dict[str, Any]:
    """Run one scenario in a fresh agent. ``on_agent`` sees that agent before the
    first turn, e.g. to log its model calls."""
    reset_sandbox()
    cfg = config_module.load_config(model=model, guards=guards)
    agent = ChatAgent(cfg, index, on_tool_call=recorder.on_tool_call,
                      on_event=recorder.on_event, on_guard=recorder.on_guard)
    if on_agent is not None:
        on_agent(agent)
    turns_out: list[dict[str, Any]] = []
    last_built: dict[str, Any] | None = None
    for spec in scenario["turns"]:
        log = recorder.new_turn()
        started = time.perf_counter()
        try:
            reply, error = agent.ask(spec["user"]), None
        except Exception as exc:  # noqa: BLE001 - a crash is a result, not an abort
            reply, error = "", f"{type(exc).__name__}: {exc}"
        seconds = time.perf_counter() - started
        for call in log.calls:
            if _status(call) == "built":   # build_case, or a change_settings rebuild
                last_built = call.result
        checks = evaluate_turn(spec, log, reply, last_built)
        if error:
            checks.insert(0, {"check": "turn completed", "ok": False, "detail": error})
        turns_out.append({
            "user": spec["user"],
            "reply": reply,
            "seconds": round(seconds, 1),
            "calls": [{"name": c.name, "args": c.args, "status": _status(c),
                       "forced": c.forced} for c in log.calls],
            "skipped": log.skipped,
            "events": log.events,
            "guards": log.guards,
            "checks": checks,
        })
    hard = [c for t in turns_out for c in t["checks"] if not c.get("soft")]
    return {
        "model": model,
        "guard_mode": guards,
        "id": scenario["id"],
        "title": scenario["title"],
        "passed": all(c["ok"] for c in hard),
        "checks_passed": sum(c["ok"] for c in hard),
        "checks_total": len(hard),
        "model_read": any(c["ok"] for t in turns_out for c in t["checks"]
                          if c["check"] == "model chose to read"),
        "seconds": round(sum(t["seconds"] for t in turns_out), 1),
        "turns": turns_out,
    }


def rescore(path: Path) -> int:
    """Re-apply the CURRENT reply checks to a stored results.json.

    For fixing a check (a regex too strict or too loose) without spending another
    hour of model time. Only reply_all / reply_none are recomputed — tool, case and
    grounding checks depend on tool results that are kept only in summary form.
    """
    specs = {s["id"]: s for suite in SUITES.values() for s in suite}
    results = json.loads(path.read_text())
    for res in results:
        spec = specs.get(res["id"])
        if spec is None:
            continue
        for turn, turn_spec in zip(res["turns"], spec["turns"]):
            reply_only = {k: turn_spec[k] for k in ("reply_all", "reply_none") if k in turn_spec}
            fresh = evaluate_turn(reply_only, TurnLog(), turn["reply"], None)
            kept = [c for c in turn["checks"] if not c["check"].startswith("reply ")]
            turn["checks"] = kept + fresh
        hard = [c for t in res["turns"] for c in t["checks"] if not c.get("soft")]
        res["passed"] = all(c["ok"] for c in hard)
        res["checks_passed"] = sum(c["ok"] for c in hard)
        res["checks_total"] = len(hard)
    path.write_text(json.dumps(results, indent=1))
    models = list(dict.fromkeys(r["model"] for r in results))
    write_report(results, models, path.parent / "report.md")
    print(f"rescored: {path.parent / 'report.md'}")
    return 0


def unload_local_model(model: str) -> None:
    """Free a local Ollama model before the next one loads.

    Ollama keeps a model resident for minutes after its last request. On a small
    machine (WSL with ~8 GB) the 3B's server (~4 GB resident) plus the 7B does
    not fit, and the eval gets killed for memory. Hosted models: nothing to do.
    """
    import urllib.request

    cfg = config_module.load_config(require_key=False, model=model)
    if not cfg.base_url or "11434" not in cfg.base_url:
        return
    root = cfg.base_url.split("/v1")[0]
    body = json.dumps({"model": cfg.chat_model, "keep_alive": 0}).encode()
    try:
        urllib.request.urlopen(urllib.request.Request(
            f"{root}/api/generate", body, {"Content-Type": "application/json"}),
            timeout=60).read()
    except OSError as exc:
        print(f"(could not unload {cfg.chat_model}: {exc})", flush=True)


def select(suite: str, only: list[str]) -> list[tuple[str, dict[str, Any]]]:
    chosen = [(name, s) for name in (["sim", "qa"] if suite == "all" else [suite])
              for s in SUITES[name]]
    if only:
        chosen = [(n, s) for n, s in chosen if any(s["id"].startswith(o) for o in only)]
    return chosen


# --------------------------------------------------------------------------- #
# Report                                                                      #
# --------------------------------------------------------------------------- #

def _guard_table(results: list[dict[str, Any]], models: list[str]) -> list[str]:
    """Guard detections per guard and model: how often the model needed help.

    "first" counts the runs where this guard was the first to fire. In shadow
    mode a first mistake is left uncorrected, so later detections in the same
    run may only be knock-on effects of it; the first one is the clean signal.
    """
    counts: dict[tuple[str, str], dict[str, list[int]]] = {}
    for r in results:
        events = [g for t in r["turns"] for g in t.get("guards", [])]
        for i, g in enumerate(events):
            cell = counts.setdefault((g["kind"], g["name"]), {}).setdefault(
                r["model"], [0, 0])
            cell[0] += 1
            cell[1] += i == 0
    if not counts:
        return []
    lines = ["## Guard detections", "",
             "| Guard | " + " | ".join(models) + " |", "|---" * (len(models) + 1) + "|"]
    for (kind, name) in sorted(counts):
        cells = [f"{c[0]} (first in {c[1]})" if (c := counts[(kind, name)].get(m))
                 else "–" for m in models]
        lines.append(f"| {kind} · {name} | " + " | ".join(cells) + " |")
    lines += ["", "Each cell: detections over all runs (runs where this guard fired "
              "first). In shadow mode a detection was only logged, not acted on.", ""]
    return lines


def write_report(results: list[dict[str, Any]], models: list[str], path: Path) -> str:
    modes = sorted({r.get("guard_mode", "on") for r in results}) or ["on"]
    lines = ["# Agent evaluation", "",
             f"Run {datetime.now():%Y-%m-%d %H:%M}. Models: {', '.join(models)}. "
             f"Guards: {', '.join(modes)}.", ""]
    for suite, title in (("sim", "A. Simulation manager"), ("qa", "B. Solver information")):
        rows = [r for r in results if r["suite"] == suite]
        if not rows:
            continue
        lines += [f"## {title}", ""]
        ids = list(dict.fromkeys(r["id"] for r in rows))
        header = "| Scenario | " + " | ".join(models) + " |"
        lines += [header, "|---" * (len(models) + 1) + "|"]
        for sid in ids:
            cells = []
            for m in models:
                runs = [r for r in rows if r["id"] == sid and r["model"] == m]
                if not runs:
                    cells.append("–")
                    continue
                npass = sum(r["passed"] for r in runs)
                c_ok = sum(r["checks_passed"] for r in runs)
                c_all = sum(r["checks_total"] for r in runs)
                mark = "✅" if npass == len(runs) else ("❌" if npass == 0 else "⚠️")
                read = ""
                if suite == "qa":
                    read = " 📖" if all(r["model_read"] for r in runs) else ""
                sec = sum(r["seconds"] for r in runs) / len(runs)
                cells.append(f"{mark} {npass}/{len(runs)} ({c_ok}/{c_all} checks, "
                             f"{sec:.0f}s){read}")
            lines.append(f"| {sid} | " + " | ".join(cells) + " |")
        lines.append("")
        for m in models:
            mine = [r for r in rows if r["model"] == m]
            if mine:
                lines.append(f"- **{m}**: {sum(r['passed'] for r in mine)}/{len(mine)} "
                             f"scenarios passed, "
                             f"{sum(r['checks_passed'] for r in mine)}/"
                             f"{sum(r['checks_total'] for r in mine)} checks")
        lines.append("")
    if any(r["suite"] == "qa" for r in results):
        lines += ["📖 = the model itself chose to open a source file (besides the "
                  "deterministic forced read). Not required to pass.", ""]

    lines += _guard_table(results, models)
    lines += ["## Failures", ""]
    for r in results:
        if r["passed"]:
            continue
        lines.append(f"### {r['model']} · {r['id']} — {r['title']}")
        for i, turn in enumerate(r["turns"], 1):
            failed = [c for c in turn["checks"] if not c["ok"] and not c.get("soft")]
            if not failed:
                continue
            lines.append(f"- turn {i}: _{turn['user']}_")
            calls = ", ".join(f"{c['name']}→{c['status']}" + (" (forced)" if c["forced"] else "")
                              for c in turn["calls"]) or "no tools"
            lines.append(f"  - tools: {calls}"
                         + (f"; guard-blocked: {', '.join(turn['skipped'])}"
                            if turn["skipped"] else ""))
            guards = turn.get("guards", [])
            if guards:
                lines.append("  - guards: " + ", ".join(
                    g["name"] + ("" if g["acted"] else " (shadow)") for g in guards))
            for c in failed:
                lines.append(f"  - ✗ `{c['check']}` {c['detail']}".rstrip())
            reply = " ".join(turn["reply"].split())
            lines.append(f"  - reply: {reply[:400]}{'…' if len(reply) > 400 else ''}")
        lines.append("")
    text = "\n".join(lines)
    path.write_text(text)
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--models", nargs="+", default=list(DEFAULT_MODELS))
    parser.add_argument("--suite", choices=("sim", "qa", "all"), default="all")
    parser.add_argument("--only", nargs="*", default=[], help="scenario id prefixes")
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--guards", choices=config_module.GUARD_MODES, default="on",
                        help="shadow: guards only report, so the model is scored alone")
    parser.add_argument("--list", action="store_true", help="list scenarios and exit")
    parser.add_argument("--rescore", type=Path, metavar="RESULTS_JSON",
                        help="re-apply the current reply checks to a stored run")
    args = parser.parse_args()
    if args.rescore:
        return rescore(args.rescore)

    chosen = select(args.suite, args.only)
    if args.list:
        for suite, s in chosen:
            print(f"{suite:4} {s['id']:32} {s['title']}  ({len(s['turns'])} turn(s))")
        return 0

    install_sandbox()
    recorder = Recorder()
    index = VectorIndex.load()
    out_dir = RESULTS_DIR / datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir.mkdir(parents=True)
    results: list[dict[str, Any]] = []
    try:
        # Model-major order: each model stays loaded in Ollama for its whole block.
        for model in args.models:
            for suite, scenario in chosen:
                for rep in range(args.repeat):
                    res = run_scenario(model, scenario, index, recorder, args.guards)
                    res["suite"], res["repeat"] = suite, rep
                    results.append(res)
                    mark = "PASS" if res["passed"] else "FAIL"
                    print(f"[{model:9}] {scenario['id']:32} {mark}  "
                          f"{res['checks_passed']}/{res['checks_total']}  {res['seconds']:.0f}s",
                          flush=True)
                    (out_dir / "results.json").write_text(json.dumps(results, indent=1))
            unload_local_model(model)
    finally:
        recorder.restore()
        shutil.rmtree(SANDBOX, ignore_errors=True)
    write_report(results, args.models, out_dir / "report.md")
    print(f"\nreport: {out_dir.relative_to(AGENT_DIR.parent)}/report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
