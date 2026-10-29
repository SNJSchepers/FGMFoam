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

"""Conversational entry point for the FGMFoam agent.

Launches a chat session that answers questions about the solver — grounded in
its own source via the RAG index — and builds the laminar slot burner case on
request. Run from the repository root:

    agent/.venv/bin/python agent/main.py                  # interactive chat
    agent/.venv/bin/python agent/main.py "set up phi 0.7 ..."   # one-shot

Configuration (models, API key, tunables) is resolved by `config.py`; the
knowledge index is auto-built on first run if missing.
"""

from __future__ import annotations

import re
import sys
import textwrap
import time
from pathlib import Path
from typing import Optional

from cases import CasePluginError

# The case plugins are loaded and validated on import (case_templates). A wrong
# declaration names its file and every problem; show that, not a traceback.
try:
    import config as config_module
    from chat import ChatAgent
    from config import Config, ConfigError
    from case_templates import describe_templates
    from providers import caveat_for, describe_presets
    from rag.index import VectorIndex, build_index, index_exists, staleness
except CasePluginError as exc:
    sys.exit(f"A case declaration has problems, so the agent did not start.\n{exc}")

try:
    from rich.console import Console
    from rich.markdown import Markdown

    _console = Console()
except ImportError:  # rich is a listed dependency; fall back if it is missing
    _console = None

# Style for the diagnostic activity lines (`· ...`). Visible but unobtrusive.
_DIAG_STYLE = "cyan"


def render_markdown(text: str) -> None:
    """Print the agent's reply, rendered as markdown with rich when available."""
    if _console is None:
        print(text)
    else:
        _console.print(Markdown(text))


def _diag(text: str) -> None:
    """Print one diagnostic activity line so the user can see what's happening."""
    line = f"· {text}"
    if _console is None:
        print(line)
    else:
        _console.print(line, style=_DIAG_STYLE, markup=False, highlight=False)


def show_tool_call(name: str, args: dict) -> None:
    """Diagnostic line for a tool the agent is invoking (e.g. read_case_file)."""
    summary = ", ".join(f"{key}={value!r}" for key, value in args.items())
    _diag(f"{name}({summary})")


def show_event(message: str) -> None:
    """Diagnostic line for non-tool activity (retrieval, forced reads)."""
    _diag(message)


def _show_timing(seconds: float) -> None:
    """Print how long the turn took (always shown; local models are slower)."""
    text = f"answered in {seconds:.1f} s"
    if _console is None:
        print(f"· {text}")
    else:
        _console.print(f"· {text}", style="green", markup=False, highlight=False)


BANNER = r"""
███████╗ ██████╗ ███╗   ███╗███████╗ ██████╗  █████╗ ███╗   ███╗
██╔════╝██╔════╝ ████╗ ████║██╔════╝██╔═══██╗██╔══██╗████╗ ████║
█████╗  ██║  ███╗██╔████╔██║█████╗  ██║   ██║███████║██╔████╔██║
██╔══╝  ██║   ██║██║╚██╔╝██║██╔══╝  ██║   ██║██╔══██║██║╚██╔╝██║
██║     ╚██████╔╝██║ ╚═╝ ██║██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║
╚═╝      ╚═════╝ ╚═╝     ╚═╝╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝
									                          						
									 ▄▄▄   ▄▄▄▄ ▄▄▄▄▄ ▄▄  ▄▄ ▄▄▄▄▄▄ 
									██▀██ ██ ▄▄ ██▄▄  ███▄██   ██   
									██▀██ ▀███▀ ██▄▄▄ ██ ▀██   ██   																
"""


HISTORY_FILE = Path.home() / ".fgmfoam_agent_history"


def _enable_line_editing() -> None:
    """Arrow keys, editing and history for the prompt.

    Plain ``input()`` only gets them when the ``readline`` module is loaded;
    without it the arrows print escape codes (``^[[D``) instead of moving the
    cursor. History (up/down) is kept across sessions in ``HISTORY_FILE``.
    ``readline`` is not available on native Windows, where this is skipped.
    """
    try:
        import atexit
        import readline
    except ImportError:
        return
    try:
        readline.read_history_file(HISTORY_FILE)
    except OSError:
        pass                                    # first run, or unreadable
    readline.set_history_length(1000)

    def save() -> None:
        try:
            readline.write_history_file(HISTORY_FILE)
        except OSError:
            pass
    atexit.register(save)


def print_banner(config: Config) -> None:
    """Print the startup logo plus a short capability / usage summary."""
    tagline = "  A code-aware CLI assistant for the FGMFoam OpenFOAM solver."
    bullets = [
        "  • Ask about the solver, its libraries, and the example cases —",
        "    answers are grounded in the source and cited.",
        "  • Set up a test case (slot burner or 1D flamelet) from a plain-language",
        "    description — /cases lists them.",
        "  • Run, monitor, stop, and delete cases (single runs or core-aware sweeps).",
    ]
    caveat = caveat_for(config.chat_model, config.provider)
    note = textwrap.fill(caveat, width=80, initial_indent="  Note: ",
                         subsequent_indent="        ") + "\n" if caveat else ""
    footer = (
        f"  Model: {config_module.describe(config)}\n"
        f"{note}"
        "  switch: /model <name>   ·   list: /models, /cases   ·   "
        "exit: 'exit' / 'quit' / Ctrl-D"
    )
    if _console is None:
        print(BANNER)
        print(tagline + "\n")
        print("\n".join(bullets))
        print("\n" + footer + "\n")
        return
    _console.print(
        BANNER, style="bold orange1", markup=False, highlight=False, soft_wrap=True
    )
    _console.print(tagline + "\n", markup=False, highlight=False)
    for bullet in bullets:
        _console.print(bullet, markup=False, highlight=False)
    _console.print("\n" + footer + "\n", style="dim", markup=False, highlight=False)


EXIT_WORDS = {"exit", "quit", ":q", "bye"}


def _ensure_index(config: Config) -> None:
    """Build the index if missing or built with a different embedding model; warn
    (but proceed) for lesser staleness such as corpus content drift."""
    if not index_exists():
        print("Building knowledge index (first run, one-time)...")
        build_index(config).save()
        print("Index ready.\n")
        return
    reasons = staleness(config.embed_model)
    if not reasons:
        return
    # A different embedding model makes the stored vectors unusable (different
    # dimensionality / not comparable), e.g. after switching provider — so rebuild
    # automatically rather than warn (and then crash on a dimension mismatch).
    if any("embedding model" in reason for reason in reasons):
        print("Embedding model changed; rebuilding the knowledge index (one-time)...")
        build_index(config).save()
        print("Index ready.\n")
        return
    print("Note: knowledge index may be out of date (" + "; ".join(reasons) + ").")
    print("      Run 'python agent/index_build.py' to refresh.\n")


def _parse_cli(argv: list[str]) -> tuple[Optional[str], list[str]]:
    """Pull a --model / --chat-model override out of argv.

    The value may be a preset name ("claude", "local"), a "preset:model-id" pair,
    or a bare model id — see :mod:`providers`.

    Returns ``(model_override, remaining_args)``; the remaining args form the
    one-shot prompt (empty -> interactive REPL).
    """
    model: Optional[str] = None
    rest: list[str] = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg in ("--model", "--chat-model") and i + 1 < len(argv):
            model = argv[i + 1]
            i += 2
            continue
        if arg.startswith(("--model=", "--chat-model=")):
            model = arg.split("=", 1)[1]
            i += 1
            continue
        rest.append(arg)
        i += 1
    return model, rest


def _handle_model_command(agent: ChatAgent, user_input: str) -> None:
    """REPL `/model [choice]`: show the active model, or switch to ``choice``.

    ``choice`` is a preset name ("claude", "mistral-large", "local-3b"), a
    "preset:model-id" pair, or a bare model id that keeps the current provider.
    Switching re-resolves the configuration, so it can move between the Mistral,
    Claude, OpenAI and local providers — the conversation is kept, and the
    knowledge index is not rebuilt because the embedding model is held steady.
    """
    parts = user_input.split(maxsplit=1)
    choice = parts[1].strip() if len(parts) > 1 else ""
    # Same implementation the conversational switch_model tool uses, so typing
    # "/model claude" and asking "switch to claude" behave identically.
    result = agent.switch_model(choice)
    if result["status"] == "switched":
        _diag(f"switched -> {result['current']}")
        for key in ("caveat", "warning"):
            if result.get(key):
                _diag(f"note: {result[key]}")
    elif result["status"] == "current":
        _diag(f"current model: {result['current']}")
        _diag("switch with /model <name>, list options with /models")
    else:
        print(f"[cannot switch to {choice!r}: {result['error']}]", file=sys.stderr)


def _handle_models_command(agent: ChatAgent) -> None:
    """REPL `/models`: list the presets, marking the active one."""
    print(describe_presets(agent.chat_model, agent.config.provider))


def _provider_error(exc: Exception) -> str:
    """Readable one-line form of a provider failure.

    Provider SDKs put the useful part — "not available in your subscription
    tier", "Rate limit exceeded" — inside a JSON body appended to the exception
    text. Surfacing that beats both a bare class name and a stack trace.
    """
    text = str(exc)
    match = re.search(r'"message"\s*:\s*"([^"]+)"', text)
    if match:
        detail = match.group(1)
        status = re.search(r"Status (\d{3})", text)
        return f"{detail} (HTTP {status.group(1)})" if status else detail
    return text.splitlines()[0] if text else type(exc).__name__


def _handle_cases_command() -> None:
    """REPL `/cases`: list the case templates that can be set up."""
    print(describe_templates())


def main() -> int:
    model_override, prompt_args = _parse_cli(sys.argv[1:])
    try:
        config = config_module.load_config(model=model_override)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 1

    _ensure_index(config)
    on_tool = show_tool_call if config.show_tools else None
    on_event = show_event if config.show_tools else None
    agent = ChatAgent(
        config, VectorIndex.load(), on_tool_call=on_tool, on_event=on_event
    )

    # One-shot mode: answer a single prompt from the command line and exit.
    if prompt_args:
        prompt = " ".join(prompt_args)
        print(f"> {prompt}\n")
        start = time.perf_counter()
        try:
            reply = agent.ask(prompt)
        except Exception as exc:  # noqa: BLE001 - report, don't traceback
            # The REPL already degrades gracefully; one-shot must too, or a
            # provider-side problem (an exhausted quota, a model the key cannot
            # reach) surfaces as a stack trace instead of the API's own message.
            print(f"\n[request failed: {_provider_error(exc)}]", file=sys.stderr)
            return 1
        render_markdown(reply)
        _show_timing(time.perf_counter() - start)
        return 0

    # Interactive REPL.
    _enable_line_editing()
    print_banner(config)
    while True:
        try:
            user_input = input("> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if not user_input:
            continue
        if user_input.lower() in EXIT_WORDS:
            break
        command = user_input.split(maxsplit=1)[0].lower()
        if command == "/cases":
            _handle_cases_command()
            continue
        if command == "/models":
            _handle_models_command(agent)
            continue
        if command == "/model":
            _handle_model_command(agent, user_input)
            continue
        try:
            start = time.perf_counter()
            reply = agent.ask(user_input)
            elapsed = time.perf_counter() - start
            print()
            render_markdown(reply)
            _show_timing(elapsed)
            print()
        except Exception as exc:  # noqa: BLE001 - keep the session alive on errors
            print(
                f"\n[request failed: {_provider_error(exc)}]\n"
                "If this is a rate limit, wait a few seconds and retry — any "
                "running simulations continue in the background regardless. "
                "You can also switch model with /model (see /models).\n",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
