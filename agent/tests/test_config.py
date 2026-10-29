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

"""Tests for preset resolution and provider/key selection in ``config.py``.

Covers the rules that make model switching predictable: presets expand to a full
provider configuration, explicit environment variables still win, each provider
looks up its own API key, and a Claude preset falls back to a usable embedding
provider (the Anthropic API has none of its own).

Runs with pytest, or standalone:

    agent/.venv/bin/python agent/tests/test_config.py
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import config as config_module  # noqa: E402
from config import ConfigError, load_config  # noqa: E402
from providers import split_choice  # noqa: E402

# Every FGM_*/provider variable the resolver consults, cleared per test so the
# developer's own shell and agent/.env cannot influence the result.
_VARS = [
    "FGM_MODEL", "FGM_PROVIDER", "FGM_BASE_URL", "FGM_CHAT_MODEL",
    "FGM_EMBED_PROVIDER", "FGM_EMBED_MODEL", "FGM_EMBED_BASE_URL",
    "FGM_API_KEY", "FGM_TOP_K", "FGM_TEMPERATURE", "FGM_SHOW_TOOLS", "FGM_GUARDS",
    "FGM_PHYSICAL_CORES", "MISTRAL_API_KEY", "ANTHROPIC_API_KEY", "OPENAI_API_KEY",
]


@contextmanager
def env(index_manifest: dict | None = None, **values: str):
    """Run with a clean environment, no .env file, and an optional fake index."""
    saved_env = {k: os.environ.get(k) for k in _VARS}
    saved_paths = (config_module.ENV_FILE, config_module.INDEX_DIR)
    for key in _VARS:
        os.environ.pop(key, None)
    os.environ.update(values)
    with tempfile.TemporaryDirectory() as tmp:
        config_module.ENV_FILE = Path(tmp) / "absent.env"
        config_module.INDEX_DIR = Path(tmp) / "index"
        if index_manifest is not None:
            config_module.INDEX_DIR.mkdir()
            (config_module.INDEX_DIR / "manifest.json").write_text(
                json.dumps(index_manifest), encoding="utf-8"
            )
        try:
            yield
        finally:
            config_module.ENV_FILE, config_module.INDEX_DIR = saved_paths
            for key, value in saved_env.items():
                os.environ.pop(key, None)
                if value is not None:
                    os.environ[key] = value


def test_default_is_mistral_large() -> None:
    """A fresh clone gets the large model: small too often answered from the
    retrieved fragments instead of the solver source."""
    with env(MISTRAL_API_KEY="k"):
        cfg = load_config()
    assert (cfg.provider, cfg.chat_model) == ("mistral", "mistral-large-4")
    assert cfg.embed_provider == "mistral" and cfg.embed_model == "mistral-embed"


def test_mistral_preset_matches_the_default() -> None:
    """`/model mistral` must give what a fresh clone runs, not a cheaper model."""
    with env(MISTRAL_API_KEY="k"):
        assert load_config(model="mistral").chat_model == "mistral-large-4"
        assert load_config(model="mistral-large-3").chat_model == "mistral-large-2512"
        assert load_config(model="mistral-small").chat_model == "mistral-small-latest"


def test_preset_aliases_resolve() -> None:
    with env(MISTRAL_API_KEY="k"):
        assert load_config(model="mistral-large").chat_model == "mistral-large-4"
        assert load_config(model="le-chonk").chat_model == "mistral-large-4"
    with env():
        assert load_config(model="ollama").chat_model == "qwen2.5-7b-16k"


def test_claude_preset_uses_opus_and_borrows_local_embeddings() -> None:
    with env(ANTHROPIC_API_KEY="k"):
        cfg = load_config(model="claude")
    assert (cfg.provider, cfg.chat_model) == ("anthropic", "claude-opus-5")
    assert cfg.api_key == "k"
    # Anthropic has no embeddings endpoint, so the preset must point elsewhere.
    assert cfg.embed_provider == "openai"
    assert cfg.embed_model == "nomic-embed-text"
    assert cfg.embed_base_url == "http://localhost:11434/v1"
    assert cfg.embed_api_key == "local"  # a local endpoint needs no real key


def test_preset_colon_model_overrides_the_model_id() -> None:
    with env(ANTHROPIC_API_KEY="k"):
        cfg = load_config(model="claude:claude-haiku-4-5")
    assert (cfg.provider, cfg.chat_model) == ("anthropic", "claude-haiku-4-5")


def test_bare_model_id_keeps_the_current_provider() -> None:
    with env(MISTRAL_API_KEY="k"):
        cfg = load_config(model="mistral-medium-latest")
    assert (cfg.provider, cfg.chat_model) == ("mistral", "mistral-medium-latest")


def test_local_preset_needs_no_key() -> None:
    with env():
        cfg = load_config(model="local-3b")
    assert (cfg.provider, cfg.chat_model) == ("openai", "qwen2.5-3b-16k")
    assert cfg.base_url == "http://localhost:11434/v1"
    assert cfg.api_key == "local"


def test_explicit_preset_beats_a_pinned_chat_model() -> None:
    """Asking for a preset by name means that preset, not a leftover .env pin.

    Someone who wants Sonnet asks for it — `/model claude-sonnet`, or
    `/model claude:claude-sonnet-5` — so a stale FGM_CHAT_MODEL must not
    silently redirect an explicit `/model claude`.
    """
    with env(ANTHROPIC_API_KEY="k", FGM_CHAT_MODEL="claude-sonnet-5"):
        cfg = load_config(model="claude")
    assert cfg.chat_model == "claude-opus-5"


def test_same_call_override_still_wins_over_a_preset() -> None:
    with env(ANTHROPIC_API_KEY="k"):
        cfg = load_config(model="claude", chat_model="claude-haiku-4-5")
    assert cfg.chat_model == "claude-haiku-4-5"


def test_existing_index_embedder_is_reused_so_switching_is_free() -> None:
    """Switching chat models must not silently invalidate the built index."""
    manifest = {"embed_model": "nomic-embed-text", "embed_provider": "openai",
                "embed_base_url": "http://localhost:11434/v1"}
    with env(manifest, MISTRAL_API_KEY="k"):
        cfg = load_config(model="mistral")
    assert cfg.embed_model == "nomic-embed-text"
    assert cfg.embed_provider == "openai"


def test_a_local_chat_model_never_reuses_a_hosted_embedder() -> None:
    """--model local must keep questions on the machine even when the index was
    built with mistral-embed and the Mistral key is still set."""
    manifest = {"embed_model": "mistral-embed", "embed_provider": "mistral"}
    with env(manifest, MISTRAL_API_KEY="k"):
        cfg = load_config(model="local-3b")
    assert (cfg.embed_provider, cfg.embed_model) == ("openai", "nomic-embed-text")
    assert cfg.embed_base_url.startswith("http://localhost")
    # Pinned explicitly, a hosted embedder is still the user's choice.
    with env(manifest, MISTRAL_API_KEY="k", FGM_EMBED_PROVIDER="mistral",
             FGM_EMBED_MODEL="mistral-embed"):
        assert load_config(model="local-3b").embed_provider == "mistral"


def test_index_embedder_is_skipped_when_its_key_is_missing() -> None:
    """A Mistral-built index must not strand a Claude user with no Mistral key."""
    manifest = {"embed_model": "mistral-embed", "embed_provider": "mistral"}
    with env(manifest, ANTHROPIC_API_KEY="k"):
        cfg = load_config(model="claude")
    assert cfg.embed_provider == "openai"
    assert cfg.embed_model == "nomic-embed-text"


def test_pinned_embed_model_wins_over_everything() -> None:
    with env(ANTHROPIC_API_KEY="k", MISTRAL_API_KEY="m",
             FGM_EMBED_PROVIDER="mistral", FGM_EMBED_MODEL="mistral-embed"):
        cfg = load_config(model="claude")
    assert (cfg.embed_provider, cfg.embed_model) == ("mistral", "mistral-embed")
    assert cfg.embed_api_key == "m"


def test_explicit_preset_beats_env_and_does_not_inherit_base_url() -> None:
    """`--model claude` must win over a .env pinned to a local provider.

    The endpoint matters as much as the provider: inheriting FGM_BASE_URL would
    point a hosted API at localhost.
    """
    with env(ANTHROPIC_API_KEY="k", FGM_PROVIDER="openai",
             FGM_CHAT_MODEL="qwen2.5-7b-16k",
             FGM_BASE_URL="http://localhost:11434/v1"):
        cfg = load_config(model="claude")
    assert (cfg.provider, cfg.chat_model) == ("anthropic", "claude-opus-5")
    assert cfg.base_url is None


def test_fgm_model_in_env_is_only_a_default() -> None:
    """A preset named by FGM_MODEL loses to explicit FGM_* variables."""
    with env(MISTRAL_API_KEY="k", FGM_MODEL="mistral",
             FGM_CHAT_MODEL="mistral-large-latest"):
        cfg = load_config()
    assert cfg.chat_model == "mistral-large-latest"


def test_local_placeholder_key_is_not_used_for_a_hosted_provider() -> None:
    """FGM_API_KEY=local must not masquerade as a Mistral key (401 otherwise)."""
    with env(FGM_API_KEY="local"):
        try:
            load_config(model="mistral")
        except ConfigError as exc:
            assert "MISTRAL_API_KEY" in str(exc)
        else:
            raise AssertionError("expected the placeholder key to be rejected")
    # ... but it is still fine for a local endpoint.
    with env(FGM_API_KEY="local"):
        assert load_config(model="local").api_key == "local"


def test_pinned_model_keeps_the_endpoint_it_was_indexed_with() -> None:
    """A pinned model name must not be paired with a preset's wrong provider.

    FGM_EMBED_MODEL=nomic-embed-text is a *local* model; switching to a hosted
    chat preset must not start sending that name to the hosted embeddings API.
    """
    manifest = {"embed_model": "nomic-embed-text", "embed_provider": "openai",
                "embed_base_url": "http://localhost:11434/v1"}
    with env(manifest, MISTRAL_API_KEY="k", FGM_EMBED_MODEL="nomic-embed-text"):
        cfg = load_config(model="mistral")
    assert cfg.embed_provider == "openai"
    assert cfg.embed_base_url == "http://localhost:11434/v1"
    assert cfg.embed_api_key == "local"


def test_missing_key_names_the_right_variable() -> None:
    with env():
        try:
            load_config(model="claude")
        except ConfigError as exc:
            assert "ANTHROPIC_API_KEY" in str(exc)
        else:
            raise AssertionError("expected a ConfigError for the missing key")


def test_unknown_provider_is_rejected() -> None:
    with env(FGM_PROVIDER="gemini", FGM_API_KEY="k"):
        try:
            load_config()
        except ConfigError as exc:
            assert "FGM_PROVIDER" in str(exc)
        else:
            raise AssertionError("expected a ConfigError for an unknown provider")


def test_status_mode_needs_no_key() -> None:
    with env():
        cfg = load_config(require_key=False, model="claude")
    assert cfg.api_key == ""


def test_split_choice_forms() -> None:
    assert split_choice("claude")[0].name == "claude"
    assert split_choice("claude:x")[1] == "x"
    assert split_choice("qwen2.5-3b-16k") == (None, "qwen2.5-3b-16k")


# --------------------------------------------------------------------------- #
# The in-conversation model switch (ChatAgent.switch_model)
# --------------------------------------------------------------------------- #


def _bare_agent(cfg):
    """A ChatAgent with only what switch_model touches — no index, no network."""
    from chat import ChatAgent

    agent = ChatAgent.__new__(ChatAgent)
    agent.config = cfg
    agent.chat_model = cfg.chat_model
    agent.client = None
    return agent


def test_dotenv_strips_inline_comments_and_quotes() -> None:
    """`FGM_MODEL=mistral  # the default` must not name a model "mistral  # ...".

    Without this the bad value survives all the way to the provider as an
    unhelpful "model not found", far from the typo that caused it.
    """
    import tempfile

    from config import _load_dotenv

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / ".env"
        path.write_text(
            "FGM_MODEL=mistral        # mistral-large-latest\n"
            "FGM_CHAT_MODEL='claude-opus-5'\n"
            'FGM_BASE_URL="http://localhost:11434/v1"\n'
            "# a whole-line comment\n"
            "FGM_API_KEY=sk-with#hash\n"
            "FGM_TOP_K=6\n"
        )
        values = _load_dotenv(path)
    assert values["FGM_MODEL"] == "mistral"
    assert values["FGM_CHAT_MODEL"] == "claude-opus-5"
    assert values["FGM_BASE_URL"] == "http://localhost:11434/v1"
    # No space before the '#', so it is part of the value, not a comment.
    assert values["FGM_API_KEY"] == "sk-with#hash"
    assert values["FGM_TOP_K"] == "6"


def test_switch_model_lists_without_switching() -> None:
    with env():
        agent = _bare_agent(load_config(model="local"))
        result = agent.switch_model()
    assert result["status"] == "current"
    assert "mistral" in [m["name"] for m in result["available"]]
    assert agent.chat_model == "qwen2.5-7b-16k", "listing must not switch"


def test_switch_model_changes_provider_and_model() -> None:
    with env(ANTHROPIC_API_KEY="k"):
        agent = _bare_agent(load_config(model="local"))
        result = agent.switch_model("claude")
    assert result["status"] == "switched"
    assert agent.config.provider == "anthropic"
    assert agent.chat_model == "claude-opus-5"
    assert agent.config.base_url is None, "a hosted preset must drop the local endpoint"
    assert result["previous"].startswith("qwen2.5-7b-16k")
    # Embeddings deliberately stay on the local embedder so the index survives.
    assert agent.config.embed_model == "nomic-embed-text"


def test_a_word_of_a_preset_name_selects_that_preset() -> None:
    from providers import split_choice
    assert split_choice("3B")[0].name == "local-3b"
    assert split_choice("sonnet")[0].name == "claude-sonnet"
    assert split_choice("small")[0].name == "mistral-small"
    # Several presets share "claude": that stays a bare model id, not a guess.
    assert split_choice("claude-opus-4-1") == (None, "claude-opus-4-1")


def test_an_impossible_model_id_is_refused_at_switch_time() -> None:
    # "/model 3B" used to be taken as a Claude model id and only failed, with a
    # 404, at the next question. Now "3B" is the local-3b preset, and a word that
    # fits no preset and no Claude id is refused straight away.
    with env(ANTHROPIC_API_KEY="k"):
        agent = _bare_agent(load_config(model="claude-sonnet"))
        result = agent.switch_model("qwen")
    assert result["status"] == "error" and "not" in result["error"]
    assert agent.chat_model == "claude-sonnet-5", "a refused switch must not change state"
    with env(ANTHROPIC_API_KEY="k"):
        agent = _bare_agent(load_config(model="claude-sonnet"))
        assert agent.switch_model("claude-haiku-4-5")["status"] == "switched"
    with env():           # a local server accepts any tag
        agent = _bare_agent(load_config(model="local"))
        assert agent.switch_model("my-finetune:latest")["status"] == "switched"


def test_switching_to_the_3b_says_to_ask_solver_questions_elsewhere() -> None:
    with env(MISTRAL_API_KEY="k"):
        agent = _bare_agent(load_config(model="local"))
        result = agent.switch_model("3B")
        assert result["status"] == "switched"
        assert "larger" in result["caveat"]
        # The large models carry no such note.
        assert "caveat" not in agent.switch_model("local")


def test_switch_model_reports_a_missing_key_without_switching() -> None:
    with env():
        agent = _bare_agent(load_config(model="local"))
        result = agent.switch_model("claude")
    assert result["status"] == "error"
    assert "ANTHROPIC_API_KEY" in result["error"]
    assert agent.chat_model == "qwen2.5-7b-16k", "a failed switch must not change state"


def test_switch_model_refuses_to_invalidate_the_index() -> None:
    """Switching must never silently strand the built index on another embedder."""
    manifest = {"embed_model": "nomic-embed-text", "embed_provider": "openai",
                "embed_base_url": "http://localhost:11434/v1"}
    with env(manifest, MISTRAL_API_KEY="k"):
        agent = _bare_agent(load_config(model="local"))
        # Pin a different embedder so the target config genuinely disagrees.
        os.environ["FGM_EMBED_PROVIDER"] = "mistral"
        os.environ["FGM_EMBED_MODEL"] = "mistral-embed"
        result = agent.switch_model("mistral")
    assert result["status"] == "error"
    assert "index_build.py --force" in result["error"]
    assert agent.chat_model == "qwen2.5-7b-16k"


def test_switch_model_warns_when_a_name_is_not_a_preset() -> None:
    """A raw model id is allowed, but must not look like a validated preset."""
    with env():
        agent = _bare_agent(load_config(model="local"))
        result = agent.switch_model("qwen2.5-3b-16k")
    assert result["status"] == "switched"
    assert agent.chat_model == "qwen2.5-3b-16k"
    assert "not a known preset" in result["warning"]
    # A real preset carries no such warning.
    with env():
        agent = _bare_agent(load_config(model="local"))
        assert "warning" not in agent.switch_model("local-3b")


def test_switch_model_is_not_runnable_by_the_stateless_tool_layer() -> None:
    import tools

    assert "switch_model" in tools.SESSION_TOOLS
    assert "must be handled by the chat agent" in (
        tools.execute_tool("switch_model", {"model": "claude"})["error"]
    )


def test_switch_model_is_not_blocked_by_the_feasibility_guard() -> None:
    """Switching models must stay callable in compound messages.

    The guard protects against fabricated parameters reaching build/run/delete.
    A model switch invents nothing and is reversible, so guarding it only broke
    "which models can I use? then switch to the local one", where the question
    half blocked the instruction half.
    """
    import chat

    assert "switch_model" not in chat._ACTION_TOOLS
    # The guard itself still covers everything that touches a case.
    for name in ("build_case", "run_case", "delete_cases", "kill_jobs"):
        assert name in chat._ACTION_TOOLS


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
