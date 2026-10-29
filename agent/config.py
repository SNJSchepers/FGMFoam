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

"""Central configuration for the FGMFoam agent.

All model names, API keys, and runtime tunables are resolved here so the rest of
the code never reads environment variables directly. This keeps the tool
configurable for a public release without editing code.

Resolution order (highest priority first):

    1. explicit keyword overrides (e.g. ``chat_model=`` from a CLI flag)
    2. a preset asked for explicitly (``--model claude`` / ``/model claude``)
    3. process environment variables
    4. an optional ``agent/.env`` file (KEY=VALUE lines)
    5. a preset named by ``FGM_MODEL`` in the environment or ``.env``
    6. built-in defaults

Real environment variables win over the ``.env`` file, matching the usual dotenv
convention. An *explicitly requested* preset outranks both, and defines the chat
side completely (provider, model, endpoint): asking for ``--model claude`` must
switch to Claude even when ``.env`` pins a local provider, and must not leave the
local ``FGM_BASE_URL`` pointing a hosted API at localhost. A preset named by
``FGM_MODEL`` is only a default, so individual ``FGM_*`` variables still win.

Chat and embeddings are configured independently (``FGM_PROVIDER`` /
``FGM_CHAT_MODEL`` vs ``FGM_EMBED_PROVIDER`` / ``FGM_EMBED_MODEL``). Two reasons:
the Anthropic API has no embeddings endpoint, so a Claude chat model must borrow
embeddings from elsewhere; and re-embedding the whole corpus every time the chat
model changes would make switching models slow for no benefit. When nothing pins
an embedding model, an existing index's model is reused so switching chat models
is free.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from providers import CHAT_PROVIDERS, EMBED_PROVIDERS, Preset, split_choice

# Stamped into every generated case's provenance record. Bump on a change that
# alters what the agent writes into a case.
AGENT_VERSION: str = "1.1.0"

AGENT_DIR: Path = Path(__file__).resolve().parent
REPO_ROOT: Path = AGENT_DIR.parent
ENV_FILE: Path = AGENT_DIR / ".env"
INDEX_DIR: Path = AGENT_DIR / "index"

# Built-in defaults, used when neither an override, the environment, the .env
# file, nor a preset supplies a value. The chat model defaults to Mistral's
# "large": noticeably better over the solver source than "small", which too often
# answered from retrieved fragments rather than the code. Switch with --model /
# /model at any time; "--model mistral-small" selects the cheaper small model.
DEFAULTS: dict[str, str] = {
    "FGM_MODEL": "",  # preset name; empty = use the FGM_PROVIDER/FGM_CHAT_MODEL pair
    "FGM_PROVIDER": "mistral",
    "FGM_BASE_URL": "",
    "FGM_CHAT_MODEL": "mistral-large-4",
    "FGM_EMBED_PROVIDER": "",  # empty = derive from the chat provider / preset
    "FGM_EMBED_MODEL": "",  # empty = reuse the existing index, else the preset
    "FGM_EMBED_BASE_URL": "",
    "FGM_TOP_K": "6",
    "FGM_TEMPERATURE": "0.2",
    "FGM_SHOW_TOOLS": "true",
    # "on": the chat-loop guards act. "shadow": they only log what they would
    # have done, so an eval measures the model's own behaviour.
    "FGM_GUARDS": "on",
    # Number of physical cores to budget for parallel runs. "auto" detects it.
    "FGM_PHYSICAL_CORES": "auto",
}

# Fallback embedding model per embedding provider, used when a chat-only preset
# (or a bare FGM_PROVIDER) leaves the embedding side unspecified.
_DEFAULT_EMBED_MODEL: dict[str, str] = {
    "mistral": "mistral-embed",
    "openai": "text-embedding-3-small",
}

# Provider-specific key variables, tried before the generic FGM_API_KEY.
_KEY_ENV: dict[str, tuple[str, ...]] = {
    "mistral": ("MISTRAL_API_KEY",),
    "anthropic": ("ANTHROPIC_API_KEY",),
    "openai": ("OPENAI_API_KEY",),
}

_TRUE_VALUES = {"1", "true", "yes", "on"}
GUARD_MODES = ("on", "shadow")


class ConfigError(RuntimeError):
    """Raised when configuration is missing or invalid."""


@dataclass(frozen=True)
class Config:
    """Resolved runtime configuration."""

    api_key: str
    provider: str  # "mistral" | "anthropic" | "openai"
    base_url: Optional[str]  # OpenAI-compatible endpoint (None for hosted APIs)
    chat_model: str
    # Embeddings are resolved separately — see the module docstring.
    embed_provider: str  # "mistral" | "openai"
    embed_model: str
    embed_base_url: Optional[str]
    embed_api_key: str
    top_k: int
    temperature: float
    show_tools: bool
    physical_cores: Optional[int]  # None = auto-detect
    preset: Optional[str] = None  # name of the preset in effect, if any
    guards: str = "on"  # "on" | "shadow" — see GUARD_MODES


def _load_dotenv(path: Path) -> dict[str, str]:
    """Parse a simple ``KEY=VALUE`` ``.env`` file (no external dependency).

    Handles the two things people actually write: surrounding quotes, and a
    trailing ``# comment``. An unquoted value ends at the first whitespace-
    preceded ``#`` — without that, ``FGM_MODEL=mistral  # the default`` silently
    becomes a model named "mistral  # the default", which fails much later as an
    unhelpful provider error. A ``#`` inside quotes, or with no space before it,
    is kept: API keys may legitimately contain one.
    """
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        value = value.strip()
        if value[:1] in ('"', "'"):
            quote = value[0]
            end = value.find(quote, 1)
            value = value[1:end] if end > 0 else value[1:]
        else:
            value = re.split(r"\s+#", value, maxsplit=1)[0].strip()
        if key:
            values[key] = value
    return values


def _index_embedding() -> tuple[str, str, str]:
    """Return ``(embed_model, embed_provider, embed_base_url)`` of the built index.

    Reusing what the index was built with is what makes switching chat models
    free: the stored vectors stay valid, so no rebuild is triggered. Returns
    empty strings when there is no readable index.
    """
    try:
        manifest = json.loads((INDEX_DIR / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return "", "", ""
    return (
        str(manifest.get("embed_model") or ""),
        # "provider" is what older manifests recorded the embedder under.
        str(manifest.get("embed_provider") or manifest.get("provider") or ""),
        str(manifest.get("embed_base_url") or ""),
    )


def _is_local(base_url: str) -> bool:
    """True for an endpoint on this machine, which needs no real API key."""
    return any(host in base_url for host in ("localhost", "127.0.0.1", "0.0.0.0", "::1"))


def _find_key(provider: str, base_url: str, env: dict[str, str]) -> str:
    """Best-effort API key for ``provider`` ("" when none is configured).

    Local OpenAI-compatible servers ignore the key entirely, so a placeholder is
    returned for them rather than treating a missing key as an error.
    """
    local_endpoint = provider == "openai" and bool(base_url) and _is_local(base_url)
    for name in (*_KEY_ENV.get(provider, ()), "FGM_API_KEY"):
        value = os.environ.get(name) or env.get(name)
        # "local" is the placeholder for a local server that ignores the key.
        # Handing it to a hosted API would turn a clear "no key set" message
        # into a confusing 401 from the provider, so it is ignored there.
        if value and not (value.strip().lower() == "local" and not local_endpoint):
            return value
    if local_endpoint:
        return "local"
    return ""


def load_config(require_key: bool = True, **overrides: Optional[str]) -> Config:
    """Resolve and validate the configuration.

    Args:
        require_key: If True (default), a missing API key raises. Set to False
            for operations that do not call the API (e.g. checking index
            status), in which case the keys are returned as empty strings.
        **overrides: Optional explicit values (e.g. ``model="claude"`` or
            ``chat_model="..."`` from a CLI flag). ``None`` values are ignored so
            callers can pass through unset flags directly.

    Returns:
        A populated :class:`Config`.

    Raises:
        ConfigError: if a required API key is missing, a provider name is
            unknown, or a numeric value is invalid.
    """
    dotenv = _load_dotenv(ENV_FILE)

    def raw(env_name: str, override_key: str) -> str:
        """Value from an override, the environment, or the .env file ("" if unset)."""
        if overrides.get(override_key) is not None:
            return str(overrides[override_key])
        return os.environ.get(env_name) or dotenv.get(env_name) or ""

    def resolve(env_name: str, override_key: str, preset_value: str = "") -> str:
        return raw(env_name, override_key) or preset_value or DEFAULTS[env_name]

    # --- Preset --------------------------------------------------------------
    # ``model`` accepts a preset name, "preset:model-id", or a bare model id.
    preset: Optional[Preset] = None
    preset_model: Optional[str] = None
    choice = raw("FGM_MODEL", "model")
    if choice:
        preset, preset_model = split_choice(choice)
        if preset is None and preset_model is None:
            raise ConfigError(f"Empty model choice: {choice!r}")

    # A preset asked for directly (--model / /model) is a deliberate choice, so
    # it outranks the environment and .env; one named by FGM_MODEL is a default.
    explicit_preset = preset if overrides.get("model") is not None else None

    def resolve_chat(env_name: str, override_key: str, preset_value: str) -> str:
        """Resolve a chat-side field, honouring an explicitly requested preset.

        When such a preset is active it *defines* the field: the environment is
        not consulted, so a hosted preset cannot inherit a local FGM_BASE_URL.
        """
        if overrides.get(override_key) is not None:
            return str(overrides[override_key])
        if explicit_preset is not None:
            return preset_value
        return (
            os.environ.get(env_name)
            or dotenv.get(env_name)
            or preset_value
            or DEFAULTS[env_name]
        )

    # --- Chat provider -------------------------------------------------------
    provider = resolve_chat(
        "FGM_PROVIDER", "provider", preset.provider if preset else ""
    ).strip().lower()
    if provider not in CHAT_PROVIDERS:
        raise ConfigError(
            f"FGM_PROVIDER must be one of {', '.join(sorted(CHAT_PROVIDERS))} "
            f"(got {provider!r})."
        )
    base_url = resolve_chat(
        "FGM_BASE_URL", "base_url", (preset.base_url if preset else "") or ""
    ).strip()

    # A bare model id from ``--model`` overrides the preset's chat model.
    chat_model = resolve_chat(
        "FGM_CHAT_MODEL",
        "chat_model",
        preset_model or (preset.chat_model if preset else ""),
    ).strip()

    api_key = overrides.get("api_key") or _find_key(provider, base_url, dotenv)
    if not api_key and require_key:
        raise ConfigError(_missing_key_message(provider, base_url))

    # --- Embeddings ----------------------------------------------------------
    embed_provider, embed_model, embed_base_url, embed_api_key = _resolve_embeddings(
        raw=raw,
        preset=preset,
        chat_provider=provider,
        chat_base_url=base_url,
        dotenv=dotenv,
        require_key=require_key,
    )

    # --- Tunables ------------------------------------------------------------
    cores_raw = resolve("FGM_PHYSICAL_CORES", "physical_cores").strip().lower()
    try:
        top_k = int(resolve("FGM_TOP_K", "top_k"))
        temperature = float(resolve("FGM_TEMPERATURE", "temperature"))
        physical_cores = None if cores_raw in ("", "auto") else int(cores_raw)
    except ValueError as exc:
        raise ConfigError(f"Invalid numeric configuration value: {exc}") from exc
    guards = resolve("FGM_GUARDS", "guards").strip().lower()
    if guards not in GUARD_MODES:
        raise ConfigError(
            f"FGM_GUARDS must be one of {', '.join(GUARD_MODES)}, not {guards!r}")

    return Config(
        api_key=str(api_key),
        provider=provider,
        base_url=base_url or None,
        chat_model=chat_model,
        embed_provider=embed_provider,
        embed_model=embed_model,
        embed_base_url=embed_base_url or None,
        embed_api_key=embed_api_key,
        top_k=top_k,
        temperature=temperature,
        show_tools=resolve("FGM_SHOW_TOOLS", "show_tools").strip().lower() in _TRUE_VALUES,
        physical_cores=physical_cores,
        preset=preset.name if preset else None,
        guards=guards,
    )


def _resolve_embeddings(
    *, raw, preset: Optional[Preset], chat_provider: str, chat_base_url: str,
    dotenv: dict[str, str], require_key: bool,
) -> tuple[str, str, str, str]:
    """Pick the embedding provider/model/endpoint/key.

    Candidates are tried in order and the first one whose API key is actually
    available wins, so switching to a chat provider you have a key for never
    strands the index on an embedder you cannot reach:

        1. whatever ``FGM_EMBED_*`` pins explicitly
        2. what the existing index was built with (avoids a needless rebuild)
        3. the preset's embedding fallback, or the chat provider itself
    """
    pinned_provider = raw("FGM_EMBED_PROVIDER", "embed_provider").strip().lower()
    pinned_model = raw("FGM_EMBED_MODEL", "embed_model").strip()
    pinned_base = raw("FGM_EMBED_BASE_URL", "embed_base_url").strip()

    # The chat provider is a natural embedding source when it has an embeddings
    # endpoint (Anthropic does not, so a Claude preset supplies its own fallback).
    if preset is not None:
        fallback = (preset.embed_provider, preset.embed_model, preset.embed_base_url or "")
    elif chat_provider in EMBED_PROVIDERS:
        fallback = (
            chat_provider,
            _DEFAULT_EMBED_MODEL[chat_provider],
            chat_base_url,
        )
    else:
        fallback = ("openai", "nomic-embed-text", "http://localhost:11434/v1")

    indexed_model, indexed_provider, indexed_base = _index_embedding()
    indexed_usable = bool(indexed_model) and indexed_provider in EMBED_PROVIDERS

    candidates: list[tuple[str, str, str]] = []
    if pinned_provider or pinned_model or pinned_base:
        # A pinned model name alone does not say *where* that model lives. If it
        # is the one the index was built with, take that endpoint too — pairing
        # a local model name like "nomic-embed-text" with a hosted preset's
        # provider would otherwise send it to an API that has never heard of it.
        if (
            pinned_model
            and not pinned_provider
            and not pinned_base
            and indexed_usable
            and pinned_model == indexed_model
        ):
            candidates.append((indexed_provider, pinned_model, indexed_base))
        # Anything pinned explicitly is used as-is; unpinned parts come from the
        # fallback so "FGM_EMBED_PROVIDER=mistral" alone works too.
        candidates.append((
            pinned_provider or fallback[0],
            pinned_model or fallback[1],
            pinned_base or fallback[2],
        ))
    else:
        if indexed_usable:
            candidates.append((indexed_provider, indexed_model, indexed_base))
        candidates.append(fallback)

    for provider, model, base in candidates:
        if provider not in EMBED_PROVIDERS:
            raise ConfigError(
                f"FGM_EMBED_PROVIDER must be one of "
                f"{', '.join(sorted(EMBED_PROVIDERS))} (got {provider!r}). "
                "The Anthropic API has no embeddings endpoint, so a Claude chat "
                "model needs one of these for the knowledge index."
            )
        key = _find_key(provider, base, dotenv)
        if key or not require_key:
            return provider, model, base, key

    provider, model, base = candidates[-1]
    raise ConfigError(
        f"No API key for the embedding provider {provider!r} (model {model!r}).\n"
        f"{_missing_key_message(provider, base)}\n"
        "Or embed locally with no key at all:\n"
        "    FGM_EMBED_PROVIDER=openai\n"
        "    FGM_EMBED_BASE_URL=http://localhost:11434/v1\n"
        "    FGM_EMBED_MODEL=nomic-embed-text"
    )


def _missing_key_message(provider: str, base_url: str) -> str:
    """A provider-specific 'how to set your key' message."""
    names = {
        "mistral": "MISTRAL_API_KEY",
        "anthropic": "ANTHROPIC_API_KEY",
        "openai": "OPENAI_API_KEY",
    }
    var = names.get(provider, "FGM_API_KEY")
    hint = (
        f"{var} is not set. Export it, e.g.\n"
        f"    export {var}=your_key_here\n"
        f"or add it to {ENV_FILE} (see agent/.env.example)."
    )
    if provider == "openai" and not base_url:
        hint += (
            "\nFor a LOCAL model instead, point it at your server and no key is "
            "needed:\n    FGM_BASE_URL=http://localhost:11434/v1"
        )
    elif provider != "openai":
        hint += "\nOr run a local model with no key at all:  --model local"
    return hint


def describe(config: Config) -> str:
    """One-line summary of where chat and embeddings are going (for the banner)."""
    chat = f"{config.chat_model} ({config.base_url or config.provider})"
    # Only mention embeddings when they come from somewhere else — which is
    # always the case for Claude, since the Anthropic API has no embeddings.
    same_source = (
        config.embed_provider == config.provider
        and config.embed_base_url == config.base_url
    )
    if same_source:
        return chat
    embed_where = config.embed_base_url or config.embed_provider
    return f"{chat}   ·   embeddings: {config.embed_model} ({embed_where})"
