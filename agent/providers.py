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

"""Provider and model presets — the friendly names used to switch models.

A *preset* is a short name (``claude``, ``mistral``, ``local``, ``openai``) that
expands into a complete provider configuration: which API to talk to, which chat
model, which endpoint, and which embedding model to pair with it. Presets are the
user-facing way to switch models:

    agent/.venv/bin/python agent/main.py --model claude
    > /model local-3b                     # live, inside the REPL
    > /models                             # list what is available

Anything a preset sets can still be overridden by the individual ``FGM_*``
environment variables (see :mod:`config`), and a raw model id is always accepted
in place of a preset name, so new models work without editing this table:

    /model claude:claude-haiku-4-5        # preset + explicit model id
    /model qwen2.5-3b-16k                  # raw id, keep the current provider

Embeddings are configured separately from chat on purpose. The Anthropic API has
no embeddings endpoint, and re-embedding the corpus every time you switch chat
models would be slow and pointless, so each preset only supplies a *fallback*
embedding configuration used when nothing else pins one.
"""

from __future__ import annotations

import difflib
import re
from dataclasses import dataclass
from typing import Optional

# Chat providers the agent can talk to.
#   "mistral"   — the hosted Mistral API (mistralai SDK)
#   "anthropic" — the Claude API (anthropic SDK)
#   "openai"    — the OpenAI API *or* any OpenAI-compatible server (Ollama, vLLM,
#                 LM Studio, ...), distinguished only by whether base_url is set
CHAT_PROVIDERS: frozenset[str] = frozenset({"mistral", "anthropic", "openai"})

# Providers that expose an embeddings endpoint. Anthropic deliberately does not,
# so a Claude chat preset must borrow embeddings from one of these.
EMBED_PROVIDERS: frozenset[str] = frozenset({"mistral", "openai"})

# The local OpenAI-compatible endpoint used by the "local" presets.
OLLAMA_BASE_URL: str = "http://localhost:11434/v1"

# The local presets point at Ollama tags built from an upstream open-weights
# model with a larger context window than Ollama's default — the agent's prompt
# plus the retrieved source chunks do not fit in 4096 tokens, and num_ctx cannot
# be set over the OpenAI-compatible endpoint, so it has to be baked into a tag.
# The tags are named after the model they wrap (Qwen2.5, by Alibaba) and the
# context size; see the README for the one-time `ollama create` step. Licences
# differ per size: the 7B is Apache-2.0, the 3B is under the Qwen Research
# licence (non-commercial).


@dataclass(frozen=True)
class Preset:
    """A named, ready-to-use provider configuration."""

    name: str
    provider: str
    chat_model: str
    summary: str
    base_url: Optional[str] = None
    # Fallback embedding configuration (used only when nothing else pins one).
    embed_provider: str = "mistral"
    embed_model: str = "mistral-embed"
    embed_base_url: Optional[str] = None
    # Environment variable holding this provider's API key ("" = no key needed).
    key_env: str = ""


# Ordered so the listing reads hosted-first, local-last.
PRESETS: tuple[Preset, ...] = (
    Preset(
        name="mistral",
        provider="mistral",
        # Pinned rather than mistral-large-latest, which still points at Large 3.
        chat_model="mistral-large-4",
        summary="Mistral API — Large 4, \"Le Chonk\" (the default)",
        key_env="MISTRAL_API_KEY",
    ),
    Preset(
        name="mistral-large-3",
        provider="mistral",
        chat_model="mistral-large-2512",
        summary="Mistral API — Large 3 (cheaper)",
        key_env="MISTRAL_API_KEY",
    ),
    Preset(
        name="mistral-small",
        provider="mistral",
        chat_model="mistral-small-latest",
        summary="Mistral API — small (cheapest)",
        key_env="MISTRAL_API_KEY",
    ),
    Preset(
        name="claude",
        provider="anthropic",
        chat_model="claude-opus-5",
        summary="Claude API — Opus 5 (most capable)",
        # Claude has no embeddings endpoint: default to the local embedder so a
        # Claude session needs exactly one API key and no second paid service.
        embed_provider="openai",
        embed_model="nomic-embed-text",
        embed_base_url=OLLAMA_BASE_URL,
        key_env="ANTHROPIC_API_KEY",
    ),
    Preset(
        name="claude-sonnet",
        provider="anthropic",
        chat_model="claude-sonnet-5",
        summary="Claude API — Sonnet 5 (faster, cheaper)",
        embed_provider="openai",
        embed_model="nomic-embed-text",
        embed_base_url=OLLAMA_BASE_URL,
        key_env="ANTHROPIC_API_KEY",
    ),
    Preset(
        name="claude-haiku",
        provider="anthropic",
        chat_model="claude-haiku-4-5",
        summary="Claude API — Haiku 4.5 (fastest, cheapest)",
        embed_provider="openai",
        embed_model="nomic-embed-text",
        embed_base_url=OLLAMA_BASE_URL,
        key_env="ANTHROPIC_API_KEY",
    ),
    Preset(
        name="openai",
        provider="openai",
        chat_model="gpt-4o",
        summary="OpenAI API — GPT-4o",
        embed_provider="openai",
        embed_model="text-embedding-3-small",
        key_env="OPENAI_API_KEY",
    ),
    Preset(
        name="openai-mini",
        provider="openai",
        chat_model="gpt-4o-mini",
        summary="OpenAI API — GPT-4o mini (cheap)",
        embed_provider="openai",
        embed_model="text-embedding-3-small",
        key_env="OPENAI_API_KEY",
    ),
    Preset(
        name="local",
        provider="openai",
        chat_model="qwen2.5-7b-16k",
        summary="Local Ollama — Qwen2.5 7B (nothing leaves the machine)",
        base_url=OLLAMA_BASE_URL,
        embed_provider="openai",
        embed_model="nomic-embed-text",
        embed_base_url=OLLAMA_BASE_URL,
    ),
    Preset(
        name="local-3b",
        provider="openai",
        chat_model="qwen2.5-3b-16k",
        summary="Local Ollama — Qwen2.5 3B (fits a small GPU)",
        base_url=OLLAMA_BASE_URL,
        embed_provider="openai",
        embed_model="nomic-embed-text",
        embed_base_url=OLLAMA_BASE_URL,
    ),
)

PRESET_BY_NAME: dict[str, Preset] = {p.name: p for p in PRESETS}

# Accepted spellings that are not listed separately in /models.
PRESET_ALIASES: dict[str, str] = {
    "mistral-large": "mistral",
    "large": "mistral",            # the default large model, not mistral-large-3
    "le-chonk": "mistral",
    "claude-opus": "claude",
    "ollama": "local",
}


def context_tokens(chat_model: str) -> Optional[int]:
    """Context window of a local re-tagged model, read from its tag.

    The local tags carry it in the name (``qwen2.5-3b-16k`` -> 16384), because
    Ollama cannot be told num_ctx over the OpenAI-compatible endpoint and silently
    drops the START of a prompt that does not fit: the system prompt and the
    question go first. Hosted models have windows far beyond this agent's
    prompts, so they return None (no budgeting needed).
    """
    match = re.search(r"-(\d+)k\b", chat_model or "")
    return int(match.group(1)) * 1024 if match else None


def find_preset(name: str) -> Optional[Preset]:
    """Return the preset called ``name`` (case-insensitive), or None."""
    key = name.strip().lower()
    return PRESET_BY_NAME.get(PRESET_ALIASES.get(key, key))


def split_choice(choice: str) -> tuple[Optional[Preset], Optional[str]]:
    """Interpret a user's model choice.

    Accepts three forms, in this order:

        "claude"                    -> (claude preset, None)
        "claude:claude-haiku-4-5"   -> (claude preset, "claude-haiku-4-5")
        "qwen2.5-3b-16k"             -> (None, "qwen2.5-3b-16k")

    The ``preset:model`` form is what keeps this table from going stale: a model
    released after this file was written is reachable without editing it.

    Returns:
        ``(preset, model)`` where either may be None. ``(None, model)`` means
        "keep the current provider, just change the model name".
    """
    text = choice.strip()
    if not text:
        return None, None
    head, sep, tail = text.partition(":")
    if sep:
        preset = find_preset(head)
        # Only treat "a:b" as preset:model when "a" really is a preset — this
        # leaves room for model ids that legitimately contain a colon.
        if preset is not None:
            return preset, (tail.strip() or None)
    preset = find_preset(text)
    if preset is not None:
        return preset, None
    preset = _preset_by_word(text)
    if preset is not None:
        return preset, None
    return None, text


def _preset_by_word(text: str) -> Optional[Preset]:
    """The one preset with ``text`` as a word of its name or model id ("3B" ->
    local-3b, "opus" -> claude), or None when no preset or several do."""
    key = text.strip().lower()
    matches = [p for p in PRESETS
               if key in p.name.split("-") + p.chat_model.lower().split("-")]
    return matches[0] if len(matches) == 1 else None


# What the hosted APIs' own model ids look like. An id that fits none of them is
# almost certainly a mistyped preset ("/model 3B" on the Claude API was accepted
# and only failed, with a 404, at the first question). OpenAI-compatible servers
# are not checked: local tags (Ollama, vLLM) are arbitrary names.
_MODEL_ID_SHAPES: dict[str, re.Pattern[str]] = {
    "anthropic": re.compile(r"^claude-[\w.-]+$", re.IGNORECASE),
    "mistral": re.compile(
        r"^(?:open-)?(?:mistral|ministral|magistral|mixtral|codestral|devstral|"
        r"pixtral|voxtral)[\w.-]*$", re.IGNORECASE),
}


def implausible_model(provider: str, model: str) -> Optional[str]:
    """Why ``model`` cannot be a model id of ``provider``, or None if it can be."""
    shape = _MODEL_ID_SHAPES.get(provider)
    if shape is None or shape.match(model.strip()):
        return None
    key = model.strip().lower()

    def likeness(preset: Preset) -> float:
        # Against the preset's name and its model id: "qwen" is local, not openai.
        return max(difflib.SequenceMatcher(None, key, text).ratio()
                   + (0.5 if key in text else 0.0)
                   for text in (preset.name, preset.chat_model.lower()))

    close = [p.name for p in sorted(PRESETS, key=likeness, reverse=True)[:2]
             if likeness(p) >= 0.4]
    hint = f"; did you mean {' or '.join(close)}?" if close else ""
    label = {"anthropic": "Claude", "mistral": "Mistral"}[provider]
    return f"'{model}' is not a preset or a {label} model id{hint} (/models lists the presets)"


def describe_presets(active_model: str = "", active_provider: str = "") -> str:
    """Render the preset table for the REPL's ``/models`` command."""
    width = max(len(p.name) for p in PRESETS)
    lines = ["Model presets  (switch with /model <name>):", ""]
    for preset in PRESETS:
        marker = (
            " *"
            if preset.chat_model == active_model and preset.provider == active_provider
            else "  "
        )
        lines.append(f"{marker}{preset.name.ljust(width)}  {preset.summary}")
    lines += [
        "",
        "  * = active.  Also accepted: /model <preset>:<model-id> (e.g. "
        "claude:claude-haiku-4-5)",
        "  or a bare model id to keep the current provider (e.g. /model "
        "mistral-large-latest).",
    ]
    return "\n".join(lines)
