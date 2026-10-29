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

"""Provider-agnostic LLM client for chat completions and embeddings.

Three chat providers, selected by ``config.provider``:

  - ``"mistral"``  : the hosted Mistral API.
  - ``"anthropic"``: the Claude API (Anthropic Messages API).
  - ``"openai"``   : the OpenAI API, or any OpenAI-compatible endpoint at
    ``config.base_url`` — e.g. a LOCAL server such as Ollama or vLLM, so no code
    or data leaves the machine.

Every adapter exposes the same small surface the agent needs:

  - ``complete(messages, model, temperature, tools, tool_choice)`` -> a response
    with the OpenAI shape (``.choices[0].message`` with ``.content`` and
    ``.tool_calls``), so ``chat.py`` parses any provider unchanged.
  - ``embed(texts, model)`` -> ``list[list[float]]``.
  - error-classification helpers (``is_retryable`` / ``is_response_validation_error``
    / ``backoff_seconds``) used by ``chat.py``'s retry loop, so all
    provider-specific exception handling lives here.

Mistral and OpenAI already share the OpenAI request/response shape. The Claude
API does not — it takes a separate ``system`` argument, tool schemas under
``input_schema``, and returns content *blocks* rather than a message with
``tool_calls`` — so :class:`_AnthropicClient` translates in both directions. The
translation lives in the module-level ``to_anthropic_*`` / ``from_anthropic_*``
helpers, which are pure functions and covered by ``agent/tests/test_llm.py``.

Embeddings are resolved independently of chat (``make_embed_client``): the Claude
API has no embeddings endpoint, so a Claude chat model borrows embeddings from
Mistral or from a local OpenAI-compatible embedder.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Optional

from config import Config, ConfigError

# HTTP statuses worth retrying (rate limit + transient server errors). 529 is
# Anthropic's "overloaded": capacity-side and short-lived, so retrying with
# backoff usually succeeds where failing the turn would not.
_RETRYABLE_STATUS: frozenset[int] = frozenset({429, 500, 502, 503, 504, 529})
_MAX_BACKOFF_S: float = 20.0

# The Anthropic API requires an explicit output cap. 16k is the documented
# default for non-streaming requests: large enough for any answer this agent
# produces, small enough to stay well inside the SDK's HTTP timeout.
ANTHROPIC_MAX_TOKENS: int = 16000


def _status_of(exc: Exception) -> Any:
    """Best-effort HTTP status code from any SDK's exception shape."""
    return (
        getattr(exc, "status_code", None)
        or getattr(getattr(exc, "raw_response", None), "status_code", None)
        or getattr(getattr(exc, "response", None), "status_code", None)
    )


def _backoff_seconds(exc: Exception, attempt: int) -> float:
    """Honour a Retry-After header if present, else exponential backoff."""
    for holder in ("raw_response", "response"):
        headers = getattr(getattr(exc, holder, None), "headers", None)
        if headers:
            try:
                value = headers.get("retry-after") or headers.get("Retry-After")
            except AttributeError:
                value = None
            if value:
                try:
                    return min(float(value), _MAX_BACKOFF_S)
                except (TypeError, ValueError):
                    pass
    return min(2.0 * (2 ** attempt), _MAX_BACKOFF_S)


# --------------------------------------------------------------------------- #
# OpenAI-shaped response objects, so one parser in chat.py serves all providers
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _Function:
    name: str
    arguments: str  # JSON-encoded, as the OpenAI/Mistral SDKs return it


@dataclass(frozen=True)
class _ToolCall:
    id: str
    function: _Function
    type: str = "function"


@dataclass(frozen=True)
class _Message:
    content: str
    tool_calls: list[_ToolCall] = field(default_factory=list)


@dataclass(frozen=True)
class _Choice:
    message: _Message


@dataclass(frozen=True)
class _Response:
    choices: list[_Choice]


# --------------------------------------------------------------------------- #
# OpenAI  <->  Anthropic translation (pure functions, unit-tested)
# --------------------------------------------------------------------------- #


def to_anthropic_tools(tools: Optional[list[dict[str, Any]]]) -> list[dict[str, Any]]:
    """Convert OpenAI ``{"type":"function","function":{...}}`` tool definitions.

    The Anthropic equivalent is flat, and the JSON Schema lives under
    ``input_schema`` rather than ``parameters``.
    """
    converted: list[dict[str, Any]] = []
    for tool in tools or []:
        fn = tool.get("function", tool)
        converted.append({
            "name": fn["name"],
            "description": fn.get("description", ""),
            "input_schema": fn.get("parameters", {"type": "object", "properties": {}}),
        })
    return converted


def to_anthropic_tool_choice(tool_choice: Any) -> Optional[dict[str, Any]]:
    """Convert an OpenAI tool_choice to the Anthropic form.

    ``"auto"`` -> ``{"type": "auto"}``; ``"any"``/``"required"`` -> forced tool
    use; a named OpenAI choice -> ``{"type": "tool", "name": ...}``.
    """
    if tool_choice is None:
        return None
    if isinstance(tool_choice, str):
        if tool_choice in ("any", "required"):
            return {"type": "any"}
        if tool_choice == "none":
            return {"type": "none"}
        return {"type": "auto"}
    if isinstance(tool_choice, dict):
        name = tool_choice.get("function", {}).get("name") or tool_choice.get("name")
        if name:
            return {"type": "tool", "name": name}
        kind = tool_choice.get("type")
        if kind in ("auto", "any", "none"):
            return {"type": kind}
    return {"type": "auto"}


def to_anthropic_messages(
    messages: list[dict[str, Any]],
) -> tuple[str, list[dict[str, Any]]]:
    """Split an OpenAI message list into ``(system_prompt, anthropic_messages)``.

    Four shape differences are handled here:

    * ``system`` messages are not part of the conversation — they are hoisted
      into the top-level ``system`` argument (joined, so a mid-conversation
      nudge still applies).
    * an assistant turn carrying ``tool_calls`` becomes content *blocks*: an
      optional ``text`` block plus one ``tool_use`` block per call.
    * ``tool`` results come back as a ``user`` message of ``tool_result`` blocks,
      and **consecutive results are merged into one message** — returning
      parallel tool results as separate messages trains the model to stop
      issuing parallel calls.
    * empty assistant text is dropped: the API rejects empty text blocks.
    """
    system_parts: list[str] = []
    out: list[dict[str, Any]] = []
    pending_results: list[dict[str, Any]] = []

    def flush_results() -> None:
        if pending_results:
            out.append({"role": "user", "content": list(pending_results)})
            pending_results.clear()

    for message in messages:
        role = message.get("role")
        content = message.get("content") or ""

        if role == "tool":
            pending_results.append({
                "type": "tool_result",
                "tool_use_id": message.get("tool_call_id", ""),
                "content": content if isinstance(content, str) else json.dumps(content),
            })
            continue

        flush_results()

        if role == "system":
            if content:
                system_parts.append(content)
        elif role == "assistant":
            blocks: list[dict[str, Any]] = []
            if isinstance(content, str) and content.strip():
                blocks.append({"type": "text", "text": content})
            for call in message.get("tool_calls") or []:
                fn = call["function"]
                try:
                    args = json.loads(fn["arguments"] or "{}")
                except (TypeError, ValueError):
                    args = {}
                blocks.append({
                    "type": "tool_use",
                    "id": call["id"],
                    "name": fn["name"],
                    "input": args,
                })
            if blocks:  # an assistant turn with no content at all is not sendable
                out.append({"role": "assistant", "content": blocks})
        else:  # "user"
            out.append({"role": "user", "content": content})

    flush_results()
    return "\n\n".join(system_parts), out


def from_anthropic_response(response: Any) -> _Response:
    """Wrap an Anthropic ``Message`` in the OpenAI response shape.

    Text blocks are concatenated into ``.content``; ``tool_use`` blocks become
    ``.tool_calls`` with JSON-encoded arguments, which is what ``chat.py``
    expects and what round-trips back through :func:`to_anthropic_messages`.
    """
    texts: list[str] = []
    calls: list[_ToolCall] = []
    for block in getattr(response, "content", []) or []:
        kind = getattr(block, "type", None)
        if kind == "text":
            texts.append(getattr(block, "text", ""))
        elif kind == "tool_use":
            calls.append(_ToolCall(
                id=getattr(block, "id", ""),
                function=_Function(
                    name=getattr(block, "name", ""),
                    arguments=json.dumps(getattr(block, "input", {}) or {}),
                ),
            ))
    return _Response(choices=[_Choice(_Message(content="".join(texts), tool_calls=calls))])


# --------------------------------------------------------------------------- #
# Adapters
# --------------------------------------------------------------------------- #


class _MistralClient:
    """Adapter for the hosted Mistral API."""

    provider = "mistral"

    def __init__(self, api_key: str, base_url: Optional[str] = None) -> None:
        try:
            from mistralai import Mistral
        except ImportError as exc:
            raise ConfigError(
                "The Mistral provider needs the 'mistralai' package:\n"
                "    agent/.venv/bin/pip install -r agent/requirements.txt"
            ) from exc
        self._sdk = Mistral(api_key=api_key)

    def complete(
        self, *, messages, model, temperature, tools=None, tool_choice=None
    ) -> Any:
        kwargs: dict[str, Any] = {
            "model": model, "messages": messages, "temperature": temperature,
        }
        if tools is not None:
            kwargs["tools"] = tools
            if tool_choice is not None:
                kwargs["tool_choice"] = tool_choice
        return self._sdk.chat.complete(**kwargs)

    def embed(self, *, texts, model) -> list[list[float]]:
        resp = self._sdk.embeddings.create(model=model, inputs=texts)
        return [item.embedding for item in resp.data]

    def is_response_validation_error(self, exc: Exception) -> bool:
        # mistral-large occasionally returns citations as structured chunks the
        # SDK cannot parse; chat.py retries those with a plain-text nudge.
        from mistralai.models import ResponseValidationError

        return isinstance(exc, ResponseValidationError)

    def is_retryable(self, exc: Exception) -> bool:
        from mistralai.models import SDKError

        if _status_of(exc) in _RETRYABLE_STATUS:
            return True
        if isinstance(exc, SDKError):
            text = str(exc).lower()
            return "rate limit" in text or "429" in text
        return False

    def backoff_seconds(self, exc: Exception, attempt: int) -> float:
        return _backoff_seconds(exc, attempt)


class _OpenAICompatibleClient:
    """Adapter for the OpenAI API or any OpenAI-compatible endpoint.

    A ``base_url`` pointing at localhost is how the agent runs fully offline
    (Ollama, vLLM, LM Studio); with no ``base_url`` this is the hosted OpenAI API.
    """

    provider = "openai"

    def __init__(self, api_key: str, base_url: Optional[str] = None) -> None:
        try:
            from openai import OpenAI
        except ImportError as exc:
            raise ConfigError(
                "The OpenAI-compatible provider needs the 'openai' package:\n"
                "    agent/.venv/bin/pip install -r agent/requirements.txt"
            ) from exc
        # Local servers ignore the key, but the SDK requires a non-empty string.
        self._sdk = OpenAI(base_url=base_url, api_key=api_key or "local")

    def complete(
        self, *, messages, model, temperature, tools=None, tool_choice=None
    ) -> Any:
        kwargs: dict[str, Any] = {
            "model": model, "messages": messages, "temperature": temperature,
        }
        if tools is not None:
            kwargs["tools"] = tools
            if tool_choice is not None:
                kwargs["tool_choice"] = tool_choice
        return self._sdk.chat.completions.create(**kwargs)

    def embed(self, *, texts, model) -> list[list[float]]:
        resp = self._sdk.embeddings.create(model=model, input=texts)
        return [item.embedding for item in resp.data]

    def is_response_validation_error(self, exc: Exception) -> bool:
        return False  # OpenAI-shaped responses parse directly.

    def is_retryable(self, exc: Exception) -> bool:
        try:
            from openai import APIConnectionError, RateLimitError
        except ImportError:
            return False
        if isinstance(exc, (APIConnectionError, RateLimitError)):
            return True
        return _status_of(exc) in _RETRYABLE_STATUS

    def backoff_seconds(self, exc: Exception, attempt: int) -> float:
        return _backoff_seconds(exc, attempt)


class _AnthropicClient:
    """Adapter for the Claude API (Anthropic Messages API).

    Translates the agent's OpenAI-shaped messages/tools into the Messages API
    shape on the way in, and the returned content blocks back on the way out, so
    the rest of the agent is provider-agnostic.
    """

    provider = "anthropic"

    def __init__(self, api_key: str, base_url: Optional[str] = None) -> None:
        try:
            import anthropic
        except ImportError as exc:
            raise ConfigError(
                "The Claude provider needs the 'anthropic' package:\n"
                "    agent/.venv/bin/pip install -r agent/requirements.txt"
            ) from exc
        self._anthropic = anthropic
        self._sdk = anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()

    def complete(
        self, *, messages, model, temperature, tools=None, tool_choice=None
    ) -> _Response:
        # `temperature` is accepted for interface parity with the other adapters
        # and intentionally unused — see the note below.
        del temperature
        system, converted = to_anthropic_messages(messages)
        kwargs: dict[str, Any] = {
            "model": model,
            "max_tokens": ANTHROPIC_MAX_TOKENS,
            "messages": converted,
        }
        if system:
            kwargs["system"] = system
        # ``temperature`` is deliberately not forwarded: the Anthropic SDK's
        # messages.create() has no sampling parameters at all (verified against
        # anthropic 1.4.0 — passing one raises TypeError before any request is
        # made). FGM_TEMPERATURE therefore applies to the other providers only.
        if tools:
            kwargs["tools"] = to_anthropic_tools(tools)
            choice = to_anthropic_tool_choice(tool_choice)
            if choice is not None:
                kwargs["tool_choice"] = choice
        return from_anthropic_response(self._sdk.messages.create(**kwargs))

    def embed(self, *, texts, model) -> list[list[float]]:
        raise ConfigError(
            "The Anthropic API has no embeddings endpoint, so Claude cannot build "
            "the knowledge index. Point the embedder somewhere else, e.g. locally "
            "with no extra key:\n"
            "    FGM_EMBED_PROVIDER=openai\n"
            "    FGM_EMBED_BASE_URL=http://localhost:11434/v1\n"
            "    FGM_EMBED_MODEL=nomic-embed-text\n"
            "or at Mistral with MISTRAL_API_KEY set:\n"
            "    FGM_EMBED_PROVIDER=mistral\n"
            "    FGM_EMBED_MODEL=mistral-embed"
        )

    def is_response_validation_error(self, exc: Exception) -> bool:
        return False  # Content blocks are parsed by the SDK into typed objects.

    def is_retryable(self, exc: Exception) -> bool:
        anthropic = self._anthropic
        if isinstance(exc, (anthropic.RateLimitError, anthropic.APIConnectionError,
                            anthropic.APITimeoutError)):
            return True
        if isinstance(exc, anthropic.InternalServerError):
            return True
        if _status_of(exc) in _RETRYABLE_STATUS:
            return True
        # 529 "overloaded" is not always mapped to a typed exception or exposed
        # as a status attribute, so fall back to the message text.
        return "overloaded" in str(exc).lower()

    def backoff_seconds(self, exc: Exception, attempt: int) -> float:
        return _backoff_seconds(exc, attempt)


_ADAPTERS = {
    "mistral": _MistralClient,
    "openai": _OpenAICompatibleClient,
    "anthropic": _AnthropicClient,
}


def make_client(config: Config):
    """Return the chat adapter for the configured chat provider."""
    return _ADAPTERS[config.provider](config.api_key, config.base_url)


def make_embed_client(config: Config):
    """Return the adapter for the configured *embedding* provider.

    Separate from :func:`make_client` because the embedding provider can differ
    from the chat provider — always so for Claude, which has no embeddings API.
    """
    return _ADAPTERS[config.embed_provider](config.embed_api_key, config.embed_base_url)
