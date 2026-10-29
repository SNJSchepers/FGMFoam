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

"""Tests for the OpenAI <-> Anthropic translation in ``llm.py``.

The agent speaks one internal dialect (the OpenAI chat-completions shape) and
``_AnthropicClient`` translates it to and from the Anthropic Messages API. That
translation is the one piece of provider code that cannot be checked by simply
running the agent against a local model, so it is covered here.

Runs with pytest, or standalone with no test dependency at all:

    agent/.venv/bin/python agent/tests/test_llm.py
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from llm import (  # noqa: E402
    from_anthropic_response,
    to_anthropic_messages,
    to_anthropic_tool_choice,
    to_anthropic_tools,
)


# --- stand-ins for the Anthropic SDK's response blocks ---------------------- #


@dataclass
class _TextBlock:
    text: str
    type: str = "text"


@dataclass
class _ToolUseBlock:
    id: str
    name: str
    input: dict[str, Any]
    type: str = "tool_use"


@dataclass
class _FakeMessage:
    content: list[Any]


def test_system_messages_are_hoisted() -> None:
    system, messages = to_anthropic_messages([
        {"role": "system", "content": "You are the FGMFoam agent."},
        {"role": "user", "content": "what is Yc?"},
    ])
    assert system == "You are the FGMFoam agent."
    assert messages == [{"role": "user", "content": "what is Yc?"}]


def test_multiple_system_messages_are_joined() -> None:
    # chat.py appends a second system message to nudge a model into plain text.
    system, messages = to_anthropic_messages([
        {"role": "system", "content": "first"},
        {"role": "user", "content": "q"},
        {"role": "system", "content": "second"},
    ])
    assert system == "first\n\nsecond"
    assert [m["role"] for m in messages] == ["user"]


def test_assistant_tool_call_becomes_tool_use_block() -> None:
    _, messages = to_anthropic_messages([
        {"role": "user", "content": "read the lookup"},
        {
            "role": "assistant",
            "content": "",
            "tool_calls": [{
                "id": "call_1",
                "type": "function",
                "function": {
                    "name": "read_case_file",
                    "arguments": json.dumps({"path": "solver/src/lookUp/FGM/FGMlib.c"}),
                },
            }],
        },
        {"role": "tool", "tool_call_id": "call_1", "name": "read_case_file",
         "content": '{"status": "ok"}'},
    ])
    assistant = messages[1]
    assert assistant["role"] == "assistant"
    # The empty assistant text must NOT become a text block (the API rejects one).
    assert [b["type"] for b in assistant["content"]] == ["tool_use"]
    block = assistant["content"][0]
    assert block["id"] == "call_1"
    assert block["name"] == "read_case_file"
    assert block["input"] == {"path": "solver/src/lookUp/FGM/FGMlib.c"}

    result = messages[2]
    assert result["role"] == "user"
    assert result["content"][0]["type"] == "tool_result"
    assert result["content"][0]["tool_use_id"] == "call_1"


def test_assistant_text_and_tool_call_keeps_both_blocks() -> None:
    _, messages = to_anthropic_messages([
        {"role": "user", "content": "go"},
        {
            "role": "assistant",
            "content": "Reading the source first.",
            "tool_calls": [{
                "id": "c1", "type": "function",
                "function": {"name": "read_case_file", "arguments": "{}"},
            }],
        },
    ])
    assert [b["type"] for b in messages[1]["content"]] == ["text", "tool_use"]


def test_parallel_tool_results_merge_into_one_user_message() -> None:
    """Separate messages per result would train the model out of parallel calls."""
    _, messages = to_anthropic_messages([
        {"role": "user", "content": "go"},
        {
            "role": "assistant",
            "content": "",
            "tool_calls": [
                {"id": "a", "type": "function",
                 "function": {"name": "run_status", "arguments": "{}"}},
                {"id": "b", "type": "function",
                 "function": {"name": "read_case_file", "arguments": "{}"}},
            ],
        },
        {"role": "tool", "tool_call_id": "a", "content": "one"},
        {"role": "tool", "tool_call_id": "b", "content": "two"},
    ])
    assert len(messages) == 3, "the two tool results must share one user message"
    blocks = messages[2]["content"]
    assert [b["tool_use_id"] for b in blocks] == ["a", "b"]


def test_malformed_tool_arguments_do_not_raise() -> None:
    _, messages = to_anthropic_messages([
        {"role": "user", "content": "go"},
        {"role": "assistant", "content": "", "tool_calls": [{
            "id": "c1", "type": "function",
            "function": {"name": "build_case", "arguments": "{not json"},
        }]},
    ])
    assert messages[1]["content"][0]["input"] == {}


def test_tools_are_converted_to_input_schema() -> None:
    schema = {"type": "object", "properties": {"path": {"type": "string"}}}
    converted = to_anthropic_tools([{
        "type": "function",
        "function": {"name": "read_case_file", "description": "read", "parameters": schema},
    }])
    assert converted == [
        {"name": "read_case_file", "description": "read", "input_schema": schema}
    ]


def test_tool_choice_conversion() -> None:
    assert to_anthropic_tool_choice(None) is None
    assert to_anthropic_tool_choice("auto") == {"type": "auto"}
    assert to_anthropic_tool_choice("any") == {"type": "any"}
    assert to_anthropic_tool_choice("required") == {"type": "any"}
    assert to_anthropic_tool_choice(
        {"type": "function", "function": {"name": "read_case_file"}}
    ) == {"type": "tool", "name": "read_case_file"}


def test_response_is_wrapped_in_the_openai_shape() -> None:
    wrapped = from_anthropic_response(_FakeMessage(content=[
        _TextBlock(text="The progress variable is "),
        _TextBlock(text="Yc = Y_CO2 / M_CO2."),
    ]))
    message = wrapped.choices[0].message
    assert message.content == "The progress variable is Yc = Y_CO2 / M_CO2."
    assert not message.tool_calls  # chat.py branches on falsiness here


def test_response_tool_calls_round_trip_back_to_anthropic() -> None:
    """A tool call must survive: Anthropic -> chat.py history -> Anthropic."""
    wrapped = from_anthropic_response(_FakeMessage(content=[
        _ToolUseBlock(id="tu_1", name="build_case", input={"T_inlet": 300.0}),
    ]))
    call = wrapped.choices[0].message.tool_calls[0]
    assert call.id == "tu_1"
    assert call.function.name == "build_case"
    assert json.loads(call.function.arguments) == {"T_inlet": 300.0}

    # Now feed it back through the history shape chat.py stores.
    _, messages = to_anthropic_messages([
        {"role": "user", "content": "build it"},
        {"role": "assistant", "content": wrapped.choices[0].message.content,
         "tool_calls": [{"id": call.id, "type": "function",
                         "function": {"name": call.function.name,
                                      "arguments": call.function.arguments}}]},
    ])
    assert messages[1]["content"][0]["input"] == {"T_inlet": 300.0}


def test_sdk_has_no_sampling_parameters() -> None:
    """Guards the reason temperature is never forwarded to Claude.

    ``messages.create()`` accepts no sampling parameters, for any model, so
    passing one raises TypeError before a request is even made. If a future SDK
    reintroduces them this test fails and the omission can be revisited.
    """
    import inspect

    import anthropic

    params = inspect.signature(anthropic.Anthropic(api_key="x").messages.create).parameters
    assert not {"temperature", "top_p", "top_k"} & set(params)


def test_request_built_for_the_anthropic_sdk() -> None:
    """The outgoing request must match what the Messages API expects.

    Builds a real ``_AnthropicClient`` and swaps in a recorder for the SDK, so
    this checks the request shape without any network call or API key.
    """
    from llm import ANTHROPIC_MAX_TOKENS, _AnthropicClient

    class _Recorder:
        def __init__(self) -> None:
            self.kwargs: dict[str, Any] = {}

        class _Messages:
            def __init__(self, outer: "_Recorder") -> None:
                self._outer = outer

            def create(self, **kwargs: Any) -> _FakeMessage:
                self._outer.kwargs = kwargs
                return _FakeMessage(content=[_TextBlock(text="ok")])

        @property
        def messages(self) -> "_Recorder._Messages":
            return _Recorder._Messages(self)

    client = _AnthropicClient.__new__(_AnthropicClient)  # skip SDK construction
    recorder = _Recorder()
    client._sdk = recorder

    reply = client.complete(
        messages=[
            {"role": "system", "content": "sys"},
            {"role": "user", "content": "hi"},
        ],
        model="claude-opus-5",
        temperature=0.2,
        tools=[{"type": "function", "function": {
            "name": "read_case_file", "description": "d",
            "parameters": {"type": "object", "properties": {}}}}],
        tool_choice="auto",
    )

    sent = recorder.kwargs
    assert sent["model"] == "claude-opus-5"
    assert sent["max_tokens"] == ANTHROPIC_MAX_TOKENS  # required by the API
    assert sent["system"] == "sys"                     # not a chat message
    assert sent["messages"] == [{"role": "user", "content": "hi"}]
    assert sent["tools"][0]["input_schema"] == {"type": "object", "properties": {}}
    assert sent["tool_choice"] == {"type": "auto"}
    # The SDK has no sampling parameters at all, so temperature must be omitted.
    assert "temperature" not in sent
    assert reply.choices[0].message.content == "ok"


def test_temperature_is_never_forwarded_for_any_model() -> None:
    from llm import _AnthropicClient

    class _Recorder:
        kwargs: dict[str, Any] = {}

        class messages:  # noqa: N801 - mimics the SDK's attribute layout
            @staticmethod
            def create(**kwargs: Any) -> _FakeMessage:
                _Recorder.kwargs = kwargs
                return _FakeMessage(content=[_TextBlock(text="ok")])

    for model in ("claude-haiku-4-5", "claude-opus-5", "claude-sonnet-5"):
        client = _AnthropicClient.__new__(_AnthropicClient)
        client._sdk = _Recorder()
        client.complete(messages=[{"role": "user", "content": "hi"}],
                        model=model, temperature=0.2)
        assert "temperature" not in _Recorder.kwargs, model


def test_overloaded_is_retryable() -> None:
    """529 "overloaded" is transient capacity — back off rather than fail the turn."""
    from llm import _RETRYABLE_STATUS, _AnthropicClient

    import anthropic

    assert 529 in _RETRYABLE_STATUS

    client = _AnthropicClient.__new__(_AnthropicClient)
    client._anthropic = anthropic

    class _Overloaded(Exception):
        status_code = 529

    assert client.is_retryable(_Overloaded())
    # Even with no status attribute, the message alone is enough.
    assert client.is_retryable(Exception("Error code: 529 - ... 'Overloaded' ..."))
    # A genuine client error must still fail fast.
    class _BadRequest(Exception):
        status_code = 400

    assert not client.is_retryable(_BadRequest())


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
