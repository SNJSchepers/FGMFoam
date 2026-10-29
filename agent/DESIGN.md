# Design: Conversational FGMFoam agent with RAG grounding

Status: **implemented and in use.** This document is the design rationale behind the agent in `agent/`; see `agent/README.md` for installation and usage. Where the two disagree, the README describes what the code does today and wins.

This document describes turning the current one-shot parameter-extraction agent
into a conversational assistant that (a) answers questions about the FGMFoam
solver grounded in its actual source and docs, and (b) still builds the OpenFOAM
case deterministically. Decisions locked in discussion: index **solver code +
repo docs**, embeddings via the **Mistral API** (`mistral-embed`).

> **Since implementation:** the agent is no longer tied to one provider. Chat runs
> on Mistral, the Claude API, the OpenAI API, or any local OpenAI-compatible
> server, selected by a named preset (`--model claude`, `/model local-3b`); see
> `agent/providers.py` and the README's "Choosing a model" section. Embeddings
> are configured independently of chat — the Anthropic API has no embeddings
> endpoint, and pinning the embedder keeps the index valid across chat-model
> switches — so `mistral-embed` is now one option among several rather than the
> fixed choice recorded above.

---

## 1. Goal

One agent that can do three things in a single conversation:

1. **Answer questions** about the solver and how to use it, grounded in the real
   code/docs (RAG), with citations — e.g. correctly answer "can it run a
   turbulent case?" by reading `turbulenceProperties` (`simulationType laminar`)
   and the solver source, instead of guessing.
2. **Gather parameters** through natural back-and-forth (no rigid question loop).
3. **Build the case** by calling the existing deterministic pipeline
   (`validator` → `cantera_tool` → `of_writer`) as a tool.

The deterministic case-building stays in Python and is unchanged; the LLM only
*decides when* to call it.

## 2. Why this shape

- The current pipeline never lets the model *respond* — it only extracts JSON
  and composes follow-up questions, so general/solver questions are dropped.
- A small chat model answers general CFD questions from training, but
  **hallucinates on solver specifics** it has never seen. RAG over the actual
  code/docs is what makes answers about *this* solver trustworthy.
- The corpus is tiny (≈140 KB), so retrieval can be a flat in-memory cosine
  search — **no vector database, no LangChain/LlamaIndex** (consistent with the
  framework-free style of the existing code).

## 3. Corpus and indexing strategy

### Sources (from repo inventory)

| Source | What | Approx size |
|--------|------|-------------|
| `solver/applications/FGMFoam/**/*.{C,H}` | the solver application | ~0.9k lines |
| `solver/src/**/*.{C,H,c,h}` | the 3 libraries (combustionModel, lookUp incl. the lowercase `FGMlib.c` lookup algorithm, lookUpBoundaryCondition) | ~3.3k lines |
| `solver/src/combustionModel/README` | combustion model notes | small |
| `README.md` | top-level project/usage docs | ~11 KB |
| `testcases/*/` (both `laminarSlotBurner` and `laminarFlamelet`) | the full case directories: `0.orig/*` fields, `system.orig/*` dicts (blockMesh, control, schemes, solution, setFields, decomposePar), `constant/*Properties`, `constant/polyMesh/boundary`, and the `run`/`clean` scripts | small (text only) |
| `agent/knowledge/*.md` | curated domain facts not stated in the code | small |
| `tables/database.fgm` (header only) | FGM table metadata: `[PRESSURE]`, `[DIMENSION]`, `[VARIABLES]` | ~35 lines |

**Excluded:** build artifacts (`*.o`, `*.dep`, `Make/`, `lnInclude/`),
`agent/.venv/`, the generated `laminarSlotBurner_agent/` case, the bulk numeric
data of the FGM tables (only each table's text header is indexed), the large
`polyMesh` geometry data
(`points`/`faces`/`owner`/`neighbour` — kept only the small `boundary` patch
list), logs (`log.*`), and editor cruft (`*.swp`, `.DS_Store`). The case
text files are included on purpose — they are how the agent answers "what
boundary conditions / models does the case use?" and how it contrasts the two
example cases.

### Chunking

- **Code (`.C`/`.H`):** sliding line windows (~40 lines, ~10 overlap). Simple and
  robust; a function-aware splitter is a possible later refinement but not needed
  for a corpus this size.
- **Markdown / README:** split on headings, with oversized sections
  re-windowed.
- **Case dicts:** one chunk per file (they are short).
- Every chunk carries metadata `{path, start_line, end_line, kind}` used for
  plain-text citations like `(solver/src/.../FGM.C:120-160)`.
- Estimated total: **~150–300 chunks** — negligible.

### Embeddings (`mistral-embed`)

- Embed all chunks at **index-build time** in batches (1024-dim vectors).
- One-time cost for the whole corpus is trivial; only changed files are
  re-embedded on rebuild (staleness by content hash).

### Vector store (local, file-based)

`agent/index/` holding:

- `embeddings.npy` — float32 matrix `(n_chunks, 1024)`
- `chunks.jsonl` — chunk text + metadata, row-aligned with the matrix
- `manifest.json` — per-source-file content hashes + model name, for staleness
  detection and reproducibility

Retrieval is **hybrid** (`VectorIndex.search_hybrid`): cosine similarity in NumPy
over the in-memory matrix, plus a small **lexical bonus** that rewards query terms
appearing in a chunk's **file path** (weight 0.35) or text (0.05). The dense-only
score misses exact-name matches that live in identifiers — e.g. "enthalpy boundary
condition" ranked `lookUpEnthalpyFvPatchScalarField.C` past rank 100 with
`nomic-embed-text`, so the forced read grabbed the unrelated `pEqn.H`; the path
bonus pulls the correctly-named source into the top-k. Both signals are computed
over the whole corpus (tiny, so this is instant). A CLI `python agent/index_build.py`
(re)builds the index; the chat entrypoint auto-builds if missing and warns if
sources changed.

## 4. Conversation + tool-calling orchestration

Single chat loop holding full message history. Per user turn:

```
user message
   │
   ├─ embed query (mistral-embed) ─▶ cosine top-k over index ─▶ context block
   │   (query = previous user turn + current message)          (chunks + [path:line])
   ▼
[ system prompt (role + grounding rules + param schema)
  + context block
  + conversation history
  + tool definitions ]                     ─▶  mistral chat model
   │
   ├─(a) grounded answer with citations            → shown to user
   ├─(b) asks for missing parameters in prose      → shown to user
   └─(c) tool call: build_case(...) / compute_*    → Python executes
                                                      → result returned to model
                                                      → model summarizes → user
```

**Tools (Mistral function-calling):**

- `build_case(equivalence_ratio, inlet_velocity, T_inlet, T_plate, sim_time,
  sim_timestep, n_procs?, sim_write_steps?, p_0?, case_name?)` — validates args,
  runs Cantera + writes the case via the existing modules, returns the summary.
  ``case_name`` chooses the output directory (auto-prefixed with
  ``laminarSlotBurner_``; default ``laminarSlotBurner_agent``). Python-side
  validation still enforces the required set; missing/invalid args come back as a
  structured error the model relays as a follow-up question.
- `run_case(case_name?, background?, timeout_s?)` — runs the solver pipeline
  (`runner.py`): blockMesh → setFields → decomposePar → mpirun FGMFoam -parallel
  → reconstructPar. By default the solver runs in the **background** (returns a
  log to tail); set `background=false` to wait. A non-zero FGMFoam exit counts as
  success if the run reached endTime (older FGMFoam builds aborted on shutdown). Build and run
  are separate steps — only run when the user explicitly asks.
- `run_status(case_name?)` — reports whether a (background) run is still going or
  finished, the latest simulated time, and the log tail.
- `change_settings(case_name?, sim_time?, sim_timestep?, sim_write_steps?)` — edits
  run-control keys in an already-built case's `system/controlDict` (`endTime`,
  `deltaT`, and `writeInterval = endTime / sim_write_steps`) WITHOUT rebuilding;
  the mesh, fields and thermochemistry are untouched
  (`of_writer.set_control_entries`). Values are validated in Python and take effect
  on the next run/restart. It never starts the solver.
- `restart_case(case_name?, timeout_s?)` — continues an already-run case from its
  latest saved time without re-meshing or re-initialising (`runner.plan_restart` +
  `runner.restart_pipeline_background`): reuses `constant/polyMesh`, runs
  `decomposePar -force -latestTime` only when the decomposition is missing or stale
  (e.g. the case was reconstructed; a case killed mid-run keeps current processors
  and is reused as-is), sets `startFrom latestTime`, then launches FGMFoam (+
  `reconstructPar`) in the background like `run_case`. Refuses if the case has not
  run (no time > 0), is currently running (stop it with `kill_jobs` first), or
  already reached `endTime` (raise it with `change_settings` first).
- `run_cases(case_names, max_parallel?)` — runs a sweep of built cases, sequential
  by default. `runner.plan_batch` enforces the physical-core budget
  (`runner.physical_core_count`): cases needing more cores than exist are refused,
  and concurrency is capped so concurrent runs never exceed the cores (the rest
  run sequentially). A detached `run_batch.py` orchestrates with that concurrency.
- `kill_jobs(case_names)` — stops/discards jobs. A running case: `runner.kill_case`
  finds the case's processes by working directory (`/proc/<pid>/cwd` under the case
  dir) and kills them individually — never the batch orchestrator or sibling cases.
  A queued case: a `.cancel` marker makes `run_batch.py` skip it when its turn comes.
- `delete_cases(case_names?, all?)` — removes agent-generated case directories
  (`of_writer.delete_case`/`generated_cases`). It can NEVER touch the template
  cases: `case_path` always prefixes the name with `laminarSlotBurner_`, so a
  request to delete `laminarSlotBurner`/`laminarFlamelet` resolves to a
  non-existent prefixed name, and `all` only enumerates `laminarSlotBurner_*` dirs.
- `read_case_file(path)` — returns the exact, verbatim contents of a solver/case/
  doc file (scoped to those roots; the FGM table returns its header only), so the
  model quotes ground truth for file-content questions instead of reconstructing
  from retrieval. The main mitigation for "describe-a-dictionary" hallucinations.
- `compute_thermochemistry(equivalence_ratio, T_inlet, p_0?)` — read-only; lets
  the agent answer "what's the equilibrium CO2 at φ=0.8?" without building a case.

Parameter gathering becomes **emergent**: the model asks for what's missing and
only calls `build_case` when ready — replacing the rigid `gather_parameters`
loop and `compose_followup`. Validation/determinism stay in Python.

**Grounding rules (system prompt):** only assert solver behaviour supported by
the retrieved context (or clearly-general OpenFOAM knowledge); cite sources as
`[path:line]`; if the context doesn't cover it, say so and point to the relevant
file rather than inventing. If top similarity is below a threshold, the context
is marked "no strong match" so the model doesn't over-anchor on irrelevant
chunks. The prompt also scopes the agent's own capabilities honestly: it must not
claim to read, analyse, plot, or interpret simulation *results* (there is no tool
for the numeric field data) — it can only report run status and the `log.FGMFoam`
tail (`run_status`) and should point users to ParaView / OpenFOAM post-processing
for results analysis.

**Forced source reads (deterministic):** prompt rules alone do not reliably make
the model open a file, so when the user asks how something is *implemented/works*
(keyword heuristic, `_wants_implementation`) **and** a source file is among the
retrieved hits, the agent **reads the most relevant source file itself, in Python,
and injects it into the conversation** before the model answers
(`ChatAgent._inject_forced_read`): it picks the highest-ranked C/C++ hit (and
prefers the sibling `.C/.c` implementation over a `.H/.h` header), calls
`read_case_file`, and records a well-formed assistant tool-call + tool-result so
the model answers from the actual source. This guarantees grounding for **any**
provider. Picking the *right* C/C++ hit depends on retrieval surfacing it: the
hybrid path-lexical bonus (§3) is what lifts an exactly-named source (e.g.
`lookUpEnthalpyFvPatchScalarField.C` for an "enthalpy boundary condition" question)
into the top-k so this forced read opens it instead of an unrelated `.H`. An
earlier version forced a *named* `tool_choice`
(`{"type":"function","function":{"name":"read_case_file"}}`) instead, but that was
**rejected**: hosted Mistral honours it, but local Ollama models treat `tool_choice`
as advisory (verified — with the full system prompt present, neither the 3B nor the
7B was forced to read), so on the local/data-safe path the read silently never
happened. Reading in Python removes that provider dependence. Tools stay available,
so the model can still read further files itself.

A purely retrieval-based gate ("force when the *top* retrieved chunk is C/C++")
was implemented and evaluated live, then rejected: the README and knowledge notes
describe the internals in prose and outrank the `.C/.H` source for these
questions, so source is rarely the top hit and the gate never fired for
"detail/explain/break down/walk me through the X"; conversely, "a source chunk
anywhere in the hits" alone over-fires, because a source chunk lands in nearly
every query's top-k (including setup/run requests, e.g. `correctPhi.H` matching
"phi 0.7"). The keyword intent signal AND-ed with "a source file is in the hits"
is what cleanly separates implementation questions from setup/run requests in this
doc-heavy corpus.

**History-aware retrieval (current-message force gate):** the index query is the
*previous user turn + the current message* (`ChatAgent._retrieval_query`), so a
context-dependent follow-up that names no topic ("I don't understand", "why?")
inherits the previous topic and still surfaces the right files instead of
retrieving noise. Only the immediately preceding user turn is prepended — more
history dilutes the embedding. The forced-read keyword gate, however, is judged on
the **current message only**, so a topic switch after an implementation question
("now set up a case at phi 0.7") is never falsely forced into a `read_case_file`,
even though its retrieval query still carries the prior "explain …" text. The
conversation sent to the model was already the full history; only retrieval and the
gate changed.

**Action-tool guards (code-enforced, provider-independent).** Two deterministic
checks in the tool loop stop a small model from acting when it should only answer
— they hold for any model because they live in Python, not the prompt:

- **Separate steps:** within one turn, once a *prepare* tool (`build_case` /
  `change_settings`) has run, a *run* tool (`run_case` / `run_cases` /
  `restart_case`) is deferred, not executed — the solver can never launch
  unprompted just because the model chained the calls. Running needs a new message.
- **Feasibility questions:** if the current message asks whether something is
  *possible / supported* (`_is_feasibility_question`, e.g. "is it possible to run
  at 10 bar?", "can the solver do X?"), **every** action tool is deferred for that
  turn and the model must answer in words. This was added after the 3B answered a
  10-bar feasibility question by calling `build_case` with *fabricated* parameters,
  building a case at the default pressure, and claiming it was 10 bar. The markers
  are info-seeking phrasings chosen NOT to fire on polite imperatives ("can you run
  it?", "can you set up a case …?"), which carry real intent and must still execute.

A second round of guards came from the eval harness (`agent/evals/`). Each one
fixes a failure the 3B actually showed:

- **Do not start:** "prepare it, but do not start it" (`_forbids_starting`)
  - defers every run tool in that turn;
  - forces `launch=false` on `run_campaign`.

  The 3B launched a "prepare only" campaign and then called `run_cases` with
  invented case names.
- **Case type must be asked for:** `build_case` / `run_campaign` are only run for a
  case type the conversation names, through its template name, a plugin
  `keywords` entry, or a parameter unique to it. Otherwise the model gets
  `missing_case_type` and asks. The 3B turned "a simulation at 350 K" into a slot
  burner on its own.
- **Claim check:** a reply that says a case was built, or a run started, when no
  tool did so this turn (`_claim_correction`, only for sentences about a case,
  simulation or campaign):
  - first gets one system nudge to actually do it, with the draft discarded;
  - if the model still claims it, a correction is appended.

  The 3B answered "rebuild it at 450 K" with only a thermochemistry call and
  replied that the case had been built.
- **Case type from the user's words:** when `build_case`/`run_campaign` arrive
  without `case_type` and the user's messages name exactly one case (current
  message first, then the conversation), it is filled in
  (`_case_type_from_user`). The 3B dropped it even after "laminar flamelet".
- **Swept variable from the user's words:** a `run_campaign` without `variable`
  gets it only from the message, in this order:
  1. the parameter whose phrase is followed by the value list ("sweeps the plate
     thickness over 2, 3 and 4 mm"), longest phrase winning;
  2. the one sweepable parameter the message names;
  3. the one whose unit the user's values carry ("300, 400 and 500 K" →
     `T_inlet`).

  Otherwise it stays unset and the tool asks.
- **Announced but not made:** a reply that says "I will call the … function" with
  no tool call this turn gets the same single nudge as a false claim.
- **Invented values:** a *required* parameter is accepted only if its value
  appears in what the user wrote, directly or after a unit conversion
  (`_invented_required`). Otherwise the call returns `missing_parameters`. After
  `build_case` asked for the temperatures, the 3B called it again in the same
  turn with `T_inlet=300` and `T_plate=300` of its own making.
- **Stated values:** a parameter phrase (`ParamSpec.phrases`, or the name with
  spaces), a few filler words, then a number with a unit that fits the
  parameter, as in "a plate thickness of 5 mm" or "the slot (hole size) 6 mm
  wide". The value is put into the `build_case` / `change_settings` call,
  replacing a value the model dropped or mangled (`_stated_parameters`). Phrases
  match case-insensitively, units case-sensitively; lists ("300, 400 and 500 K")
  are left to campaigns. Told "rebuild with a plate thickness of 5 mm", the 3B
  called `change_settings` with only the old run controls, nudges included.
- **CPU-plan answers:** when the previous turn ended in `missing_cpu_plan` and
  the message gives a plan ("2 processors per case, 2 at the same time"), the
  chat loop makes the pending `run_campaign` call with it, injected like the
  forced read, and the model only reports. `launch=false` still applies after
  "do not start". The 3B answered the plan with no call at all, 3 out of 3 times.
- **Requested values:** when an action tool ran but none of its arguments carries a
  value the user gave with a unit (`450 K`, `6 mm`, converted to SI), the model
  gets one nudge to apply it (`_unapplied_values`); if the value is still missing
  afterwards, the reply says it was not applied. The 3B dropped "450 K" by
  calling `change_settings` with the old run controls, and turned "6 mm" into
  `0.0006`.
- **Excerpts without retrieved ranges** (a file chosen through a note) show the
  run of lines densest in the question's words (`_densest_window`), not the top
  of the file. Otherwise the 3B was shown the includes and constructor of
  `FGM_PM_UL_NA.C` and answered about `update()` without having seen it.
- **When a small model repeats a wrong fact with the right file read, check the
  source first.** The 3B called the lookUpEnthalpy table "a neural network" in
  3 of 3 runs. The cause was the class description in
  `lookUpEnthalpyFvPatchScalarField.H` itself ("values obtained from a neural
  network", a leftover comment), not the model. The comment now describes the
  table inversion. The knowledge notes, which used to say "it is NOT a neural
  network", now state only what the lookup is.
- **Note-steered source reads:** if the best-ranked knowledge note that names a
  source file ranks above every code chunk, the forced read opens that file (if
  the note names several, the one sharing most words with the question). Code
  that is hard for a small embedder to match ("table lookup out of range" against
  bilinear interpolation code) is reached through a note that names it. That is
  how `FGMlib.c` is found for the out-of-range question.
- **Context budget for local models:**
  - The window is read from the tag (`providers.context_tokens`, e.g. `-16k`).
  - File reads are cut to an excerpt around the retrieved lines (`_fit_file_read`).
  - The oldest whole turns are dropped before the prompt overflows.

  Ollama silently drops the *start* of an over-long prompt, which is the system
  prompt. The agent's system prompt plus tool schemas alone measure about 10.3k
  tokens, so the old 8k `local-3b` tag lost its instructions on every turn. It
  is now `qwen2.5-3b-16k`.

In the tool layer (`tools.py`):

- Single-case tools default to the **case built last in this conversation**
  rather than `<template>_agent`. With no case named and none built, they ask.
- `build_case` reads the case type from a `<template>_` case-name prefix.
- A rebuild **inherits** the parameters the user chose for the previous case of
  that type and didn't restate. They are labelled `inherited` in the provenance
  record. A campaign's own builds neither inherit nor become the last case.
- A length above 0.1 m is refused as a probable mm-for-metres mistake, unless
  the user's own text gives that number in mm or cm. In that case the chat loop
  converts it first (`_convert_stated_units`, including "2, 3 and 4 mm").
- `build_case` refuses parameter values it cannot read. `coerce_types` drops
  those, which is right for a required value (it gets asked for again) but turned
  an optional one into its silent default: `plate_thickness=[2, 3, 4]` built a case
  with the default plate.
- `change_settings` given a condition or geometry parameter rebuilds the case
  (keeping the user's other choices from its provenance record) instead of
  ignoring the value and answering "updated".
- `run_campaign` never defaults its case type and never guesses the swept
  variable. An earlier "the flamelet has one sweepable parameter" rule turned
  "equivalence ratios 0.6, 0.7, 0.8" into a T_inlet sweep at 0.6–0.8 K.
- Temperatures outside 200–3000 K are refused as a wrong unit or a value meant
  for another parameter.

## 5. Module layout

**New:**

| File | Responsibility |
|------|----------------|
| `agent/config.py` | source of truth for models/keys/tunables (defaults + env + CLI) |
| `agent/runner.py` | runs the OpenFOAM solver pipeline for a built case; restart-from-latest-time planning + launch; physical-core detection + batch planning |
| `agent/run_batch.py` | detached orchestrator running several cases with bounded concurrency |
| `agent/.env.example` | documented template of every configuration variable |
| `agent/rag/corpus.py` | enumerate + read source files (globs + ignore rules) |
| `agent/rag/chunker.py` | split files into metadata-tagged chunks |
| `agent/rag/embeddings.py` | `mistral-embed` wrapper (batching, retry) |
| `agent/rag/index.py` | build / save / load index; hybrid (cosine + lexical) top-k search |
| `agent/index_build.py` | CLI to (re)build the index |
| `agent/tools.py` | tool schemas + implementations wrapping existing modules |
| `agent/chat.py` | conversational loop (RAG retrieval + tool calling) |

**Reused unchanged:** `validator.py`, `cantera_tool.py`, `of_writer.py` (wrapped
by `tools.py`).

**Changed:** `main.py` is now the conversational entrypoint (delegates to
`chat.py`); `extractor.py` has been removed — its extraction/follow-up role is
absorbed by the chat loop and the `build_case` tool. `requirements.txt` gains an
explicit `numpy` (already present transitively via Cantera).

**Dependencies:** `numpy` (vector index) and `rich` (terminal markdown
rendering); embeddings use the already-installed `mistralai` SDK. No vector DB,
no agent framework.

## 6. Models and configuration

This is intended for a **public GitHub release**, so every model and runtime
choice is user-configurable with friendly defaults — no code editing required.

### Models

- **Chat + tool orchestration:** default `mistral-large-latest`, which is markedly
  better at reasoning over the C++ (measured in the eval: README "Measured
  results"); `mistral-small` is the cheaper preset and `local-3b` the free,
  fully local one.
- **Embeddings:** `mistral-embed`.

### Configuration mechanism

A single **`agent/config.py`** is the source of truth. Precedence:
**CLI flag > environment variable > `agent/.env` > built-in default.** The local
provider adds the `openai` SDK (used only when `FGM_PROVIDER=openai`).

| Variable | Default | Purpose |
|----------|---------|---------|
| `MISTRAL_API_KEY` | *(required)* | API key |
| `FGM_CHAT_MODEL` | `mistral-small-latest` | answering / tool orchestration |
| `FGM_EMBED_MODEL` | `mistral-embed` | embeddings (recorded in index manifest) |
| `FGM_TOP_K` | `6` | retrieved chunks per query |
| `FGM_TEMPERATURE` | `0.2` | chat sampling temperature |
| `FGM_SHOW_TOOLS` | `true` | print a dim diagnostic line per tool call |
| `FGM_GUARDS` | `on` | `shadow` = guards detect and log but do not act (evals: `run_eval.py --guards shadow`) |
| `FGM_PHYSICAL_CORES` | `auto` | physical cores to budget for batch runs |

- CLI flags (`--chat-model`, `--embed-model`, …) override env vars for one-off
  runs.
- A committed **`.env.example`** and a README "Configuration" table document
  every variable.
- **Embedding-model safety:** the index manifest stores `FGM_EMBED_MODEL`; if it
  no longer matches, the tool warns and offers to rebuild (vectors are not
  comparable across embedding models).
- **Validation:** missing key or an unknown/blank model name fails at startup
  with a message naming the exact variable to set.

### Provider scope

All API access is funnelled through `config.py` + a thin client wrapper
(`agent/llm.py`), so the provider is a localized choice rather than a rewrite. Two
providers are supported (`FGM_PROVIDER`):

- **`mistral`** (default) — the hosted Mistral API (`mistralai` SDK).
- **`openai`** — any OpenAI-compatible endpoint at `FGM_BASE_URL` (`openai` SDK),
  i.e. a **local** server such as Ollama or vLLM, so neither code nor data leaves
  the machine. Both the chat model and the embeddings run there.

`llm.make_client(config)` returns the adapter for the configured provider; each
exposes the same `complete` / `embed` surface plus provider-specific
error-classification and the "force some tool" token (`any` for Mistral,
`required` for OpenAI), used by the forced-read degrade chain (named tool_choice →
force-some-tool → plain auto) so engines with limited `tool_choice` support still
work. Switching the embedding model (e.g. when changing provider) invalidates the
stored vectors, so `main._ensure_index` rebuilds the index automatically. The
deterministic case build/run pipeline is provider-independent — the model only
routes and explains, so a smaller local model is viable.

## 7. Risks / edge cases

- **Tool-call reliability:** always validate args in Python; never trust the
  model to enforce the required set.
- **Stale index:** manifest hash check warns when solver files changed since the
  last build.
- **Hallucination:** mitigated by grounding rules + citations + the larger chat
  model; the similarity threshold avoids anchoring on weak matches.
- **Network:** embeddings and chat both need the API (already true today); index
  build is a one-time call.
- **Context size:** corpus and top-k are small — no context-window pressure.

## 8. Verification plan

1. **Index build:** chunk count > 0, embedding matrix shape matches, reloads
   identically.
2. **Retrieval unit checks:** "turbulence model" returns the
   `turbulenceProperties` chunk; "progress variable" returns FGM source chunks.
3. **Grounded Q&A (live):** "can it run a turbulent simulation?" → a
   laminar-only answer citing `turbulenceProperties`.
4. **Tool call (live):** a full setup prompt triggers `build_case` and produces
   `testcases/laminarSlotBurner_agent/` with no leftover placeholders (reuses the
   existing checks).
5. **Regression:** the deterministic validator/Cantera/writer behaviour is
   unchanged.

## 9. Phasing

1. **RAG core** — *done.* `config.py`, `rag/{corpus,chunker,embeddings,index}.py`,
   `index_build.py`, `.env`/`.env.example`; index builds (68 files → 249 chunks)
   and retrieval verified; `agent/index/` and `.env` git-ignored.
2. **Conversation** — *done.* `tools.py` (`build_case`, `compute_thermochemistry`)
   + `chat.py` (RAG retrieval + grounded tool-calling loop); `main.py` rewritten
   as the conversational entrypoint using `config.py` and the index. `extractor.py`
   removed. Grounded Q&A, the thermochemistry tool, and case building verified live.
3. **Polish** — *done.* README rewritten for the conversational interface; live
   end-to-end tests pass.

## 10. Decisions

1. **Chat model — pending.** Default `mistral-large-latest` for answers
   (configurable via env var), with `mistral-small-latest` available to minimise
   cost / stay within free-tier rate limits. Large is paid (~$0.50/$1.50 per M
   tokens) but ~a fraction of a cent per turn here; a free experimentation tier
   also exists. Awaiting final choice of default.
2. **Confirmed — include the testcase directories.** Both `laminarSlotBurner`
   and `laminarFlamelet`, text files only (see §3 exclusions).
3. **Confirmed — add the read-only `compute_thermochemistry` tool.**
4. **Confirmed — git-ignore `agent/index/` and auto-build on first run.**
