# Style — documentation and comments

What this file governs: every prose surface of the repository — rustdoc, `//`
comments, test names, the `docs/` tree, README, the record, commit messages. What
it does not govern: code formatting, which [rustfmt.toml](../rustfmt.toml) pins
to the defaults, deliberately, so that hand-laid tables inside comments survive.

It is a rulebook, not a description of current practice. Where the two
disagree, the gap is listed in [What this changes](#what-this-changes) at the
bottom, with the reason.

## The readers

Most of this library was written by coding agents and is maintained the same
way. That fact shapes the audience more than anything else: the most frequent
reader of every prose surface is an agent that opens the repository with no
memory of any previous session. Documentation here is not a courtesy to
outsiders — it is the project's working memory, and the interface between one
session and the next.

1. **The researcher deciding whether to build on the library.** Knows more
   about the mathematics than the docs do; lives in Sage; is evaluating an
   instrument, and asks two questions in order. First, capability: *what does
   this let me compute that I currently cannot?* — a table past Sage's wall,
   a single coefficient whose full product no machine holds, an open
   conjecture checked one degree higher. That is why they came; if the docs
   do not surface it, they never look further. Second, trust: *which* `H̃` is
   this, and how do you know? — they have been burned by a normalization
   mismatch that produced plausible wrong answers, and results they publish
   will carry their name, not the crate's. Capability recruits this reader;
   trust keeps them. They read the README for the first question, docs.rs
   and the validation scripts for the second — and if they find a convention
   left unstated, they leave, correctly.
2. **The working agent.** Arrives with full command of the mathematics and of
   Rust, and none of this repository: every session starts from nothing but
   what the tree says. It reads by search — a grep hit and a window around it,
   not the module top to bottom — and it verifies rather than trusts: it will
   run the doctest, the pinning test, the harness. Two failure modes dominate,
   and much of this guide exists to block them. First, **confident
   pattern-matching to the wrong convention**: an agent has seen every
   `K_{λμ}` variant in training and will supply one fluently; only the doc
   stating *which* one, with a runnable pin, stops the wrong-by-a-twist
   answer. Second, **trusting stale text**: a human might smell that a
   paragraph is out of date; a fresh context reads it as current instructions.
   This reader is also the old "maintainer, five years out" — except the wait
   is not five years; the maintainer with zero context arrives at the next
   session.
3. **The human directing the agents.** Steers by reading the record and the
   commit messages, and reviews the diffs. The narrative genres —
   `docs/record/` and commit bodies — are primarily this conversation: the agent
   reports what happened and what it learned; the human decides what happens
   next.
4. **The auditor** — a Sage packager, a distro reviewer, someone reading the
   licensing story. Reads every claim looking for the check behind it.

The Rust engineer who wants symmetric functions as a component exists, but is
rare and is served for free: the precision reader 1 demands and the contracts,
pins, and examples reader 2 needs are everything that reader ever wanted.

## Three rules that govern everything

1. **Front-load for the partial reader.** Summary sentence, then contract,
   then mathematics, references last. Agents land on a search hit and read a
   window; humans skim listings. Both must meet the facts a caller cannot
   work without — what comes back, what is required, what silently goes
   wrong — before the treatise, and the treatise stays findable below it.
2. **A claim carries its evidence.** A formula names the equation it
   implements. A convention names the test that pins it. A number names the
   harness that produced it. An adjective — "fast", "small", "safe" — is
   deleted or replaced by its number.
3. **Genres do not mix.** *Reference* (rustdoc, README) is present tense and
   describes what is. *Record* (`docs/record/`) is past tense,
   anchored to commits, and describes what happened, including what failed.
   *Spec* is definitional — and has the shortest lifespan of any genre: a
   spec exists to be implemented against, and merges into the subsystem's
   record file when the implementation lands, leaving no separate document
   behind (the Specs section below gives the disposition; the clean-room spec
   is the exception, frozen outside the record as permanent evidence). A fourth
   genre is mortal by design: the **working paper** — the capability survey,
   the packaging audits, the release checklist — scoped to one question and
   expected to be absorbed or deleted when the question closes. The genres are roles, not
   paths: `docs/` may be reorganized or thinned at release, and a file may
   die — provided every role keeps a home and no durable surface is left
   citing a corpse. Before a working paper is deleted, whatever the reference
   or the record leans on moves out of it first, and the inbound links move
   with it (grep for them); a dangling evidence pointer silently converts a
   backed claim into an unbacked one. Time-indexed material placed in
   the reference rots there: [lib.rs](../src/lib.rs) opened the crate's front
   page with a "Roadmap (the marked seams)" listing `bignum`, `python`,
   plethysm, and an optimized LR backend as future work, all four of which had
   shipped long before anyone removed the list — and the last two remnants of
   it outlived even that, to be found by an audit two rewrites later. A human
   reader finds that embarrassing; a fresh-context agent finds it
   *directive* — stale reference does not just age, it redirects the next
   session toward plans already executed or retired. The README stayed correct
   over the same period because it is treated as reference and maintained; the
   lesson is the rule. That front page is reference throughout today, and
   the exhibit stays anyway: this rule is backed by a cost the project
   actually paid, not by a hypothetical.

## Rustdoc — the reference

### The first sentence

One complete sentence, ending in a period, that stands alone: rustdoc uses it
as the summary line in listings and search results, so it is read out of
context far more often than in context. It names the mathematical object and
what the item does with it, by the object's precise name — the precise name
is the one both the researcher and the agent search for.

Two in-tree models:

> The full Schur expansion of the skew Schur function `s_{outer/inner}`.
> — [skew_lr.rs](../src/skew_lr.rs), `expand_skew`

> A partition λ = (λ₁ ≥ λ₂ ≥ … ≥ λ_k > 0), stored as its nonzero parts.
> — [partition.rs](../src/partition.rs), `Partition`

Not a fragment, not a restatement of the name, not "Computes the thing this
function computes."

### The division of labor: module docs own the model, item docs own the contract

A module doc is the treatise: the mathematical objects, the conventions chosen
among the ones circulating, the engines and when each wins, the range and what
it opens, the references. [llt.rs](../src/llt.rs) is the model of the form.

An item doc is the contract, and stays lean because the module doc exists —
but it must stand on its own, because the item is where a search lands: the
contract cannot require having read the treatise. In order:

1. **Summary sentence** (above).
2. **Contract.** What comes back and in what state (sorted? zero-free?
   deduplicated?); what the arguments must satisfy; behavior at the degenerate
   inputs; `# Panics` / `# Errors`; cost *in shape terms*.
3. **Mathematics** — only what is specific to this item; link the module doc
   for the model.
4. **`# Examples`** — a doctest (next section).
5. **References** — `[KEY] Def 3.2` pointers into the module's reference list.

Degenerate inputs are convention choices, not corner cases, so they are
documented rather than left to be discovered: the empty partition, degree 0,
`k = 1`, equal shapes in a skew pair. `expand_skew` documenting "the empty
vector when `inner ⊄ outer`, and `[(∅, 1)]` when the shapes are equal" is the
standard; `z(∅) = 1` and `partitions_of(0)` likewise.

Cost is stated machine-independently: "cost is `#SYT(ν)` for all monomials at
once" ([llt.rs](../src/llt.rs), route R1) survives every hardware generation;
"takes 0.3s" is false somewhere already.

**A `foo` / `foo_in` pair states its range once**, on the one a caller lands
on, and the other points at it. The crate is full of these — `character` /
`character_in`, `character_table` / `character_table_in`,
`hall_littlewood` and its table — and a wall written out twice drifts twice:
`character_table_in` argued the memory wall in TB while `character_table`
argued it in GB, neither wrong and both maintained by hand. The generic form
says what its ring parameter is *for*, which is the question its own caller
arrived with.

### Sentence discipline in item docs

An item doc is read at a search hit, by someone who needs the contract and is
not reading an argument — the transactional register the [Voice](#voice)
section asks for, made checkable. Two rules, and they bind `///` only: a module
doc is the treatise, and a `//!` sentence that takes forty words to separate
two conventions is doing its job.

**One sentence carries one fact, and stops at 25 words.** The count is a proxy
for what actually costs the reader — how many facts must be held at once — so a
30-word list of degenerate inputs passes on sight and a 20-word chain of "and …
so … which" does not. `expand_skew`'s `# Panics` was a verbless fragment
followed by one 49-word sentence carrying three more facts: the retry in the
wider type, the refusal above it, and the shape that shows how far off the wall
is. It is now three sentences, and a reader who needs only the first stops
after it.

**Elide the subject, never the verb.** `Returns every ν with c^outer_{inner,ν}
≠ 0` and `Panics if a coefficient exceeds u128` are the forms the Rust
ecosystem writes and this guide keeps: the elided subject is the item the doc
sits on, and no reader resolves it wrongly. A missing verb is a different
thing. All 54 verbless `# Panics` sections in `src/` opened like `character`'s
— "If the value exceeds `i128`, around n ≈ 58." — which leaves the reader to
supply the relationship, and which `std` does not write either. The summary
sentence is the one place a verbless line is house style, because rustdoc
renders it as a title.

`scripts/check_doc_sentences.py` holds both, and `scripts/preflight.sh` runs
it. The verb rule is absolute. The 25-word cap is reported rather than gated —
685 of 3048 item-doc sentences are over it — and what fails the build is a
ratchet set just under the longest sentence that survived the last pass, so the
tail cannot grow back while it is being cut down.

**This is two rules of ASD-STE100 and deliberately not the rest.** The
standard also bans `-ing` forms and every word outside a 900-word dictionary;
those stay rejected, because the coined vocabulary below is already
one-word-one-meaning for a domain that dictionary does not cover.

⚠️ **The metaphor ban was rejected here and is now adopted** — see [No
aphorisms, no metaphors](#no-aphorisms-no-metaphors), which governs. The
rejection argued that "the convention minefield" is the heading that makes an
agent stop and read, and that the coined vocabulary covers the rest. Both
halves failed on contact. The vocabulary argument was backwards: `seam`,
`layer` and `peel` are not evidence that metaphor is safe, they are the eight
cases where a metaphor was *converted into a term* by being pinned to one
meaning — the undefined ones are exactly what the ban is for. And the
stop-and-read argument valued a heading's effect on attention over its
content, which is the trade this file refuses everywhere else. What was taken
from the standard is what serves a reader arriving mid-file with no context:
short sentences, one fact each, no grammatical gap to fill, and no figure to
decode.

### Examples are convention pins

Every public entry-point family carries at least one doctest whose value is

- **small enough to verify by hand**, and
- **chosen to reveal the convention** — the value that *distinguishes* this
  normalization from its rivals, not the one every convention agrees on.

`s_2 · s_1 = s_3 + s_{21}` pins a product. A `K̃_{λμ}(q,t)` value where the
(q,t) order matters pins which modified Macdonald convention shipped; a charge
value that differs from cocharge pins that choice. The researcher pastes the
example into Sage and compares — that one paste answers their trust question —
and the agent executes it, which turns the stated convention from prose it
must trust into a fact it has checked. Doctests are also the only
documentation CI executes, so they are the only documentation that cannot rot.

### Conventions get their own section

Any family where the literature circulates more than one normalization gets a
module-doc section stating each one, which public function returns which, the
dictionary between them, and the trap. "The conventions in circulation" in
[llt.rs](../src/llt.rs) is the model, and its two properties are the
requirement:

- each entry says **what silently goes wrong** — "returns a wrong-by-a-twist
  answer rather than an error" — because that is the failure the reader must
  recognize, and because an agent fluent in every circulating convention will
  confidently apply the wrong one unless the difference is stated where it
  bites, and
- each entry names the **test that pins it**, because a stated convention with
  no pin is a promise; a pinned one is a fact any session can re-verify by
  running it.

State the Sage equivalent by name — `kfpoly`, the `Ht` basis — or state
explicitly that none exists. That single line is what reader 1 greps for, and
its absence is the most expensive omission a module can have.

### Range and performance in rustdoc

Properties **of the code** belong: asymptotics, allocation counts, byte sizes,
scaling shapes. "One allocation per term", "cost independent of `p(n)`",
"14.1 MB per call on `[8,7,6,5,4,3]²`" are deterministic facts a reader can
reproduce exactly.

**Range belongs too, and is not optional.** The researcher planning a
computation campaign needs the practical limits before anything else: which
degrees are routine, where the fixed-width wall sits and what crosses it, what
has actually been computed. State range in reproducible terms — "the
fixed-width wall is at total degree 24 and is `z_γ` rather than the answers,
which are 16 bits; `bignum` carries it to 32" is a fact about `i128`, not
about a laptop — and point at the record's tables for the largest runs.

Properties **of one machine and one rival's version** do not belong: seconds,
RSS, ×Sage and ×Symmetrica ratios, and **profile shares** — "83% of the
profile", "53.8% allocator" — which are a sampling run on one binary and read
as durable facts. Those live in the record with their harness and context, and
the rustdoc names the record file that owns them — as a backticked path, the
docs.rs-safe form.

Filing the number off does not convert a profile share into a durable fact.
"The single largest thing in the profile", "spent essentially all of its
runtime here", "the broadest win" are the same sampling run with the one
checkable part removed — still one binary, now unfalsifiable. What survives a
recompile is the **mechanism** that produced the share: "a doomed division
still runs the whole elimination", "the same list would otherwise be rebuilt
thousands of times". Say that, and let the record hold the percentage. A
superlative is worth writing only where it is structural — "the cheapest
product in the crate" is true of multiset union by construction.

⚠️ A percentage that describes the **mathematics** is not a profile share and
stays: "72% of its trial divisions fail" is a property of the atom family,
true on any machine, and it is the reason the necessary-condition pre-pass
exists. The test is whether recompiling could change it. The distinction matters because
rustdoc ships inside the crate forever, while the record is dated and expected
to age.

A ratio against another project is the worst case, because **it is a claim
about two codebases and the reference can only version one of them**. When
Symmetrica or lrcalc improves, every `2.5x` in this crate silently becomes a
lie that no test fails on, no reader can date, and nobody thinks to re-run —
and it is exactly the number a researcher deciding between the two would act
on. What survives the rival's next release is the *shape* of the difference:
"Sage prices this like a full expansion", "the solve needed `O(p(n)²)` Kostka
numbers to read `p(n)` of them", "out of Sage's range at `J[3,2,1]²`" — the
algorithmic fact that produced the ratio, which stays true when the constant
factor moves. State that, cite the record for the measurement, and let the
dated document carry the number.

### Say what it opens

Where the library answers something no other package can — a
single-coefficient Schubert query on a product no machine can materialize, an
`st` product past the wall Sage dies on — the module doc says so, plainly, at
the top. That sentence is the capability the researcher came looking for, and
burying it under implementation notes is the most expensive modesty a module
can have. Lead with the newly askable question; the speedup is the means, not
the point — "aimed at a gap rather than at parity" is the README's framing,
and each gap-aimed module carries it down.

Honesty runs both directions, in the house manner: where the capability exists
elsewhere, say so and state the actual differentiator — scale, exactness, or
query shape. Behind every such claim stands a measured survey of the
incumbents — dated, versioned against them, methodology caveats attached,
with its own "what *does* exist, to be fair" section. Today that survey is
[research-gaps.md](research-gaps.md), a working paper from the library's
creation; if it does not survive to release, its measured walls move into the
record — they are dated measurements, which is exactly what the record holds
— and the claims re-point. A capability claim points at the survey, wherever
it lives, rather than restating it.

Internal docs (`pub(crate)` and below) may keep measured numbers where they
justify a design — the `Key` packing comment in
[skew_lr.rs](../src/skew_lr.rs) citing "a measured 38% of wall time" is the
model — provided the harness that produced them is named.

### What rustdoc must not contain

- **Roadmaps and future work.** The record owns the future. (The lib.rs
  exhibit, above.)
- **History.** "This used to be…" and "changed in…" belong to git and the
  record; the reference describes the present. The tempting case is the
  cautionary tale: [character.rs](../src/character.rs) closed its range
  statement with "this module previously returned `i64` and wrapped
  silently — χ^λ(1³⁶) came back *negative*", three lines that answer no
  question a caller has. The incident is in the commit that fixed it
  (`c8cf34b`), where a reader looking for it will be looking.
- **Reviewer-talk.** "Now correct", "simplified", "cleaned up" — statements
  addressed to a diff reviewer are noise the moment the commit merges.
- **A number that measures the design which lost.** `character_table`'s only
  figure was "recursing per entry instead was worth 0.66–0.77x against
  Symmetrica's `chartafel` — behind": a ×-ratio for the *rejected* approach,
  which a skimming reader attaches to the function it sits on. It was also
  the pre-fix number the record carries as "was 0.6x" — the shipped sweep is
  2.9× at degree 12 — so the reference advertised a deficit the change had
  already closed. The alternative that lost belongs to the record with the
  one that won, and a comparison the reference keeps must be about the code
  a caller is about to run.
- **Explanations resting on items the reader cannot open.** The same doc
  credited `p_expand`, which is private — and is not the function this path
  uses (`p_expand_shared` is). A private helper can be named in a `//`
  comment beside the code; leaning on one in rustdoc gives the reader a
  dead end and the writer a place to be wrong unnoticed.

## Mathematical notation

**Bracketed math goes inside backticks; the rest needs nothing.** The hazard
is exactly one character pair: rustdoc parses a bare `[q,t]` as an intra-doc
link, so `ℚ[q,t]`, `ℤ[t]`, `ℕ[α]` and `f[g]` must be `` `ℚ[q,t]` `` and so
on. Bare `Σ_ν`, `⟨H̃_μ, h_ν⟩`, `μ ⊢ n`, `⊗`, `λ'` are safe and are left
alone — backticking them is noise that buys nothing. The distinction is not
stylistic: every one of the 54 bracket sites fixed in the sweep to zero
`cargo doc` warnings was a bracket, and no bare glyph ever warned.

**In `//` comments nothing needs backticks at all**, because rustdoc never
parses them. Reserve the markup for `///` and `//!`, where it does work.

Display formulas go in ` ```text ` blocks, as in
[skew_lr.rs](../src/skew_lr.rs) and [llt.rs](../src/llt.rs).

**No LaTeX rendering** (KaTeX header injection or similar). The crate's docs
build with zero dependencies like the crate itself; Unicode covers the notation
this subject needs; and source, terminal, and docs.rs all render the same
glyphs. This is a decision, not an accident.

**Prototyped and rejected, 2026-08-01**, on `deltaop.rs` with both KaTeX
0.16 and MathJax 3 injected through `--html-in-header`. Recorded so it is not
re-explored at full price; the display math genuinely does render better, so
the reason for rejecting it is not that it fails to work.

- **Both libraries behave identically**, because both post-process the
  rendered DOM and the damage is upstream, in rustdoc's markdown pass.
  Choosing between them is a question of size and accessibility, not of
  whether this works.
- **Markdown eats the LaTeX first.** `_` is subscript in TeX and *emphasis*
  in markdown, so `\tilde{H}_\mu … \rangle_*` becomes `<em>` tags and the
  `$$…$$` delimiters end up split across elements, at which point nothing
  renders. Display math survives only inside a raw `<div>`, whose contents
  CommonMark passes through untouched; inline math survives only if every
  underscore is escaped `\_`. An inline `<span>` does **not** work — only
  block-level HTML skips markdown.
- **Smart quotes break primes.** `a'` becomes `a’`, and this crate's notation
  is full of `a'`, `l'`, `λ'`, `H'_λ`, `Δ'`. Each one needs `\prime`.
- **The failure is silent, which is what decides it.** A formula with an
  unescaped `_` or `'` emits *no* `cargo doc` warning, passes
  `scripts/preflight.sh`, and renders as raw LaTeX in the published page.
  Every other documentation hazard in this tree is caught by a check; this
  one would be caught by a human noticing.
- Converting only display math leaves two notations on one screen — the
  prototype had `⟨F, H̃_μ⟩_*` typeset and `⟨F,s_κ⟩_*` monospace three lines
  apart — and converting the inline math means touching every formula in
  every prose line, each carrying the silent-failure risk above.

Revisit if rustdoc gains native math support, so the source stops being
markdown-mangled, or if a check can validate the escaping — the objection is
the silence, not the syntax.

**ASCII in identifiers, Unicode in prose.** `lambda`, never `λ`, as a variable
name. Use the literature's letter when the object has no better role name;
use the role name when it has one; never both names for one object in one
module.

| notation | meaning | in identifiers |
|---|---|---|
| `λ, μ, ν` | partitions | `lambda`, `mu`, `nu` — or role names: `outer`/`inner` for a skew pair, `shape`, `content` |
| `λ'` | conjugate partition | `conjugate` |
| `ℓ(λ)`, `\|λ\|` | length, size | `len`, `size` |
| `z_λ` | centralizer order | `z` |
| `s, h, e, m, p, f` | the six classical bases | `Schur`, `Homogeneous`, `Elementary`, `Monomial`, `PowerSum`, `Forgotten` |
| `q, t, α` | parameters | `q`, `t`, `alpha` |
| `ℚ[q,t]` vs `ℚ(q,t)` | polynomial ring vs fraction field | `QtPoly` vs `Frac` — never blur the two in prose |
| `H̃`, `Q'` | modified Macdonald, Hall–Littlewood Q′ | `Ht`/`htilde`, `hall_littlewood` |
| `c^λ_{μν}`, `K_{λμ}(…)` | LR coefficient, the Kostka family | spell out which Kostka: `K`, `K(t)`, `K(q,t)`, `K̃(q,t)` are four different objects |

### One word, one meaning across the tree

A word that carries two senses is a normalization trap in prose: the reader
resolves it wrongly and never learns they did. The senses do not have to
collide inside one file to cost something — a term whose meaning depends on
which directory you are in is a term nobody can grep.

Three words are spoken for, and a fourth is retired. The coined vocabulary
that follows them is fixed to one sense each; a new coinage earns its place
the same way — defined once at the thing it names, then used, never
re-explained.

Spoken for:

- **`layer`** is the live state set of a step-indexed traversal — the data
  structure, nothing else: the merged partial fillings in
  [skew_lr.rs](../src/skew_lr.rs), the β-mask maps in
  [convert.rs](../src/convert.rs), the chain states in
  [kostka.rs](../src/kostka.rs). Qualify it (`row layer`, `skew-LR layer`,
  `DP layers`) wherever the Python boundary's *three layers*
  ([python.md](policies/python.md)) is close enough to be misread.
- **The algorithm around a layer is named, not called "the layer".** Write
  `SkewLr`, `StripLr`, `AutoLr`, or "the layer DP" — "counting beats
  `SkewLr`", never "counting beats the layer". A comparison names the two
  things being compared.
- **`level`** is a mathematical parameter and is not available: the ribbon
  level `k` in [llt.rs](../src/llt.rs), the `(perm, level, stufe)` memo key in
  [schubert.rs](../src/schubert.rs), and `level_arg` in
  [failure.md](policies/failure.md).
- **`frontier` is retired.** It was doing all three jobs at once — data
  structure, engine, and the limit of what is computable — and the third sense
  has its own words already: a hard limit is a **wall** (the fixed-width wall,
  Sage's plethysm wall, "the computational wall"), and what a change to the
  checks does is **widen coverage**.

The coined terms, each fixed to one sense:

- **`seam`** — a `Ring` method where swapping the coefficient type changes
  behavior with no call site edited, in the sense of Feathers, *Working
  Effectively with Legacy Code*, ch. 4; what is swapped here is exactness, not
  testability. Defined at the `Ring` trait in [coeff.rs](../src/coeff.rs), and
  that is the only place it is explained.
- **`atom`** — a factor `q^a − t^b` in the `(q,t)` fraction fields
  ([bh.rs](../src/bh.rs), [frac.rs](../src/frac.rs),
  [deltaop.rs](../src/deltaop.rs)). ⚠️ The literature's *Demazure atom* and
  *Lascoux–Schützenberger atom* are unrelated objects; nothing in the tree
  computes them today, so the word is free, but a change that introduces them
  must rename one of the two rather than let both stand.
- **`route`** — a dispatch path through the engines to an answer (the product
  route, the LLT route, the two slow routes). It describes *how the answer is
  reached*, never the answer and never the engine: name the engine when
  comparing engines, exactly as with `layer`.
- **`ladder`** — a benchmark sweep over an increasing parameter, usually
  degree, whose rows are one workload at successive sizes (the degree ladder,
  the mains ladder, the measured ladder). It is a *measurement* word. The
  straightening and offset ladders in [llt.rs](../src/llt.rs) are algorithmic
  structures and predate the measurement sense; qualify those at every use so
  a reader never has to guess which is meant.
- **`range`** — how far a family is exact and usable: the fixed-width wall,
  the degrees that finish, what a `bignum` build changes. It is a property of
  the **entry point**, not of the family
  ([failure-and-overflow.md](record/failure-and-overflow.md)) — a whole-degree
  table and a single shape have different ranges. The `## Range` section of a
  module doc is where it lives, and R9 in
  [failure.md](policies/failure.md) requires one per fixed-width public
  family. ⚠️ It was called **reach** until this entry; `reach` is now only the
  ordinary verb ("values a caller can reach", "the reachable state space"),
  and the identifier `abacus_reach` is the one survivor, deliberately.
- **`ceiling`** — the largest value something holds: `3.40e38` for `u128`, 32
  for the `u32` orientation mask, ~1.15× for a speedup that cannot be beaten.
  A **wall** is the other half of the same limit — the point in the *input*
  space where a family stops, stated in n or in degree — and `range` is the
  description of what lies below it. So `i128` is a ceiling and n ≈ 58 is the
  wall it produces, and a sentence names whichever one it is actually stating.
  ⚠️ The two were interchangeable until this entry, and
  [hl.rs](../src/hl.rs) is the exhibit: "the character ceiling near n ≈ 58"
  and "no wall at all" sat two sentences apart, saying the same kind of thing
  in two words, in the paragraph whose whole point is that the limit is not
  a function of n. Converged on touch, as with the spelling rule, not swept.
- **`lex-monic`** — leading coefficient ±1 under the lexicographic order on
  exponent pairs that [qt.rs](../src/qt.rs) already sorts terms by; for
  `q^a − t^b` the leading term is `q^a` when `a > 0` and `−t^b` when `a = 0`.
  It is what keeps exact division inside `ℤ[q,t]` — elimination divides only
  by the leading coefficient, and a unit needs no field — so "up to sign" is
  the whole content of the claim, not a hedge. Defined in the arithmetic
  section of [bh.rs](../src/bh.rs).
- **`peel`** — the part-removal recursion (Morris, and the merged peel DAG in
  [gjmod.rs](../src/gjmod.rs)), not a general "strip one off" verb.
- **`wall clock` is always two words, and never bare `wall`.** A bare `wall`
  is the capacity boundary above. `min-of-3` and `interleaved` are the words
  for how a timing was taken; "wall" alone never means elapsed time.
- **`battery`** is a suite of checks a family must carry
  ([validation.md](policies/validation.md)). The power source is written *on
  battery* and only ever appears inside a ⚠️ caveat on a measurement
  ([record/README.md](record/README.md), "Power state"), so the two never
  collide in the same sentence.

## Citations

Bracketed keys — `[LLT]`, `[HHL]`, `[KMS]` — defined once per module in a
`## References` block. Each entry carries authors, title, an arXiv link, and
**which results are used, by equation number**. "See [HHL]" sends the reader
on an expedition; "[HHL] Def 3.2" is a grid reference.

**The block goes at the end of the module doc**, with the link definitions
after it. Rustdoc renders the module doc verbatim at the top of the module's
page, so a bibliography placed early is the first screen a reader gets: the
`## References` blocks in this tree all sat 2–5% into their module docs until
2026-08-01, which put 25 lines of arXiv links between `llt.rs`'s summary
sentence and its first word about what an LLT polynomial is. The model comes
first; the bibliography is what you consult after it.

The block in [llt.rs](../src/llt.rs) is the model, down to its most valuable
line: recording that [KMS] is the *normative* source for the straightening
rules because [LLT] §7's printing of the same rules carries two misprints. A
sentence like that is a day of someone's life, saved. When sources disagree,
say which one this crate follows and why.

Make the keys resolve. A bare `[KMS]` in prose is an unresolved link warning.
Reference definitions at the bottom of the module doc fix it:

```text
//! [KMS]: https://arxiv.org/abs/q-alg/9508006
//! [HHL]: https://arxiv.org/abs/math/0409538
```

**They fix the module doc only.** Rustdoc scopes link definitions to the doc
comment they appear in, so a definition in `//!` does nothing for a `[KMS]`
inside a `///` on an item — in [llt.rs](../src/llt.rs) the definitions
resolved 21 module-doc mentions and left all 29 item-doc ones warning. Item
docs therefore **escape** instead: `\[KMS\]`, which renders as `[KMS]` and
sends the reader to the module's `## References` block, where the link lives.
Repeating the definitions per item would resolve them too, at the cost of a
block of URLs on every documented function.

A key with no URL — a journal-only reference like `[GJ]` or `[BH]`, or an
unpublished preprint like `[GH]` — is escaped everywhere, module doc
included. There is nothing to link to, and the `## References` entry is what
carries the citation.

## Comments in the code

A `//` comment earns its line by stating something the code cannot. Four
things qualify:

1. **The invariant relied on right here.** "Strict edges run u < v, so the
   constraint on the vertex being assigned is always against an
   already-coloured endpoint" ([hopf.rs](../src/hopf.rs)).
2. **The mathematical fact that makes the step valid.** "Column counts are
   automatically weakly decreasing and positive."
   ([partition.rs](../src/partition.rs)) — one sentence, at the line that
   would otherwise look unjustified.
3. **The measured reason for a shape that looks wrong.** The `num-bigint`
   block in [Cargo.toml](../Cargo.toml) is the model: the decision, the
   measurement behind it (`examples/coeff_sizes.rs`), the licensing constraint
   — everything a future "why not GMP?" needs. A comment citing a measurement
   names the harness that produced it.
4. **The trap.** "Pinned to 0.4 … a 0.5 here would compile and then produce a
   *second, incompatible* BigInt that pyo3 cannot convert."

Nothing else does. No narrating mechanics, no restating the line below, no
talking to the reviewer.

**When the decision has a story, the comment points at it.** The comment
itself carries the conclusion and the premise — enough that the code can be
maintained safely if nothing else were readable — and then names where the
narrative lives: the harness, the record file, or both. The two jobs split by
reader: the conclusion-in-place serves whoever must *not break* the code; the
pointer serves whoever wants to *challenge or extend* the decision, and pays
for the full story only when they ask for it. A pointer never substitutes for
the conclusion — "see the record" with no stated reason is the outsourcing
this section exists to forbid. The trigger for a pointer is a story worth the
trip: a rejected alternative, a ranking that inverted, a measurement that
surprised; a pointer on a plain invariant is noise. Write pointers as
backticked repo paths to **files** — `docs/record/llt.md`, the form
[llt.rs](../src/llt.rs) uses: greppable, so a reorganization's link sweep
finds them; inert on docs.rs, where a relative markdown link would 404; and
verifiable, since a file either exists or does not. Never a coordinate inside
the file — see [Specs, and how they end](#specs-and-how-they-end) for what
section numbers cost when they drift.

- **Assert messages are documentation printed at the worst moment.** Write
  them as sentences stating the violated requirement in the problem's terms:
  `"an abacus needs at least ℓ(λ) beads to hold λ"`
  ([partition.rs](../src/partition.rs)), not `"bad rows"`.
- **No TODO / FIXME / commented-out code.** The tree has zero today; keep it
  at zero. Future work goes in the record, where it has context and gets
  triaged. Rejected code goes nowhere — git keeps the bytes, and the record's
  "measured and rejected" notes keep the reason.
- Section banners (`// ===== the tuple model`) are welcome in long modules.

## Tests

- **Test names are propositions** — declarative sentences of the law being
  checked: `conjugate_is_involutive_and_correct`,
  `k_quotient_is_independent_of_the_padding`,
  `to_schur_is_a_ring_homomorphism`. The test list should read as a table of
  contents of what is known to hold.
- When the name cannot hold the whole law, a doc comment states it: "The
  Littlewood decomposition: `|λ| = |k-core| + k · Σ |quotient|`."
- **Every assertion in a sweep names its counterexample**:
  `assert_eq!(…, "s→h→s at {lam}")` — a red test that does not say *which
  input* failed wastes the sweep that found it.
- **Fixtures state their provenance** in the file, and the generating script
  is committed (`scripts/gen_sage_oracle.sage`), so an auditor can regenerate
  rather than trust.
- Prefer law-shaped tests — two computations sharing no code, agreeing — over
  value pins, where both exist. When a value pin exists to fix a *convention*,
  its comment says so; changing one of those changes what the library means,
  which an ordinary regression value does not.

## Research drivers — examples/

The examples directory is where the library is an instrument rather than a
dependency: conjecture checks, coefficient hunts, calibration probes. The
genre has three obligations, and
[delta_conjecture.rs](../examples/delta_conjecture.rs) models the first two:

- **Open with the question**, the run lines, and what the output means.
- **State the interpretation contract** — which disagreement is a bug and
  which is a discovery: "The **rise** version is a theorem, so a mismatch
  there is a bug in this crate. The **valley** version is open — a mismatch
  there … would be a counterexample and is printed as one rather than
  swallowed." A driver that checks an open conjecture is instrumented to
  *notice*: the one output it must never mangle is the interesting one.
- **Record propose and dispose together** where the driver embodies a search
  strategy ([find_nonzero.rs](../examples/find_nonzero.rs) is the model): why
  the approach should find something, and the known reason its raw output
  cannot be trusted without the exact re-check it feeds.

Bench and profile drivers are instruments of the record instead: their numbers
land in `docs/record/` with the harness named, per the record's rules.

## The record — docs/record/

The record is where narrative belongs, and its discipline is the two-clause
form the commit titles use: what was done, and what was learned.

One file per subsystem, plus `README.md` as the index. The five written
specifications merged into those files rather than surviving beside them. The
directory was `docs/roadmap/` until 2026-07-31 — a roadmap when nothing was
written, which became the record as the plans were executed. That is the
natural life of a plan under execution and the right outcome, but the label
had come to point the wrong way: these files are the memory, not the plan,
and a name is a signal a fresh-context reader takes at face value. An agent
sent to "check the roadmap" landed in the past; an agent hunting for what had
already been tried did not think to look under "roadmap" at all. What remains
genuinely forward-looking is thin — [release-readiness.md](release-readiness.md)
and each file's open tail (below) — and no longer shares a name with the
history.

- **Negative results are first-class.** "Measured and rejected"
  ([memory.md](record/memory.md)), "three designs measured, three dead",
  "it is 0.5x, and why". A dead end recorded with its measurement stays dead;
  one recorded nowhere gets re-explored at full price, and not in some distant
  year — the next session that wanders near it has no memory of the last one.
  The record is the working agent's only long-term memory.
- **Numbers carry their context**: the harness, the input, the build flags,
  and — until CI exists — the standing caveat that every number is from one
  machine ([release-readiness](release-readiness.md) says it plainly; keep
  saying it until it stops being true).
- **One owner per number.** The subsystem file under `docs/record/` owns its
  benchmarks. README may quote headline numbers as the shop window, but every
  quoted number points at the record entry that owns it, so an update has one
  place to land and staleness is detectable.
- **Every entry stands alone.** The reader of the record — human or agent —
  does not have the conversation that produced it. The entry itself names the
  workload, the question it was answering, and the commit; the tree is the
  only context that survives the session that wrote it.
- **Lessons get stated as findings**, one bold sentence naming the case, so
  they can be found again and checked: "the node-count cost model ranked E3
  above E2 because it omitted the size of the running element". ⚠️ This rule
  asked for the maxim — "a cost model with a factor missing will rank engines
  confidently and wrongly" — until the ban in [No aphorisms, no
  metaphors](#no-aphorisms-no-metaphors). The maxim is not more portable, only
  less falsifiable: the finding carries the engines and the omitted factor, so
  a later session can tell whether its own case is the same one.

### How a learning ages

The reference answers *what is true of the code today* and stays bounded — it
grows with the code, not with the calendar. The record answers *how we came to
know it* and grows monotonically. The two stay compatible because a learning
moves through fixed stages:

1. **Capture** — the commit body, written at the moment, states what was
   built, what was measured, what was learned. It is anchored for free —
   date, tree state, diff — and immutable.
2. **Consolidate** — the subsystem file absorbs the durable version: the
   table, the premise, the portable lesson. Append-mostly, not append-only:
   when a later result supersedes an entry, annotate the entry rather than
   rewriting it — the ⚠️ on the record's Phase 5, marking the GMP item as
   describing "what was built, not what ships", is the model. The record's
   whole value is that it can be trusted backwards.
3. **Promote** — when a learning hardens into a fact about the code, its
   *conclusion* moves into the reference at the point of use — the dependency
   comment in Cargo.toml, a conventions entry, an item's contract — compressed
   to the decision plus a pointer back to the derivation. The reference
   absorbs conclusions; it never absorbs journeys.
4. **Demote** — when a promoted conclusion is reversed, it leaves the
   reference entirely and becomes one more chapter of the record. Nothing is
   deleted; it is re-filed as history.

Three rules keep an ever-growing record usable rather than merely large:

- **A rejection records its premise.** "Measured and rejected" binds only
  while its premise holds — the workload shape, the allocator's behavior, the
  rival's version, the representation of the day. State the condition, not
  just the verdict — "churn matters when the sizes are diverse or the
  buffers are retained" ([memory.md](record/memory.md)) — so a later
  session can tell a dead end from a door someone has since unlocked. A
  verdict with no premise is a permanent wall no agent will ever re-test; a
  verdict with its premise is re-opened exactly when the premise falls.
- **State at the head, history in the body, open questions at the tail.**
  Each subsystem file opens with a short present-tense digest — where the
  subsystem stands, with pointers — followed by the chapters in order, and
  closes with what is genuinely open: routes not yet tried, conjectures not
  yet checked, walls not yet pushed. The digest is held to reference
  discipline: it is the one part of the record that must not go stale. The
  body is grepped, not loaded. The tail is the one home future work has —
  rustdoc is forbidden from carrying it, so an open question not written
  here evaporates with the session that noticed it — and it sits next to the
  history that makes it intelligible: "X and Y were measured and rejected, so
  Z is the open route" is a tail entry whose premises are the chapters above
  it. When a tail item is executed, it moves up into a chapter; when it is
  abandoned, it moves up with its disposal reason. The tail shrinks by
  promotion, never by silent deletion.
- **A measured lesson gets one home; everywhere else cites it.** A finding
  that governs how numbers across the record are read belongs in exactly one
  block, with its evidence; every file that depends on it carries the marker
  and a pointer, not a retelling. The ~1.8× battery drift was independently
  re-derived in four files before it was consolidated into "Power state"
  ([record/README.md](record/README.md)) — four narrations of one fact, each
  citing its own evidence, none of them wrong and none of them the place to
  update. The failure mode is not length: it is that a reader cannot tell
  which copy is authoritative, and a corrected number leaves the other three
  standing. Restating a lesson is how the record acquires contradictions;
  citing it is how the record stays trustworthy backwards.

## Specs, and how they end

A spec is a working document with a lifespan: it exists so an implementation
can be built against something fixed, and its planning function expires the
day the implementation lands. While it is alive, it is definitional and
self-contained: notation defined before use, definitions cited to textbooks
(Macdonald; Fulton) rather than to this repository, readable by an
implementer holding nothing else. Where it records formula-by-formula
verification, that ledger is part of the spec: a formula listed without its
verification status is a claim, not a spec line.

When the implementation lands, the spec stops being a plan. What remains is
evidence of how the implementation was derived and checked — record material —
so it **merges into the subsystem's record file** and stops existing as a
separate document. Its parts have three fates:

- **Measurements, dead ends, methodology, verification ledgers** move into the
  record file. This is the bulk of what a landed spec is actually made of, and
  it is exactly what the record holds.
- **Conventions, references and traps** are copied into the module doc, where
  a reader lands. Copied, not moved: the module doc may delegate *depth* to
  the record, never the convention itself.
- **Plan scaffolding dies** — proposed API signatures the code has superseded,
  crate-fit tables of `file.rs:NNN` pointers that drifted the day they were
  written, correctness checklists now realized as named tests, performance
  targets already scored against measurements.

**This reverses a rule that stood in this file for one day, and the reversal
is the more instructive half.** The 2026-07-31 reorganization retired the five
specs *intact* as `docs/record/<subsystem>-spec.md`, on the argument that they
were cited by section number from source — `§1.3(a)`, `§3.1`, `§7.1` — so that
dissolving them would break every citation at once. A spec, the reasoning went,
is a stable numbered document, and that is what makes it citable.

Both halves were checked the next day and both were false. The citations are
**provenance, not retrieval**: at nearly every site the citing comment already
states the fact it needs in full, and the total unique content behind all
37 citations was on the order of tens of lines, most of it already written in
`src/`. And the numbering was not stable — it had already drifted. `llt-spec`'s
open-questions list ran 1,2,3,4,5,**7,6**; `schubert-spec`'s ran
1,**2,2**,4,5,6,7,**9,10**, and its four "§7 Q6" references all meant item
**7**. A citation format that is already wrong in four places is not a stable
target; it is an unverified claim wearing the costume of a precise one.

The general lesson, which is why this is written down rather than quietly
fixed: **a cross-reference is only as good as the thing that checks it.** A
doctest is checked by CI; a test name is checked by the compiler; a section
number is checked by nobody, so it decays silently and takes the reader's
trust with it. Prefer a pointer to something executable — a test name, a
function — over a pointer into prose, and prefer a fact stated in place over
either. Where the record genuinely needs to be named, name the *file*, whose
existence a link check can verify, and not a coordinate inside it.

Clean-room specs are the exception, and
[cleanroom-spec-skew-lr.md](cleanroom-spec-skew-lr.md) is the model twice
over. While alive it carries two extra obligations: state at the top *why the
document exists* and what it deliberately excludes; and contain only
mathematics, functional requirements, performance requirements, and interface
— **no implementation technique**, because the document's legal function is
to prove the implementer needed none. And it does not retire into the record
at all: its function is legal, not historical, so it stays at
`docs/cleanroom-spec-skew-lr.md` where NOTICE.md and the README cite it, and
survives as long as the licensing story does. Filing it under "record" would
subordinate evidence to narrative; it is the one document whose location is
part of its argument.

## README

The shop window, read by the researcher and the auditor, and the standard the
other reference surfaces should meet — it is the one document that has stayed
current.

- **Every claim ships with its check, in the same breath**: "4678 computations
  driven by Sage itself (`scripts/check_backend.py`)".
- **Caveats travel with the claim, not behind a footnote**: "1.84x
  like-for-like, or 4.37x with a cache Symmetrica's wrapper does not have and
  could equally adopt"; "including the ones that went the wrong way". This
  habit is what makes the rest of the numbers believable. Guard it — the day a
  number appears without its caveat is the day the rest stop being believed.
- **A test count is not a validation claim, and does not belong on a durable
  surface.** "444 unit and integration tests" says nothing about what is
  covered, and it is wrong again on the next commit that adds one: the README
  said 202 and the record's digest said 407 while the true figure was 444, and
  all three had been written by someone who checked at the time. Name the
  suites and what each one holds — the law suites, the committed fixtures, the
  non-field coefficient ring — which is the claim a reader wanted and which
  stays true. Counts of *checked values* are different and stay, because they
  are the size of the evidence: "8647 computations driven by Sage itself".

## Commit messages

- **The title states what is now true that was not**, with the number when
  there is one — and a failed experiment gets the same prominence as a win:
  "Measure the LLM route: 29x slower, and the solve is 99.99% of it".
- **The body is the record entry in miniature**: what was built, what was
  measured, what inverted, what was learned — the Schubert commit (`3bd3c00`)
  is the model.
- **Corrections to earlier claims get their own paragraph** ("Also corrected:
  …"), never a silent fix. The record's value is exactly that it can be
  trusted backwards.
- **A commit body is read cold** — in `git log`, by someone who has not opened
  this file — which costs three things. A rule restated from here takes the
  instruction form rather than the heading: "do not call a technique a trick;
  name the mechanism", not "nothing here is a trick". The headings below are
  written to be memorable to a reader who already holds the rule, and the
  reader of a message is not that reader; `237d6ac` landed four prose rules
  and restated all four as headings. Coined vocabulary is spelled out or
  avoided — `layer`, `range`, `profile share`, `the durable half` all resolve
  here and nowhere the message's reader is standing. And metaphor meets the
  test in [Voice](#voice).
- **Round a ratio toward the claim you can defend.** `b3f328d`'s title read
  "2.0-4.0x cumulative" against a measured 2.03-3.95x: both ends rounded
  outward, so the title claims a wider win than the harness produced. Round
  inward or quote the measurement. This is the README's caveats-travel-with-
  the-claim habit applied to arithmetic, and rounding is where it slips
  unnoticed, because no one reads a rounded number as a claim.
- **A message is immutable, so it is written in the vocabulary that will
  outlive it.** `b3f328d` titled a data structure "the Pieri frontier"
  18 minutes before `5c09a33` retired the word; [convert.rs](../src/convert.rs)
  was swept and reads 25 `layer` to zero `frontier` today, and the message
  cannot be. The title also named no type, so "the Pieri frontier" invites the
  leading-edge reading — the third sense that got the word retired — where
  "the Pieri layer is a β-mask map, not a `Schur<C>`" cannot be misread.
  Prefer the word the tree will still use, name the type, and expect retired
  vocabulary to survive in the log: messages are the one prose surface a sweep
  cannot reach.

## Voice

- **Impersonal in the reference** ("Returns every ν with…"); "we" is welcome
  in the record and README, where there is a narrator.
- **Bold the term being defined**, once, at its definition — the house
  pattern in the specs and module docs.
- **Density is a genre property.** The record may be dense; the reference must
  be scannable. In an item doc, every sentence past the contract must survive
  the question "does the caller need this to use the function correctly?" —
  what fails moves up to the module doc or out to the record.
- **No hedging on what a test pins** ("should be", "probably") — say what
  holds and name the pin. **No hype** — a superlative is replaced by its
  number, which is more impressive anyway.
- **Name the consequence, not the category.** "X is load-bearing, not
  cosmetic" asserts that X matters without saying how it fails, and nobody
  says it out loud; the phrase was struck from all 28 sites it had reached.
  Write what breaks instead — "the order matters: `G_ν` is not symmetric in
  its components", "skip the primitive part and the answers leave the ring".
  Same length, and it tells the reader something. The same applies to
  "crucial", "critical", "essential" and "key" standing alone — and to
  "the fallible form of X", which names the return type's shape where the
  reader needs its *meaning*: `try_character`'s `None` is overflow and never
  an input error, which is the one thing "fallible" leaves open.
- **Define a subset by its shape, not by who uses it.** `frac.rs` and
  `afrac.rs` both opened with "ℚ(q,t) — but only the part of it Macdonald
  actually inhabits", which leaves a reader unable to tell whether their
  expression fits. They now name the restriction — "restricted to
  denominators that are products of binomials `1 − qᵃtᵇ`" — and the family
  that motivated it becomes the *reason*, which the paragraph below always
  went on to give anyway. Watch for "the useful part", "what X actually
  needs", and any verb standing in for a definition.
- **The reference is transactional.** An item doc answers what the caller
  must do, in the fewest sentences that answer it — the researcher is
  deciding whether to trust a number, not reading an argument. A `# Panics`
  says what panics, when, and what to call instead; the reasoning behind the
  choice belongs to the policy that made it. `character`'s used to end
  "panicking is the deliberate default: the alternative this replaced was a
  silently wrong answer", which defends a decision against an alternative the
  caller will never meet — [failure.md](policies/failure.md) owns that
  argument, in the row that pairs `try_character` with `character_in`.
- **An ordinal names its siblings.** `bh.rs` opened "this is a third route to
  the same numbers, and the reason to prefer it is measured rather than
  aesthetic" — which asserts the *kind* of reason it has instead of giving
  one, and leaves the reader to guess routes one and two. Say how many there
  are, name them where they live, and give the reason itself. Same for "the
  second engine", "the other implementation", "unlike the older path".
- **Nothing here is a trick.** "The same trick `Frac` plays for ℚ(q,t)",
  "β-number bit tricks", "the bit-trick containment test" — the word tells a
  reader who does not already know the mechanism that there *is* one and
  declines to say what. Name it: the representation, the β-number bit
  arithmetic, the abacus test. Reserve the register for what it fits, which
  in this tree is nothing.
- **"The whole point of X" is emphasis standing in for a sentence.** It had
  reached nineteen sites, one of them stuttering ("the whole point of the
  whole-degree entry points"). The plain forms say the same thing and stop
  shouting: "that is why the bound is `QAlgebra`", "that sharing is what the
  whole-degree entry points are for". Same for "the whole reason".
- **The reference does not narrate itself.** "The reason is worth recording",
  "as one question: did anything that must not happen happen?", "note that" —
  the doc is not a voice with opinions about its own contents. State the
  reason; state the predicate. `laws_hold` now says which two fields it tests
  and that both must be empty, which is what a caller reading it needs.
- **American English** (`normalize`, `memoize`, `summarized`), matching the
  Rust ecosystem's own API vocabulary — and because one spelling is one grep:
  a reader searching `normalize` must not miss `normalise`. The tree currently
  mixes in a few British forms; converge on touch, don't sweep.
- Prose wraps at 80 columns, in `.md` and in doc comments alike.

### No aphorisms, no metaphors

An **aphorism** is a general maxim stated as if self-evident: "a cost model
with a factor missing will rank engines confidently and wrongly", "a verdict
with no premise is a permanent wall". It takes one measured case and states
the universal it suggests. A **metaphor** names the thing as something it is
not: "filing the number off", "a hole the size of a rewrite", "the durable
half", "the shop window". Both compress by assuming the reader supplies what
was left out, and neither is allowed in this tree.

Three costs, and they compound:

- **Neither can be checked.** The case behind the maxim had a workload, a
  measurement, a version; the maxim has none, so nothing fails when it stops
  being true. This is the defect the profile-share rule already names, one
  level up: a claim with the falsifiable part removed.
- **Both read as settled.** A universal sentence is read by a fresh-context
  agent as a rule of the project and applied to cases the original never
  covered. A metaphor is read as a rule it cannot state.
- **A metaphor needs the answer to decode.** "Filing the number off does not
  convert a profile share into a durable fact" works only for a reader who
  already knows the number was a percentage and that deleting it was the
  evasion. "Deleting the percentage does not make the claim portable" says it
  once, to anyone.

**The replacement is the specific claim the figure generalized**, and it is
usually the same length. Not "a cost model with a factor missing will rank
engines confidently and wrongly" but "the node-count model ranked E3 above E2
because it omitted the size of the running element, which E2 never inflates".
Not "the record is the working agent's only long-term memory" but "a dead end
recorded nowhere is re-explored by the next session at full price".

**The one carve-out is the defined term.** `seam`, `atom`, `layer`, `peel`,
`ladder`, `wall`, `ceiling`, `lex-monic` began as metaphors and stopped being
metaphors when each was pinned to a single meaning at a single place — that
pinning is the mechanism this file already uses, and [One word, one
meaning](#one-word-one-meaning-across-the-tree) is where a new one earns its
place. An undefined figure is not a term; it is a metaphor with ambitions.
Headings get no carve-out: "the convention minefield" was a heading and a
metaphor, and is now "The conventions in circulation".

## What this changes

Current practice already embodies most of this guide; these are the deltas,
each deliberate:

1. **Doctests everywhere that matters.** The crate has exactly one doctest
   today (the record's test census). Every public entry-point family gets a
   convention-pinning `# Examples` block. This is the largest ask, and the
   highest-value one: it is the only documentation CI executes, and it is the
   researcher's trust check and the agent's runnable convention pin in one
   block.
2. ~~**`# Panics` / `# Errors` sections** wherever a public function can.~~
   **Done** — 46 added, and every `pub fn` in `src/` outside `python.rs` that
   can panic now carries one; it was four when this delta was written. The
   sweep found a public accessor returning a wrong answer rather than
   panicking at all. Chapter in
   [record/failure-and-overflow.md](record/failure-and-overflow.md).
3. **Backticked math and resolving citation keys**, which retires ~160 of the
   173 `cargo doc` warnings and turns every `[KEY]` into a working link
   (release-readiness Phase 1 counted them; this guide makes the fix the
   standing rule, not a one-time cleanup).
4. ~~**lib.rs is rewritten as reference.**~~ **Done** — the crate front page
   describes the present and links the record; the embedded roadmap, stale by
   four shipped features, is gone. Rule 3 keeps the exhibit, because the cost
   it records was paid.
5. **Seconds and ×-ratios leave public rustdoc** for the record, which owns
   them with harness and caveats. Deterministic counts (allocations, bytes,
   asymptotics) stay.
6. ~~**One spelling** (American), where today there are two.~~ **Done** —
   138 British forms across 44 files, and `scripts/check_spelling.py` keeps
   the count at zero. The delta was written as a small tidy; it was the
   largest single prose defect in the tree, and it grew for as long as it had
   no gate. A hand sweep of the same rule, one session earlier, found four.
7. **Range and capability move to the front of module docs.** The README and
   the best modules already lead with what a family opens and how far it
   runs; make that uniform, with a measured incumbent survey — today
   research-gaps.md; at release, wherever its walls are re-homed — behind
   every "no other package" claim.
8. **Item-doc sentences get a cap and a gate.** Every `# Panics` in `src/`
   opened with a verbless fragment; all 54 now lead with the verb, and
   `scripts/check_doc_sentences.py` keeps it that way. The 25-word cap is the
   standing rule with 685 sentences still over it, held from behind by a
   ratchet that only comes down. See [Sentence discipline in item
   docs](#sentence-discipline-in-item-docs) for what was taken from ASD-STE100
   and what was left.

Enforcement is cheap and already planned: `cargo doc --no-deps --all-features`
gated at `-D warnings` in CI (release-readiness Phase 0/1) covers rules the
compiler can see; the checklist below covers the rest at review time.

## The checklist, before a `pub` item ships

- [ ] First sentence is complete, stands alone in a listing, and names the
      object by its precise, searchable name.
- [ ] Contract states result order/zero-freeness, requirements on arguments,
      and behavior at the degenerate inputs (∅, degree 0, equal shapes…).
- [ ] `# Panics` if it can; `# Errors` if it returns `Result` — each opening
      with its verb, not with `If`.
- [ ] No sentence carries two facts, and none runs past 25 words
      (`scripts/check_doc_sentences.py`).
- [ ] A doctest pins the convention with a hand-checkable value.
- [ ] Every formula cites `[KEY] (eq)`, and every `[KEY]` resolves.
- [ ] The Sage equivalent is named, or its absence stated.
- [ ] Range is stated in reproducible terms; a "no other package" claim, where
      true, is made — and points at the measured survey that backs it.
- [ ] No seconds, no ×-ratios, no future work — those link to the record.
