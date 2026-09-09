<!-- The commit message is the report: what was wrong, what changed, how it
     was measured, what was rejected, what is open. Say here only what the
     reviewer needs beyond it. -->

## Checks

The list is CLAUDE.md's definition of done; delete the lines that do not apply.

- [ ] `scripts/preflight.sh` is green.
- [ ] `scripts/preflight_python.sh` is green (anything under `src/python.rs`,
      `python/symfn/` or `docsite/` changed).
- [ ] A new `pub` item meets the checklist at the end of `docs/style.md`:
      summary sentence, contract, degenerate inputs, `# Panics`, a
      convention-pinning doctest, resolving citations.
- [ ] A new failure path fits a row of the mechanism table in
      `docs/policies/failure.md`.
- [ ] A new or changed `#[pyfunction]` sits in a row of the home table in
      `docs/policies/python.md`, and `symfn.pyi` follows it.
- [ ] A new family or engine carries the evidence its row of the table in
      `docs/policies/validation.md` requires.
- [ ] Every measurement is in the record file that owns the subsystem, with
      its harness and power state, and none is in rustdoc.
- [ ] Anything this closes or opens in `docs/todo-1.0.md` or a record file's
      open tail is changed in the same pull request.
