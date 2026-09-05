# Security

symfn computes symmetric functions. It does no I/O, opens no network
connection, reads no environment variable, and runs no code it did not ship
with. The default build has no dependencies; the `python` feature adds four
crates, every one checked by `cargo deny` on each push. Two `unsafe` blocks
exist, both reviewed and each stating its obligations beside it
(`docs/record/failure-and-overflow.md`, "The unsafe review").

## Supported versions

The latest release on the `0.9` line. Release candidates and anything built
from a branch are not supported.

## Reporting

Report a vulnerability privately, through the repository's Security tab on
GitHub ("Report a vulnerability"), so that a fix can ship before the report is
public. There is one maintainer; you will get an answer, and it may take
longer than a week.

## What is and is not a vulnerability here

In scope: anything that lets a caller corrupt memory, read memory it did not
write, or execute code — through the Python boundary above all, since that is
where every consumer reaches the library. A value that crosses the boundary
wrong in a way the caller cannot detect is also in scope.

Not in scope, and welcome as ordinary issues:

- **A wrong value.** Every value that leaves the library is exact or the call
  fails loudly (`docs/policies/failure.md`). A wrong value is a correctness
  bug, the most serious kind this library can have, and it goes to the issue
  tracker with the "A value is wrong" template.
- **A panic on malformed input.** A partition that is not a partition, a
  permutation that is not a permutation: these are contract violations and
  the documented `# Panics` sections say so. Through the Python boundary they
  arrive as typed exceptions, never as a crash of the interpreter — if one
  crashes the interpreter, that is in scope above.
- **Time or memory exhausted by a large input.** The library computes what it
  is asked to; the caches are bounded by a budget the caller sets
  (`cache_budget`, `set_cache_budget`), and the walls of each family are
  stated where the family lives. Nothing here is written for adversarial
  input, and the hashers are not cryptographic.
