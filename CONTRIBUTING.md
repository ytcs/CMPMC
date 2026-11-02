# Contributing

## Build and test

```bash
cmake -B build -DCMPMC_WARNINGS_AS_ERRORS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

`ctest -L <label>` runs a subset: `unit`, `validation`, `phonon`, `charge`,
`superconductor`, `materials`, `transport`, `external-data`.

## Code layout and conventions

- The library is header-only C++20 under `include/cmpmc/`. Module boundaries and
  dependency direction are described in [docs/architecture.md](docs/architecture.md).
- Internal units are SI, with energies in eV. Write literals with `cmpmc::units`
  constants (`3.0 * units::GPa`).
- Physics kernels are free functions or plain structs without engine dependencies. The
  transport loop lives only in `transport/engine.hpp`.
- Every header must compile on its own; the build checks this (`cmpmc_header_check`).
- Hot paths must not allocate, throw or use virtual dispatch.
- Comments explain the physics and cite the source in the form
  `Author, arXiv:XXXX.XXXXX, Eq. (n)` or a journal reference. Comments do not
  describe development history; that belongs in commit messages and `CHANGELOG.md`.

## Physics changes

Any change to a model, default or constant needs all of the following:

1. **A source.** A published reference, a closed-form derivation written out in the
   code comment, or a measured value with its origin.
2. **A test** in `tests/` that compares the implementation with that reference and has
   an explicit tolerance. Register it in `tests/CMakeLists.txt` with labels.
3. **An entry** in [docs/validation.md](docs/validation.md). If the model is not yet
   validated, guard it with `require_validated` and list it in
   [docs/experimental.md](docs/experimental.md) and [docs/limitations.md](docs/limitations.md).
4. **A changelog line** in [CHANGELOG.md](CHANGELOG.md) when the change alters results.

## Tests

- Use the helpers in `tests/test_support.hpp`: `check(name, value, reference, rel_tol,
  abs_tol)`, `info(name, value)`, `result()`, and `skip(reason)` for tests that need
  external data.
- Use fixed seeds. Keep each test under a minute in a Release build where possible.

## Commits and pull requests

- Keep each commit to one logical change, with a message that says what changed and why.
- Update the documentation in the same pull request as the code.
