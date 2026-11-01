# Canon tables

The Rust crate [bigcurve](https://github.com/hypertidy/bigcurve-rs) is
the canonical bigcurve engine; `src/bigcurve.h` (frozen, bug fixes only)
follows it. `canon-tables.txt` is written by bigcurve-rs
(`tests/canon.rs`, see its header for the format) and copied here
unchanged. Each case has its inputs, settings, every forward call the
Rust engine made (analytic maps on the sphere, so no PROJ is needed) and
the output.

`check.cpp` replays those calls through the header, driven as
`src/densify.cpp` drives it, and requires the same output: bitwise on
Linux with glibc, elsewhere the same vertex counts and segments with
coordinates to 1e-9. Cases with Planar edges or ellipsoid placement are
Rust-only for now (the R package gets them when it moves onto the
crate, #16) and are counted, not checked.

    c++ -std=c++17 -O2 -o canon/check canon/check.cpp && ./canon/check

The header turns floating-point contraction off with a pragma: with
fused multiply-adds (GCC with `-march=x86-64-v3`, clang on arm64), 88 of
the 164 cases differed. An explicit `-ffp-contract=fast` overrides the
pragma in clang; don't build the engine with it.

After a deliberate change in bigcurve-rs: regenerate there
(`CANON_WRITE=1 cargo test --all-features --test canon`), copy
`tests/canon/canon-tables.txt` here, and fix the header until the check
is green. This directory is not part of the R package build.
