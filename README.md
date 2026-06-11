# PHANTOM_oct
Profile &amp; Halo Analysis for Numerous Theoretical Dark Matter


## Installation (Octave)

**Requirements:** Octave >= 4.5.0

Install directly from GitHub:
```octave
pkg install "https://github.com/matc-thaher/phantom-octave/archive/refs/heads/main.zip"
pkg load phantom-octave
```

For the fitting functions, also install the `optim` package:
```octave
pkg install -forge optim
pkg load optim
```

For `colossus_query`, Octave >= 7.0 is required (needs `jsonencode`).


## Documentation

Full documentation for all functions is available on the
[PHANTOM Wiki](https://github.com/matc-thaher/PHANTOM/wiki).

This Octave package shares the same API and function signatures
as the MATLAB version. All usage examples and parameter descriptions
apply directly.
