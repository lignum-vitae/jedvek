# Jedvek

[![crates.io](https://img.shields.io/crates/v/jedvek.svg)](https://crates.io/crates/jedvek)
[![docs.rs](https://docs.rs/jedvek/badge.svg)](https://docs.rs/jedvek)
[![Build Status](https://github.com/lignum-vitae/jedvek/workflows/Build%20and%20test/badge.svg)](https://github.com/lignum-vitae/jedvek/actions?workflow=Build%20and%20test)
[![Build Status](https://github.com/lignum-vitae/jedvek/workflows/Clippy%20check%20-%20lint/badge.svg)](https://github.com/lignum-vitae/jedvek/actions?workflow=Clippy%20check%20-%20lint)
[![GitHub stars](https://img.shields.io/github/stars/lignum-vitae/jedvek?style=social)](https://github.com/lignum-vitae/jedvek/stargazers)

Jedvek /jɛdvɛk/ [derived from the Czech jednorozměrný (one-dimensional) + vektor (vector)] is a
Rust library for multidimensional data structures backed by a contiguous 1D memory footprint.

## Installation

Add Jedvek as a dependency in your Cargo.toml:

### Via Crates.io (Recommended)

To use the latest stable release, add jedvek as a dependency:

```toml
[dependencies]
jedvek = "0.X.X" # Always use the latest version available on crates.io
```

Or, use Cargo:

```nginx
cargo add jedvek
```

### Via Git Repository (Latest Development Build)

```toml
[dependencies]
jedvek = { git = "https://github.com/lignum-vitae/jedvek.git" }
```

Then run:

`cargo build`

## Matrices

The Matrix2D struct provides a convenient way to use a two-dimensional matrix that is
a one-dimensional vector under the hood.
This allows for efficient memory usage while enabling standard matrix operations.

Matrix2D implements a wide range of methods and traits, including:

- Linear Algebra: dot product, inverse, transpose.
- Arithmetic: Implementations of multiplication by matrices, vectors, and scalars and division by scalars
- Manipulation: shape, size, full, reshape, map
- Conversion: from_flat, From, TryFrom
- Utility: new, max, min, is_empty

## Project layout

Within these crates, the following modules `jedvek::<module name>` are provided

| Module          | Description                                                                                            |
| --------------- | ------------------------------------------------------------------------------------------------------ |
|        -        | `Matrix2D` struct and implementations                                                                  |
| `substitution`  | Backward and forward substitution algorithms                                                           |
| `decomposition` | Decomposition algorithms including LU decomposition and LU decomposition with partial pivoting         |

### Running Examples

Working examples of the available algorithms as well as a full list of available algorithms can be found in the
[`examples/`](https://github.com/lignum-vitae/jedvek/tree/main/examples) directory.

Run any example with the following command:

`cargo run --example <example_name>`

Do not include `.rs` when running examples.

## Contributing

We welcome contributions! Please read our:

- [Code of Conduct](https://github.com/lignum-vitae/jedvek/blob/main/docs/CODE_OF_CONDUCT.md)
- [Contribution Guidelines](https://github.com/lignum-vitae/jedvek/blob/main/docs/CONTRIBUTING.md)

> [!NOTE]
> Before submitting a PR, install [just](https://github.com/casey/just) and run `just check`
> to pull the latest changes from the main branch as well as to format, test, and lint your code.
> Just can be installed using `cargo install just`, curl, or your favourite package manager.

Stay connected via our **[Discord Server](https://discord.gg/PdVZCtcgaH)**

## Stability

This project is in the alpha stage. APIs may change without warning until version
1.0.0.

## Special Mentions

A special thank you to [Mittei](https://github.com/mittei) for his contributions
to the [Spindalis project](https://github.com/lignum-vitae/spindalis).
The `Arr2D` implementation he authored serves as the structural backbone for
Jedvek’s `Matrix2D` struct.
