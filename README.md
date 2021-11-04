# Causal Set MCMC

This code performs MCMC simulations with $d=2$ causal sets representing Minkowski space, with the aim of reproducing the results of [https://arxiv.org/abs/1110.6244](1110.6244).
The code is not thoroughly tested, so use at your own risk.

## Installation

Clone this repository and run `cargo build --release`. Notice that the speed difference between the development and the release version is significant, so make sure to always use the release version for any significant computation.

## Usage

Compiling produces a binary `causal_sets`. Run this binary with a data directory, e.g. `causal_sets data/some_run`. The code expects to find a file `parameters.yml` in this directory. For an example see `resources/parameters.yml`.
The current version of the phase diagram computed from this file can be seen below.

![Phase Diagram](/resources/phase_diagram.png)