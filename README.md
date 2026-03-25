# Tempological Control
Julia source files and scripts for the paper [Tempological Control of Network Dynamics](https://arxiv.org/abs/2510.10926) (Zhang, Rock & Cornelius, 2026).

## System requirements
You will need an installation of the [Julia](https://julialang.org/) programming language. This code was tested with Julia 1.10.9 on a 2019 Intel Core i9 iMac computer running macOS Sequoia (15.7.2).

All Julia dependences are listed in the `Project.toml` file.

## Installation guide
After cloning this repository, navigate to the corresponding folder and start julia in project mode:
```bash
julia --project
```
Then, within julia, execute the following to install all dependencies:
```julia
using Pkg
Pkg.Instantiate()
```
Installation may take approximately 10 minutes.

## Demo
In your terminal, navigate to the `scripts` folder in this repository and start julia in project mode:
```bash
julia --project --threads=auto
```
Within julia, run the `demo.jl` script by including it:
```julia
include("demo.jl")
```
Running this demo will take approximately 3 minutes using Julia with four (4) threads.

### Expected output:
Thsi script replicates Fig. 2 of the main text. I.e. the tempological control success rate as a function of $m$, for $n = 10$, $q = 0.44$, and various values of $p$.

## Instructions for use
The  `demo.jl` script demonstrates the main functionality of the package. I.e. Kuramoto snapshot generation, switched system creation, and the tempological control algorithm itself.
