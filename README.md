# local-fdr-interpretation-paper

Reproduction code for the paper:

> "Interpretation of local false discovery rates under the zero-assumption", Daniel Xiang, Nikolaos Ignatiadis, Peter McCullagh

The main code for the reproduction is available in the file `interpretation_sim.jl`. The code is written in Julia, however, several R package are called through the `RCall.jl` package. Precise specifications of all packages used within Julia are provided in the `Project.toml` and `Manifest.toml` files, and for R in the `R_session_info.txt` file.