#!/bin/bash
export PATH="$HOME/.juliaup/bin:$PATH"
julia --project=. test/test_vulnerability_regime_outputs.jl > test_out.log 2>&1
