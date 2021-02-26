# PhotometricFilters

[![Build Status](https://github.com/juliaastro/PhotometricFilters.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/PhotometricFilters.jl/actions)
[![Coverage](https://codecov.io/gh/juliaastro/PhotometricFilters.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/PhotometricFilters.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/PhotometricFilters.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/PhotometricFilters.jl/dev)

**WIP** @mileslucas

## Installation

## Usage

```julia
using PhotometricFilters
using PhotometricFilters: SDSS_u, SDSS_g, SDSS_r, SDSS_i, SDSS_z

filts = [SDSS_u(), SDSS_g(), SDSS_r(), SDSS_i(), SDSS_z()]
combined = sum(filts)
```

plotting works out of the box

```julia
using Plots
plot(wave(combined), throughput(combined) .+ 0.1, label="")
plot!(filts)
```

![](sdss.png)

Using [Unitful.jl](https://github.com/painterqubits/Unitful.jl) is built in to all functionality

```julia

julia> filt_units = SDSS_u(units=true);

julia> fwhm(filt_units)
600.0 Å
```

## Citations
