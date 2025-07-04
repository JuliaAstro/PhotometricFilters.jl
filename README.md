# PhotometricFilters

[![Build Status](https://github.com/juliaastro/PhotometricFilters.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/PhotometricFilters.jl/actions)
[![Coverage](https://codecov.io/gh/juliaastro/PhotometricFilters.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/PhotometricFilters.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.org/PhotometricFilters/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/PhotometricFilters.jl/dev)

## Installation

## Usage

```julia
using PhotometricFilters
using PhotometricFilters: SDSS_u, SDSS_g, SDSS_r, SDSS_i, SDSS_z, fwhm

filts = [SDSS_u(), SDSS_g(), SDSS_r(), SDSS_i(), SDSS_z()]
```

plotting works out of the box

```julia
using Plots, ColorSchemes
plot(filts, palette=palette(:magma, 6), fill=(0, 0.2))
```

![](sdss.png)

Using [Unitful.jl](https://github.com/painterqubits/Unitful.jl) is built in to all functionality

```julia

julia> filt_units = SDSS_u(units=true);

julia> fwhm(filt_units)
600.0 Å
```

## Citations
