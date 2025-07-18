```@meta
CurrentModule = PhotometricFilters
```

# PhotometricFilters

This package provides access to, and operations on, photometric filter curves. Such filter curves are defined by a filter's transmission as a function of wavelength. Transmission and wavelength vectors are therefore the foundation of a filter curve, but it is also important to note whether the filter is used for photon counter or energy counter detectors, as the integrals used to calculate statistics over a filter curve are different between these two types of detectors.

## Types

All photometric filter types should be subtypes of the [`AbstractFilter`](@ref PhotometricFilters.AbstractFilter) type. We define a minimal API that should be followed so that new types can make use of the generic filter operations we define.

The simplest concrete filter type we provide to represent photometric filters is [`PhotometricFilter`](@ref).

```@docs
PhotometricFilter
```

Users can construct their own filter curvers from raw data using this type.

## Accessing Filter Curves

We provide a modest collection of filter curves through a data dependency. The available filter curves are accessible via the `FILTER_NAMES` module constant,

```@example 1
using PhotometricFilters
PhotometricFilters.FILTER_NAMES |> println
```

These included filter curves can be accessed like so,
```@example 1
using PhotometricFilters: SDSS_u, SDSS_g, SDSS_r, SDSS_i, SDSS_z, fwhm
filts = [SDSS_u(), SDSS_g(), SDSS_r(), SDSS_i(), SDSS_z()]
```

**NOTE THAT THESE INCLUDED FILTER CURVES ARE NOT GUARANTEED TO BE UP-TO-DATE.** If you are using a filter/instrument that may have recently had its filter curves updated (e.g., JWST/NIRCAM), you should use our SVO query interface to make sure you get the most up-to-date filter curves. SVO also provides additional metadata that is useful for some applications (e.g., filter zeropoints). If you know the SVO-designated name of the filter you want, you can use [`get_filter`](@ref) to retrieve its transmission data, which returns an instance of [`SVOFilter`](@ref PhotometricFilters.SVOFilter).

```@docs
get_filter
PhotometricFilters.SVOFilter
```

If you'd like to perform a search on the filters available through the SVO filter service, you can use [`query_filters`](@ref).

```@docs
query_filters
```

## Supported Operations
We include functions for performing many common operations on photometric filters, summarized below.

### Accessors
```@docs
name
wave
throughput
detector_type
```

### Applying Filter Curves to Spectra

```@docs
apply
apply!
integrate
F_nu
F_lambda
```

### Statistics
```@docs
PhotometricFilters.reference_wavelength
PhotometricFilters.mean_wavelength
PhotometricFilters.central_wavelength
PhotometricFilters.effective_wavelength
PhotometricFilters.pivot_wavelength
PhotometricFilters.min_wave
PhotometricFilters.max_wave
PhotometricFilters.fwhm
PhotometricFilters.width
```

## Internals
```@docs
PhotometricFilters.AbstractFilter
PhotometricFilters.get_metadata
```

## Index

```@index
```
