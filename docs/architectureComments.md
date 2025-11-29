# General comments on the code architecture

## Repeated code and lack of common Julia features

This code was developed for people who are not familiar with software development, so I purposely avoided using this like modules, abstract types, and the other peculiar facets of Julia that are not easily understandable (although I made an exception with the use of `Union` in `GmshVtkIO.jl`). As a side effect, there's a good amount of repeated code here (`src/2D`, `src/3D`, and `src/3DMembrane`) that is far less than ideal. 

Unfortunately, by avoiding modules and using the `include()` functionality seems to make code completion difficult to configure, and I never managed to get it working well with this codebase. If you manage to get it working, I'll greatly appreciate it if you post a GitHub issue for documentation.

## Inefficient `MaterialPoint` data structure

You may notice that each `MaterialPoint` and `SurfacePoint` keeps a copy of the vertex list and `updateMaterialVertices!()` and `updateSurfaceVertices!()` are rather inefficient. It is more efficient to keep a separate `vertices` array and instead let each `Point` data structure keep a much smaller array that simply indexes into the `vertices` array.

This inefficiency is simply a relic from [Sinai's MPM implementation](https://github.com/vinhphunguyen/MPM-Julia) that we referenced early into our research, and I never got around to refactoring it to be more memory efficient. Although, for very large simulations, there *could* be some caching behavior that makes it faster at the cost of the larger memory footprint, since the relevant vertex data is spatially closer to each `Point` instead of being in a separate data structure. No clue, and I have no desire to test this.

## StaticArrays

Julia is an interesting language with powerful and fast capabilities. However, at on version 1.9+, it's not trivial to control where data is stored in memory. You can read more about it in the [Julia documentation](https://docs.julialang.org/en/v1/manual/performance-tips), but basically, you want to prevent as many unnecessary heap allocations as possible. Hence, we have `Util.jl` containing our custom fixed-sized vector/matrix definitions.
