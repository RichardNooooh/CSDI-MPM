using Base.Threads
using Folds
using StaticArrays
using LinearAlgebra
using WriteVTK
using ProgressMeter
using Dates
using Logging
# using BenchmarkTools

include("../.venv/lib64/gmsh.jl")

# Small cutoff value near zero for small-mass or floating point comparisons
const EPSILON = 1e-12

# Custom defined SVector data types
const Vec2{T} = SVector{2, T}
const Vec3{T} = SVector{3, T}
const Vec4{T} = SVector{4, T}
const Vec5{T} = SVector{5, T}
const Vec6{T} = SVector{6, T}
const Vec7{T} = SVector{7, T}
const Vec8{T} = SVector{8, T}
const Vec16{T} = SVector{16, T}
const Vec32{T} = SVector{32, T}

# Custom defined SMatrix data types (last number required for stack allocation)
const Mat22{T} = SMatrix{2, 2, T, 4}
const Mat33{T} = SMatrix{3, 3, T, 9}

# NOTE: these are immutable.
# Commonly used vectors and matrices
const ZERO_VEC2 = Vec2{Float64}(0.0, 0.0)
const ZERO_VEC3 = Vec3{Float64}(0.0, 0.0, 0.0)
const IDENTITY_MAT22 = Mat22{Float64}(1.0, 0.0, 
                                      0.0, 1.0)
const ZERO_MAT22     = Mat22{Float64}(0.0, 0.0, 
                                      0.0, 0.0)
const IDENTITY_MAT33 = Mat33{Float64}(1.0, 0.0, 0.0, 
                                      0.0, 1.0, 0.0, 
                                      0.0, 0.0, 1.0)
const ZERO_MAT33     = Mat33{Float64}(0.0, 0.0, 0.0, 
                                      0.0, 0.0, 0.0, 
                                      0.0, 0.0, 0.0)

"""
A lookup table for the VTK index and its corresponding position in space
"""
struct VTKLookup
    ordered :: Vector{Int64}   # VTK expects ascending order
    to_pos  :: Vector{Int64}   # fast Int→Int map from vtk_idx to lookup index, 0 means "tag not present"
    x       :: Vector{Float64}       # coordinates – overwritten every output step
    y       :: Vector{Float64}
    z       :: Vector{Float64}
end

function VTKLookup(tags::Vector{Int64})::VTKLookup
    ordered = sort!(collect(tags))
    to_pos  = zeros(Int, maximum(ordered))
    for (i, tag) in enumerate(ordered)
        to_pos[tag] = i
    end
    n = length(ordered)
    VTKLookup(ordered, to_pos, zeros(n), zeros(n), zeros(n))
end

# Gets the unique values of a sorted array, `a`, and places them at the beginning of the array in order.
# Assumes `a` is sorted. Returns the last index that was written to. Indices after that is considered invalid.
function arrangeUniqueInPlace!(a::Vector{UInt32})::UInt8
    write_index::UInt8 = 1
    read_index::UInt8 = 2

    while read_index <= length(a)
        if a[read_index] > a[write_index]
            write_index += 1
            a[write_index] = a[read_index]
        end
        read_index += 1
    end
    return write_index
end

# only appends to a file
function recordDataOverTime(file_name::String, t::Float64, data_points::AbstractVector{Float64})::Nothing
    open(file_name, "a") do f
        data::String = ""
        for data_point::Float64 in data_points
            data *= " $data_point"
        end
        println(f, t, data)
    end
    return nothing
end
