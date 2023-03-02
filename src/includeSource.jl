
using DrWatson
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using Printf
using Convex
import SCS

const MOI = Convex.MOI

include(srcdir("Enumerations.jl"))
include(srcdir("SCvxParams.jl"))
include(srcdir("ScaleParams.jl"))
include(srcdir("generateGuess.jl"))
include(srcdir("computeDiscreteMatricies.jl"))
include(srcdir("flowMap.jl"))
include(srcdir("defectCost.jl"))
include(srcdir("printStatus.jl"))