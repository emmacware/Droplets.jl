using Droplets
using Test

@testset "Droplets.jl" begin
    FT = Float64
    Ns = 2
    scale = Ns * (Ns - 1) / 2 / (Ns / 2)
    coagsettings = coag_settings{FT}(Ns=Ns,scale=scale,Δt=1.0)

    ξ = [Int(2),Int(3)]
    R = [FT(1.0),FT(2.0)]
    X = 3pi/4 .*R.^3
    drops = droplet_attributes{FT}(ξ, X)
    coag_data = coagulation_run{FT}(Ns)

    @test golovin(drops, (1,2), coagsettings) == 1500*(X[1]+X[2])

    pair_Ps!(1, (1,2), drops,coag_data, golovin,coagsettings)
    @test coag_data.pαdt[1] == 3*1500*(X[1]+X[2])

    coag_data.pαdt[1]=1
    coag_data.ϕ[1]=0.5
    lowest_zero = Ref(false)
    sdm_update!((1,2),1, drops,coag_data)

    @test drops.ξ == [2,1]

end

# @testset "Golovin Test" begin
#     include("golovintest.jl")
# end

@testset "Coalescence" begin
    include("coalescence_tests.jl")
end

@testset "Condensation" begin
    include("condensation_tests.jl")
end

@testset "Turbulence" begin
    include("turbulence_tests.jl")
end

@testset "MPDATA" begin
    include("mpdata_tests.jl")
end

@testset "SCM" begin
    include("scm_tests.jl")
end