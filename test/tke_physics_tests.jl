using Test
using Droplets


FT = Float64

function make_grid(nz, dz; θ_profile=nothing, qv=nothing, e=nothing, u=nothing, v=nothing)
    drops = droplet_attributes_1d(Int[], FT[], FT[], FT[], Int[], FT[], fill(1:0, nz), Int[])
    spatialset = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz * dz))
    grid = create_scm_grids(nz, dz, drops, spatial=spatialset, output=nothing)

    P_ref = FT(91500.0)
    grid.states.P .= P_ref
    grid.states.qv .= something(qv, fill(FT(8e-3), nz))
    grid.states.θ .= something(θ_profile, fill(FT(289.0), nz))
    grid.states.e .= something(e, fill(FT(0.1), nz))
    grid.states.ql_tmp .= zeros(FT, nz)
    grid.states.ρ .= ρ_calc_θ.(P_ref, grid.states.θ, grid.states.qv, constants)
    grid.wind.u .= something(u, zeros(FT, nz))
    grid.wind.v .= something(v, zeros(FT, nz))

    T_from_theta!(grid.states.T_tmp, grid.states.θ, grid.states.P, grid.states.qv, constants)
    grid.states.θl_tmp .= θl.(grid.states.P, grid.states.T_tmp, grid.states.ql_tmp, grid.states.qv, constants)
    return grid
end

@testset "saturation_specific_humidity" begin
    c = Constants{FT}()
    P = FT(95000.0)
    qv = FT(8e-3)
    θ = FT(290.0)
    T = T_from_theta(θ, P, qv, c)

    qsat = Droplets.saturation_specific_humidity(θ, P, qv, c)
    es = esat(T)

    # sat() hardcodes 0.622/0.378 while saturation_specific_humidity uses
    # constants.ϵ = Rd/Rv = 0.622559... 
    @test_broken sat(qsat, P) ≈ es rtol=1e-10

    # warmer air holds more vapor at saturation
    qsat_cold = Droplets.saturation_specific_humidity(FT(280.0), P, qv, c)
    qsat_warm = Droplets.saturation_specific_humidity(FT(300.0), P, qv, c)
    @test qsat_warm > qsat_cold
end

@testset "my25_stability_functions: neutral case reduces to the bare A-coefficients" begin
    nz, dz = 10, 50.0
    grid = make_grid(nz, dz)  # uniform θ, zero wind -> N2=0, S2=0 -> GH=GM=0
    l = fill(FT(30.0), nz)
    K_h = zeros(FT, nz); K_m = zeros(FT, nz); K_e = zeros(FT, nz)
    GH = zeros(FT, nz); GM = zeros(FT, nz); SM = zeros(FT, nz); SH = zeros(FT, nz)
    tke = tke_settings{FT}()

    k = 5
    my25_stability_functions(l, K_h, K_m, K_e, GH, GM, SM, SH, k, tke, grid, constants)

    @test GH[k] ≈ 0 atol=1e-8
    @test GM[k] ≈ 0 atol=1e-8
    # with GH=GM=0: A10 = 1, SH = A1, SM = A4 (the level-2.5 neutral limit)
    @test SH[k] ≈ 0.74 rtol=1e-6
    @test SM[k] ≈ 0.699 rtol=1e-6

    q = sqrt(2 * max(grid.states.e[k], tke.e_min))
    @test K_m[k] ≈ l[k] * q * SM[k] rtol=1e-10
    @test K_h[k] ≈ l[k] * q * SH[k] rtol=1e-10
    @test K_e[k] ≈ 0.2 * l[k] * q rtol=1e-10
    @test K_m[k] > 0 && K_h[k] > 0 && K_e[k] > 0
end

@testset "my25_stability_functions: GH is always clamped within tke.GH_lims" begin
    nz, dz = 10, 50.0
    tke = tke_settings{FT}()
    θ_unstable = FT[300.0 - 5.0 * k for k in 1:nz]
    grid = make_grid(nz, dz; θ_profile=θ_unstable, e=fill(FT(1e-4), nz))
    l = fill(FT(50.0), nz)
    K_h = zeros(FT, nz); K_m = zeros(FT, nz); K_e = zeros(FT, nz)
    GH = zeros(FT, nz); GM = zeros(FT, nz); SM = zeros(FT, nz); SH = zeros(FT, nz)

    for k in 2:nz-1
        my25_stability_functions(l, K_h, K_m, K_e, GH, GM, SM, SH, k, tke, grid, constants)
        @test tke.GH_lims[1] <= GH[k] <= tke.GH_lims[2]
        @test isfinite(K_h[k]) && isfinite(K_m[k])
    end
end

@testset "tke_update!: zero-production decay matches the governing ODE's analytic solution" begin
    # check semi-implicit...
    nz, dz = 1, 40.0
    e0 = FT(1.0)
    grid = make_grid(nz, dz; e=fill(e0, nz))
    tke = tke_settings{FT}()
    turbdata = turbulence_data(FT, nz)
    l = fill(FT(20.0), nz)
    zeros_nz = zeros(FT, nz)

    T = FT(50.0)
    dt = FT(0.01)  # small relative to the decay timescale for the semi-implicit scheme to track the ODE closely
    for _ in 1:round(Int, T / dt)
        tke_update!(l, zeros_nz, zeros_nz, zeros_nz, zeros_nz, grid, tke, dt, constants, turbdata)
    end

    k = 2 * sqrt(FT(2)) * tke.my_diss / l[1]
    e_analytic = (e0^(-FT(0.5)) + (k / 2) * T)^(-2)
    @test grid.states.e[1] ≈ e_analytic rtol=1e-4
end

@testset "tke_update!: zero production decays TKE monotonically toward zero" begin
    nz, dz = 4, 40.0
    grid = make_grid(nz, dz; e=fill(FT(1.0), nz))
    tke = tke_settings{FT}()
    turbdata = turbulence_data(FT, nz)
    l = fill(FT(20.0), nz)
    zeros_nz = zeros(FT, nz)
    dt = FT(10.0)

    e_prev = copy(grid.states.e)
    for _ in 1:10
        tke_update!(l, zeros_nz, zeros_nz, zeros_nz, zeros_nz, grid, tke, dt, constants, turbdata)
        @test all(grid.states.e .<= e_prev .+ 1e-12)  # pure dissipation: monotone decreasing
        e_prev = copy(grid.states.e)
    end
    @test all(grid.states.e .>= 0)
end

@testset "bott1997term: zero for a vertically uniform profile" begin
    nz, dz = 8, 50.0
    grid = make_grid(nz, dz)  # uniform θ, qv, ql=0 everywhere
    for k in 2:nz-1
        @test Droplets.bott1997term(grid, k, constants) ≈ 0 atol=1e-10
    end
end

@testset "bott_mixing_length!" begin
    nz, dz = 12, 50.0
    tke = tke_settings{FT}()
    # qv drops below 0.008 above k=6 -> inversion at k=6
    qv = FT[k <= 6 ? 0.010 : 0.003 for k in 1:nz]
    grid = make_grid(nz, dz; qv=qv)
    l = zeros(FT, nz)
    Droplets.bott_mixing_length!(l, grid, tke, constants)

    @test all(isfinite, l)
    @test all(l .> 0)
    # above the inversion, l0 = dz -- caps the mixing length close to dz there
    @test l[nz] <= dz + 1e-6
    # k==1 isn't floored at dz (per the source's `k == 1 ? lf : max(lf, dz)`);
    # every other level is floored at dz
    @test all(l[2:end] .>= dz - 1e-9)
end

@testset "edmfx_mixing_length!: l_W (Blackadar neutral-scaling) formula, no-TKE/no-stratification limit" begin
    nz, dz = 10, 50.0
    tke = tke_settings{FT}()
    # zero TKE and zero wind/buoyancy forcing: l_TKE forcing term a_pd = c_m*(S2-N2)*sqrt(e) = 0
    # (since e=0), so l[k] must reduce to min(l_W, min(dz,z)) via the l_N=z fallback
    grid = make_grid(nz, dz; e=zeros(FT, nz))
    l = zeros(FT, nz)
    Droplets.edmfx_mixing_length!(l, grid, tke, constants)

    vk = tke.vk
    for k in 1:nz
        z = grid.centers_z[k]
        l_W = vk * z * tke.l_inf / (vk * z + tke.l_inf)
        @test l[k] ≈ max(min(l_W, l_W, z), min(dz, z)) rtol=1e-10
        @test l[k] > 0
    end
    @test all(l .<= tke.l_inf + 1e-9)
end
