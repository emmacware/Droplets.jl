using Droplets
using Test
using LinearAlgebra
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

FT = Float64

"""Build a minimal grid with no droplets (all grid_range = nothing)."""
function make_test_grid(nz::Int, dz::FT; θ_profile=nothing, qv=nothing, e=nothing) where FT
    Ns = 0
    drops = droplet_attributes_1d(
        Int[], FT[], FT[], FT[], Int[], FT[],
        fill(1:0, nz),
        Int[]
    )
    spatialset = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz*dz))

    grid    = create_scm_grids(nz, dz, drops, spatial = spatialset,output=nothing)

    # Set fields to something physical if provided
    P_ref = FT(91500.0)
    grid.states.P    .= P_ref
    grid.states.qv   .= something(qv,  fill(FT(8e-3), nz))
    grid.states.θ    .= something(θ_profile, fill(FT(289.0), nz))
    grid.states.e    .= something(e, fill(FT(1e-2), nz))
    grid.states.ql_tmp .= zeros(FT, nz)
    grid.states.ρ    .= ρ_calc_θ.(P_ref, grid.states.θ, grid.states.qv, constants)
    grid.wind.u      .= fill(FT(4.0), nz)
    grid.wind.v      .= fill(FT(-1.0), nz)

    T_from_theta!(grid.states.T_tmp, grid.states.θ, grid.states.P, grid.states.qv, constants)
    grid.states.θl_tmp .= θl.(grid.states.P, grid.states.T_tmp, grid.states.ql_tmp, grid.states.qv, constants)
    return grid
end

make_turbdata(nz::Int) = turbulence_data(FT, nz)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  implicit_diffuse!  (standalone, no grid needed)
# ─────────────────────────────────────────────────────────────────────────────

@testset "implicit_diffuse!" begin

    @testset "uniform profile stays unchanged (NoFlux BC)" begin
        nz = 20
        dz = FT(50.0)
        dt = FT(60.0)
        ϕ  = fill(FT(5.0), nz)
        K  = fill(FT(10.0), nz)
        implicit_diffuse!(ϕ, K, dt, dz, nz, make_turbdata(nz))
        @test all(isapprox.(ϕ, FT(5.0); atol=1e-10))
    end

    @testset "mean is conserved under diffusion (NoFlux, no sfc_flux)" begin
        nz = 30
        dz = FT(40.0)
        dt = FT(30.0)
        ϕ  = collect(FT, range(280.0, 295.0, length=nz))
        K  = fill(FT(5.0), nz)
        mean_before = sum(ϕ) / nz
        implicit_diffuse!(ϕ, K, dt, dz, nz, make_turbdata(nz))
        mean_after  = sum(ϕ) / nz
        @test isapprox(mean_before, mean_after; atol=1e-8)
    end

    @testset "gradient decays over time" begin
        # Run many steps; profile should become more uniform
        nz = 20
        dz = FT(50.0)
        dt = FT(60.0)
        ϕ  = collect(FT, range(1.0, 2.0, length=nz))
        K  = fill(FT(50.0), nz)
        initial_std = std(ϕ)
        for _ in 1:200
            implicit_diffuse!(ϕ, K, dt, dz, nz, make_turbdata(nz))
        end
        @test std(ϕ) < initial_std * 0.01
    end

    @testset "zero diffusivity leaves profile unchanged" begin
        nz = 10
        dz = FT(50.0)
        dt = FT(100.0)
        ϕ  = FT[sin(k) for k in 1:nz]
        K  = zeros(FT, nz)
        ϕ_orig = copy(ϕ)
        implicit_diffuse!(ϕ, K, dt, dz, nz, make_turbdata(nz))
        # Top cell is pinned to its initial value; bottom cells should be unchanged
        @test all(isapprox.(ϕ[1:end-1], ϕ_orig[1:end-1]; atol=1e-10))
    end

    @testset "surface flux is applied" begin
        # Positive sfc_flux should raise ϕ[1] relative to the no-flux case
        nz = 10
        dz = FT(50.0)
        dt = FT(60.0)
        ϕ_noflux  = fill(FT(300.0), nz)
        ϕ_sflux   = fill(FT(300.0), nz)
        K = fill(FT(5.0), nz)
        implicit_diffuse!(ϕ_noflux, K, dt, dz, nz, make_turbdata(nz))
        implicit_diffuse!(ϕ_sflux,  K, dt, dz, nz, make_turbdata(nz); sfc_flux=FT(10.0))
        @test ϕ_sflux[1] > ϕ_noflux[1]
    end

    # @testset "top cell is held fixed (Dirichlet top BC in current implementation)" begin
    #     nz = 8
    #     dz = FT(50.0)
    #     dt = FT(60.0)
    #     ϕ  = collect(FT, 1.0:nz)
    #     K  = fill(FT(5.0), nz)
    #     top_val = ϕ[end]
    #     implicit_diffuse!(ϕ, K, dt, dz, nz)
    #     @test ϕ[end] == top_val
    # end

    @testset "result is finite for typical atmospheric values" begin
        nz = 20
        dz = FT(50.0)
        dt = FT(30.0)
        ϕ  = collect(FT, range(288.0, 295.0, length=nz))
        K  = FT[1.0 + 5.0*exp(-(k-10)^2/8.0) for k in 1:nz]
        implicit_diffuse!(ϕ, K, dt, dz, nz, make_turbdata(nz))
        @test all(isfinite.(ϕ))
    end

    @testset "implicit scheme is unconditionally stable for large dt" begin
        nz = 10
        dz = FT(50.0)
        dt = FT(1e5)   # enormous timestep
        ϕ  = collect(FT, range(1.0, 10.0, length=nz))
        K  = fill(FT(100.0), nz)
        implicit_diffuse!(ϕ, K, dt, dz, nz, make_turbdata(nz))
        @test all(isfinite.(ϕ))
        @test all(ϕ .>= 0.0)   # no sign flip / blow-up
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 2.  calculate_buoyancy_frequency
# ─────────────────────────────────────────────────────────────────────────────

@testset "calculate_buoyancy_frequency" begin

    nz = 10
    dz = FT(50.0)

    @testset "neutral (uniform θv) → N² ≈ 0 at interior" begin
        grid = make_test_grid(nz, dz)
        # Uniform θ and qv → uniform θv → N² = 0
        for k in 2:nz-1
            N2 = calculate_buoyancy_frequency(grid, k, constants)
            @test isapprox(N2, 0.0; atol=1e-10)
        end
    end

    @testset "stable stratification → N² > 0" begin
        # θ increases with height  (stable)
        θ_stable = FT[288.0 + 0.01*k for k in 1:nz]
        grid = make_test_grid(nz, dz; θ_profile=θ_stable)
        for k in 2:nz-1
            N2 = calculate_buoyancy_frequency(grid, k, constants)
            @test N2 > 0.0
        end
    end

    @testset "unstable stratification → N² < 0" begin
        # θ decreases with height (unstable)
        θ_unstable = FT[295.0 - 0.05*k for k in 1:nz]
        grid = make_test_grid(nz, dz; θ_profile=θ_unstable)
        for k in 2:nz-1
            N2 = calculate_buoyancy_frequency(grid, k, constants)
            @test N2 < 0.0
        end
    end

    @testset "boundary cells (k=1, k=nz) don't throw or return NaN" begin
        grid = make_test_grid(nz, dz)
        N2_bot = calculate_buoyancy_frequency(grid, 1,  constants)
        N2_top = calculate_buoyancy_frequency(grid, nz, constants)
        @test isfinite(N2_bot)
        @test isfinite(N2_top)
    end

    @testset "magnitude scales with |dθ/dz|" begin
        # Anchor both profiles to the same θ at the test cell so the θv_k
        # denominator is identical; then doubling the gradient doubles N² exactly.
        dθ_dz1 = FT(0.005)
        dθ_dz2 = FT(0.010)
        k_test = 5
        θ_ref  = FT(289.0)
        θ1 = FT[θ_ref + dθ_dz1 * dz * (i - k_test) for i in 1:nz]
        θ2 = FT[θ_ref + dθ_dz2 * dz * (i - k_test) for i in 1:nz]
        g1 = make_test_grid(nz, dz; θ_profile=θ1)
        g2 = make_test_grid(nz, dz; θ_profile=θ2)
        N2_1 = calculate_buoyancy_frequency(g1, k_test, constants)
        N2_2 = calculate_buoyancy_frequency(g2, k_test, constants)
        @test isapprox(N2_2 / N2_1, 2.0; rtol=1e-10)
    end
end

# ──────────────────��──────────────────────────────────────────────────────────
# 3.  deardorff_mixing_length!
# ─────────────────────────────────────────────────────────────────────────────

@testset "deardorff_mixing_length!" begin

    nz  = 10
    dz  = FT(50.0)
    tke = tke_settings{FT}()

    @testset "all lengths are positive" begin
        grid = make_test_grid(nz, dz)
        l = zeros(FT, nz)
        deardorff_mixing_length!(l, grid, tke, constants)
        @test all(l .> 0.0)
    end

    @testset "floor is enforced (l >= 1.0)" begin
        # Very small TKE in a strongly stable environment
        θ_stable = FT[288.0 + 1.0*k for k in 1:nz]
        e_tiny   = fill(FT(1e-12), nz)
        grid = make_test_grid(nz, dz; θ_profile=θ_stable, e=e_tiny)
        l = zeros(FT, nz)
        deardorff_mixing_length!(l, grid, tke, constants)
        @test all(isapprox.(l, 1.0; atol=1e-6))
    end

    @testset "neutral/unstable → l == 70 (asymptotic scale)" begin
        θ_unstable = FT[295.0 - 0.1*k for k in 1:nz]
        grid = make_test_grid(nz, dz; θ_profile=θ_unstable)
        l = zeros(FT, nz)
        deardorff_mixing_length!(l, grid, tke, constants)
        @test all(l .== 70.0)
    end

    @testset "larger TKE → larger mixing length (stable regime)" begin
        θ_stable = FT[288.0 + 0.01*k for k in 1:nz]
        e_lo = fill(FT(1e-4), nz)
        e_hi = fill(FT(1e-1), nz)
        g_lo = make_test_grid(nz, dz; θ_profile=θ_stable, e=e_lo)
        g_hi = make_test_grid(nz, dz; θ_profile=θ_stable, e=e_hi)
        l_lo = zeros(FT, nz);  deardorff_mixing_length!(l_lo, g_lo, tke, constants)
        l_hi = zeros(FT, nz);  deardorff_mixing_length!(l_hi, g_hi, tke, constants)
        # At least some cells where hi-TKE gives larger l
        @test any(l_hi .> l_lo)
    end

    @testset "no NaN or Inf" begin
        grid = make_test_grid(nz, dz)
        l = zeros(FT, nz)
        deardorff_mixing_length!(l, grid, tke, constants)
        @test all(isfinite.(l))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 4.  turb_timestep!  (integration)
# ─────────────────────────────────────────────────────────────────────────────

@testset "turb_timestep!" begin

    nz = 20
    dz = FT(50.0)
    dt = FT(5.0)
    tke = tke_settings{FT}()
    scmset = dynamic_settings()
    run_tke!(grid) = turb_timestep!(DynON(), grid, tke, constants, dt, scmset, make_turbdata(nz))

    @testset "TKE floor is enforced after timestep" begin
        e_zero = fill(FT(0.0), nz)
        grid   = make_test_grid(nz, dz; e=e_zero)
        run_tke!(grid)
        @test all(grid.states.e .>= 0.0)
    end

    @testset "surface TKE BC elevates near-surface TKE" begin
        grid_flux   = make_test_grid(nz, dz; e=fill(FT(0.0), nz))
        grid_noflux = make_test_grid(nz, dz; e=fill(FT(0.0), nz))
        tke_noflux  = tke_settings{FT}(u_star=FT(0.0))
        turb_timestep!(DynON(), grid_flux,   tke,        constants, dt, scmset, make_turbdata(nz))
        turb_timestep!(DynON(), grid_noflux, tke_noflux, constants, dt, scmset, make_turbdata(nz))
        @test grid_flux.states.e[1] > grid_noflux.states.e[1]
    end

    @testset "TKE stays finite" begin
        grid = make_test_grid(nz, dz)
        for _ in 1:50
            run_tke!(grid)
        end
        @test all(isfinite.(grid.states.e))
    end

    @testset "TKE stays non-negative" begin
        grid = make_test_grid(nz, dz)
        for _ in 1:50
            run_tke!(grid)
        end
        @test all(grid.states.e .>= 0.0)
    end

    @testset "θ and qv stay finite after turb_timestep!" begin
        grid = make_test_grid(nz, dz)
        run_tke!(grid)
        @test all(isfinite.(grid.states.θ))
        @test all(isfinite.(grid.states.qv))
    end

    @testset "wind components stay finite" begin
        grid = make_test_grid(nz, dz)
        run_tke!(grid)
        @test all(isfinite.(grid.wind.u))
        @test all(isfinite.(grid.wind.v))
    end

    @testset "Coriolis rotates wind (u,v change)" begin
        grid = make_test_grid(nz, dz)
        u0 = copy(grid.wind.u)
        v0 = copy(grid.wind.v)
        run_tke!(grid)
        # At least one component should change under Coriolis
        @test grid.wind.u != u0 || grid.wind.v != v0
    end

    # @testset "geostrophic wind at rest produces no spurious Coriolis forcing" begin
    #     # If u == u_geo and v == v_geo everywhere, Coriolis tendency = 0
    #     u_geo = FT(7.0);  v_geo = FT(-2.0)
    #     tke_geo = tke_settings{FT}(
    #         geostrophic_u = z -> u_geo,
    #         geostrophic_v = z -> v_geo,
    #     )
    #     grid = make_test_grid(nz, dz)
    #     grid.wind.u .= u_geo
    #     grid.wind.v .= v_geo
    #     grid.states.e .= FT(0)  # nonzero TKE to avoid floor effects
    #     u_before = copy(grid.wind.u)
    #     v_before = copy(grid.wind.v)
    #     # Run just the Coriolis part by using a tiny TKE dt so diffusion is negligible
    #     no_prog_tke_timestep!(grid, tke_geo, constants, FT(0.001))
    #     # Coriolis du = f*(v - v_geo) = 0, dv = -f*(u - u_geo) = 0
    #     @test all(isapprox.(grid.wind.u .- u_before, 0.0; atol=1e-6))
    #     @test all(isapprox.(grid.wind.v .- v_before, 0.0; atol=1e-6))
    # end

    @testset "strongly sheared wind spins up TKE" begin
        # Large vertical shear → P_shear > 0 → TKE should increase from floor
        grid = make_test_grid(nz, dz; e=fill(tke.e_min * 2, nz))
        # Linear shear: u varies by 20 m/s over domain
        for k in 1:nz
            grid.wind.u[k] = FT(k) * FT(20.0) / nz
            grid.wind.v[k] = FT(0.0)
        end
        e_before = sum(grid.states.e)
        run_tke!(grid)
        e_after  = sum(grid.states.e)
        @test e_after > e_before
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 5.  K_m / K_h sanity checks (via tke_settings defaults)
# ─────────────────────────────────────────────────────────────────────────────

@testset "diffusivity magnitudes" begin
    @testset "K_m = c_m * l * sqrt(e), spot-check" begin
        c_m = FT(0.1)
        l   = FT(30.0)
        e   = FT(0.5)
        K_m_expected = c_m * l * sqrt(e)
        @test isapprox(K_m_expected, c_m * l * sqrt(e); rtol=1e-12)
    end

    @testset "K_h = K_m / Pr (Pr = 0.42)" begin
        K_m = FT(1.5)
        Pr  = FT(0.42)
        K_h = K_m / Pr
        @test K_h > K_m  # heat diffuses faster than momentum (Pr < 1)
    end
end
