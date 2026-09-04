using Test
using Droplets
using Distributions

# Wiring smoke test for `run_scm!`. Builds a grid/droplets directly (via
# `init_droplets_scm`/`create_scm_grids`) instead of `initialize_scm_environment`,
# since the latter's `hydrostatic_pysdm=true` path needs OrdinaryDiffEq, which isn't
# a package dependency. This exercises the same `run_scm!` API surface the Examples
# driver scripts use (default `prescribed_w` threading, `step_forcing!` overrides),
# to catch API drift between `run_scm!`'s signature and the scripts that call it.

function make_smoke_env(FT, nz, dz, prescribed_w)
    Ns = 200
    spatialsettings = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz*dz), dt=FT(1.0), t_max=5,
        dt_output=FT(5.0), area_per_grid=FT(1.0))
    coagsettings = coag_settings{FT}(Ns=Ns, ΔV=dz*spatialsettings.area_per_grid, n0=FT(1e8))
    condensationsettings = condensation_settings{FT}(kappa=FT(0.61), Δt=FT(1.0))
    # NoFlux, not the mpdata_settings_1d default of Periodic(): Periodic hits a
    # documented pre-existing bug (test/mpdata_tests.jl:247, antiosc!'s GCz_tmp[nz+1]
    # is never written) that isn't what this test is targeting.
    mpdatasettings = mpdata_settings_1d(nz, vertical_boundary_condition=NoFlux())
    tkesettings = tke_settings{FT}()
    diagnosticsettings = diagnostic_settings{FT}()

    dist = LogNormal(log(0.05e-6), log(1.5))
    qv_profile = fill(FT(0.01), nz)
    droplets = init_droplets_scm(dist, coagsettings, spatialsettings, qv_profile;
        deterministic_multiplicity=true)

    grid = create_scm_grids(nz, dz, droplets; spatial=spatialsettings, output=nothing)
    grid.states.P .= FT(90000.0)
    grid.states.P_faces .= FT(90000.0)
    grid.states.θ .= FT(290.0)
    grid.states.qv .= FT(0.01)
    T_from_theta!(grid.states.T_tmp, grid.states.θ, grid.states.P, grid.states.qv, constants)
    grid.states.ρ .= ρ_ideal_gas.(grid.states.P, grid.states.T_tmp, grid.states.qv, Ref(constants))
    # mpdata_scm! mass-weights by ρ_dry (its g_factor) regardless of density_feedback;
    # left at create_scm_grids' zero default this divides by zero
    grid.states.ρ_dry .= calc_ρ_dry_from_ρ.(grid.states.ρ, grid.states.qv)
    # matches how initialize_scm_environment seeds grid.wind.w from prescribed_w
    grid.wind.w .= prescribed_w.(grid.faces_z)

    # condensation/turbulence stay off: this test targets run_scm!'s wiring
    # (prescribed_w/step_forcing! threading), not thermodynamic correctness
    dyn = dynamic_settings(turbulence=DynOFF(), motion=DynON(), advection=DynON(),
        radiation=DynOFF(), condensation=DynOFF(), turbulent_droplet_diffusion_on=DynOFF(),
        keep_grid_filled=DynOFF(), recycling=DynOFF(), settling=DynOFF(),
        spinupsaturation=DynOFF(), coalescence=DynOFF(), top_escape=DynOFF(),
        thermo_feedback=DynOFF(), density_feedback=DynON())

    scmsettings = scm_settings{FT}(spatial=spatialsettings, coagsettings=coagsettings,
        condsettings=condensationsettings, tkesettings=tkesettings,
        mpdatasettings=mpdatasettings, diagnosticsettings=diagnosticsettings, dynamics=dyn,
        Δt=FT(1.0), n_cond=1, n_coag=1, spinup_time=FT(0.0), n_rad=1)

    scmdata = scm_data(grid, droplets, scmsettings)

    return grid, droplets, scmdata, scmsettings
end

@testset "run_scm!: default prescribed_w is actually used" begin
    FT = Float64
    prescribed_w(z) = FT(-0.5)  # constant downward motion, no step_forcing! override
    grid, droplets, scmdata, scmsettings = make_smoke_env(FT, 10, FT(50.0), prescribed_w)

    z_loc_before = copy(droplets.z_loc)

    run_scm!(grid, droplets, scmdata, constants, scmsettings, prescribed_w)

    @test all(isfinite, grid.states.qv)
    @test all(isfinite, droplets.X)
    # every droplet should have moved down by exactly w*t_max = -0.5*5 = -2.5 m,
    # except any that started close enough to the domain floor to get clamped by
    # update_droplet_positions!'s z_loc >= -1 floor (init_droplets_scm's z_loc is
    # randomly placed within its cell regardless of deterministic_multiplicity, so
    # which droplets those are isn't reproducible across platforms/RNG streams)
    unclamped = z_loc_before .> 2.5
    @test any(unclamped)
    @test all(droplets.z_loc[unclamped] .- z_loc_before[unclamped] .≈ -2.5)
end

@testset "run_scm!: step_forcing! overrides prescribed_w per step" begin
    FT = Float64
    # small nonzero so grid.wind.w isn't exactly zero (keeps mpdata well-behaved);
    # step_forcing!'s much larger override is what should actually move the droplets
    prescribed_w(z) = FT(-0.001)
    grid, droplets, scmdata, scmsettings = make_smoke_env(FT, 10, FT(50.0), prescribed_w)

    z_loc_before = copy(droplets.z_loc)

    run_scm!(grid, droplets, scmdata, constants, scmsettings, prescribed_w;
        step_forcing! = (grid, droplets, i, dt) -> (z -> FT(-0.5)))

    @test all(isfinite, droplets.X)
    # step_forcing! always returns a nonzero override, so prescribed_w's zero must
    # never be used: droplets should move exactly as in the previous test (modulo
    # the domain-floor clamp — see the comment in the test above)
    unclamped = z_loc_before .> 2.5
    @test any(unclamped)
    @test all(droplets.z_loc[unclamped] .- z_loc_before[unclamped] .≈ -2.5)
end
