using Test
using Droplets
using Random
using Distributions

# scm_advance_tests.jl and condensation_bisection_tests.jl/coalescence_spatial_tests.jl
# test each sub-function's (condensation_time_step_spatial!, coalescence_timestep!,
# turb_timestep!, recycle_precipitation!, ...) own DynON/DynOFF dispatch directly

FT = Float64

function make_wiring_env(nz, dz; Ns=200, dynamics, qv_val=FT(0.01), θ_val=FT(290.0),
        P_val=FT(90000.0), area_per_grid=FT(1.0), w_val=FT(0.0), dist=nothing)
    spatial = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz * dz), dt=FT(1.0), t_max=1,
        dt_output=FT(1.0), area_per_grid=area_per_grid)
    coagsettings = coag_settings{FT}(Ns=Ns, ΔV=dz * area_per_grid, n0=FT(1e8))
    condsettings = condensation_settings{FT}(kappa=FT(0.61), Δt=FT(1.0))
    mpdatasettings = mpdata_settings_1d(nz, vertical_boundary_condition=NoFlux())
    tkesettings = tke_settings{FT}()
    diagset = diagnostic_settings{FT}()

    d = something(dist, LogNormal(log(0.05e-6), log(1.5)))
    droplets = init_droplets_scm(d, coagsettings, spatial, fill(qv_val, nz); deterministic_multiplicity=true)

    grid = create_scm_grids(nz, dz, droplets; spatial=spatial, output=nothing)
    grid.states.P .= P_val
    grid.states.P_faces .= P_val
    grid.states.θ .= θ_val
    grid.states.qv .= qv_val
    T_from_theta!(grid.states.T_tmp, grid.states.θ, grid.states.P, grid.states.qv, constants)
    grid.states.ρ .= ρ_ideal_gas.(grid.states.P, grid.states.T_tmp, grid.states.qv, Ref(constants))
    grid.states.ρ_dry .= calc_ρ_dry_from_ρ.(grid.states.ρ, grid.states.qv)
    grid.states.P_dry .= calc_P_dry_from_P.(grid.states.P, grid.states.qv, Ref(constants))
    grid.wind.w .= w_val

    scmsettings = scm_settings{FT}(spatial=spatial, coagsettings=coagsettings,
        condsettings=condsettings, tkesettings=tkesettings, mpdatasettings=mpdatasettings,
        diagnosticsettings=diagset, dynamics=dynamics, Δt=FT(1.0), n_cond=1, n_coag=1,
        spinup_time=FT(0.0), n_rad=1)
    scmdata = scm_data(grid, droplets, scmsettings)

    return grid, droplets, scmdata, scmsettings
end

all_off(; kwargs...) = dynamic_settings(; turbulence=DynOFF(), condensation=DynOFF(),
    REM=DynOFF(), motion=DynOFF(), spinupsaturation=DynOFF(), radiation=DynOFF(),
    coalescence=DynOFF(), settling=DynOFF(), advection=DynOFF(), recycling=DynOFF(),
    top_escape=DynOFF(), thermo_feedback=DynOFF(), turbulent_droplet_diffusion_on=DynOFF(),
    keep_grid_filled=DynOFF(), density_feedback=DynOFF(), kwargs...)

step!(grid, droplets, scmdata, scmsettings; w=z -> FT(0.0)) =
    single_column_timestep(grid, scmsettings.Δt, droplets, scmdata.coagdata, scmdata.conddata,
        scmdata.raddata, scmdata.turbdata, w, scmdata.mpdatatmp, constants, scmsettings, 1)

@testset "single_column_timestep dynamics wiring" begin

    @testset "condensation: dynamics.condensation gates droplet growth" begin
        qv_super = FT(1.01) * esat(FT(285.0)) * FT(0.622) / FT(90000.0)
        for (label, dyn_cond, should_change) in
                [("off", DynOFF(), false), ("on", DynON(), true)]
            grid, droplets, scmdata, scmsettings = make_wiring_env(3, FT(50.0);
                dynamics=all_off(condensation=dyn_cond), qv_val=qv_super, θ_val=theta_from_T(FT(285.0), FT(90000.0), qv_super, Constants{FT}()))
            X_before = copy(droplets.X)
            step!(grid, droplets, scmdata, scmsettings)
            changed = !all(droplets.X .== X_before)
            @test changed == should_change
        end
    end

    @testset "coalescence: dynamics.coalescence gates multiplicity loss" begin
        Random.seed!(11)
        for (dyn_coag, should_change) in [(DynOFF(), false), (DynON(), true)]
            grid, droplets, scmdata, scmsettings = make_wiring_env(1, FT(50.0); Ns=300,
                dynamics=all_off(coalescence=dyn_coag), area_per_grid=FT(0.01))
            # init_droplets_scm seeds dry-aerosol-sized volumes (~0.05µm radius) with a
            # negligible golovin kernel; bump to cloud-droplet size so collisions are likely
            droplets.X .= radius_to_volume.(rand(FT, length(droplets.X)) .* 10e-6 .+ 10e-6)
            ξ_before = sum(droplets.ξ)
            for _ in 1:20
                step!(grid, droplets, scmdata, scmsettings)
            end
            changed = sum(droplets.ξ) != ξ_before
            @test changed == should_change
        end
    end

    @testset "turbulence: dynamics.turbulence gates TKE evolution" begin
        for (dyn_turb, should_change) in [(DynOFF(), false), (DynON(), true)]
            grid, droplets, scmdata, scmsettings = make_wiring_env(5, FT(50.0);
                dynamics=all_off(turbulence=dyn_turb))
            grid.states.e .= 0.0
            e_before = copy(grid.states.e)
            step!(grid, droplets, scmdata, scmsettings)
            changed = !all(grid.states.e .== e_before)
            @test changed == should_change
        end
    end

    @testset "advection: dynamics.advection gates MPDATA transport" begin
        for (dyn_adv, should_change) in [(DynOFF(), false), (DynON(), true)]
            grid, droplets, scmdata, scmsettings = make_wiring_env(6, FT(50.0);
                dynamics=all_off(advection=dyn_adv, thermo_feedback=dyn_adv), w_val=FT(0.3))
            grid.states.qv .= range(FT(0.005), FT(0.015), length=6)  # gradient for MPDATA to act on
            qv_before = copy(grid.states.qv)
            step!(grid, droplets, scmdata, scmsettings; w=z -> FT(0.3))
            changed = !all(grid.states.qv .== qv_before)
            @test changed == should_change
        end
    end

    @testset "motion: dynamics.motion gates droplet vertical displacement" begin
        for (dyn_mot, should_change) in [(DynOFF(), false), (DynON(), true)]
            grid, droplets, scmdata, scmsettings = make_wiring_env(5, FT(50.0);
                dynamics=all_off(motion=dyn_mot, advection=DynOFF(), settling=DynON()))
            z_before = copy(droplets.z_loc)
            step!(grid, droplets, scmdata, scmsettings)
            changed = !all(droplets.z_loc .== z_before)
            @test changed == should_change
        end
    end

    @testset "density_feedback: dynamics.density_feedback picks update_density!'s branch" begin
        cst = Constants{FT}()
        for (dyn_df,) in [(DynON(),), (DynOFF(),)]
            grid, droplets, scmdata, scmsettings = make_wiring_env(4, FT(50.0);
                dynamics=all_off(density_feedback=dyn_df))
            step!(grid, droplets, scmdata, scmsettings)
            if dyn_df isa DynON
                @test grid.states.ρ ≈ ρ_calc_θ.(grid.states.P, grid.states.θ, grid.states.qv, cst)
            else
                @test grid.states.ρ ≈ calc_ρ_from_ρ_dry.(grid.states.ρ_dry, grid.states.qv)
                @test grid.states.P ≈ calc_P_from_P_dry.(grid.states.P_dry, grid.states.qv, Ref(cst))
            end
        end
    end

    @testset "turbulence: qt_tmp/θl_tmp are re-seeded from qv/ql before diffusion (advance.jl:69)" begin
        # single_column_timestep re-derives qt_tmp = qv + ql_tmp (and θl_tmp) from the
        # current state just before turb_timestep! diffuses them; if that seeding were
        # ever dropped or moved after the diffusion call, the column would diffuse a
        # stale/uninitialized qt_tmp instead of the real state. With diffusion otherwise
        # conservative (no-flux top, flux-form interior), total column qt should change
        # by exactly the surface LHF input and nothing else.
        nz, dz = 5, 50.0
        grid, droplets, scmdata, scmsettings = make_wiring_env(nz, dz;
            dynamics=all_off(turbulence=DynON()))

        Droplets.compute_ql_all_cells!(grid.states, constants)
        qt_before_total = sum(grid.states.qv .+ grid.states.ql_tmp)
        ρ1 = grid.states.ρ[1]

        step!(grid, droplets, scmdata, scmsettings)

        Droplets.compute_ql_all_cells!(grid.states, constants)
        qt_after_total = sum(grid.states.qv .+ grid.states.ql_tmp)

        tke = scmsettings.tkesettings
        expected_input = scmsettings.Δt * tke.LHF / (constants.L * ρ1) / dz
        @test (qt_after_total - qt_before_total) ≈ expected_input rtol=1e-6
    end

    @testset "recycling: dynamics.recycling gates surface_precipitation bookkeeping" begin
        for (dyn_rec, should_change) in [(DynOFF(), false), (DynON(), true)]
            grid, droplets, scmdata, scmsettings = make_wiring_env(10, FT(30.0); Ns=50,
                dynamics=all_off(recycling=dyn_rec))
            # force a raindrop below the domain floor
            droplets.z_loc[droplets.I[1]] = FT(-0.5)
            droplets.X[droplets.I[1]] = radius_to_volume(FT(50e-6))
            droplets.ξ[droplets.I[1]] = 10
            Droplets.fill_grid_ranges!(droplets)
            step!(grid, droplets, scmdata, scmsettings)
            changed = any(!=(0.0), grid.output.surface_precipitation)
            @test changed == should_change
        end
    end
end
