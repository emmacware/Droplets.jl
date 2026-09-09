using Test
using Droplets

#  check bisection algorithm itself directly against the implicit equation it is solving

FT = Float64

function make_one_droplet(FT, R, dry_r)
    droplet_attributes_1d{FT}(
        [100], [radius_to_volume(R)], [dry_r^3], [FT(0.0)], [1], [FT(0.0)], [1:1], [1],
    )
end

""" Computed independently here (not by calling into the implementation) so tests
can deliberately place Senv on either side of it."""
kohler_Scrit(T, kappa, dry_r3) = 1 + sqrt(4 * Droplets.akk(T)^3 / (27 * kappa * dry_r3))

"""Ground truth for the implicit growth equation
    r^2 = R_old^2 + Δt*c*(Senv - 1 - A/r + b/r^3)
independant"""
function reference_root(R_old, kappa, T, Senv, c, dry_r3, Δt)
    A = Droplets.akk(T)
    b = kappa * dry_r3
    S1 = Senv - 1
    g(r) = r^2 - R_old^2 - Δt * c * (S1 - A / r + b / r^3)
    lo, hi = cbrt(dry_r3), max(R_old, cbrt(dry_r3)) * FT(3.0)
    @assert g(lo) < 0 && g(hi) > 0 "reference bracket doesn't contain a root -- test setup is off"
    for _ in 1:200
        mid = (lo + hi) / 2
        g(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

@testset "dXkappakohler_bisection (DynOFF): matches an independent reference solve" begin
    kappa = FT(0.61)
    dry_r = FT(0.05e-6)
    dry_r3 = dry_r^3
    T = FT(290.0)
    c = FT(2.0) / (FK(T, Constants{FT}()) + FD(T, Constants{FT}()))
    cst = Constants{FT}()
    Scrit = kohler_Scrit(T, kappa, dry_r3)
    @test Scrit > 1  

    @testset "supersaturated growth (below Scrit, ordinary bisection branch)" begin
        R0 = FT(5e-6)
        # Senv must exceed not just 1 but this droplet's own Köhler-equilibrium S at its
        # current radius (S_eq ≈ 1 + A/R0 ≈ 1.00023 here) to actually grow
        Senv = FT(1.0008)
        @assert Senv < Scrit  # must land in the ordinary branch, not the closed-form one
        droplets = make_one_droplet(FT, R0, dry_r)
        dXkappakohler_bisection(DynOFF(), droplets, 1, kappa, T, Senv, c, FT(0.0), cst, nothing, FT(1.0), 20)
        r_new = volume_to_radius(droplets.X[1])
        @test r_new > R0  # grew
        r_ref = reference_root(R0, kappa, T, Senv, c, dry_r3, FT(1.0))
        @test r_new ≈ r_ref rtol=1e-3  # bisection's own break criterion is ~1e-4 relative in r^2
    end

    @testset "subsaturated evaporation (not fully to dry radius)" begin
        R0 = FT(8e-6)
        Senv = FT(0.99)
        @assert Senv < Scrit
        droplets = make_one_droplet(FT, R0, dry_r)
        dXkappakohler_bisection(DynOFF(), droplets, 1, kappa, T, Senv, c, FT(0.0), cst, nothing, FT(1.0), 20)
        r_new = volume_to_radius(droplets.X[1])
        @test r_new < R0
        @test r_new > dry_r
        r_ref = reference_root(R0, kappa, T, Senv, c, dry_r3, FT(1.0))
        @test r_new ≈ r_ref rtol=1e-3
    end

    @testset "too few iterations leaves a detectably worse (but still bracketed) answer" begin
        R0 = FT(3e-6)
        Senv = FT(1.0002)
        @assert Senv < Scrit
        r_ref = reference_root(R0, kappa, T, Senv, c, dry_r3, FT(2.0))

        d_1iter = make_one_droplet(FT, R0, dry_r)
        dXkappakohler_bisection(DynOFF(), d_1iter, 1, kappa, T, Senv, c, FT(0.0), cst, nothing, FT(2.0), 1)
        d_default = make_one_droplet(FT, R0, dry_r)
        dXkappakohler_bisection(DynOFF(), d_default, 1, kappa, T, Senv, c, FT(0.0), cst, nothing, FT(2.0), 20)

        err_1iter = abs(volume_to_radius(d_1iter.X[1]) - r_ref)
        err_default = abs(volume_to_radius(d_default.X[1]) - r_ref)
        @test err_default < err_1iter
        @test err_default / r_ref < 1e-3
    end

    @testset "cross-check against find_equilibrium_radius: bisection near equilibrium barely moves the droplet" begin
        # find_equilibrium_radius solves the *steady-state* Köhler equation
        # (Sm1 = A/r - b/r^3) via an independent bisection
        S_env = FT(0.98)
        droplets = make_one_droplet(FT, FT(1e-7), dry_r)
        find_equilibrium_radius(droplets, 1, kappa, T, S_env, cst)
        R_eq = volume_to_radius(droplets.X[1])

        dXkappakohler_bisection(DynOFF(), droplets, 1, kappa, T, S_env, c, FT(0.0), cst, nothing, FT(1.0), 20)
        R_after = volume_to_radius(droplets.X[1])
        @test R_after ≈ R_eq rtol=1e-4
    end

    @testset "Senv > Scrit: unconditionally-activated closed-form branch matches its own formula" begin
        # once haze droplets are past the critical (max) Köhler supersaturation there
        # is no equilibrium root; verify
        dry_r_big = FT(0.2e-6)
        A = Droplets.akk(T)
        b = kappa * dry_r_big^3
        Scrit = 1 + sqrt(4 * A^3 / (27 * b))
        Senv = Scrit + FT(0.05)  # past the hump

        R0 = FT(0.3e-6)
        droplets = make_one_droplet(FT, R0, dry_r_big)
        Δt = FT(1.0)
        dXkappakohler_bisection(DynOFF(), droplets, 1, kappa, T, Senv, c, FT(0.0), cst, nothing, Δt, 20)

        r_new = volume_to_radius(droplets.X[1])
        expected = sqrt(R0^2 + 2 * c * (Senv - 1) * Δt)
        @test r_new ≈ expected rtol=1e-10
    end

    @testset "strongly subsaturated: droplet clamps at the dry radius, never below it" begin
        droplets = make_one_droplet(FT, FT(0.5e-6), dry_r)
        dXkappakohler_bisection(DynOFF(), droplets, 1, kappa, T, FT(0.5), c, FT(0.0), cst, nothing, FT(50.0), 20)
        r_new = volume_to_radius(droplets.X[1])
        @test r_new >= dry_r * (1 - 1e-9)
    end

    @testset "growth rate is monotone increasing in Senv" begin
        R0 = FT(5e-6)
        Xs = map(FT[0.98, 1.0, 1.005, 1.01]) do S
            d = make_one_droplet(FT, R0, dry_r)
            dXkappakohler_bisection(DynOFF(), d, 1, kappa, T, S, c, FT(0.0), cst, nothing, FT(1.0), 20)
            d.X[1]
        end
        @test issorted(Xs)
    end
end

@testset "limitsat" begin
    @test limitsat(DynON(), 1.5) == 1.01
    @test limitsat(DynON(), 1.005) == 1.005
    @test limitsat(DynON(), 0.9) == 0.9
    @test limitsat(DynOFF(), 1.5) == 1.5
end

@testset "apply_thermo_feedback!" begin
    cst = Constants{FT}()
    P = FT(95000.0)
    qv = FT(8e-3)
    T0 = FT(285.0)
    state = (T_tmp=[T0], θ=[theta_from_T(T0, P, qv, cst)], P=[P], qv=[qv], ρ=[ρ_ideal_gas(P, T0, qv, cst)])

    dqv = FT(-1e-4)  # condensation: qv decreases, so T/θ/ρ get updated with warming
    apply_thermo_feedback!(DynON(), state, 1, dqv, cst)
    @test state.T_tmp[1] ≈ T0 - dqv * cst.L / cst.Cp_air
    @test state.θ[1] ≈ theta_from_T(state.T_tmp[1], P, qv, cst)
    @test state.ρ[1] ≈ ρ_ideal_gas(P, state.T_tmp[1], qv, cst)

    state2 = (T_tmp=[T0], θ=[theta_from_T(T0, P, qv, cst)], P=[P], qv=[qv], ρ=[ρ_ideal_gas(P, T0, qv, cst)])
    apply_thermo_feedback!(DynOFF(), state2, 1, dqv, cst)
    @test state2.T_tmp[1] == T0  # strict
end

@testset "condensation_time_step_spatial!(DynON, ...): qt and θl are conserved across a step" begin
    # both qt and θl should be invariant through condensation
    nz, nsd = 1, 4
    cst = Constants{FT}()
    condsettings = condensation_settings{FT}(kappa=FT(0.61))
    spatial = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(100.0))
    dynamics = dynamic_settings()  # thermo_feedback = DynON() by default

    Rs = FT[3e-6, 5e-6, 7e-6, 9e-6]
    dry_r = FT(0.05e-6)
    droplets = droplet_attributes_1d{FT}(
        fill(100, nsd), radius_to_volume.(Rs), fill(dry_r^3, nsd),
        zeros(FT, nsd), ones(Int, nsd), zeros(FT, nsd), [1:nsd], collect(1:nsd),
    )

    T0 = FT(285.0); P0 = FT(97000.0); qv0 = FT(0.012)
    state = (
        T_tmp=[T0], θ=[theta_from_T(T0, P0, qv0, cst)], P=[P0], qv=[qv0],
        ρ=[ρ_ideal_gas(P0, T0, qv0, cst)], droplets=droplets, spatial=spatial, ql_tmp=[FT(0.0)],
    )

    compute_ql_at_cell!(state, 1, cst)
    qt_before = state.qv[1] + state.ql_tmp[1]
    θl_before = θl(state.P[1], state.T_tmp[1], state.ql_tmp[1], state.qv[1], cst)
    qv_before = state.qv[1]

    conddata = condensation_data(FT, nz, nsd)
    condensation_time_step_spatial!(DynON(), droplets, state, nz, FT(1.0), conddata, cst,
        condsettings, spatial, nothing, dynamics, FT(0.0))

    compute_ql_at_cell!(state, 1, cst)
    qt_after = state.qv[1] + state.ql_tmp[1]
    θl_after = θl(state.P[1], state.T_tmp[1], state.ql_tmp[1], state.qv[1], cst)

    @test state.qv[1] != qv_before  # sanity: the step actually did something
    @test qt_after ≈ qt_before rtol=1e-10
    @test θl_after ≈ θl_before rtol=1e-10
end

@testset "condensation_time_step_spatial!(DynOFF, ...): strict no-op" begin
    nz, nsd = 1, 2
    cst = Constants{FT}()
    condsettings = condensation_settings{FT}(kappa=FT(0.61))
    spatial = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(100.0))
    dynamics = dynamic_settings()
    droplets = droplet_attributes_1d{FT}(
        fill(100, nsd), fill(radius_to_volume(FT(5e-6)), nsd), fill(FT(0.05e-6)^3, nsd),
        zeros(FT, nsd), ones(Int, nsd), zeros(FT, nsd), [1:nsd], collect(1:nsd),
    )
    T0 = FT(285.0); P0 = FT(97000.0); qv0 = FT(0.01)
    state = (
        T_tmp=[T0], θ=[theta_from_T(T0, P0, qv0, cst)], P=[P0], qv=[qv0],
        ρ=[ρ_ideal_gas(P0, T0, qv0, cst)], droplets=droplets, spatial=spatial,
    )
    conddata = condensation_data(FT, nz, nsd)
    X_before = copy(droplets.X)
    qv_before = copy(state.qv)

    condensation_time_step_spatial!(DynOFF(), droplets, state, nz, FT(1.0), conddata, cst,
        condsettings, spatial, nothing, dynamics, FT(0.0))

    @test droplets.X == X_before
    @test state.qv == qv_before
end
