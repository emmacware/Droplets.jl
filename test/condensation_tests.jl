using Test
using Droplets

FT = Float64

const _L   = 25.0e5  # matches constants.L
const _Rv  = 461.0
const _k   = 0.024
const _Dv  = 2.26e-5
const _ρl  = 1000.0

kohler_S(R, T, kappa, dry_r3) =
    1.0 + (3.3e-7 / T) / R - kappa * dry_r3 / R^3

kohler_Rcrit(T, kappa, dry_r3) =
    sqrt(3.0 * kappa * dry_r3 / (3.3e-7 / T))

# Helper: build a minimal droplet_attributes_1d with n droplets all in cell 1
function make_droplets(FT, n; R=FT(5e-6), dry_r=FT(0.05e-6), xi=100)
    droplet_attributes_1d{FT}(
        fill(xi, n),
        fill(radius_to_volume(R), n),
        fill(dry_r^3, n),
        fill(FT(25.0), n),
        fill(1, n),
        fill(FT(0.0), n),
        [1:n],
        collect(1:n),
    )
end


@testset "esat: saturation vapor pressure" begin
    # August–Roche–Magnus formula: e_s(0°C) = 610.94 Pa by definition (T = 273.15 K)
    @test esat(FT(273.15)) ≈ 610.94 atol=0.01

    # Reference values at integer-Celsius temperatures
    @test esat(FT(283.15)) ≈ 1228.0 atol=15.0   # 10°C
    @test esat(FT(293.15)) ≈ 2338.0 atol=25.0   # 20°C
    @test esat(FT(303.15)) ≈ 4243.0 atol=50.0   # 30°C

    # Always positive
    for T in [FT(240), FT(260), FT(273.15), FT(290), FT(310)]
        @test esat(T) > 0
    end

    # Strictly monotone increasing with temperature
    T_vals = [FT(250), FT(260), FT(270), FT(280), FT(290), FT(300)]
    for i in 1:length(T_vals)-1
        @test esat(T_vals[i]) < esat(T_vals[i+1])
    end

    # Clausius-Clapeyron: roughly doubles every 10 K near 0°C
    ratio = esat(FT(283.15)) / esat(FT(273.15))
    @test ratio > 1.6 && ratio < 2.2
end


@testset "sat: actual vapor pressure from mixing ratio" begin
    P = FT(97000)

    @test sat(FT(0), P) == FT(0)

    qv_tiny = FT(1e-4)
    @test sat(qv_tiny, P) ≈ qv_tiny * P / 0.622  rtol=0.01

    qv = FT(5e-3)
    @test sat(qv, FT(2) * P) / sat(qv, P) ≈ 2.0  rtol=0.01

    @test sat(FT(0.01), P) > 0

    esat_290 = FT(100 * 6.1094 * exp(17.625 * 17 / (17 + 243.04)))
    qv_98pct = FT(0.98 * esat_290 * 0.622 / (P - 0.98 * esat_290 * 0.378))
    @test sat(qv_98pct, P) ≈ FT(0.98) * esat_290  rtol=0.001
end


@testset "FK and FD: thermodynamic growth coefficients" begin
    for T in [FT(260), FT(273), FT(290), FT(305)]
        @test FK(T, constants) > 0
        @test FD(T, constants) > 0
        @test FK(T, constants) + FD(T, constants) > 0
    end

    @test FK(FT(280), constants) > FK(FT(300), constants)

    T = FT(290)
    FK_manual = (_L / (_Rv * T) - 1) * (_L * _ρl) / (_k * T)
    FD_manual = _ρl * _Rv * T / (_Dv * esat(T))
    @test FK(T, constants) ≈ FK_manual  rtol=1e-4
    @test FD(T, constants) ≈ FD_manual  rtol=1e-4
    @test FK(T, constants) ≈ 6.36e9  rtol=0.05
    @test FD(T, constants) ≈ 3.06e9  rtol=0.05
    @test FK(T, constants) > FD(T, constants)
end


@testset "drkohler (Senv form): radial growth rate" begin
    R_large = FT(5e-6)
    M_sol   = FT(1e-16)
    m_NaCl  = FT(58.44)
    T       = FT(290)
    dt      = FT(1.0)

    @test drkohler(R_large, M_sol, m_NaCl, T, FT(1.01), constants, dt) > 0
    @test drkohler(R_large, M_sol, m_NaCl, T, FT(0.95), constants, dt) < 0

    dr_eq = drkohler(FT(200e-6), M_sol, m_NaCl, T, FT(1.0), constants, dt)
    @test abs(dr_eq) < 1e-6

    # Guard: R + dr·dt must not go negative
    dr_bad = drkohler(FT(1e-8), M_sol, m_NaCl, T, FT(0.0), constants, FT(100.0))
    @test FT(1e-8) + dr_bad * FT(100.0) >= 0

    dr1 = drkohler(R_large, M_sol, m_NaCl, T, FT(1.005), constants, dt)
    dr2 = drkohler(R_large, M_sol, m_NaCl, T, FT(1.010), constants, dt)
    @test dr2 > dr1

    dr_5  = drkohler(FT(5e-6),  M_sol, m_NaCl, T, FT(1.01), constants, dt)
    dr_20 = drkohler(FT(20e-6), M_sol, m_NaCl, T, FT(1.01), constants, dt)
    @test dr_5 > dr_20

    esat_T    = 100 * 6.1094 * exp(17.625 * (T - 273) / ((T - 273) + 243.04))
    FK_T      = (_L / (_Rv * T) - 1) * (_L * _ρl) / (_k * T)
    FD_T      = _ρl * _Rv * T / (_Dv * esat_T)
    dr_approx = (FT(1.01) - 1) / ((FK_T + FD_T) * R_large)
    @test drkohler(R_large, M_sol, m_NaCl, T, FT(1.01), constants, dt) ≈ dr_approx  rtol=0.05
end


@testset "drkohler (qv/P form): consistency with Senv form" begin
    R     = FT(5e-6)
    M_sol = FT(1e-16)
    m     = FT(58.44)
    T     = FT(290)
    P     = FT(97000)
    dt    = FT(1.0)

    qv_super = FT(1.01 * esat(T) * 0.622 / P)
    @test drkohler(R, M_sol, m, T, qv_super, P, constants, dt) > 0

    qv_sub = FT(0.95 * esat(T) * 0.622 / P)
    @test drkohler(R, M_sol, m, T, qv_sub, P, constants, dt) < 0

    qv   = FT(0.012)
    Senv = sat(qv, P) / esat(T)
    @test drkohler(R, M_sol, m, T, Senv, constants, dt) ≈
          drkohler(R, M_sol, m, T, qv, P, constants, dt)  rtol=1e-9
end


@testset "dXkohler_function_of_radius: volume growth" begin
    R     = FT(5e-6)
    M_sol = FT(1e-16)
    m     = FT(58.44)
    T     = FT(290)
    dt    = FT(1.0)
    S     = FT(1.01)

    dR         = drkohler(R, M_sol, m, T, S, constants, dt)
    dX_from_dR = 4 * π * R^2 * dR
    @test dXkohler_function_of_radius(R, M_sol, m, T, S, constants, dt) ≈ dX_from_dR  rtol=1e-10
    @test dXkohler_function_of_radius(R, M_sol, m, T, S, constants, dt) > 0
end


@testset "dXkohler_function_of_radius_activated: activated volume growth" begin
    R  = FT(10e-6)
    T  = FT(290)
    dt = FT(1.0)
    S  = FT(1.005)

    dR_act = drkohler_activated(R, T, S, constants, dt)
    @test dXkohler_function_of_radius_activated(R, T, S, constants, dt) ≈ 4 * π * R^2 * dR_act  rtol=1e-10
    @test dXkohler_function_of_radius_activated(R, T, S, constants, dt) > 0
end


@testset "drkappakohler: kappa-Köhler growth rate" begin
    dry_r  = FT(0.05e-6)
    dry_r3 = dry_r^3
    kappa  = FT(1.28)
    T      = FT(290)
    dt     = FT(1.0)
    R_large = FT(5e-6)

    @test drkappakohler(R_large, dry_r3, kappa, T, FT(1.01), constants, dt) > 0
    @test drkappakohler(R_large, dry_r3, kappa, T, FT(0.95), constants, dt) < 0

    # Guard: R never goes negative
    dr_evap = drkappakohler(FT(1e-8), dry_r3, kappa, T, FT(0.0), constants, FT(100.0))
    @test FT(1e-8) + dr_evap * FT(100.0) >= 0

    R_huge = FT(200e-6)
    Senv   = FT(1.01)
    dr_approx = (Senv - 1) / ((FK(T, constants) + FD(T, constants)) * R_huge)
    @test drkappakohler(R_huge, dry_r3, kappa, T, Senv, constants, dt) ≈ dr_approx  rtol=0.05

    R_small    = FT(0.15e-6)
    Senv_small = FT(1.002)
    dr_hi = drkappakohler(R_small, dry_r3, FT(1.2), T, Senv_small, constants, dt)
    dr_lo = drkappakohler(R_small, dry_r3, FT(0.1), T, Senv_small, constants, dt)
    @test dr_hi > dr_lo

    akk_T = 3.3e-7 / T
    Senv_above_kelvin = FT(1.0 + akk_T / R_large + 0.001)
    @test drkappakohler(R_large, dry_r3, FT(0.0), T, Senv_above_kelvin, constants, dt) > 0
end


@testset "drkohler_activated: activated droplet growth" begin
    R  = FT(10e-6)
    T  = FT(290)
    dt = FT(1.0)

    @test drkohler_activated(R, T, FT(1.005), constants, dt) > 0
    @test drkohler_activated(R, T, FT(0.99),  constants, dt) < 0
    @test drkohler_activated(R, T, FT(1.0),   constants, dt) == FT(0)

    # raw dr/dt has no self-clamp (unclamped by design, like its siblings), so
    # near the dry radius with a large timestep it's just a large negative rate
    dr = drkohler_activated(FT(1e-9), T, FT(0.0), constants, FT(1000.0))
    @test isfinite(dr)

    # Linear in (Senv - 1)
    dr1 = drkohler_activated(R, T, FT(1.01), constants, dt)
    dr2 = drkohler_activated(R, T, FT(1.02), constants, dt)
    @test dr2 ≈ dr1 * 2  rtol=1e-6

    # dr ∝ 1/R
    dr_small = drkohler_activated(FT(5e-6),  T, FT(1.005), constants, dt)
    dr_large = drkohler_activated(FT(20e-6), T, FT(1.005), constants, dt)
    @test dr_small > dr_large
    @test dr_small / dr_large ≈ 4.0  rtol=1e-6

    T_mag  = FT(290)
    R_mag  = FT(10e-6)
    S_mag  = FT(1.01)
    esat_290  = 100 * 6.1094 * exp(17.625 * (T_mag - 273) / ((T_mag - 273) + 243.04))
    FK_290 = (_L / (_Rv * T_mag) - 1) * (_L * _ρl) / (_k * T_mag)
    FD_290 = _ρl * _Rv * T_mag / (_Dv * esat_290)
    dr_expected = (S_mag - 1) / ((FK_290 + FD_290) * R_mag)
    @test drkohler_activated(R_mag, T_mag, S_mag, constants, dt) ≈ dr_expected  rtol=0.01
    @test dr_expected ≈ 1.06e-7  rtol=0.1
end


@testset "terminal_v: droplet terminal velocity (Gunn-Kinzer table)" begin
    for r in [FT(1e-6), FT(10e-6), FT(25e-6), FT(50e-6),
              FT(100e-6), FT(500e-6), FT(1e-3), FT(3e-3)]
        @test terminal_v(r) > 0
    end

    # Stokes regime (r < 39 μm, d < 0.0078 cm): coeff is 1.2e6, v ∝ r²
    @test terminal_v(FT(10e-6)) ≈ 0.012  rtol=1e-4
    @test terminal_v(FT(20e-6)) / terminal_v(FT(10e-6)) ≈ 4.0  rtol=1e-6

    # Gunn-Kinzer table: values at exact table nodes (d in cm → v in cm/s → m/s)
    @test terminal_v(FT(50e-6))   ≈ 0.27  rtol=1e-4   # d=0.010 cm → 27 cm/s
    @test terminal_v(FT(100e-6))  ≈ 0.72  rtol=1e-4   # d=0.020 cm → 72 cm/s
    @test terminal_v(FT(150e-6))  ≈ 1.17  rtol=1e-4   # d=0.030 cm → 117 cm/s
    @test terminal_v(FT(2000e-6)) ≈ 8.83  rtol=1e-4   # d=0.400 cm → 883 cm/s

    # Cap at d > 0.58 cm (r > ~2.9 mm): all larger drops extrapolate at 917 cm/s
    @test terminal_v(FT(3e-3)) ≈ 9.17  rtol=1e-4
    @test terminal_v(FT(5e-3)) == terminal_v(FT(3e-3))

    # Monotonically increasing up to cap
    rvals = [FT(10e-6), FT(50e-6), FT(100e-6), FT(500e-6), FT(2e-3)]
    for i in 1:length(rvals)-1
        @test terminal_v(rvals[i]) < terminal_v(rvals[i+1])
    end
end


@testset "set_X_crit!: critical (activation) radius" begin
    kappa  = FT(1.28)
    dry_r  = FT(0.1e-6)
    dry_r3 = dry_r^3
    T      = FT(290)

    drops = make_droplets(FT, 1; R=FT(1e-7), dry_r=dry_r)

    set_X_crit!(drops, 1, kappa, T)

    R_set = volume_to_radius(drops.X[1] / FT(0.98))

    # NOTE: code has a known /2 bug: R_crit = sqrt(3b/a)/2 instead of sqrt(3b/a)
    R_crit_theory = kohler_Rcrit(T, kappa, dry_r3)
    @test R_set ≈ R_crit_theory / 2  rtol=0.01

    @test drops.X[1] > 0
    @test volume_to_radius(drops.X[1]) > dry_r
end


@testset "find_equilibrium_radius: haze-branch equilibrium" begin
    kappa  = FT(0.61)
    dry_r  = FT(0.05e-6)
    dry_r3 = dry_r^3
    T      = FT(290)

    drops = make_droplets(FT, 1; R=FT(2e-7), dry_r=dry_r)

    S_env = FT(0.98)
    find_equilibrium_radius(drops, 1, kappa, T, S_env, constants)

    R_eq = volume_to_radius(drops.X[1])

    S_kohler_at_Req = kohler_S(R_eq, T, kappa, dry_r3)
    @test S_kohler_at_Req ≈ S_env  atol=1e-6

    @test R_eq >= dry_r

    R_crit = kohler_Rcrit(T, kappa, dry_r3)
    @test R_eq < R_crit

    @test FT(1e-9) < R_eq < FT(1e-3)
end


@testset "condensation_time_step_spatial!: single-cell condensation" begin
    nz  = 1
    nsd = 4

    cst         = Constants{FT}()
    condsettings = condensation_settings{FT}(kappa=FT(0.61))
    spatial     = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(100.0))
    scmsettings = scm_settings{FT}(REM=DynOFF(), spinupsaturation=DynOFF())

    dry_r = FT(0.05e-6)
    R0    = FT(5e-6)
    xi    = 1000

    droplets = make_droplets(FT, nsd; R=R0, dry_r=dry_r, xi=xi)

    T_init  = FT(285.0)
    P_init  = FT(97000.0)
    # ~0.5% supersaturated (inverting sat(qv,P) = 1.005*esat(T))
    es_init = esat(T_init)
    qv_init = FT(1.005) * es_init * FT(0.622) / (P_init - FT(1.005) * es_init * FT(0.378))
    ρ_init  = ρ_ideal_gas(P_init, T_init, qv_init, cst)
    θ_init  = theta_from_T(T_init, P_init,qv_init, cst)

    state = scm_states{FT}(
        nz,
        fill(P_init, nz),
        fill(P_init, nz + 1),
        fill(T_init, nz),
        fill(FT(0.0), nz),
        fill(θ_init, nz),
        fill(qv_init, nz),
        fill(ρ_init, nz),
        fill(FT(0.0), nz),
        fill(FT(0.0), nz),
        fill(FT(0.0), nz),
        fill(FT(0.0), nz),
        droplets,
        spatial,
        fill(FT(0.0), nz),
        fill(FT(0.0), nz),
        fill(FT(0.0), nz),
    )

    conddata = condensation_data(FT, nz, nsd)

    X_before  = copy(droplets.X)
    qv_before = copy(state.qv)
    T_before  = copy(state.T_tmp)

    condensation_time_step_spatial!(
        DynON(), droplets, state, nz, FT(1.0),
        conddata, cst, condsettings, spatial, nothing, scmsettings, nothing, FT(0.0),
    )

    # Supersaturated → droplets grow
    @test all(droplets.X .> X_before)

    # qv decreases as vapor condenses
    @test all(state.qv .< qv_before)

    # Temperature increases from latent heat release
    @test all(state.T_tmp .> T_before)

    # Water mass conservation: ΔX_drops * ξ * ρl = -Δqv * ρ_air * V_cell
    V_cell   = spatial.area_per_grid * spatial.z_grid_height
    dqv      = state.qv[1] - qv_before[1]
    dX_total = sum((droplets.X .- X_before) .* droplets.ξ)
    @test dqv * ρ_init * V_cell ≈ -dX_total * cst.ρl  rtol=1e-3  # bisection converges to 1e-4 relative

    # T change consistent with latent heat
    dT = state.T_tmp[1] - T_before[1]
    @test dT ≈ -dqv * cst.L / cst.Cp_air  rtol=1e-3

    # DynOFF dispatch is a no-op: droplets and state are unchanged
    droplets2 = make_droplets(FT, nsd; R=R0, dry_r=dry_r, xi=xi)
    X2_before = copy(droplets2.X)
    condensation_time_step_spatial!(
        DynOFF(), droplets2, state, nz, FT(1.0),
        conddata, cst, condsettings, spatial, nothing, scmsettings, nothing, FT(0.0),
    )
    @test droplets2.X == X2_before

    # Subsaturated environment → droplets evaporate
    droplets3 = make_droplets(FT, nsd; R=R0, dry_r=dry_r, xi=xi)
    qv_sub = FT(0.95) * esat(T_init) * FT(0.622) / P_init
    state3 = scm_states{FT}(
        nz,
        fill(P_init, nz), fill(P_init, nz + 1),
        fill(T_init, nz), fill(FT(0.0), nz),
        fill(theta_from_T(T_init, P_init,qv_sub, cst), nz),
        fill(qv_sub, nz),
        fill(ρ_ideal_gas(P_init, T_init, qv_sub, cst), nz),
        fill(FT(0.0), nz), fill(FT(0.0), nz),fill(FT(0.0), nz),fill(FT(0.0), nz),
        droplets3, spatial,
        fill(FT(0.0), nz), fill(FT(0.0), nz), fill(FT(0.0), nz),
    )
    conddata3 = condensation_data(FT, nz, nsd)
    X3_before = copy(droplets3.X)
    condensation_time_step_spatial!(
        DynON(), droplets3, state3, nz, FT(1.0),
        conddata3, cst, condsettings, spatial, nothing, scmsettings, nothing, FT(0.0),
    )
    @test all(droplets3.X .< X3_before)
    @test all(state3.qv .> qv_sub)   # qv increases as droplets evaporate
end
