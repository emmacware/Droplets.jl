using Test
using Droplets

FT = Float64

# -----------------------------------------------------------------------
# Physical constants used independently of the code's `constants` struct
# -----------------------------------------------------------------------
const _L   = 2.26e6          # J/kg, latent heat of vaporization
const _Rv  = 461.0           # J/(kg·K), gas constant for water vapor
const _k   = 0.024           # W/(m·K), thermal conductivity of air
const _Dv  = 2.26e-5         # m²/s, diffusivity of water vapor in air
const _ρl  = 1000.0          # kg/m³, liquid water density

# -----------------------------------------------------------------------
# Helper: analytic Köhler supersaturation for the kappa formulation
#   S_kohler(R) = 1 + akk(T)/R - kappa*dry_r3/R³
# At equilibrium, S_env == S_kohler(R_eq).
# -----------------------------------------------------------------------
kohler_S(R, T, kappa, dry_r3) =
    1.0 + (3.3e-7 / T) / R - kappa * dry_r3 / R^3

# Analytic critical radius (dS/dR = 0):  R_crit = sqrt(3·κ·r³/akk)
kohler_Rcrit(T, kappa, dry_r3) =
    sqrt(3.0 * kappa * dry_r3 / (3.3e-7 / T))

# -----------------------------------------------------------------------
@testset "esat: saturation vapor pressure" begin
    # August–Roche–Magnus formula: e_s(0°C) = 610.94 Pa exactly
    @test esat(FT(273.0)) ≈ 610.94 atol=0.01

    # Reference values from the Magnus formula, computed independently:
    #   e_s = 100 · 6.1094 · exp(17.625·(T−273) / ((T−273)+243.04))
    #   10°C (283 K): 610.94·exp(17.625·10/253.04)  ≈ 1225 Pa
    #   20°C (293 K): 610.94·exp(17.625·20/263.04)  ≈ 2338 Pa
    #   30°C (303 K): 610.94·exp(17.625·30/273.04)  ≈ 4228 Pa
    @test esat(FT(283)) ≈ 1225.0 atol=10.0
    @test esat(FT(293)) ≈ 2338.0 atol=20.0
    @test esat(FT(303)) ≈ 4228.0 atol=40.0

    # Always positive
    for T in [FT(240), FT(260), FT(273), FT(290), FT(310)]
        @test esat(T) > 0
    end

    # Strictly monotone increasing with temperature
    T_vals = [FT(250), FT(260), FT(270), FT(280), FT(290), FT(300)]
    for i in 1:length(T_vals)-1
        @test esat(T_vals[i]) < esat(T_vals[i+1])
    end

    # Clausius-Clapeyron: roughly doubles every 10 K near 0°C
    ratio = esat(FT(283)) / esat(FT(273))
    @test ratio > 1.6 && ratio < 2.2
end


@testset "sat: actual vapor pressure from mixing ratio" begin
    P = FT(97000)

    # Zero mixing ratio → zero vapor pressure
    @test sat(FT(0), P) == FT(0)

    # For very small qv: sat ≈ qv·P / 0.622 (the 0.378·qv term is negligible)
    qv_tiny = FT(1e-4)
    @test sat(qv_tiny, P) ≈ qv_tiny * P / 0.622  rtol=0.01

    # Linearly proportional to P (for fixed qv)
    qv = FT(5e-3)
    @test sat(qv, FT(2) * P) / sat(qv, P) ≈ 2.0  rtol=0.01

    # Positive for positive qv
    @test sat(FT(0.01), P) > 0

    # Absolute magnitude: near-saturated marine boundary layer
    #   T=290K → esat≈1934 Pa; for 98% RH: sat ≈ 0.98·1934 ≈ 1895 Pa
    #   From sat formula: qv·P/(0.378·qv+0.622) = 1895 →  qv ≈ 1895·0.622/(P−1895·0.378)
    #   ≈ 1178.9/(97000−716.7) ≈ 1178.9/96283 ≈ 0.01224 kg/kg
    esat_290 = FT(100 * 6.1094 * exp(17.625 * 17 / (17 + 243.04)))  # T=290K → (T-273)=17°C
    qv_98pct = FT(0.98 * esat_290 * 0.622 / (P - 0.98 * esat_290 * 0.378))
    @test sat(qv_98pct, P) ≈ FT(0.98) * esat_290  rtol=0.001
end


@testset "FK and FD: thermodynamic growth coefficients" begin
    # FK = (L/RvT - 1) · L·ρl / (k·T)  — should be large and positive at cloud T
    # FD = ρl·Rv·T / (Dv·esat(T))       — same

    for T in [FT(260), FT(273), FT(290), FT(305)]
        @test FK(T) > 0
        @test FD(T) > 0
        @test FK(T) + FD(T) > 0     # denominator of growth equation is positive
    end

    # FK should decrease as T increases (L/RvT dominates and decreases)
    @test FK(FT(280)) > FK(FT(300))

    # Independent magnitude check at T = 290 K.
    # From physical constants:
    #   FK = (L/(Rv·T) - 1) · L·ρl / (k·T)
    #      = (2.26e6/461/290 − 1) · 2.26e6·1000 / (0.024·290)
    #      = 15.91 · 3.247e8  ≈  5.16e9  s/m²
    #
    #   FD = ρl·Rv·T / (Dv·esat(T))
    #      = 1000·461·290 / (2.26e-5 · 1934)
    #      ≈  3.06e9  s/m²
    T = FT(290)
    FK_manual = (_L / (_Rv * T) - 1) * (_L * _ρl) / (_k * T)
    FD_manual = _ρl * _Rv * T / (_Dv * esat(T))
    @test FK(T) ≈ FK_manual  rtol=1e-4
    @test FD(T) ≈ FD_manual  rtol=1e-4

    @test FK(T) ≈ 5.16e9  rtol=0.05   # ~5.2 × 10⁹ s/m²
    @test FD(T) ≈ 3.06e9  rtol=0.05   # ~3.1 × 10⁹ s/m²

    # FK > FD at typical cloud temperatures (latent-heat term dominates diffusion term)
    @test FK(T) > FD(T)
end


@testset "drkohler (Senv form): radial growth rate" begin
    R_large = FT(5e-6)   # 5 μm — curvature and solute terms are negligible
    M_sol   = FT(1e-16)  # kg solute mass
    m_NaCl  = FT(58.44)  # g/mol NaCl
    T       = FT(290)
    dt      = FT(1.0)

    # Supersaturated environment → positive growth
    @test drkohler(R_large, M_sol, m_NaCl, T, FT(1.01), dt) > 0

    # Subsaturated environment → evaporation (negative dr)
    @test drkohler(R_large, M_sol, m_NaCl, T, FT(0.95), dt) < 0

    # Very large R, S = 1.0 exactly: dr ≈ 0 (Kelvin and solute both vanish)
    dr_eq = drkohler(FT(200e-6), M_sol, m_NaCl, T, FT(1.0), dt)
    @test abs(dr_eq) < 1e-6   # m/s, effectively zero

    # Guard: R + dr·dt must never be negative
    dr_bad = drkohler(FT(1e-8), M_sol, m_NaCl, T, FT(0.0), FT(100.0))
    @test FT(1e-8) + dr_bad * FT(100.0) >= 0

    # Higher supersaturation → faster growth
    dr1 = drkohler(R_large, M_sol, m_NaCl, T, FT(1.005), dt)
    dr2 = drkohler(R_large, M_sol, m_NaCl, T, FT(1.010), dt)
    @test dr2 > dr1

    # Larger R → slower radial growth rate (dr ∝ 1/R for large drops)
    dr_5  = drkohler(FT(5e-6),  M_sol, m_NaCl, T, FT(1.01), dt)
    dr_20 = drkohler(FT(20e-6), M_sol, m_NaCl, T, FT(1.01), dt)
    @test dr_5 > dr_20

    # Absolute magnitude at R=5μm, S=1.01, T=290K.
    # At this radius the Kelvin term (akk/R = 3.3e-7/290/5e-6 ≈ 2.3e-4) and
    # solute term are both small vs (S-1)=0.01, so dr ≈ (S-1)/((FK+FD)·R):
    esat_T    = 100 * 6.1094 * exp(17.625 * (T - 273) / ((T - 273) + 243.04))
    FK_T      = (_L / (_Rv * T) - 1) * (_L * _ρl) / (_k * T)
    FD_T      = _ρl * _Rv * T / (_Dv * esat_T)
    dr_approx = (FT(1.01) - 1) / ((FK_T + FD_T) * R_large)
    @test drkohler(R_large, M_sol, m_NaCl, T, FT(1.01), dt) ≈ dr_approx  rtol=0.05
end


@testset "drkohler (qv/P form): consistency with Senv form" begin
    R     = FT(5e-6)
    M_sol = FT(1e-16)
    m     = FT(58.44)
    T     = FT(290)
    P     = FT(97000)
    dt    = FT(1.0)

    # Construct qv such that sat(qv,P)/esat(T) = 1.01
    # For small qv: sat ≈ qv·P/0.622, so qv ≈ 1.01·esat(T)·0.622/P
    qv_super = FT(1.01 * esat(T) * 0.622 / P)
    @test drkohler(R, M_sol, m, T, qv_super, P, dt) > 0   # supersaturated

    qv_sub = FT(0.95 * esat(T) * 0.622 / P)
    @test drkohler(R, M_sol, m, T, qv_sub, P, dt) < 0     # subsaturated

    # Both overloads must give the same answer when Senv = sat(qv,P)/esat(T)
    qv   = FT(0.012)
    Senv = sat(qv, P) / esat(T)
    @test drkohler(R, M_sol, m, T, Senv, dt) ≈
          drkohler(R, M_sol, m, T, qv, P, dt)  rtol=1e-9
end


@testset "dXkohler_function_of_radius: volume growth" begin
    R     = FT(5e-6)
    M_sol = FT(1e-16)
    m     = FT(58.44)
    T     = FT(290)
    dt    = FT(1.0)
    S     = FT(1.01)

    # dX = 4π R² · dR  (surface-area weighted)
    dR       = drkohler(R, M_sol, m, T, S, dt)
    dX_from_dR = 4 * π * R^2 * dR
    @test dXkohler_function_of_radius(R, M_sol, m, T, S, dt) ≈ dX_from_dR  rtol=1e-10

    # Sign follows dR: supersaturated → positive volume change
    @test dXkohler_function_of_radius(R, M_sol, m, T, S, dt) > 0

    # BUG: the qv/P overload calls drkohler without the timestep argument.
    # This should throw a MethodError (missing timestep).
    P = FT(97000)
    @test_throws MethodError dXkohler_function_of_radius(R, M_sol, m, T, FT(0.01), P, dt)
end


@testset "dXkohler_function_of_radius_activated: activated volume growth" begin
    R  = FT(10e-6)
    T  = FT(290)
    dt = FT(1.0)
    S  = FT(1.005)

    dR_act = drkohler_activated(R, T, S, dt)
    @test dXkohler_function_of_radius_activated(R, T, S, dt) ≈ 4 * π * R^2 * dR_act  rtol=1e-10

    # Sign: positive for supersaturated environment
    @test dXkohler_function_of_radius_activated(R, T, S, dt) > 0
end


@testset "drkappakohler: kappa-Köhler growth rate" begin
    dry_r  = FT(0.05e-6)     # 50 nm dry radius
    dry_r3 = dry_r^3
    kappa  = FT(1.28)        # NaCl kappa
    T      = FT(290)
    dt     = FT(1.0)
    R_large = FT(5e-6)

    # Supersaturated → grow
    @test drkappakohler(R_large, dry_r3, kappa, T, FT(1.01), dt) > 0

    # Subsaturated → shrink
    @test drkappakohler(R_large, dry_r3, kappa, T, FT(0.95), dt) < 0

    # Guard: R never goes negative
    dr_evap = drkappakohler(FT(1e-8), dry_r3, kappa, T, FT(0.0), FT(100.0))
    @test FT(1e-8) + dr_evap * FT(100.0) >= 0

    # For very large R (both Kelvin and kappa terms → 0):
    #   dr ≈ (Senv - 1) / (denom · R)
    R_huge = FT(200e-6)
    Senv   = FT(1.01)
    dr_approx = (Senv - 1) / ((FK(T) + FD(T)) * R_huge)
    @test drkappakohler(R_huge, dry_r3, kappa, T, Senv, dt) ≈ dr_approx  rtol=0.05

    # Higher kappa (stronger solute effect) lowers the effective supersaturation
    # needed to grow a small droplet — at fixed S_env and small R, higher kappa → more growth
    R_small = FT(0.15e-6)
    Senv_small = FT(1.002)
    dr_hi_kappa  = drkappakohler(R_small, dry_r3, FT(1.2),  T, Senv_small, dt)
    dr_lo_kappa  = drkappakohler(R_small, dry_r3, FT(0.1),  T, Senv_small, dt)
    @test dr_hi_kappa > dr_lo_kappa

    # kappa = 0 (pure water droplet): only Kelvin term remains
    # dr = (Senv - 1 - akk/R) / (denom · R), positive only when Senv > 1 + akk/R
    akk_T = 3.3e-7 / T
    Senv_above_kelvin = FT(1.0 + akk_T / R_large + 0.001)
    @test drkappakohler(R_large, dry_r3, FT(0.0), T, Senv_above_kelvin, dt) > 0
end


@testset "drkohler_activated: activated droplet growth" begin
    R  = FT(10e-6)
    T  = FT(290)
    dt = FT(1.0)

    # Supersaturated → grow
    @test drkohler_activated(R, T, FT(1.005), dt) > 0

    # Subsaturated → shrink
    @test drkohler_activated(R, T, FT(0.99),  dt) < 0

    # Exactly saturated (Senv = 1): dr = 0 exactly (no Kelvin or solute terms)
    @test drkohler_activated(R, T, FT(1.0), dt) == FT(0)

    # Guard: R never goes negative
    dr = drkohler_activated(FT(1e-9), T, FT(0.0), FT(1000.0))
    @test FT(1e-9) + dr * FT(1000.0) >= 0

    # dr = (Senv-1)/(denom·R) → linear in (Senv-1)
    dr1 = drkohler_activated(R, T, FT(1.01), dt)
    dr2 = drkohler_activated(R, T, FT(1.02), dt)
    @test dr2 ≈ dr1 * 2  rtol=1e-6   # exactly linear

    # dr ∝ 1/R: smaller drop grows faster in radius at same S
    dr_small = drkohler_activated(FT(5e-6),  T, FT(1.005), dt)
    dr_large = drkohler_activated(FT(20e-6), T, FT(1.005), dt)
    @test dr_small > dr_large
    @test dr_small / dr_large ≈ 4.0  rtol=1e-6   # inversely proportional to R

    # -----------------------------------------------------------------------
    # Absolute magnitude: classic diffusional growth result.
    # At R = 10 μm, S = 1.01, T = 290 K:
    #   dr/dt = (S−1) / ((FK+FD) · R)
    #         = 0.01 / (8.22e9 · 10e-6)
    #         ≈ 1.22e-7  m/s  ≈ 0.12 μm/s
    # (well-known result from Pruppacher & Klett / Seinfeld & Pandis)
    #
    # Compute denom independently from physical constants to avoid circularity:
    T_mag  = FT(290)
    R_mag  = FT(10e-6)
    S_mag  = FT(1.01)
    esat_290  = 100 * 6.1094 * exp(17.625 * (T_mag - 273) / ((T_mag - 273) + 243.04))
    FK_290 = (_L / (_Rv * T_mag) - 1) * (_L * _ρl) / (_k * T_mag)
    FD_290 = _ρl * _Rv * T_mag / (_Dv * esat_290)
    dr_expected = (S_mag - 1) / ((FK_290 + FD_290) * R_mag)

    @test drkohler_activated(R_mag, T_mag, S_mag, dt) ≈ dr_expected  rtol=0.01
    @test dr_expected ≈ 1.2e-7  rtol=0.1   # order-of-magnitude sanity check
end


@testset "terminal_v: droplet terminal velocity" begin
    # All terminal velocities must be positive
    for r in [FT(1e-6), FT(10e-6), FT(25e-6), FT(35e-6), FT(55e-6),
              FT(100e-6), FT(500e-6), FT(1e-3), FT(3e-3)]
        @test terminal_v(r) > 0
    end

    # Stokes regime (r < 30 μm): v ∝ r² → ratio test
    # v(20μm)/v(10μm) = (20/10)² = 4.0  exactly (same formula, same regime)
    @test terminal_v(FT(20e-6)) / terminal_v(FT(10e-6)) ≈ 4.0  rtol=1e-6

    # -----------------------------------------------------------------------
    # Absolute magnitudes derived from the parameterization formulas directly.
    # (All formulas in cm-space; final result in m/s = cm/s × 1e-2.)
    #
    # Stokes regime (r < 30 μm): vt [cm/s] = 1.19e6 · r_cm²
    #   r = 10 μm → r_cm = 1e-3 cm → vt = 1.19e6·(1e-3)² = 1.19 cm/s = 0.0119 m/s
    @test terminal_v(FT(10e-6)) ≈ 0.0119  rtol=1e-4

    # Stokes law cross-check (independent): v = 2·ρl·r²·g / (9·η)
    #   with ρl=1000, g=9.8, η≈1.8e-5 Pa·s → 0.0121 m/s  (within 2%)
    @test terminal_v(FT(10e-6)) ≈ 2 * 1000 * (10e-6)^2 * 9.8 / (9 * 1.8e-5)  rtol=0.03

    # Intermediate regime (30 μm < r < 60 μm): vt [cm/s] = 8e3 · r_cm
    #   r = 40 μm → r_cm = 4e-3 cm → vt = 8e3·4e-3 = 32 cm/s = 0.32 m/s
    @test terminal_v(FT(40e-6)) ≈ 0.32  rtol=1e-4

    # Intermediate regime (30 μm < r < 60 μm): v ∝ r → ratio test
    # v(50μm)/v(40μm) = 50/40 = 1.25  exactly
    @test terminal_v(FT(50e-6)) / terminal_v(FT(40e-6)) ≈ 1.25  rtol=1e-6

    # Large-drop regime (60 μm < r < 2 mm): vt [cm/s] = 2.01e3 · sqrt(r_cm)
    #   r = 100 μm → r_cm = 0.01 cm → vt = 2.01e3·√0.01 = 201 cm/s = 2.01 m/s
    @test terminal_v(FT(100e-6)) ≈ 2.01  rtol=1e-4

    # Large-drop regime: v ∝ r^0.5 → ratio test
    # v(400μm)/v(100μm) = (400/100)^0.5 = 2.0
    @test terminal_v(FT(400e-6)) / terminal_v(FT(100e-6)) ≈ 2.0  rtol=1e-6

    # Terminal velocity caps at 2 mm radius — anything larger returns the same value
    @test terminal_v(FT(3e-3)) == terminal_v(FT(2e-3))
    @test terminal_v(FT(5e-3)) == terminal_v(FT(2e-3))

    # Monotonically increasing from 10 μm to 2 mm
    rvals = [FT(10e-6), FT(25e-6), FT(45e-6), FT(100e-6), FT(500e-6), FT(2e-3)]
    for i in 1:length(rvals)-1
        @test terminal_v(rvals[i]) < terminal_v(rvals[i+1])
    end
end


@testset "set_X_crit!: critical (activation) radius" begin
    kappa  = FT(1.28)         # NaCl
    dry_r  = FT(0.1e-6)       # 100 nm dry radius
    dry_r3 = dry_r^3
    T      = FT(290)

    drops = droplet_attributes_1d{FT}(
        Int[1],
        FT[radius_to_volume(FT(1e-7))],   # initial X (will be overwritten)
        FT[dry_r3],
        FT[0.0],
        Int[1],
        Vector{Union{UnitRange{Int},Nothing}}([1:1]),
        Int[1],
    )

    set_X_crit!(drops, 1, kappa, T)

    # Recover the radius that was set (undoing the 0.98 factor applied in the code)
    R_set = volume_to_radius(drops.X[1] / FT(0.98))

    # Theoretical critical radius: R_crit = sqrt(3·κ·r_d³ / akk(T))
    # (from dS_kohler/dR = 0 → akk/R² = 3κr_d³/R⁴ → R² = 3κr_d³/akk)
    R_crit_theory = kohler_Rcrit(T, kappa, dry_r3)

    # This test FAILS if the code has an extra /2 factor on R_crit
    @test R_set ≈ R_crit_theory  rtol=0.01

    # At the true R_crit, the Köhler curve has a local maximum, i.e. dS/dR = 0:
    #   -akk/R² + 3κr_d³/R⁴ ≈ 0
    akk_T = 3.3e-7 / T
    dSdR = -akk_T / R_crit_theory^2 + 3 * kappa * dry_r3 / R_crit_theory^4
    @test abs(dSdR) < 1e3   # essentially zero (numerically, not in code)

    # Sanity: X must be positive and non-tiny
    @test drops.X[1] > 0
    @test volume_to_radius(drops.X[1]) > dry_r   # activated radius > dry radius
end


@testset "find_equilibrium_radius: haze-branch equilibrium" begin
    kappa  = FT(0.61)        # ammonium sulfate
    dry_r  = FT(0.05e-6)     # 50 nm dry radius
    dry_r3 = dry_r^3
    T      = FT(290)

    drops = droplet_attributes_1d{FT}(
        Int[1],
        FT[radius_to_volume(FT(2e-7))],   # initial guess
        FT[dry_r3],
        FT[0.0],
        Int[1],
        Vector{Union{UnitRange{Int},Nothing}}([1:1]),
        Int[1],
    )

    S_env = FT(0.98)   # subsaturated — stable haze equilibrium exists on left Köhler branch
    find_equilibrium_radius(drops, 1, kappa, T, S_env)

    R_eq = volume_to_radius(drops.X[1])

    # The equilibrium radius must satisfy the Köhler equation:
    #   S_env = 1 + akk(T)/R_eq - kappa·dry_r3/R_eq³
    # This test FAILS if the Newton step has a sign error and fails to converge
    S_kohler_at_Req = kohler_S(R_eq, T, kappa, dry_r3)
    @test S_kohler_at_Req ≈ S_env  atol=1e-6

    # R_eq must be larger than the dry radius
    @test R_eq >= dry_r

    # R_eq must be below the critical radius (haze branch is sub-critical)
    R_crit = kohler_Rcrit(T, kappa, dry_r3)
    @test R_eq < R_crit

    # R_eq must be physically reasonable
    @test FT(1e-9) < R_eq < FT(1e-3)
end
