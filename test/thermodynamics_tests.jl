using Test
using Droplets


FT = Float64

@testset "mixing_ratio / specific_humidity: inverse of each other" begin
    for q in FT[0.0, 1e-4, 5e-3, 0.02]
        @test specific_humidity(mixing_ratio(q)) ≈ q rtol=1e-12
    end
    for r in FT[0.0, 1e-4, 5e-3, 0.02]
        @test mixing_ratio(specific_humidity(r)) ≈ r rtol=1e-12
    end
end

@testset "calc_P_dry_from_P / calc_P_from_P_dry: round trip" begin
    c = Constants{FT}()
    P = FT(97000.0)
    for qv in FT[0.0, 1e-3, 8e-3, 0.02]
        P_dry = calc_P_dry_from_P(P, qv, c)
        @test calc_P_from_P_dry(P_dry, qv, c) ≈ P rtol=1e-12
        @test P_dry <= P  # dry partial pressure never exceeds total pressure
    end
end

@testset "calc_ρ_dry_from_ρ / calc_ρ_from_ρ_dry: round trip" begin
    ρ = FT(1.1)
    for qv in FT[0.0, 1e-3, 8e-3, 0.02]
        ρ_dry = calc_ρ_dry_from_ρ(ρ, qv)
        @test calc_ρ_from_ρ_dry(ρ_dry, qv) ≈ ρ rtol=1e-12
        @test ρ_dry <= ρ
    end
end

@testset "T_virtual / θ_virtual: reduce to identity at qv=0, exceed T at qv>0" begin
    c = Constants{FT}()
    T = FT(290.0)
    θ = FT(300.0)
    @test T_virtual(T, FT(0.0), c) ≈ T
    @test θ_virtual(θ, FT(0.0), c) ≈ θ
    @test T_virtual(T, FT(0.01), c) > T   # water vapor is lighter than dry air
    @test θ_from_θv(θ_virtual(θ, FT(0.01), c), FT(0.01), c) ≈ θ rtol=1e-12
end

@testset "theta_from_T / T_from_theta: round trip, both dry_theta_convention branches" begin
    for dry_conv in (false, true)
        c = Constants{FT}(dry_theta_convention=dry_conv)
        T = FT(285.0)
        P = FT(95000.0)
        for qv in FT[0.0, 2e-3, 1e-2]
            θ = theta_from_T(T, P, qv, c)
            @test T_from_theta(θ, P, qv, c) ≈ T rtol=1e-10

            θvec = zeros(FT, 3)
            theta_from_T!(θvec, fill(T, 3), fill(P, 3), fill(qv, 3), c)
            @test all(θvec .≈ θ)

            Tvec = zeros(FT, 3)
            T_from_theta!(Tvec, fill(θ, 3), fill(P, 3), fill(qv, 3), c)
            @test all(Tvec .≈ T)
        end
    end
end

@testset "ρ_ideal_gas: satisfies the ideal gas law it was derived from" begin
    P = FT(96000.0)
    T = FT(288.0)
    qv = FT(8e-3)

    c_moist = Constants{FT}(dry_theta_convention=false)
    ρ = ρ_ideal_gas(P, T, qv, c_moist)
    @test ρ * c_moist.Rd * T_virtual(T, qv, c_moist) ≈ P rtol=1e-10

    c_dry = Constants{FT}(dry_theta_convention=true)
    ρd = ρ_ideal_gas(P, T, qv, c_dry)
    r = mixing_ratio(qv)
    Rm = (c_dry.Rd + r * c_dry.Rv) / (1 + r)
    @test ρd * Rm * T ≈ P rtol=1e-10

    ρvec = zeros(FT, 4)
    ρ_ideal_gas!(ρvec, fill(P, 4), fill(T, 4), fill(qv, 4), c_moist)
    @test all(ρvec .≈ ρ)
end

@testset "ρ_calc_θ: consistent with theta_from_T/ρ_ideal_gas composed directly" begin
    c = Constants{FT}()
    P = FT(94000.0)
    θ = FT(291.0)
    qv = FT(6e-3)
    ρ = ρ_calc_θ(P, θ, qv, c)
    T = T_from_theta(θ, P, qv, c)
    @test ρ ≈ ρ_ideal_gas(P, T, qv, c) rtol=1e-12

    ρvec = zeros(FT, 3)
    ρ_calc_θ!(ρvec, fill(P, 3), fill(θ, 3), fill(qv, 3), c)
    @test all(ρvec .≈ ρ)
end

@testset "θl / θ_from_θl: round trip via T (not θ directly, per θl's own definition)" begin
    for dry_conv in (false, true)
        c = Constants{FT}(dry_theta_convention=dry_conv)
        P = FT(93000.0)
        T = FT(283.0)
        qv = FT(7e-3)
        ql = FT(2e-4)

        θl_val = θl(P, T, ql, qv, c)
        # θ_from_θl inverts θl's own (P0/P)^.. * (T - L ql/Cp) construction, so composing
        # them should recover theta_from_T(T,...) 
        @test θ_from_θl(P, θl_val, ql, qv, c) ≈ theta_from_T(T, P, qv, c) rtol=1e-10

        # ql=0 (no liquid): θl reduces to ordinary potential temperature
        @test θl(P, T, FT(0.0), qv, c) ≈ theta_from_T(T, P, qv, c) rtol=1e-10
    end
end

@testset "compute_ql_at_cell! / compute_ql_all_cells!: agree and match the manual sum" begin
    FT_ = Float64
    nz = 4
    ξ = [50, 50, 20, 20, 20]
    R = FT_[5e-6, 8e-6, 3e-6, 3e-6, 3e-6]
    X = radius_to_volume.(R)
    cell_id = [1, 1, 2, 3, 3]
    droplets = droplet_attributes_1d{FT_}(
        ξ, X, zeros(FT_, 5), zeros(FT_, 5), cell_id, zeros(FT_, 5),
        [1:2, 3:3, 4:5, 1:0], collect(1:5),
    )
    spatial = spatial_settings_1d{FT_}(Nz=nz, Z_max=FT_(4 * 50.0), area_per_grid=FT_(10.0))
    ρ = fill(FT_(1.1), nz)
    c = Constants{FT_}()

    state = (droplets=droplets, ρ=ρ, spatial=spatial, ql_tmp=zeros(FT_, nz))

    for k in 1:nz
        compute_ql_at_cell!(state, k, c)
    end
    ql_grid_range = copy(state.ql_tmp)

    state2 = (droplets=droplets, ρ=ρ, spatial=spatial, ql_tmp=zeros(FT_, nz))
    compute_ql_all_cells!(state2, c)

    @test ql_grid_range ≈ state2.ql_tmp rtol=1e-12

    V = spatial.area_per_grid * spatial.z_grid_height
    for k in 1:nz
        idx = findall(==(k), cell_id)
        expected = isempty(idx) ? zero(FT_) :
            sum(X[idx] .* ξ[idx]) * c.ρl / (ρ[k] * V)
        @test ql_grid_range[k] ≈ expected rtol=1e-10
    end
    @test ql_grid_range[4] == 0.0  # empty cell 4
end
