using Droplets
using Test

# Helper: allocate tmp for nz-cell 1D grid (GCz has nz+1 faces)
make1d_tmp(nz) = mpdata_tmp_1d(zeros(nz), zeros(nz + 1))

@testset "MPDATA 1D" begin

    # ===================================================================
    @testset "limit — Periodic" begin
        N = 5
        @test limit(Periodic(), 1,  N) == 1      # in-range
        @test limit(Periodic(), N,  N) == N      # in-range
        @test limit(Periodic(), 0,  N) == N      # underflow wraps to end
        @test limit(Periodic(), N+1,N) == 1      # overflow wraps to start
        @test limit(Periodic(), -1, N) == N - 1  # double underflow
        @test limit(Periodic(), N+2,N) == 2      # double overflow
    end

    @testset "limit — NoFlux" begin
        N = 5
        @test limit(NoFlux(), 1,   N) == 1   # in-range
        @test limit(NoFlux(), N,   N) == N   # in-range
        @test limit(NoFlux(), N+1, N) == N   # one above: stays at N
        @test limit(NoFlux(), N+2, N) == N-1 # two above: reflects inward
        @test limit(NoFlux(), 0,   N) == 1   # one below: reflects to 1
    end

    # ===================================================================
    @testset "flux — standard (not infinite_gauge)" begin
        # positive GC → donor is left cell
        @test flux(2.0, 5.0,  0.3) ≈  0.3 * 2.0
        # negative GC → donor is right cell
        @test flux(2.0, 5.0, -0.3) ≈ -0.3 * 5.0
        # zero GC → zero flux
        @test flux(2.0, 5.0,  0.0) == 0.0
    end

    @testset "flux — infinite_gauge returns GC regardless of ϕ" begin
        @test flux(2.0, 5.0,  0.3; infinite_gauge=true) ==  0.3
        @test flux(2.0, 5.0, -0.3; infinite_gauge=true) == -0.3
        @test flux(0.0, 0.0,  0.5; infinite_gauge=true) ==  0.5
    end

    # ===================================================================
    @testset "donor_cell_pass! 1D — mass conservation, Periodic" begin
        nz = 10
        ϕ  = rand(nz)
        GCz = fill(0.3, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        ϕ0  = sum(ϕ)
        donor_cell_pass!(ϕ, tmp, Periodic())
        @test sum(ϕ) ≈ ϕ0 atol=1e-12
    end

    @testset "donor_cell_pass! 1D — mass conservation, NoFlux + zero boundary velocity" begin
        # With GCz[1]=GCz[end]=0 and NoFlux, no mass crosses boundaries
        nz = 10
        ϕ  = rand(nz)
        GCz = fill(0.3, nz + 1); GCz[1] = 0.0; GCz[end] = 0.0
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        ϕ0  = sum(ϕ)
        donor_cell_pass!(ϕ, tmp, NoFlux(), topcellreservoir=false)
        @test sum(ϕ) ≈ ϕ0 atol=1e-12
    end

    @testset "donor_cell_pass! 1D — topcellreservoir makes top cell non-conservative" begin
        nz = 6
        GCz_internal = fill(0.2, nz + 1); GCz_internal[1] = 0.0; GCz_internal[end] = 0.0

        ϕ_conservative = rand(nz)
        tmp1 = make1d_tmp(nz); tmp1.GCz_step .= GCz_internal
        ϕ_conservative_copy = copy(ϕ_conservative)
        donor_cell_pass!(ϕ_conservative_copy, tmp1, NoFlux(), topcellreservoir=false)

        ϕ_reservoir = copy(ϕ_conservative)
        tmp2 = make1d_tmp(nz); tmp2.GCz_step .= GCz_internal
        donor_cell_pass!(ϕ_reservoir, tmp2, NoFlux(), topcellreservoir=true)

        # Conservative version conserves; reservoir version does not (top cell frozen,
        # but mass still leaves other cells toward it)
        @test sum(ϕ_conservative_copy) ≈ sum(ϕ_conservative) atol=1e-12
        @test abs(sum(ϕ_reservoir) - sum(ϕ_conservative)) > 1e-12
    end

    @testset "donor_cell_pass! 1D — uniform field is unchanged, Periodic" begin
        nz = 8
        ϕ  = fill(3.0, nz)
        GCz = fill(0.4, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        donor_cell_pass!(ϕ, tmp, Periodic())
        @test all(ϕ .≈ 3.0)
    end

    @testset "donor_cell_pass! 1D — positive velocity moves mass upward" begin
        nz = 8
        ϕ  = zeros(nz); ϕ[3] = 1.0    # single spike at cell 3
        ϕ0 = copy(ϕ)
        GCz = fill(0.4, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        donor_cell_pass!(ϕ, tmp, Periodic())
        # Cell 3 should lose mass; cell 4 (upward) should gain
        @test ϕ[3] < ϕ0[3]
        @test ϕ[4] > ϕ0[4]
    end

    @testset "donor_cell_pass! 1D — negative velocity moves mass downward" begin
        nz = 8
        ϕ  = zeros(nz); ϕ[5] = 1.0    # single spike at cell 5
        ϕ0 = copy(ϕ)
        GCz = fill(-0.4, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        donor_cell_pass!(ϕ, tmp, Periodic())
        @test ϕ[5] < ϕ0[5]
        @test ϕ[4] > ϕ0[4]
    end

    @testset "donor_cell_pass! 1D — infinite_gauge uses GC directly as flux" begin
        nz = 6
        ϕ  = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        GCz = fill(0.3, nz + 1)

        ϕ_std = copy(ϕ)
        tmp1 = make1d_tmp(nz); tmp1.GCz_step .= GCz
        donor_cell_pass!(ϕ_std, tmp1, Periodic(), infinite_gauge=false)

        ϕ_ig  = copy(ϕ)
        tmp2 = make1d_tmp(nz); tmp2.GCz_step .= GCz
        donor_cell_pass!(ϕ_ig, tmp2, Periodic(), infinite_gauge=true)

        # With infinite_gauge each face contributes ±GC regardless of ϕ;
        # the results differ from the standard upwind case
        @test ϕ_std ≠ ϕ_ig
        # Both should conserve mass (Periodic)
        @test sum(ϕ_std) ≈ sum(ϕ) atol=1e-12
        @test sum(ϕ_ig)  ≈ sum(ϕ) atol=1e-12
    end

    # ===================================================================
    @testset "find_extrema! 1D — interior cell" begin
        nz = 5
        ϕ  = [1.0, 5.0, 3.0, 2.0, 4.0]
        tmp = make1d_tmp(nz)
        find_extrema!(ϕ, tmp, vbc=Periodic())
        # Cell 2 (value 5): neighbors are cells 1 (1.0) and 3 (3.0)
        @test tmp.minmax.localmax[2] == 5.0
        @test tmp.minmax.localmin[2] == 1.0
        # Cell 3 (value 3): neighbors 5.0 and 2.0
        @test tmp.minmax.localmax[3] == 5.0
        @test tmp.minmax.localmin[3] == 2.0
    end

    @testset "find_extrema! 1D — periodic boundary cell 1" begin
        nz = 5
        ϕ  = [1.0, 5.0, 3.0, 2.0, 8.0]
        tmp = make1d_tmp(nz)
        find_extrema!(ϕ, tmp, vbc=Periodic())
        # Cell 1: km = limit(Periodic,0,5) = 5 (value 8), kp = 2 (value 5)
        @test tmp.minmax.localmax[1] == 8.0
        @test tmp.minmax.localmin[1] == 1.0
    end

    @testset "find_extrema! 1D — NoFlux boundary cell 1" begin
        nz = 5
        ϕ  = [1.0, 5.0, 3.0, 2.0, 8.0]
        tmp = make1d_tmp(nz)
        find_extrema!(ϕ, tmp, vbc=NoFlux())
        # Cell 1: km = limit(NoFlux,0,5) = 1 (self-reflection), kp = 2
        @test tmp.minmax.localmax[1] == max(ϕ[1], ϕ[1], ϕ[2])  # 5.0
        @test tmp.minmax.localmin[1] == min(ϕ[1], ϕ[1], ϕ[2])  # 1.0
    end

    # ===================================================================
    @testset "beta_function_1d — values clipped to [0,1]" begin
        nz = 6
        ϕ  = [1.0, 3.0, 0.5, 4.0, 2.0, 1.5]
        tmp = make1d_tmp(nz)
        find_extrema!(ϕ, tmp, vbc=Periodic())
        for k in 1:nz
            km1 = limit(Periodic(), k - 1, nz)
            kp1 = limit(Periodic(), k + 1, nz)
            for ig in (false, true)
                β_in, β_out = beta_function_1d(ϕ, tmp, k, km1, kp1,
                    0.3, 0.3, ig)
                @test 0.0 ≤ β_in  ≤ 1.0
                @test 0.0 ≤ β_out ≤ 1.0
            end
        end
    end

    @testset "beta_function_1d — local max: β_in=0 (can't bring more in)" begin
        # Cell k is the local maximum → β_in numerator = 0
        nz = 5
        ϕ  = [1.0, 5.0, 1.0, 1.0, 1.0]   # cell 2 is a spike
        tmp = make1d_tmp(nz)
        find_extrema!(ϕ, tmp, vbc=Periodic())
        # k=2, km1=1, kp1=3
        β_in, β_out = beta_function_1d(ϕ, tmp, 2, 1, 3, 0.4, 0.4, false)
        # localmax[2] = 5 = ϕ[2] → numerator of β_in = 0
        @test β_in == 0.0
        # β_out should be positive (cell has mass to give)
        @test β_out > 0.0
    end

    @testset "beta_function_1d — local min: β_out=0 (can't send more out)" begin
        nz = 5
        ϕ  = [5.0, 0.5, 5.0, 5.0, 5.0]   # cell 2 is a valley
        tmp = make1d_tmp(nz)
        find_extrema!(ϕ, tmp, vbc=Periodic())
        β_in, β_out = beta_function_1d(ϕ, tmp, 2, 1, 3, 0.4, 0.4, false)
        # localmin[2] = 0.5 = ϕ[2] → numerator of β_out = 0
        @test β_out == 0.0
        @test β_in  > 0.0
    end

    # ===================================================================
    @testset "antiosc! 1D — scaled fluxes ≤ original in magnitude" begin
        nz = 8
        ϕ   = [1.0, 3.0, 2.0, 5.0, 4.0, 1.0, 2.0, 3.0]
        GCz = fill(0.4, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        settings = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                      nonoscillatory=true)
        find_extrema!(ϕ, tmp, vbc=Periodic())
        antiosc!(ϕ, tmp, settings, vbc=Periodic())
        @test all(abs.(tmp.GCz_step) .≤ abs.(GCz) .+ 1e-14)
    end

    @testset "antiosc! 1D — Periodic: face 1 is consistent with periodic wrap (bug check)" begin
        # With Periodic BC, face 1 and face nz+1 should be the same wrap-around face.
        # Currently antiosc! sets tmp.GCz_tmp[1] = tmp.GCz_tmp[end] where end=nz+1,
        # but GCz_tmp[nz+1] is never written in the loop — this should be 0, which is wrong.
        # This test documents the expected (correct) behavior so the bug is visible.
        nz = 6
        ϕ   = [1.0, 4.0, 2.0, 5.0, 3.0, 0.5]
        GCz = fill(0.3, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        settings = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                      nonoscillatory=true)
        find_extrema!(ϕ, tmp, vbc=Periodic())
        antiosc!(ϕ, tmp, settings, vbc=Periodic())
        # Face 1 = the periodic wrap face between cell nz and cell 1.
        # The loop writes faces 2..nz; tmp.GCz_step[1] should match face nz
        # (the other end of the wrap), but due to the bug it is 0.
        # Correct assertion: tmp.GCz_step[1] == tmp.GCz_step[nz] (periodic symmetry)
        # If this fails, the bug is present.
        @test_broken tmp.GCz_step[1] == tmp.GCz_step[nz]
    end

    @testset "antiosc! 1D — NoFlux: interior faces bounded" begin
        nz = 6
        ϕ   = [1.0, 4.0, 2.0, 5.0, 3.0, 1.0]
        GCz = fill(0.3, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        settings = mpdata_settings_1d(nz, vertical_boundary_condition=NoFlux(),
                                      nonoscillatory=true)
        find_extrema!(ϕ, tmp, vbc=NoFlux())
        antiosc!(ϕ, tmp, settings, vbc=NoFlux())
        @test all(abs.(tmp.GCz_step[2:nz]) .≤ abs.(GCz[2:nz]) .+ 1e-14)
    end

    # ===================================================================
    @testset "compute_antidiffusive_velocity! 1D — zero for uniform field" begin
        nz = 8
        ϕ   = fill(3.0, nz)
        GCz = fill(0.4, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        settings = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic())
        compute_antidiffusive_velocity!(ϕ, tmp, settings, vbc=Periodic())
        # Uniform field → zero gradient → zero antidiffusive correction
        @test all(tmp.GCz_step .≈ 0.0)
    end

    @testset "compute_antidiffusive_velocity! 1D — nonzero for non-uniform field" begin
        nz = 6
        ϕ   = [1.0, 3.0, 2.0, 4.0, 1.0, 2.0]
        GCz = fill(0.3, nz + 1)
        tmp = make1d_tmp(nz); tmp.GCz_step .= GCz
        settings = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic())
        compute_antidiffusive_velocity!(ϕ, tmp, settings, vbc=Periodic())
        # Should produce at least some nonzero antidiffusive velocities
        @test any(abs.(tmp.GCz_step) .> 1e-14)
    end

    @testset "compute_antidiffusive_velocity! 1D — infinite_gauge vs standard differ" begin
        nz = 6
        ϕ   = [1.0, 4.0, 2.0, 5.0, 3.0, 2.0]
        GCz = fill(0.35, nz + 1)

        tmp_n = make1d_tmp(nz); tmp_n.GCz_step .= GCz
        s_n   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                   infinite_gauge=false)
        compute_antidiffusive_velocity!(ϕ, tmp_n, s_n, vbc=Periodic())

        tmp_g = make1d_tmp(nz); tmp_g.GCz_step .= GCz
        s_g   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                   infinite_gauge=true)
        compute_antidiffusive_velocity!(ϕ, tmp_g, s_g, vbc=Periodic())

        # infinite_gauge uses denominator=2; standard uses ψL+ψR+eps → must differ
        @test tmp_n.GCz_step ≠ tmp_g.GCz_step
    end

    # ===================================================================
    @testset "mpdata_step! 1D — mass conservation, Periodic, single step" begin
        nz = 20
        ϕ   = rand(nz)
        GCz = fill(0.3, nz + 1)
        tmp = make1d_tmp(nz)
        s   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=2)
        ϕ0  = sum(ϕ)
        mpdata_step!(ϕ, GCz, tmp, s)
        @test sum(ϕ) ≈ ϕ0 atol=1e-12
    end

    @testset "mpdata_step! 1D — mass conservation, Periodic, 20 steps" begin
        nz = 20
        ϕ   = rand(nz)
        GCz = fill(0.25, nz + 1)
        tmp = make1d_tmp(nz)
        s   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=2)
        ϕ0  = sum(ϕ)
        for _ in 1:20; mpdata_step!(ϕ, GCz, tmp, s); end
        @test sum(ϕ) ≈ ϕ0 atol=1e-10
    end

    @testset "mpdata_step! 1D — uniform field unchanged" begin
        nz = 12
        ϕ   = fill(4.0, nz)
        GCz = fill(0.3, nz + 1)
        tmp = make1d_tmp(nz)
        s   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=2)
        mpdata_step!(ϕ, GCz, tmp, s)
        @test all(ϕ .≈ 4.0)
    end

    @testset "mpdata_step! 1D — Gaussian advects correct distance, Periodic" begin
        nz = 60
        x0 = 20; σ = 5.0
        x   = collect(1.0:nz)
        ϕ   = exp.(-(x .- x0).^2 ./ (2σ^2))
        v   = 1.0; dt = 0.5; dz = 1.0
        GCz = fill(v * dt / dz, nz + 1)   # CFL = 0.5
        tmp = make1d_tmp(nz)
        s   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=2)
        nsteps = 10
        for _ in 1:nsteps; mpdata_step!(ϕ, GCz, tmp, s); end
        expected_shift = round(Int, v * dt * nsteps / dz)
        peak_idx = argmax(ϕ)
        @test abs(peak_idx - (x0 + expected_shift)) ≤ 2
    end

    @testset "mpdata_step! 1D — positive field stays non-negative" begin
        nz = 20
        ϕ   = abs.(randn(nz))
        GCz = fill(0.3, nz + 1)
        tmp = make1d_tmp(nz)
        s   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=2)
        for _ in 1:10; mpdata_step!(ϕ, GCz, tmp, s); end
        @test all(ϕ .≥ -1e-12)
    end

    @testset "mpdata_step! 1D — nonoscillatory: step function stays in [0,1]" begin
        nz = 30
        ϕ   = [k ≤ 15 ? 0.0 : 1.0 for k in 1:nz]
        GCz = fill(0.4, nz + 1)

        ϕ_noo = copy(ϕ)
        tmp   = make1d_tmp(nz)
        s_noo = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                   n_corr=2, nonoscillatory=true)
        for _ in 1:8; mpdata_step!(ϕ_noo, GCz, tmp, s_noo); end
        @test all(-1e-12 .≤ ϕ_noo .≤ 1.0 + 1e-12)
    end

    @testset "mpdata_step! 1D — nonoscillatory reduces overshoot vs standard" begin
        nz = 30
        ϕ   = [k ≤ 15 ? 0.0 : 1.0 for k in 1:nz]
        GCz = fill(0.4, nz + 1)

        ϕ_std = copy(ϕ)
        tmp1  = make1d_tmp(nz)
        s_std = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                   n_corr=2, nonoscillatory=false)
        for _ in 1:8; mpdata_step!(ϕ_std, GCz, tmp1, s_std); end

        ϕ_noo = copy(ϕ)
        tmp2  = make1d_tmp(nz)
        s_noo = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                   n_corr=2, nonoscillatory=true)
        for _ in 1:8; mpdata_step!(ϕ_noo, GCz, tmp2, s_noo); end

        # Nonoscillatory must not overshoot
        @test maximum(ϕ_noo) ≤ 1.0 + 1e-12
        @test minimum(ϕ_noo) ≥ 0.0 - 1e-12
        # Standard MPDATA on a step function at CFL=0.4 typically overshoots
        over_std = maximum(ϕ_std) > 1.0 + 1e-10 || minimum(ϕ_std) < -1e-10
        over_noo = maximum(ϕ_noo) > 1.0 + 1e-10 || minimum(ϕ_noo) < -1e-10
        # If standard overshoots, nonoscillatory must not
        if over_std
            @test !over_noo
        end
    end

    @testset "mpdata_step! 1D — n_corr=2 more accurate than n_corr=1 for Gaussian" begin
        nz  = 50; x0 = 25; σ = 5.0
        x   = collect(1.0:nz)
        ϕ0  = exp.(-(x .- x0).^2 ./ (2σ^2))
        GCz = fill(0.3, nz + 1)
        nsteps = 20
        shift  = round(Int, 0.3 * nsteps)
        ϕ_exact = circshift(ϕ0, shift)   # exact solution on periodic grid

        ϕ1  = copy(ϕ0)
        tmp1 = make1d_tmp(nz)
        s1  = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=1)
        for _ in 1:nsteps; mpdata_step!(ϕ1, GCz, tmp1, s1); end

        ϕ2  = copy(ϕ0)
        tmp2 = make1d_tmp(nz)
        s2  = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=2)
        for _ in 1:nsteps; mpdata_step!(ϕ2, GCz, tmp2, s2); end

        err1 = sum((ϕ1 .- ϕ_exact).^2)
        err2 = sum((ϕ2 .- ϕ_exact).^2)
        @test err2 < err1
    end

    @testset "mpdata_step! 1D — infinite_gauge: both variants conserve mass" begin
        nz  = 16
        ϕ   = rand(nz)
        GCz = fill(0.3, nz + 1)
        ϕ0  = sum(ϕ)

        ϕ_n  = copy(ϕ)
        tmp_n = make1d_tmp(nz)
        s_n  = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                  n_corr=2, infinite_gauge=false)
        mpdata_step!(ϕ_n, GCz, tmp_n, s_n)

        ϕ_g  = copy(ϕ)
        tmp_g = make1d_tmp(nz)
        s_g  = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                  n_corr=2, infinite_gauge=true)
        mpdata_step!(ϕ_g, GCz, tmp_g, s_g)

        @test sum(ϕ_n) ≈ ϕ0 atol=1e-12
        @test sum(ϕ_g) ≈ ϕ0 atol=1e-12
        # Results differ between the two modes
        @test !(ϕ_n ≈ ϕ_g)
    end

    @testset "mpdata_step! 1D — infinite_gauge + nonoscillatory: step function" begin
        nz  = 30
        ϕ   = [k ≤ 15 ? 0.0 : 1.0 for k in 1:nz]
        GCz = fill(0.35, nz + 1)
        tmp = make1d_tmp(nz)
        s   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(),
                                 n_corr=2, infinite_gauge=true, nonoscillatory=true)
        for _ in 1:8; mpdata_step!(ϕ, GCz, tmp, s); end
        @test all(-1e-12 .≤ ϕ .≤ 1.0 + 1e-12)
        @test sum(ϕ) ≈ 15.0 atol=1e-10
    end

    @testset "mpdata_step! 1D — CFL stability check (CFL < 1)" begin
        # MPDATA is stable for CFL ∈ (0, 1). Verify no NaN/Inf after many steps.
        for cfl in (0.1, 0.5, 0.9)
            nz  = 20
            ϕ   = rand(nz)
            GCz = fill(cfl, nz + 1)
            tmp = make1d_tmp(nz)
            s   = mpdata_settings_1d(nz, vertical_boundary_condition=Periodic(), n_corr=2)
            for _ in 1:20; mpdata_step!(ϕ, GCz, tmp, s); end
            @test all(isfinite.(ϕ))
        end
    end

end
