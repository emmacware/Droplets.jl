using Droplets
using Test


FT = Float64

# GCx is (nx+1, ny) sized (face-centered in x), GCy is (nx, ny+1) (face-centered in y),
# matching the 1D convention where GCz has nz+1 faces for nz cells.
make2d_tmp(nx, ny) = mpdata_tmp(zeros(nx, ny), zeros(nx + 1, ny), zeros(nx, ny + 1))

@testset "MPDATA 2D" begin

    @testset "flux_horizontal!/flux_vertical!: mass conservation, Periodic" begin
        nx, ny = 6, 5
        ϕ = rand(nx, ny)
        ϕ0 = sum(ϕ)
        GCx = fill(0.3, nx + 1, ny)
        GCy = fill(0.2, nx, ny + 1)
        ϕ_tmp = zeros(nx, ny)
        flux_horizontal!(ϕ, ϕ_tmp, GCx, bc=Periodic())
        flux_vertical!(ϕ, ϕ_tmp, GCy, bc=Periodic())
        ϕ .+= ϕ_tmp
        @test sum(ϕ) ≈ ϕ0 atol=1e-10
    end

    @testset "donor_cell_pass! 2D — mass conservation, Periodic" begin
        nx, ny = 8, 6
        ϕ = rand(nx, ny)
        ϕ0 = sum(ϕ)
        GCx = fill(0.25, nx + 1, ny)
        GCy = fill(0.15, nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        tmp.GCx_step .= GCx
        tmp.GCy_step .= GCy
        donor_cell_pass!(ϕ, tmp, hbc=Periodic(), vbc=Periodic())
        @test sum(ϕ) ≈ ϕ0 atol=1e-10
    end

    @testset "donor_cell_pass! 2D — uniform field unchanged, Periodic" begin
        nx, ny = 6, 6
        ϕ = fill(3.0, nx, ny)
        GCx = fill(0.3, nx + 1, ny)
        GCy = fill(0.2, nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        tmp.GCx_step .= GCx
        tmp.GCy_step .= GCy
        donor_cell_pass!(ϕ, tmp, hbc=Periodic(), vbc=Periodic())
        @test all(ϕ .≈ 3.0)
    end

    @testset "donor_cell_pass! 2D — pure x-advection matches 1D donor cell along each row" begin
        nx, ny = 10, 3
        ϕ = rand(nx, ny)
        ϕ0 = copy(ϕ)
        GCx = fill(0.3, nx + 1, ny)
        GCy = zeros(nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        tmp.GCx_step .= GCx
        tmp.GCy_step .= GCy
        donor_cell_pass!(ϕ, tmp, hbc=Periodic(), vbc=Periodic())

        # with zero y-velocity, each column of ny should evolve independently and
        # identically to the corresponding 1D donor-cell pass along x
        for j in 1:ny
            row = ϕ0[:, j]
            tmp1d = mpdata_tmp_1d(zeros(nx), zeros(nx + 1))
            tmp1d.GCz_step .= 0.3
            donor_cell_pass!(row, tmp1d, Periodic())
            @test ϕ[:, j] ≈ row atol=1e-12
        end
    end

    @testset "find_extrema! 2D — interior cell neighbors (4-point stencil)" begin
        nx, ny = 4, 4
        ϕ = FT[1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
        tmp = make2d_tmp(nx, ny)
        find_extrema!(ϕ, tmp, hbc=Periodic(), vbc=Periodic())
        # cell (2,2)=6 has neighbors (1,2)=2,(3,2)=10,(2,1)=5,(2,3)=7
        @test tmp.minmax.localmax[2, 2] == max(6, 2, 10, 5, 7)
        @test tmp.minmax.localmin[2, 2] == min(6, 2, 10, 5, 7)
    end

    @testset "mpdata_step! 2D — mass conservation, Periodic, single step" begin
        nx, ny = 12, 10
        ϕ = rand(nx, ny)
        ϕ0 = sum(ϕ)
        GCx = fill(0.2, nx + 1, ny)
        GCy = fill(0.15, nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        settings = mpdata_settings((nx, ny), horizontal_boundary_condition=Periodic(),
            vertical_boundary_condition=Periodic(), n_corr=2)
        mpdata_step!(ϕ, GCx, GCy, tmp, settings)
        @test sum(ϕ) ≈ ϕ0 atol=1e-8
    end

    @testset "mpdata_step! 2D — mass conservation over multiple steps" begin
        nx, ny = 10, 8
        ϕ = rand(nx, ny)
        ϕ0 = sum(ϕ)
        GCx = fill(0.25, nx + 1, ny)
        GCy = fill(0.1, nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        settings = mpdata_settings((nx, ny), horizontal_boundary_condition=Periodic(),
            vertical_boundary_condition=Periodic(), n_corr=2)
        for _ in 1:15
            mpdata_step!(ϕ, GCx, GCy, tmp, settings)
        end
        @test sum(ϕ) ≈ ϕ0 atol=1e-6
    end

    @testset "mpdata_step! 2D — uniform field stays uniform" begin
        nx, ny = 8, 8
        ϕ = fill(2.0, nx, ny)
        GCx = fill(0.3, nx + 1, ny)
        GCy = fill(0.2, nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        settings = mpdata_settings((nx, ny), horizontal_boundary_condition=Periodic(),
            vertical_boundary_condition=Periodic(), n_corr=2)
        mpdata_step!(ϕ, GCx, GCy, tmp, settings)
        @test all(ϕ .≈ 2.0)
    end

    @testset "mpdata_step! 2D — 2D Gaussian blob advects the correct distance" begin
        nx, ny = 30, 30
        x0, y0, σ = 10, 15, 3.0
        ϕ = [exp(-((i - x0)^2 + (j - y0)^2) / (2σ^2)) for i in 1:nx, j in 1:ny]
        vx, vy, dt = 1.0, 0.0, 0.5
        GCx = fill(vx * dt, nx + 1, ny)
        GCy = fill(vy * dt, nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        settings = mpdata_settings((nx, ny), horizontal_boundary_condition=Periodic(),
            vertical_boundary_condition=Periodic(), n_corr=2)
        nsteps = 8
        for _ in 1:nsteps
            mpdata_step!(ϕ, GCx, GCy, tmp, settings)
        end
        expected_shift = round(Int, vx * dt * nsteps)
        peak = argmax(ϕ)
        @test abs(peak[1] - (x0 + expected_shift)) <= 2
        @test peak[2] == y0
    end

    @testset "mpdata_step! 2D — nonoscillatory keeps a step function within bounds" begin
        nx, ny = 20, 6
        ϕ = [i <= 10 ? 0.0 : 1.0 for i in 1:nx, j in 1:ny]
        GCx = fill(0.35, nx + 1, ny)
        GCy = zeros(nx, ny + 1)
        tmp = make2d_tmp(nx, ny)
        settings = mpdata_settings((nx, ny), horizontal_boundary_condition=Periodic(),
            vertical_boundary_condition=Periodic(), n_corr=2, nonoscillatory=true)
        for _ in 1:6
            mpdata_step!(ϕ, GCx, GCy, tmp, settings)
        end
        @test all(-1e-10 .<= ϕ .<= 1.0 + 1e-10)
    end

    @testset "mpdata_step! 2D — CFL stability check" begin
        for cfl in (0.1, 0.4, 0.8)
            nx, ny = 10, 10
            ϕ = rand(nx, ny)
            GCx = fill(cfl, nx + 1, ny)
            GCy = fill(cfl / 2, nx, ny + 1)
            tmp = make2d_tmp(nx, ny)
            settings = mpdata_settings((nx, ny), horizontal_boundary_condition=Periodic(),
                vertical_boundary_condition=Periodic(), n_corr=2)
            for _ in 1:15
                mpdata_step!(ϕ, GCx, GCy, tmp, settings)
            end
            @test all(isfinite.(ϕ))
        end
    end
end
