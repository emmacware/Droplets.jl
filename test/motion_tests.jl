using Test
using Droplets
using Random


FT = Float64

function make_motion_droplets(z_loc::Vector{FT}, R::Vector{FT})
    n = length(z_loc)
    droplet_attributes_1d{FT}(
        fill(100, n), radius_to_volume.(R), zeros(FT, n), copy(z_loc),
        ones(Int, n), zeros(FT, n), [1:n], collect(1:n),
    )
end

@testset "update_droplet_positions!" begin
    spatial = spatial_settings_1d{FT}(Nz=10, Z_max=FT(500.0))
    dz = spatial.z_grid_height
    Δt = FT(2.0)
    w_zero(z) = FT(0.0)

    @testset "DynON/DynON/DynON: advection + settling both apply" begin
        z0 = FT[50.0, 200.0]
        R = FT[20e-6, 20e-6]  # sizable settling velocity
        droplets = make_motion_droplets(z0, R)
        w_const(z) = FT(1.0)
        update_droplet_positions!(DynON(), DynON(), DynON(), droplets, w_const, Δt, spatial, nothing, 1)

        vt = terminal_v_X(droplets.X[1])
        expected = z0 .+ w_const.(z0) .* Δt .- vt * Δt
        @test droplets.z_loc ≈ clamp.(expected, -1, spatial.Z_max + dz)
        @test droplets.cell_id == floor.(Int, droplets.z_loc ./ dz) .+ 1
    end

    @testset "DynON/DynON/DynOFF: advection only, no settling" begin
        z0 = FT[50.0, 200.0]
        R = FT[20e-6, 20e-6]
        droplets = make_motion_droplets(z0, R)
        w_const(z) = FT(1.5)
        update_droplet_positions!(DynON(), DynON(), DynOFF(), droplets, w_const, Δt, spatial, nothing, 1)
        @test droplets.z_loc ≈ z0 .+ w_const.(z0) .* Δt
    end

    @testset "DynON/DynOFF/DynON: settling only, no advection" begin
        z0 = FT[50.0, 200.0]
        R = FT[30e-6, 10e-6]
        droplets = make_motion_droplets(z0, R)
        w_ignored(z) = FT(999.0)  # must be ignored when advection is DynOFF
        update_droplet_positions!(DynON(), DynOFF(), DynON(), droplets, w_ignored, Δt, spatial, nothing, 1)
        vt = terminal_v_X.(radius_to_volume.(R))
        @test droplets.z_loc ≈ z0 .- vt .* Δt
    end

    @testset "DynOFF motion: strict no-op regardless of advection/settling" begin
        z0 = FT[50.0, 200.0]
        droplets = make_motion_droplets(z0, FT[10e-6, 10e-6])
        z_before = copy(droplets.z_loc)
        cid_before = copy(droplets.cell_id)
        update_droplet_positions!(DynOFF(), DynON(), DynON(), droplets, z -> FT(5.0), Δt, spatial, nothing, 1)
        @test droplets.z_loc == z_before
        @test droplets.cell_id == cid_before
    end

    @testset "DynON/DynOFF/DynOFF: no-op (neither advection nor settling)" begin
        z0 = FT[50.0, 200.0]
        droplets = make_motion_droplets(z0, FT[10e-6, 10e-6])
        z_before = copy(droplets.z_loc)
        update_droplet_positions!(DynON(), DynOFF(), DynOFF(), droplets, z -> FT(5.0), Δt, spatial, nothing, 1)
        @test droplets.z_loc == z_before
    end

    @testset "position is clamped to the domain, never goes negative past the floor" begin
        z0 = FT[0.5]
        R = FT[500e-6]  # large droplet, fast settling, would go well below zero
        droplets = make_motion_droplets(z0, R)
        update_droplet_positions!(DynON(), DynON(), DynON(), droplets, w_zero, FT(100.0), spatial, nothing, 1)
        @test droplets.z_loc[1] == -1.0

        z0_top = FT[499.0]
        droplets_top = make_motion_droplets(z0_top, FT[1e-6])
        update_droplet_positions!(DynON(), DynON(), DynON(), droplets_top, z -> FT(100.0), FT(100.0), spatial, nothing, 1)
        @test droplets_top.z_loc[1] == spatial.Z_max + dz
    end
end

@testset "fill_grid_ranges!" begin
    n = 10
    cell_id = [3, 1, 2, 1, 3, 2, 2, 1, 3, 1]
    droplets = droplet_attributes_1d{FT}(
        fill(1, n), fill(FT(1e-15), n), zeros(FT, n), zeros(FT, n),
        copy(cell_id), zeros(FT, n), fill(1:0, 3), collect(1:n),
    )

    Droplets.fill_grid_ranges!(droplets)

    # I sorted by cell_id, grid_range[g] indexes exactly the droplets with cell_id==g
    for g in 1:3
        idxs = droplets.I[droplets.grid_range[g]]
        @test all(i -> cell_id[i] == g, idxs)
        @test sort(idxs) == sort(findall(==(g), cell_id))
    end
    @test length(droplets.I) == n
    @test issorted(droplets.cell_id[droplets.I])
end

@testset "nondivergent_vector_field_2d" begin
    grid = (16, 16)
    dx = dy = 1.0 / 16
    dt = 0.01
    omega = 0.05
    sf = make_stream_function(grid, dx, dy, omega, 0.5, 0.5)
    u, v, GCx, GCy = nondivergent_vector_field_2d(dx, dy, grid, dt, sf)

    @test all(isfinite, u)
    @test all(isfinite, v)
    @test all(abs.(GCx) .< 1.0)
    @test all(abs.(GCy) .< 1.0)

    # a solid body rotation stream function should give a genuinely nonzero flow
    @test any(!=(0.0), u) || any(!=(0.0), v)

    # CFL violation must actually throw (dt way too large for this grid)
    @test_throws AssertionError nondivergent_vector_field_2d(dx, dy, grid, 100.0, sf)
end

@testset "turbulent droplet diffusion schemes: finite, bounded, no crashes" begin
    Random.seed!(42)
    nz = 12
    dz = FT(50.0)
    droplets = make_motion_droplets(collect(range(FT(10.0), FT(590.0), length=24)), fill(FT(10e-6), 24))
    grid_range = map(1:nz) do g
        idx = findall(z -> floor(Int, z / dz) + 1 == g, droplets.z_loc)
        isempty(idx) ? (1:0) : minimum(idx):maximum(idx)
    end
    droplets = droplet_attributes_1d{FT}(droplets.ξ, droplets.X, droplets.dry_r3, droplets.z_loc,
        floor.(Int, droplets.z_loc ./ dz) .+ 1, droplets.w_prime, grid_range, droplets.I)

    grid = (nz=nz, dz=dz, centers_z=collect(dz/2:dz:(nz-0.5)*dz),
        states=(qv=fill(FT(0.01), nz), e=fill(FT(0.1), nz)))
    tke = tke_settings{FT}()
    dt = FT(1.0)
    l = fill(FT(30.0), nz)

    d1 = deepcopy(droplets)
    turbulent_droplet_diffusion!(DynON(), l, d1, grid, tke, dt)
    @test all(isfinite, d1.z_loc)
    @test all(isfinite, d1.w_prime)

    d2 = deepcopy(droplets)
    turbulent_droplet_diffusion_wellmixed!(DynON(), l, d2, grid, tke, dt)
    @test all(isfinite, d2.z_loc)
    @test all(0 .<= d2.z_loc .<= nz * dz)

    d3 = deepcopy(droplets)
    K_centers = fill(FT(5.0), nz)
    turbulent_droplet_diffusion_visser!(DynON(), d3, grid, K_centers, dt)
    @test all(isfinite, d3.z_loc)
    @test all(0 .<= d3.z_loc .<= nz * dz)

    # DynOFF variants dont deploy
    d4 = deepcopy(droplets)
    turbulent_droplet_diffusion!(DynOFF(), l, d4, grid, tke, dt)
    @test d4.z_loc == droplets.z_loc
    d5 = deepcopy(droplets)
    turbulent_droplet_diffusion_wellmixed!(DynOFF(), l, d5, grid, tke, dt)
    @test d5.z_loc == droplets.z_loc
    d6 = deepcopy(droplets)
    turbulent_droplet_diffusion_visser!(DynOFF(), d6, grid, K_centers, dt)
    @test d6.z_loc == droplets.z_loc
end

@testset "stochastic_jump_diffusion!: only cloud droplets (<40µm) jump, at most one cell" begin
    Random.seed!(7)
    nz = 5
    dz = FT(50.0)
    R = FT[5e-6, 60e-6, 5e-6, 60e-6]  # alternate cloud/large droplets
    z0 = fill(FT(2.5 * dz), 4)  # all start in cell 3
    droplets = droplet_attributes_1d{FT}(
        fill(1, 4), radius_to_volume.(R), zeros(FT, 4), copy(z0),
        fill(3, 4), zeros(FT, 4), [1:0, 1:0, 1:4, 1:0, 1:0], collect(1:4),
    )
    grid = nothing
    K_centers = fill(FT(20.0), nz)

    Droplets.stochastic_jump_diffusion!(DynON(), grid, droplets, K_centers, FT(1.0), dz, nz)

    # large (>=40µm) droplets never move
    @test droplets.cell_id[2] == 3
    @test droplets.cell_id[4] == 3
    # small droplets move at most one cell away
    @test all(abs(droplets.cell_id[i] - 3) <= 1 for i in (1, 3))
end
