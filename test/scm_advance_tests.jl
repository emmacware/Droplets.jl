using Test
using Droplets

FT = Float64

function make_bookkeeping_droplets(ξ, X, cell_id, z_loc, nz)
    n = length(ξ)
    dropidx = collect(1:n)
    grid_range = map(1:nz) do g
        idx = findall(==(g), cell_id)
        isempty(idx) ? (1:0) : minimum(idx):maximum(idx)
    end
    # cell_id/dropidx must agree with grid_range's contiguous-block assumption:
    # sort dropidx by cell_id first, like fill_grid_ranges! would leave it
    sort!(dropidx, by=i -> cell_id[i])
    grid_range = map(1:nz) do g
        s = findfirst(i -> cell_id[i] == g, dropidx)
        s === nothing ? (1:0) : s:findlast(i -> cell_id[i] == g, dropidx)
    end
    droplet_attributes_1d{FT}(copy(ξ), copy(X), zeros(FT, n), copy(z_loc), copy(cell_id),
        zeros(FT, n), grid_range, dropidx)
end

@testset "recycle_precipitation!" begin
    nz = 30
    dz = FT(30.0)
    cst = Constants{FT}()
    diagset = diagnostic_settings{FT}()
    spatial = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz * dz), area_per_grid=FT(10.0), dt_output=FT(60.0))
    grid = (states=(qv=fill(FT(0.01), nz),), centers_z=collect(dz/2:dz:(nz - 0.5) * dz),
        output=(surface_precipitation=zeros(FT, 1),))
    grid.states.qv[26:end] .= FT(0.003)

    @testset "raindrop below the domain: mass recorded, droplet recycled as a clone" begin
    
        donor_ξ = 8
        ξ = [10, donor_ξ]
        X = [radius_to_volume(FT(50e-6)), radius_to_volume(FT(15e-6))]  # precip raindrop, donor cloud drop
        cell_id = [20, 20]
        z_loc = [FT(-0.5), FT(20 * dz - 1.0)]  # first droplet has fallen below the domain
        droplets = make_bookkeeping_droplets(ξ, X, cell_id, z_loc, nz)

        mass_before = sum(droplets.X .* droplets.ξ)
        precip_X = X[1]; precip_ξ = ξ[1]

        recycle_precipitation!(DynON(), droplets, grid, spatial, diagset, nothing, cst, 1)

        expected_precip_mass = precip_X * precip_ξ * cst.ρl
        @test grid.output.surface_precipitation[1] ≈
              expected_precip_mass / (spatial.area_per_grid * spatial.dt_output) rtol=1e-10

        # ξ conservation: the precipitated slot and its donor together keep donor_ξ total
        @test droplets.ξ[1] + droplets.ξ[2] == donor_ξ
        @test droplets.ξ[1] == floor(Int, donor_ξ / 2)
        @test droplets.X[1] == droplets.X[2]  # cloned from the donor's (unsplit) volume

        # total system mass decreased by exactly the recorded precipitated mass
        mass_after = sum(droplets.X .* droplets.ξ)
        @test mass_before - mass_after ≈ expected_precip_mass / cst.ρl rtol=1e-8
    end

    @testset "sub-cutoff droplet below the domain: recycled in place, no mass loss (not real precip)" begin
        aero_X = diagset.aerosol_cloud_cuttoff / 2  # comfortably below the cutoff
        ξ = [aero_X > 0 ? 5 : 5]
        X = [aero_X]
        cell_id = [1]
        z_loc = [FT(-0.2)]
        droplets = make_bookkeeping_droplets(ξ, X, cell_id, z_loc, nz)
        X_before = copy(droplets.X); ξ_before = copy(droplets.ξ)
        grid2 = (states=(qv=fill(FT(0.01), nz),), centers_z=grid.centers_z,
            output=(surface_precipitation=zeros(FT, 1),))
        grid2.states.qv[26:end] .= FT(0.003)

        recycle_precipitation!(DynON(), droplets, grid2, spatial, diagset, nothing, cst, 1)

        @test droplets.X == X_before
        @test droplets.ξ == ξ_before
        @test 0 <= droplets.z_loc[1] < spatial.z_grid_height
        @test droplets.cell_id[1] == 1
        @test grid2.output.surface_precipitation[1] == 0.0  # no real precip recorded
    end

    @testset "DynOFF: strict no-op" begin
        droplets = make_bookkeeping_droplets([10], [radius_to_volume(FT(50e-6))], [10], [FT(-1.0)], nz)
        X_before = copy(droplets.X)
        out = (surface_precipitation=zeros(FT, 1),)
        recycle_precipitation!(DynOFF(), droplets, grid, spatial, diagset, nothing, cst, 1)
        @test droplets.X == X_before
    end
end

@testset "recycle_top_escape!" begin
    nz = 10
    dz = FT(50.0)
    cst = Constants{FT}()
    spatial = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz * dz))
    P0 = FT(95000.0); T0 = FT(288.0); qv0 = FT(6e-3)  # subsaturated near the surface
    grid = (states=(θ=[theta_from_T(T0, P0, qv0, cst)], P=[P0], qv=[qv0]),)

    kappa = FT(0.61)
    dry_r = FT(0.05e-6)
    droplets = droplet_attributes_1d{FT}([100], [radius_to_volume(FT(3e-6))], [dry_r^3],
        [spatial.Z_max + FT(1.0)], [nz], [FT(0.0)], [1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:1], [1])

    T = T_from_theta(grid.states.θ[1], grid.states.P[1], grid.states.qv[1], cst)
    S_env = min(sat(grid.states.qv[1], grid.states.P[1]) / esat(T), FT(0.95))

    recycle_top_escape!(DynON(), droplets, spatial, grid, cst)

    @test droplets.z_loc[1] == FT(0.1)
    @test droplets.cell_id[1] == 1

    # BUG: the source calls
    # `find_equilibrium_radius.(droplets, k, 0.95, T, S_env, constants)`
    # positionally passing the literal 0.95 (meant as the S_env cap, per the
    # `S_env > 0.95` clamp two lines above) as the `kappa` argument, instead of the
    # droplet's real hygroscopicity. 
    ref_droplets = droplet_attributes_1d{FT}([100], [radius_to_volume(FT(3e-6))], [dry_r^3],
        [FT(0.0)], [1], [FT(0.0)], [1:1], [1])
    find_equilibrium_radius(ref_droplets, 1, kappa, T, S_env, cst)  # the droplet's real kappa
    @test_broken droplets.X[1] ≈ ref_droplets.X[1] rtol=1e-10

    @testset "DynOFF: strict no-op" begin
        d2 = droplet_attributes_1d{FT}([100], [radius_to_volume(FT(3e-6))], [dry_r^3],
            [spatial.Z_max + FT(1.0)], [nz], [FT(0.0)], [1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:0, 1:1], [1])
        z_before = copy(d2.z_loc)
        recycle_top_escape!(DynOFF(), d2, spatial, grid, cst)
        @test d2.z_loc == z_before
    end
end

@testset "keep_layer_filled!: total multiplicity is conserved" begin
    nz = 10
    dz = FT(50.0)
    cst = Constants{FT}()
    diagset = diagnostic_settings{FT}()
    spatial = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz * dz))
    grid = (states=(qv=fill(FT(0.01), nz),), centers_z=collect(dz/2:dz:(nz - 0.5) * dz))
    grid.states.qv[9:end] .= FT(0.003)  # inversion near the top -> z_inv_idx ~ 8

    # cell 5 (in the [idx_400m, z_inv_idx] range) is starved (2 droplets); cell 2
    # (a donor candidate below idx_400m) is well-stocked
    ξ_donor = fill(20, 60)
    X_donor = radius_to_volume.(range(FT(3e-6), FT(3.5e-6), length=60))
    ξ_starved = [15, 22]
    X_starved = radius_to_volume.(FT[10e-6, 11e-6])

    ξ = vcat(ξ_donor, ξ_starved)
    X = vcat(X_donor, X_starved)
    cell_id = vcat(fill(2, 60), fill(5, 2))
    z_loc = zeros(FT, length(ξ))

    droplets = make_bookkeeping_droplets(ξ, X, cell_id, z_loc, nz)
    ξ_total_before = sum(droplets.ξ)

    keep_layer_filled!(DynON(), droplets, grid, spatial, diagset, cst; min_count=10)

    @test sum(droplets.ξ) == ξ_total_before  # borrowing/cloning only ever splits ξ, never creates or destroys it
    @test all(droplets.ξ .>= 0)
    @test all(isfinite, droplets.X)
    @test all(droplets.X .> 0)

    @testset "DynOFF: strict no-op" begin
        d2 = make_bookkeeping_droplets(ξ, X, cell_id, z_loc, nz)
        ξ_before = copy(d2.ξ)
        keep_layer_filled!(DynOFF(), d2, grid, spatial, diagset, cst; min_count=10)
        @test d2.ξ == ξ_before
    end
end

@testset "fill_grid_ranges! is idempotent when called again by these functions" begin
    # recycle_precipitation!/keep_layer_filled! both call fill_grid_ranges! internally
    # at the start
    nz = 5
    cell_id = [3, 1, 2, 1, 3]
    droplets = droplet_attributes_1d{FT}(fill(1, 5), fill(FT(1e-15), 5), zeros(FT, 5), zeros(FT, 5),
        copy(cell_id), zeros(FT, 5), fill(1:0, nz), collect(1:5))
    Droplets.fill_grid_ranges!(droplets)
    ranges1 = copy(droplets.grid_range)
    I1 = copy(droplets.I)
    Droplets.fill_grid_ranges!(droplets)
    @test droplets.grid_range == ranges1
    @test droplets.I == I1
end
