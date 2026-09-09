using Test
using Droplets
using Random

# The spatially-resolved SDM path used by the actual SCM

FT = Float64

"""Build droplets_1d spread across `counts_per_cell` cells, all starting at the same
volume/multiplicity so mass conservation math is easy to check by hand."""
function make_spatial_droplets(counts_per_cell::Vector{Int}; R=FT(1.5e-6), xi=100)
    Ns = sum(counts_per_cell)
    ξ = fill(xi, Ns)
    X = fill(radius_to_volume(R), Ns)
    dry_r3 = fill(FT(0.0), Ns)
    z_loc = fill(FT(0.0), Ns)
    cell_id = Vector{Int}(undef, Ns)
    grid_range = Vector{UnitRange{Int}}(undef, length(counts_per_cell))
    offset = 0
    for (g, n) in enumerate(counts_per_cell)
        cell_id[offset+1:offset+n] .= g
        grid_range[g] = (offset+1):(offset+n)
        offset += n
    end
    w_prime = zeros(FT, Ns)
    I = collect(1:Ns)
    droplets = droplet_attributes_1d{FT}(ξ, X, dry_r3, z_loc, cell_id, w_prime, grid_range, I)
    return droplets
end

make_spatial_coagsettings(Ns; ΔV=FT(1e6)) =
    coag_settings{FT}(Ns=Ns, ΔV=ΔV, Δt=FT(1.0), golovin_kernel_coeff=FT(1500.0), kernel=golovin)

@testset "rebuild_coalescence_pairing!" begin
    counts = [4, 6, 0, 3]  # includes an empty cell
    droplets = make_spatial_droplets(counts)
    Ns = sum(counts)
    coagsettings = make_spatial_coagsettings(Ns)
    coag_data = coagulation_run_spatial{FT}(length(counts), Ns, droplets)

    Droplets.rebuild_coalescence_pairing!(droplets, coag_data)

    # number of pairs == sum of floor(N_c/2) over nonempty cells
    @test length(coag_data.pair_starts) == div(4, 2) + div(6, 2) + div(3, 2)

    # per-cell scale matches the box formula Ns_c*(Ns_c-1)/2/div(Ns_c,2)
    @test coag_data.scale[1] ≈ 4 * 3 / 2 / 2
    @test coag_data.scale[2] ≈ 6 * 5 / 2 / 3
    @test coag_data.scale[4] ≈ 3 * 2 / 2 / 1
    @test coag_data.scale[3] == 0.0  # empty cell untouched

    # every pair start index falls within its own cell's grid_range
    for g in eachindex(droplets.grid_range)
        r = droplets.grid_range[g]
        isempty(r) && continue
        starts_in_cell = filter(i -> i in r, coag_data.pair_starts)
        @test all(i -> i + 1 in r, starts_in_cell)
    end
end

@testset "coalescence_reshuffle_and_step!: per-cell mass conservation, Serial" begin
    Random.seed!(1)
    counts = [10, 14, 8]
    droplets = make_spatial_droplets(counts)
    Ns = sum(counts)
    coagsettings = make_spatial_coagsettings(Ns)
    coag_data = coagulation_run_spatial{FT}(length(counts), Ns, droplets)
    Droplets.rebuild_coalescence_pairing!(droplets, coag_data)

    mass_before = [sum(droplets.X[droplets.I[r]] .* droplets.ξ[droplets.I[r]]) for r in droplets.grid_range]

    for _ in 1:30
        Droplets.coalescence_reshuffle_and_step!(Serial(), droplets, coag_data, coagsettings.kernel, coagsettings)
    end

    mass_after = [sum(droplets.X[droplets.I[r]] .* droplets.ξ[droplets.I[r]]) for r in droplets.grid_range]

    # mass conserved independently within each cell (no cross-cell transport happens
    # during coalescence, this would fail if pairing ever crossed a grid_range boundary)
    for g in eachindex(counts)
        @test mass_after[g] ≈ mass_before[g] rtol=1e-10
    end
    @test all(droplets.ξ .>= 0)
    @test all(droplets.X .> 0)
    @test all(isfinite, droplets.dry_r3)
end

@testset "coalescence_reshuffle_and_step!: per-cell mass conservation, Parallel" begin
    Random.seed!(2)
    counts = [12, 20, 6, 16]
    droplets = make_spatial_droplets(counts)
    Ns = sum(counts)
    coagsettings = make_spatial_coagsettings(Ns)
    coag_data = coagulation_run_spatial{FT}(length(counts), Ns, droplets)
    Droplets.rebuild_coalescence_pairing!(droplets, coag_data)

    mass_before = [sum(droplets.X[droplets.I[r]] .* droplets.ξ[droplets.I[r]]) for r in droplets.grid_range]

    for _ in 1:30
        Droplets.coalescence_reshuffle_and_step!(Parallel(), droplets, coag_data, coagsettings.kernel, coagsettings)
    end

    mass_after = [sum(droplets.X[droplets.I[r]] .* droplets.ξ[droplets.I[r]]) for r in droplets.grid_range]
    for g in eachindex(counts)
        @test mass_after[g] ≈ mass_before[g] rtol=1e-10
    end
    @test all(droplets.ξ .>= 0)
end

@testset "coalescence_reshuffle_and_step!: cells never exchange droplets" begin
    # tag dry_r3 with a per-cell marker; since sdm_update! combines dry_r3 the same
    # way it combines X (weighted by γ), any post-coalescence dry_r3 that isn't an
    # integer multiple of its own cell's marker proves cross-cell mixing occurred.
    Random.seed!(3)
    counts = [8, 8]
    droplets = make_spatial_droplets(counts)
    droplets.dry_r3[droplets.grid_range[1]] .= 1.0
    droplets.dry_r3[droplets.grid_range[2]] .= 1000.0
    Ns = sum(counts)
    coagsettings = make_spatial_coagsettings(Ns)
    coag_data = coagulation_run_spatial{FT}(length(counts), Ns, droplets)
    Droplets.rebuild_coalescence_pairing!(droplets, coag_data)

    for _ in 1:20
        Droplets.coalescence_reshuffle_and_step!(Serial(), droplets, coag_data, coagsettings.kernel, coagsettings)
    end

    @test all(x -> abs(x - round(x)) < 1e-8, droplets.dry_r3[droplets.grid_range[1]])
    @test all(x -> x >= 999.0, droplets.dry_r3[droplets.grid_range[2]])
end

@testset "coalescence_timestep!(DynON/DynOFF, ...): full 1D wrapper" begin
    Random.seed!(4)
    counts = [10, 10]
    droplets = make_spatial_droplets(counts)
    Ns = sum(counts)
    coagsettings = make_spatial_coagsettings(Ns)
    coag_data = coagulation_run_spatial{FT}(length(counts), Ns, droplets)

    mass_before = sum(droplets.X .* droplets.ξ)
    ξ_sum_before = sum(droplets.ξ)

    # DynOFF is a strict no-op
    droplets_off = deepcopy(droplets)
    coag_data_off = coagulation_run_spatial{FT}(length(counts), Ns, droplets_off)
    coalescence_timestep!(DynOFF(), Serial(), none(), droplets_off, coag_data_off, coagsettings)
    @test droplets_off.ξ == droplets.ξ
    @test droplets_off.X == droplets.X

    # DynON actually pairs and updates droplets, conserving total mass
    for _ in 1:20
        coalescence_timestep!(DynON(), Serial(), none(), droplets, coag_data, coagsettings)
    end
    @test sum(droplets.X .* droplets.ξ) ≈ mass_before rtol=1e-10
    @test sum(droplets.ξ) <= ξ_sum_before  # total multiplicity is non-increasing under coalescence

    # n_substeps variant: same wrapper, run n_substeps reshuffle passes per call
    droplets2 = make_spatial_droplets(counts)
    coag_data2 = coagulation_run_spatial{FT}(length(counts), Ns, droplets2)
    mass_before2 = sum(droplets2.X .* droplets2.ξ)
    coalescence_timestep!(DynON(), Serial(), none(), droplets2, coag_data2, coagsettings, 5)
    @test sum(droplets2.X .* droplets2.ξ) ≈ mass_before2 rtol=1e-10

    droplets3 = deepcopy(droplets2)
    coag_data3 = coagulation_run_spatial{FT}(length(counts), Ns, droplets3)
    coalescence_timestep!(DynOFF(), Serial(), none(), droplets3, coag_data3, coagsettings, 5)
    @test droplets3.ξ == droplets2.ξ
end

@testset "reduce_collision_rate!: accumulated collision mass matches per-pair record" begin
    Random.seed!(5)
    counts = [16, 16]
    droplets = make_spatial_droplets(counts)
    Ns = sum(counts)
    coagsettings = make_spatial_coagsettings(Ns)
    coag_data = coagulation_run_spatial{FT}(length(counts), Ns, droplets)
    Droplets.rebuild_coalescence_pairing!(droplets, coag_data)
    Droplets.coalescence_reshuffle_and_step!(Serial(), droplets, coag_data, coagsettings.kernel, coagsettings)

    # collision_rate[cell] must equal the sum of collision_rate_pair over that cell's pairs
    for g in eachindex(counts)
        pairs_in_cell = filter(i -> i in droplets.grid_range[g], coag_data.pair_starts)
        expected = sum(coag_data.collision_rate_pair[i] for i in pairs_in_cell; init=zero(FT))
        @test coag_data.collision_rate[g] ≈ expected
    end
    @test all(coag_data.collision_rate .>= 0)
end

@testset "sdm_update! (spatial, per-pair ϕ): single event and lowest_zero flag" begin
    counts = [2]
    droplets = make_spatial_droplets(counts; xi=1)
    droplets.ξ .= [6, 2]
    droplets.X .= 1.5e-15 * ones(2)
    X0 = copy(droplets.X)
    Ns = 2
    coag_data = coagulation_run_spatial{FT}(1, Ns, droplets)

    # deterministic single collision event: pα=2.0, ϕ=0.5 -> γ=2
    sdm_update!((1, 2), FT(2.0), FT(0.5), droplets, coag_data, 1)

    @test droplets.ξ == [2, 2]
    @test droplets.X[2] ≈ 2 * X0[1] + X0[2]
    @test droplets.X[1] ≈ X0[1]
    @test coag_data.collision_rate_pair[1] ≈ 4.0  # γ*ξk = 2*2

    # lowest_zero flag set when the pair collapses and the split leaves ξ==0
    droplets2 = make_spatial_droplets([2]; xi=1)
    droplets2.ξ .= [1, 1]
    droplets2.X .= 1.5e-15 * ones(2)
    coag_data2 = coagulation_run_spatial{FT}(1, 2, droplets2)
    sdm_update!((1, 2), FT(1.0), FT(0.0), droplets2, coag_data2, 1)
    @test droplets2.ξ == [0, 1]
    @test coag_data2.lowest_zero[] == true
end

@testset "split_highest_multiplicity! (droplet_attributes_1d): copies all per-droplet fields" begin
    ξ = [0, 1, 2, 4]
    R = [FT(1.5e-6), FT(2.5e-6), FT(3.5e-6), FT(4.5e-6)]
    X = radius_to_volume.(R)
    droplets = droplet_attributes_1d{FT}(
        copy(ξ), copy(X), FT[10.0, 20.0, 30.0, 40.0], FT[100.0, 200.0, 300.0, 400.0],
        [1, 1, 1, 1], zeros(FT, 4), [1:4], collect(1:4),
    )

    split_highest_multiplicity!(droplets)

    @test droplets.ξ == [2, 1, 2, 2]
    @test droplets.X[1] == droplets.X[4] == X[4]
    # dry_r3/cell_id/z_loc/w_prime must follow X when the donor slot is overwritten
    @test droplets.dry_r3[1] == droplets.dry_r3[4]
    @test droplets.z_loc[1] == droplets.z_loc[4]
end
