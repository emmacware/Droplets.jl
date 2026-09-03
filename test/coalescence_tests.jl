using Test
using Droplets
using Combinatorics
using Random

FT = Float64

function empty_drops(Ns;FT=FT)
    coagsettings = coag_settings{FT}(Ns = Ns)
    coagdata = coagulation_run{FT}(Ns)
    drops = droplet_attributes{FT}(Int.(zeros(Ns)),zeros(Ns))
    return drops,coagsettings,coagdata
end


@testset "Initialization" begin
    coagsettings = coag_settings{FT}(Ns = 2^8)
    @test coagsettings.scale == coagsettings.Ns * (coagsettings.Ns - 1) / 2 / (coagsettings.Ns / 2)

    function test_init_method(init_func, coagsettings)
        init_method = init_func(coagsettings)
        @test sum(init_method.ξ) ≈ coagsettings.n0*coagsettings.ΔV rtol = 3e-2

        effective_vol = sum(init_method.X.*init_method.ξ) / (sum(init_method.ξ))
        effective_radius = volume_to_radius(effective_vol)
        @test effective_radius ≈ coagsettings.R0 rtol = 1e-1
        if init_func == init_ξ_const
            @test all(==(init_method.ξ[1]), init_method.ξ)
        end
    end

    test_init_method(init_ξ_const,coagsettings)
    test_init_method(init_logarithmic,coagsettings)
    test_init_method(init_uniform_sd,coagsettings)
    test_init_method(init_monodisperse,coagsettings)
end

@testset "Kernels" begin
    coagsettings = coag_settings{FT}(Ns = 2^8)
    ξ_const = init_ξ_const(coagsettings)

    droplet_probability = golovin(ξ_const, (1,2), coagsettings)
    @test droplet_probability == coagsettings.golovin_kernel_coeff * (ξ_const.X[1] + ξ_const.X[2])

    big_golovin_coefficient_settings = coag_settings{FT}(Ns = 2^15, golovin_kernel_coeff = 1700.0)
    small_golovin_coefficient_settings = coag_settings{FT}(Ns = 2^15, golovin_kernel_coeff = 1300.0)
    @test golovin(ξ_const, (1,2), big_golovin_coefficient_settings) > golovin(ξ_const, (1,2), small_golovin_coefficient_settings)

    #same-size drops cannot collide
    mono = init_monodisperse(coagsettings)
    @test hydrodynamic(mono, (1,2), coagsettings) == 0.0
    # set test for different sized drops
    # @test hydrodynamic(uniform, (1,2), coagsettings) == ?
end

@testset "Schemes" begin
    coagsettings = coag_settings{FT}(Ns = 2^8)
    coag_data = coagulation_run{FT}(coagsettings.Ns)

    function test_schemes(scheme, coag_data, coagsettings)
        Random.seed!()
        ξ_const_droplets = init_ξ_const(coagsettings)
        ctimestep = deepcopy(ξ_const_droplets)
        for _ in 1:50
            coalescence_timestep!(Serial(),scheme, ctimestep,coag_data, coagsettings)
        end
        # coag only
        @test sum(ctimestep.ξ) <= sum(ξ_const_droplets.ξ)
        # mass conservation:
        @test sum(ctimestep.X.*ctimestep.ξ) ≈ sum(ξ_const_droplets.X.*ξ_const_droplets.ξ) rtol = 1e-12
        @test minimum(ctimestep.ξ) > 0
        @test minimum(ctimestep.X) > 0
    end

    test_schemes(none(),coag_data, coagsettings)
    test_schemes(Adaptive(),coag_data, coagsettings)

    #test that deficit for Adaptive is zero
    begin
        coagsettings = coag_settings{FT}(Ns = 2^8)
        deficit_drops = init_ξ_const(coagsettings)
        coag_data = coagulation_run{FT}(coagsettings.Ns)    
        for _ in 1:10
            coalescence_timestep!(Serial(),Adaptive(), deficit_drops,coag_data, coagsettings)
        end
        @test coag_data.deficit[] == 0
    end

end

@testset "Collision Updates" begin

    function multiple_collisions(two_drops, two_drops_coag_data)
        #arrange
        two_drops.ξ .= [2,1]
        two_drops.X .= 1.5e-15 * ones(2)
        X = two_drops.X .+0
        two_drops_coag_data.ϕ[1] = 0.5
        two_drops_coag_data.pαdt[1] = 2
        @test two_drops_coag_data.lowest_zero[] == false
        
        #act
        test_pairs!(Serial(),[(1,2)],two_drops,two_drops_coag_data)

        #assert
        @test two_drops_coag_data.lowest_zero[] == true
        @test two_drops.ξ == [0,1]
        @test two_drops.X[1] == two_drops.X[2] == 2*X[1]+X[2]
    end

    function test_pairs_ξjminus_γtildeξk(two_drops, two_drops_coag_data)
        #arrange
        two_drops.ξ .= [6,2]
        two_drops.X .= 1.5e-15 * ones(2)
        X = two_drops.X .+0
        two_drops_coag_data.ϕ[1] = 0.5
        two_drops_coag_data.pαdt[1] = 2
        #act
        test_pairs!(Serial(),[(1,2)],two_drops,two_drops_coag_data)

        #assert
        @test two_drops.ξ == [2,2]
        @test two_drops.X[2] == 2*X[1]+X[2]
        @test two_drops.X[1] == X[1]
    end

    # test deficit
    function test_deficit(two_drops, two_drops_coag_data)
        # arrange
        two_drops.ξ .= [8,4]
        two_drops.X .= 1.5e-15 * ones(2)
        two_drops_coag_data.ϕ[1] = 0.5
        two_drops_coag_data.pαdt[1] = 3
        two_drops_coag_data.deficit[] = 0

        # act
        test_pairs!(Serial(),[(1,2)],two_drops,two_drops_coag_data)

        # assert
        @test two_drops_coag_data.deficit[] == 4
    end

    function test_split_highest_multiplicity(four_drops)
        #arrange
        ξ = [0,1,2,4]
        R = [FT(1.5e-6), FT(2.5e-6), FT(3.5e-6), FT(4.5e-6)]
        X = radius_to_volume.(R)
        four_drops.ξ .= ξ .+0
        four_drops.X .= X .+0

        #act
        split_highest_multiplicity!(four_drops)

        @test four_drops.ξ == [2,1,2,2]
        @test four_drops.X[1] == four_drops.X[4] == X[4]
    end

    two_drops,two_drops_coag_settings,two_drops_coag_data = empty_drops(2)
    four_drops,four_drops_coag_settings,four_drops_coag_data = empty_drops(4)

    multiple_collisions(two_drops,two_drops_coag_data)
    test_pairs_ξjminus_γtildeξk(two_drops,two_drops_coag_data)
    test_pairs_ξjminus_γtildeξk(two_drops,two_drops_coag_data)
    test_deficit(two_drops,two_drops_coag_data)
    test_split_highest_multiplicity(four_drops)

end

@testset "Probabilities" begin

    # pair_Ps
    begin
        #arrange
        L = [(1,2),(3,4)]
        ξ = [2,3,4,5]
        R = [1.5e-6,1.5e-6,1.5e-6,1.5e-6]
        X = radius_to_volume.(R)
        drops = droplet_attributes{FT}(ξ.+0,X.+0)
        coag_data = coagulation_run{FT}(4)
        coagsettings = coag_settings{FT}(Ns = 4)

        #act
        compute_pαdt!(L,drops,coag_data,golovin,coagsettings)

        #assert
        @test coag_data.pαdt[1] == 3*coagsettings.golovin_kernel_coeff*(X[1]+X[2])*(coagsettings.scale * coagsettings.Δt / coagsettings.ΔV)
        @test coag_data.pαdt[2] == 5*coagsettings.golovin_kernel_coeff*(X[3]+X[4])*(coagsettings.scale * coagsettings.Δt / coagsettings.ΔV)
    end

    function test_adaptive_limits_tstep(sixdrops,coag_data,coagsettings)
        #arrange
        L = [(1,2),(3,4),(5,6)]
        sixdrops.ξ[1:6] .= [2e9,3e10,4e9,5e10,(7e13-1),7e13]
        R = [3e-6,3e-6,3e-6,3e-6,3e-5,3e-5]
        sixdrops.X .= radius_to_volume.(R)
        coag_data.ϕ[1:3] .= [0.5,1e-9,0.2]
        t_start = 100.0
        tlim_pair2 = div(sixdrops.ξ[4], sixdrops.ξ[3])/(sixdrops.ξ[4] * golovin(sixdrops, (3,4), coagsettings)*(coagsettings.scale / coagsettings.ΔV))
        tlim_pair3 = div(sixdrops.ξ[6], sixdrops.ξ[5])/(sixdrops.ξ[6] * golovin(sixdrops, (5,6), coagsettings)*(coagsettings.scale / coagsettings.ΔV))
        expected_tlims = [t_start,tlim_pair2,tlim_pair3]


        function test_tmax(sixdrops,coag_data,coagsettings,L,t_start,expected_tlims)
            t_left = Ref(t_start .+0)
            t_max = Vector{FT}(undef, length(L))
            #act
            for (α,pair) in enumerate(L)
                pair_Ps_adaptive!(α,pair, sixdrops,coag_data,t_max,t_left,golovin,coagsettings)
            end
            #assert
            @test t_max ≈ expected_tlims
        end

        function test_tleft(sixdrops,coag_data,coagsettings,L,t_start,expected_tlims)
            t_left = Ref(t_start .+0)
            #act
            adaptive_pαdt!(L,sixdrops,coag_data,t_left,golovin,coagsettings)

            #assert
            @test coag_data.pαdt[2] ≈ div(sixdrops.ξ[4], sixdrops.ξ[3])/expected_tlims[2] * expected_tlims[3]
            @test coag_data.pαdt[3] ≈ div(sixdrops.ξ[6], sixdrops.ξ[5])#tmp_pair3 * expected_tlims[3]
            @test t_left[] ≈ 100 - expected_tlims[3]
        end

        test_tmax(sixdrops,coag_data,coagsettings,L,t_start,expected_tlims)
        test_tleft(sixdrops,coag_data,coagsettings,L,t_start,expected_tlims)
    end

    six_drops,six_drops_coag_settings,six_drops_coag_data = empty_drops(6)
    test_adaptive_limits_tstep(six_drops,six_drops_coag_data,six_drops_coag_settings)

end

