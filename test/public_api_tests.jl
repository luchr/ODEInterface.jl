@testset "Public IVP wrapper API" begin
    public_names = (
        :OptionsODE,
        :RHS_CALL_INSITU,
        :OUTPUTFCN_CALL_REASON,
        :OUTPUTFCN_CALL_STEP,
        :OUTPUTFCN_WODENSE,
        :OUTPUTFCN_DENSE,
        :OUTPUTFCN_RET_CONTINUE,
        :OUTPUTFCN_RET_CONTINUE_XCHANGED,
        :dopri5,
        :dop853,
        :odex,
        :seulex,
        :radau,
        :radau5,
        :rodas,
        :ddeabm,
        :ddebdf,
    )

    @test all(name -> isdefined(ODEInterface, name), public_names)
    @test all(name -> haskey(Docs.meta(ODEInterface), Docs.Binding(ODEInterface, name)), public_names)

    # `public` only exists from Julia 1.11 on; see the declaration in
    # src/ODEInterface.jl.
    @static if VERSION ≥ v"1.11"
        @test all(name -> Base.ispublic(ODEInterface, name), public_names)
        @test !any(name -> Base.isexported(ODEInterface, name), public_names)
    end
end
