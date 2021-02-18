let # preventing global scope
    # ======================================================
    # Testing continuous-time homogeneous system (no input)
    # ======================================================
    A = sparse([1, 1, 2, 3, 4], [1, 2, 2, 4, 3], [1., 2., 3., 4., 5.], 4, 4)
    X0 = BallInf(zeros(4), 0.1)
    cont_sys_homog = IVP(LCS(A), X0)

    # Check if the input is constant
    @test inputdim(cont_sys_homog) == 0
    # Check if the input is empty
    #@test isa(next_set(cont_sys_homog.U), ZeroSet) # there is no input set
    # Check data fields
    @test cont_sys_homog.s.A == A
    @test cont_sys_homog.x0.center == zeros(4) && cont_sys_homog.x0.radius == 0.1

    # ===================================================
    # Testing continuous-time system with constant input
    # ===================================================
    U = ConstantInput(Ball2(ones(4), 0.5))
    cont_sys = IVP(CLCCS(A, Matrix(1.0I, 4, 4), nothing, U), X0)

    # check initial state
    @test cont_sys.x0.center ≈ zeros(4) && cont_sys.x0.radius ≈ 0.1

    # recover input
    inputs = next_set(cont_sys.s.U) #collect(nextinput(cont_sys.s.U))[1]

    @test inputs.center ≈ ones(4) && inputs.radius ≈ 0.5

    # ========================================================
    # Testing continuous-time system with time-varying input
    # ========================================================
    U = VaryingInput([Ball2(0.01*i*ones(4), i*0.2) for i in 1:3])
    cont_sys = cont_sys = IVP(CLCCS(A, Matrix(1.0I, 4, 4), nothing, U), X0)
    for (i, inputs) in enumerate(cont_sys.s.U)
        if i == 1
            @test inputs.center ≈ 0.01*ones(4)
            @test inputs.radius ≈ 0.2
        elseif i == 2
            @test inputs.center ≈ 0.02*ones(4)
            @test inputs.radius ≈ 0.4
        else
            @test inputs.center ≈ 0.03*ones(4)
            @test inputs.radius ≈ 0.6
        end
    end

    # =========================================
    # Testing discrete-time homogeneous system
    # =========================================
    δ = 0.01
    discr_sys_homog = IVP(LDS(A), X0)

    # Check if the input is constant
    #@test isa(discr_sys_homog.s.U, ConstantInput) # no field U
    # Check if the input is empty
    #@test isa(next_set(discr_sys_homog.s.U), ZeroSet)
    # Check data fields
    @test discr_sys_homog.s.A == A
    @test discr_sys_homog.x0.center == zeros(4) && discr_sys_homog.x0.radius == 0.1

    # =================================================
    # Testing discrete-time system with constant input
    # =================================================
    U = ConstantInput(Ball2(ones(4), 0.5))
    discr_sys = IVP(CLCDS(A, Matrix(1.0I, 4, 4), nothing, U), X0)

    # check initial state
    @test discr_sys.x0.center ≈ zeros(4) && discr_sys.x0.radius ≈ 0.1

    # recover input
    inputs = collect(nextinput(discr_sys.s.U))[1]

    @test inputs.center ≈ ones(4) && inputs.radius ≈ 0.5

    # =====================================================
    # Testing discrete-time system with time-varying input
    # =====================================================
    U = [Ball2(0.01*i*ones(4), i*0.2) for i in 1:3]
    discr_sys = IVP(CLCDS(A, Matrix(1.0I, 4, 4), nothing, U), X0)
    Uset = inputset(discr_sys.s)
    for (i, inputs) in enumerate(cont_sys.s.U)
        if i == 1
            @test inputs.center ≈ 0.01*ones(4)
            @test inputs.radius ≈ 0.2
        elseif i == 2
            @test inputs.center ≈ 0.02*ones(4)
            @test inputs.radius ≈ 0.4
        else
            @test inputs.center ≈ 0.03*ones(4)
            @test inputs.radius ≈ 0.6
        end
    end
end
