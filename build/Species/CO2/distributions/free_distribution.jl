################################################################################
#                             NOT IMPLEMENTED YET
################################################################################


#=

mutable struct VibrationalStateDistribution <: CO2Distribution
    Q     ::Float64
    T_rot ::Float64
    fvib  ::Dict{VibrationalQuantumNumbers, Float64}
    fvib_extrapolation ::Dict{VibrationalQuantumNumbers, Float64}
    iso_ID::Symbol

    function VibrationalStateDistribution(T_rot, f_vib::Dict{VibrationalQuantumNumbers, Float64}; iso_ID = :O16C12O16)
        # creates VibrationalStateDistribution object.
        # for all states included in the fit, the state with v₃+1 is extrapolated to
        # as a closing boundary condition 
        # then calculates the partition sum (using only the values included in the fit atm)
        f_vib_extrapolate = Dict{VibrationalQuantumNumbers, Float64}()

        # boundary condition for ground state
        f_vib[VibrationalQuantumNumbers(0,0,0,0,1)] = 1.0

        # if needed, extrapolate to v₃+1 since all absorption lines have an increase in v₃   
        for vib_qn in keys(f_vib)
            @unpack v₁, v₂, l₂, v₃, r = vib_qn
            qn′ = VibrationalQuantumNumbers(v₁, v₂, l₂, v₃+1, r)

            if haskey(f_vib, qn′)
                continue
            else
                # extrapolate for closing boundary condition
                qn_pure_v₃       = VibrationalQuantumNumbers(0, 0, 0, v₃+1, 1)
                qn_pure_v₃_lower = VibrationalQuantumNumbers(0, 0, 0, v₃,   1)

                haskeys_for_T3_extrapolation = haskey(f_vib, qn_pure_v₃_lower) && haskey(f_vib, qn_pure_v₃)
                if haskeys_for_T3_extrapolation
                    f_vib_extrapolate[qn′] = f_vib[vib_qn] * f_vib[qn_pure_v₃]/f_vib[qn_pure_v₃_lower]
                end
            end
        end

        df   = new(1.0, T_rot, f_vib, f_vib_extrapolate, iso_ID)
        df.Q = partition_sum(df)
        return df
    end

    function VibrationalStateDistribution(df::MultiTemperatureDistribution)
        iso_ID = df.iso_ID
        qns = mode_population_dict_for_fitting()
        
        for qn in keys(qns)
            vib_state = VibrationalState(qn, get_CO₂_metadata(iso_ID))
            qns[qn] = fvib(df, vib_state)
        end

        df = new(df.T_rot, qns, iso_ID)
        return df
    end
end

function is_state_included_in_fit(df::VibrationalStateDistribution, vib_qn::VibrationalQuantumNumbers)
    return haskey(df.fvib, vib_qn)
end

function is_state_included_in_extrapolation(df::VibrationalStateDistribution, vib_qn::VibrationalQuantumNumbers)
    return haskey(df.fvib_extrapolation, vib_qn)
end

function fvib(df::VibrationalStateDistribution, state::Union{State, VibrationalState})
    v₁, v₂, l₂, v₃, r = state.qn.v₁, state.qn.v₂, state.qn.l₂, state.qn.v₃, state.qn.r
    vib_qn = VibrationalQuantumNumbers(v₁, v₂, l₂, v₃, r)
    fvib(df, vib_qn)
end

function fvib(df::VibrationalStateDistribution, vib_qn::VibrationalQuantumNumbers)
    if is_state_included_in_fit(df, vib_qn)
        # return the population set by the fitting algorithm
        df.fvib[vib_qn]
    elseif is_state_included_in_extrapolation(df, vib_qn)
        # all upper states in the absorption have v₃+1 to which is
        # being extrapolated. If an extrapolated value exist, return it.
        df.fvib_extrapolation[vib_qn]
    else
        return 0.0
    end
end
=#