struct InstrumentConfiguration
    profile::Union{DeltaProfile,Spectrum}

    function InstrumentConfiguration(;profile::Union{DeltaProfile,Spectrum} = DeltaProfile())
        new(profile)
    end
end