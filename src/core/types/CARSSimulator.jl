mutable struct CARSSimulator
    species   ::Vector{CARSSpecies}
    conditions::GasConditions
    lasers    ::LaserConfiguration
    instrument::InstrumentConfiguration
    grid      ::Union{UniformGrid, AdaptiveGrid}
    vertical_shift::Float64
end

function CARSSimulator(;
    species   ::Vector{T},
    conditions::GasConditions,
    lasers    ::LaserConfiguration,
    instrument::InstrumentConfiguration,
    grid_type ::Symbol = :adaptive,
    vertical_shift::Float64 = 0.0
    ) where {T<:CARSSpecies}

    if grid_type == :adaptive
        grid = AdaptiveGrid(species = species, lasers = lasers, conditions = conditions)
    elseif grid_type == :uniform
        grid = UniformGrid(species = species, lasers = lasers, conditions = conditions)
    end

    CARSSimulator(species, conditions, lasers, instrument, grid, vertical_shift)
end