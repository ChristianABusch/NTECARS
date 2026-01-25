"""
    GasConditions(;pressure::AbstractFloat, T_gas::AbstractFloat) 

Contains the basis information of the pressure in Pa and translational temperature in K.

The pressure and translational temperature are used for calculating linewidths. When fitting,
it should be remembered to update `T_gas` accordingly.
"""
mutable struct GasConditions
    pressure::AbstractFloat
    T_gas   ::AbstractFloat
end

function GasConditions(;pressure::AbstractFloat, T_gas::AbstractFloat) 
    GasConditions(pressure, T_gas)
end