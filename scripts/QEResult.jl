mutable struct QEResult
    PLAN::String
    PROG::Cmd
    CONF::Dict
    WORKSPACE::String
    INP::Vector{String}
    INP_F::String
    SUCC::Bool
    STAT::Symbol
    OUTP::Vector{String}
    OUTP_F::String
    RES::Dict
    UPDATER::Function
    WD::NamedTuple
    SCRATCH::Bool
end


#TODO  making QEResult from execute_serial() into a blockchain 
#+ by intruducing stamps 

import Base.hash

function hash(R::QEResult)::UInt64
    return UInt64(0)
end
