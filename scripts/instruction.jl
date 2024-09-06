INSTRKEY = Union{Symbol,String}

mutable struct INSTRUCTION
    #* common part of the jobs, need to complete all necessary keys
    SEED::Dict{Cmd, Dict{INSTRKEY,Any}}
    #* [(plan_name1, program1, modify_dict1), (plan_name2, program2, modify_dict2), ...] #! plan_name is used as prefix 
    MODIFY::Vector{Tuple}
end


update_strategy(S::INSTRUCTION, d::Dict{INSTRKEY,Any}) = INSTRUCTION(S.SEED, [(n, p, m ⬱ d) for (n,p,m) in S.MODIFY])

purge_strategy(S::INSTRUCTION, kd::Vector{INSTRKEY}) = INSTRUCTION(purge(S.SEED[p],kd), [(n, p, m ↓ kd) for (n,p,m) in S.MODIFY])

↑(S::INSTRUCTION,d::Dict{INSTRKEY,Any}) = update_strategy(S,d)

↑(S::INSTRUCTION,p::Pair{INSTRKEY,Any}) = update_strategy(S,Dict{INSTRKEY,Any}(p))

↓(S::INSTRUCTION, kd::Vector{INSTRKEY}) = purge_strategy(S,kd)

expand(S::INSTRUCTION) = [(plan_name, prog, (S.SEED[prog] ⬱ modify_dict) ) for (plan_name, prog, modify_dict) in S.MODIFY]
