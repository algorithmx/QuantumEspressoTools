function user_pseudo(env)
    if "USER_PSEUDO" ∉ keys(env)
        @warn "ENV[\"USER_PSEUDO\"] is not specified." 
    elseif !isdir(env["USER_PSEUDO"])
        @warn "ENV[\"USER_PSEUDO\"] is not a folder." 
    else
        if length(readdir(env["USER_PSEUDO"]))==0
            @warn "ENV[\"USER_PSEUDO\"] = $(env["USER_PSEUDO"]) folder empty." 
        end
    end
    dic = ( ("USER_PSEUDO" in keys(env)) && isdir(env["USER_PSEUDO"]) 
            ? Dict(  Find_element_name("$(env["USER_PSEUDO"])/$f")=>extract_upf("USER_PSEUDO","$(env["USER_PSEUDO"])/$f")
                    for f ∈ readdir(env["USER_PSEUDO"])  ) 
            : Dict() )
    return dic
end

