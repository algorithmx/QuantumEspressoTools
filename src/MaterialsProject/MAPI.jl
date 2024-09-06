# https://materialsproject.org/dashboard


function get_vasp(mp)
    url = "https://www.materialsproject.org/rest/v2/materials/mp-$(mp)/vasp"
    r = nothing
    try
        r = HTTP.request( 
                    "GET", 
                    url, 
                    ["X-API-KEY" => __MP_API_KEY__],
        )
    catch _e_
        println(_e_)
        return Dict()
    end
    return JSON.parse(String(r.body))
end
