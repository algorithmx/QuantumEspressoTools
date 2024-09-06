
conv_xml_dic_old(d) = ( (d isa XMLDict.XMLDictElement) 
                        ? Dict(k=>conv_xml_dic(v) for (k,v) in d) 
                        : ((eltype(d) <: XMLDict.XMLDictElement) ? conv_xml_dic.(d) : strtonum_fort(d)) )



conv_xml_dic(d) = ( (d isa OrderedCollections.OrderedDict) 
                        ? Dict(k=>conv_xml_dic(v) for (k,v) in d) 
                        : ((eltype(d) <: OrderedCollections.OrderedDict) ? conv_xml_dic.(d) : strtonum_fort(d)) )


load_xml(fn::String) = xml_dict(join(readlines(fn),"\n")) 


function conv_xml(fn::String)
    xml = load_xml(fn)
    dic = conv_xml_dic(xml)
    if "espresso" in keys(dic)
        return dic["espresso"]
    else
        throw(error("Key espresso not found."))
        return Dict()
    end
end


