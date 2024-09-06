import Base: >>

function >>(s::Vector{S}, fn::AbstractString) where {S<:AbstractString}
    # valid path ?
	open(fn,"w") do file
		write(file, join(s,"\n"))
	end
end

function >>(s::AbstractString, fn::AbstractString)
    split(s,"\n") >> fn
end


# robost mkdir()
# try create folde fdn, if failed, use default def. If def failed, use pwd().
function try_mkdir(fdn, def, msg)
	if !isdir(def)
		@info  "Default folder ( $(def) ) is not a folder. Use pwd ($(pwd())) as def."
		def=pwd()
	end
	if isdir(fdn)
		@info  msg * "\n$(fdn) exists. Use $(fdn)."
		return fdn
	end
	cnt = 0 
	while cnt < 10
		try
			mkdir(fdn)
		catch
			@info  msg * "\n$(fdn) not created. Try again $(cnt+1)."
		end
		if isdir(fdn)  return fdn   end
		sleep(rand(1:10))
		cnt += 1
	end
	if isdir(fdn)
		return fdn
	else
		@info  msg * "\n$(fdn) not created. Use $(def)."
		return def
	end
end


try_mkdir(fdn) = try_mkdir(fdn, "./", "try_mkdir($(fdn))")
