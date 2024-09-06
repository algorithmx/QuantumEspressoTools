@inline divisors(x) = [j for j=2:x-1  if x%j==0]

# yields good to very good scaling, especially if the number of processors in a pool is a
# divisor of N 3 and N r3 (the dimensions along the z-axis of the FFT grids, nr3 and nr3s,
# which coincide for NCPPs);

# Quick estimate of parallelization parameters You need to know
# • the number of k-points, Nk
# • the third dimension of the (smooth) FFT grid, N3
# • the number of Kohn-Sham states, M

# TODO
function determine_npools_1(np::Int, Nkp::Int, N3::Int)
    div_np = divisors(np)
    # if the number of processors in a pool is a divisor of N3
    npools_sugg_1 = np .÷ intersect(div_np, divisors(N3))
    # if the number of k-points is a multiple of the number of pools
    npools_sugg_2 = intersect(div_np, divisors(Nkp))
    if length(npools_sugg_1)==0
        if length(npools_sugg_2)==0
            return last([i for i in div_np if i<sqrt(np)])
        else
            return maximum(npools_sugg_2)
        end
    else
        if length(npools_sugg_2)==0
            return maximum(npools_sugg_1)
        else
            intsct = intersect(npools_sugg_1, npools_sugg_2)
            if length(intsct)==0
                return maximum(union(npools_sugg_1, npools_sugg_2))
            else
                return maximum(intsct)
            end
        end
    end
end


# TODO
function determine_npools(np::Int, Nkp::Int, N3::Int)
    div_np = divisors(np)
    # if the number of processors in a pool is a divisor of N3
    npools_sugg_1 = np .÷ intersect(div_np, divisors(N3))
    # if the number of k-points is a multiple of the number of pools
    npools_sugg_2 = intersect(div_np, divisors(Nkp))
    @show npools_sugg_1
    @show npools_sugg_2
    u = union(npools_sugg_1, npools_sugg_2)
    if length(u)!=0
        return maximum(u)
    else
        return last([i for i in div_np if i<=sqrt(np)])
    end
end