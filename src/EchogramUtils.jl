module EchogramUtils

using Statistics

export db2pow, pow2db, mag2db, db2mag, hysthresh, bwselect, vertically_smooth, vertically_bin, INmask,INfilter, meandb

"""
    db2pow(ydb)

Convert decibels to power
"""
db2pow(ydb) = 10^(ydb/10)

"""
    pow2db(y)

Convert power to decibels.
"""
pow2db(y) = 10log10(y)

"""
    mag2db(y)

Convert magnitude to decibels.
"""
mag2db(y) = 20log10(y)

"""
    db2mag(ydb)

Convert decibels to magnitude.
"""
db2mag(ydb) = 10^(ydb/20)


"""
    meandb(ydb)

Remember, the Average of the Log is Not the Log of the Average,
so averaging dB requires conversion to and from the linear domain.
"""
function meandb(ydb)
    a = db2pow.(ydb)
    pow2db(mean(a))
end

"""
    vertically_smooth(A, r, thickness)

Taking an echogram array `A` and a range vector `r`, smooth values
vertically by taking the mean over succesive `thickness` bins.

N.B. if your data is in dB, consider converting to linear (perhaps
using `dB2pow`).

"""
function vertically_smooth(A, r; thickness=1.43)
    b = zeros(size(A))
    m,n = size(A)
    for j =1:n
        for i in 0:thickness:r[end,j]
            mask = (r[:,j] .>= i)  .& (r[:,j] .< i+thickness)
            b[mask,j] .= mean(A[mask,j])
        end
    end
    return b
end

"""
    vertically_bin(A, r, thickness)

Taking an echogram array `A` and a range vector `r`, bin values
vertically by taking the mean over succesive `thickness` sized bins.

N.B. if your data is in dB, consider converting to linear (perhaps
using `dB2pow`).

"""
function vertically_bin(A, r; thickness=1.43)
    b = []
    m,n = size(A)
    for j =1:n
        for i in 0:thickness:r[end,j]
            mask = (r[:,j] .>= i)  .& (r[:,j] .< i+thickness)
            push!(b,mean(A[mask,j]))
        end
    end
    return vcat(b...)
end

"""
    INmask(A, delta)

Impulse noise filter based based on the two-sided comparison method
described by Anderson et al. (2005) and further described in Ryan et
al. (2015).

It is often desirable to linearly average to a coarser vertical
resolution (perhaps using `vertically_smooth`) before calling IN.

`A` is an echogram array and `delta` is the threshold. A sample is
marked as true if its value is more than `delta` greater than samples
on either side.

Returns a `BitArray` with the same dimensions as `A`.

"""
function INmask(A; delta=10)
    m,n = size(A)

    a = A[:,1:end-2]
    b = A[:,2:end-1]
    c = A[:,3:end]
    m1 = (b.-a) .> delta
    m2 = (b.-c) .> delta

    hcat(falses(m), m1 .& m2, falses(m))
end

function mymedian(a,b)
    (a + b) /2
end

function INfilter(A; delta=10)
    mask = INmask(A,delta=delta)

    r = A[1:end,2:end]
    l = A[1:end,1:end-1]
    x = copy(A)

    x[mask] .= mymedian.(l[mask[1:end,2:end]], r[mask[1:end,1:end-1]])
    return x
end


"""
    hysthresh(im, T1, T2)

Hysteresis thresholding

A simple port of Peter Kovesi's MATLAB method to Julia.

See http://www.peterkovesi.com/matlabfns/Spatial/hysthresh.m

Usage: bw = hysthresh(im, T1, T2)

Arguments:
             im  - image to be thresholded (assumed to be non-negative)
             T1  - upper threshold value
             T2  - lower threshold value
                   (T1 and T2 can be entered in any order, the larger of the
                   two values is used as the upper threshold)
Returns:
             bw  - the thresholded image (containing values 0 or 1)

Function performs hysteresis thresholding of an image.
All pixels with values above threshold T1 are marked as edges
All pixels that are connected to points that have been marked as edges
and with values above threshold T2 are also marked as edges.

Copyright (c) 1996-2005 Peter Kovesi
School of Computer Science & Software Engineering
The University of Western Australia
http://www.csse.uwa.edu.au/
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.


"""
function hysthresh(im, T1, T2)

    if T1 < T2  # T1 and T2 reversed - swap values 
        tmp = T1;
        T1 = T2; 
        T2 = tmp;
    end
    
    # Edge points above lower threshold. 
    aboveT2 = im .> T2 
    
    # Row and colum coords of points above upper threshold.
    aboveT1r, aboveT1c = findn(im .> T1);  
 
    # Obtain all connected regions in aboveT2 that include a point that has a
    # value above T1 
    
    bw = bwselect(aboveT2, aboveT1c, aboveT1r);
    
end


"""
    bwselect(BW, c , r)

Select objects in a binary image.

Similar to the MATLAB function of the same name.
"""
function bwselect(BW, c, r)
    # constants
    north = CartesianIndex(-1,  0)
    south = CartesianIndex( 1,  0)
    east  = CartesianIndex( 0,  1)
    west  = CartesianIndex( 0, -1)
    
    queue = CartesianIndex.(r,c)
    
    m,n = size(BW)
    out = falses(m,n)
    
    while !isempty(queue)
        node = pop!(queue)
        
        if BW[node]
            wnode = node
            enode = node + east
        
             # Move west until node is false
            while checkbounds(Bool, BW, wnode) && BW[wnode] 
                out[wnode] = true
                if checkbounds(Bool, BW, wnode + north) && !out[wnode + north]
                    push!(queue, wnode + north)
                end
                if checkbounds(Bool, BW, wnode + south) && !out[wnode + south]
                    push!(queue, wnode + south)
                end
                wnode += west
            end
        
            # Move east until node is false
            while checkbounds(Bool, BW, enode) && BW[enode]
                out[enode] = true
                if checkbounds(Bool, BW, enode + north) && !out[enode + north]
                    push!(queue, enode + north)
                end
                if checkbounds(Bool, BW, enode + south) && !out[enode + south]
                    push!(queue, enode + south)
                end
                enode += east
            end
        end
        
    end
    return out
end

end # module
