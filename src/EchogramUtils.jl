module EchogramUtils

export hysthresh, bwselect

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
