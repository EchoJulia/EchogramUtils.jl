using EchogramUtils
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@test pow2db(1000) == 30
@test db2pow(30) == 1000

@test mag2db(10) == 20
@test db2mag(20) == 10
