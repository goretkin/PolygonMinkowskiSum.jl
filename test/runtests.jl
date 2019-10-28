using PolygonMinkowskiSum
using Test

using GeometryTypes: Point
using PolygonMinkowskiSum: SimplePolygon, cross, find_ccw_violations

sweepme = SimplePolygon{Array{Point{2,Float64},1}}(Point{2,Float64}[[-0.35, -0.35], [0.35, -0.35], [0.35, 0.35], [-0.35, 0.35]])
alongme = SimplePolygon{Array{Point{2,Float64},1}}(Point{2,Float64}[[0.19857737682165164, 0.0], [0.19857737682165164, 0.14427490936067594], [0.06452170095561843, 0.19857737682165164], [-0.06452170095561849, 0.19857737682165164], [-0.19857737682165166, 0.1442749093606759], [-0.19857737682165166, 0.0], [-0.19857737682165166, -0.14427490936067594], [-0.06452170095561849, -0.19857737682165166], [0.06452170095561843, -0.19857737682165166], [0.19857737682165164, -0.144274909360676]])

@test find_ccw_violations(sweepme) == []
@test find_ccw_violations(alongme) == []

result = minkowski_sum(sweepme, alongme)

@test find_ccw_violations(result) == []
