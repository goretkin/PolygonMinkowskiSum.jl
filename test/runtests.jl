using PolygonMinkowskiSum
using Test

import Unitful: °
using GeometryTypes: Point
using PolygonMinkowskiSum: SimplePolygon

function make_regular_ngon(n; offset=0°, radius=1)
    points = [radius .* Point(cos(θ+offset), sin(θ+offset))
        for θ = range(0°, 360°, length=(n+1))][1:end-1]
    SimplePolygon(points)
end

using Plots


diamond = make_regular_ngon(4; offset=0°)
square = make_regular_ngon(4; offset=45°)

big_square = minkowski_sum(square, square)
plot(big_square, markershape=:cross)


for offset = 0°:45°:360°
    @show offset
    diamond = make_regular_ngon(4; offset=0°)
    square = make_regular_ngon(4; offset=offset#=+45°=#, radius=0.5)


    octogon = minkowski_sum(square, diamond)
    p = plot(aspect_ratio=1.0)
    plot!(p ,diamond, fillalpha=0.2)
    plot!(p, square, fillalpha=0.2)
    plot!(p, octogon, fillalpha=0.2)
    display(p)
end
