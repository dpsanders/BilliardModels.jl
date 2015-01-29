using PyPlot
using PyCall
using BilliardModels

@pyimport matplotlib.patches as patches

import PyPlot.draw
export draw

function draw(d::Disc, subplot)

    circ = patches.Circle((d.centre.x, d.centre.y), d.radius, alpha=0.5)
    subplot[:add_patch](circ)
end

function draw(p::Plane, subplot)
end

function draw(table::BilliardTable, subplot)

    xmax, xmin = -Inf, Inf
    ymax, ymin = -Inf, Inf

    for obstacle in table.obstacles
        if typeof(obstacle) == Plane
            xmin = min(xmin, obstacle.c.x)
            xmax = max(xmax, obstacle.c.x)
            ymin = min(ymin, obstacle.c.y)
            ymax = max(ymax, obstacle.c.y)
        else
            draw(obstacle, subplot)
        end
    end

    subplot[:plot]([xmin,xmax,xmax,xmin,xmin],
                 [ymin,ymin,ymax,ymax,ymin],"-k",lw=2)
end

function draw{T}(xs::Array{Vector2D{T},1}, subplot)

    x = [pt.x for pt in xs]
    y = [pt.y for pt in xs]
    subplot[:plot](x, y, "-", alpha=0.5)
    axes[:axis]("image")

end
