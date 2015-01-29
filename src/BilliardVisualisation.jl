using PyPlot
using PyCall
using BilliardModels

@pyimport matplotlib.patches as patches

import PyPlot.draw
export draw

function draw(d::Disc, subplot)  # default is the current axis

    circ = patches.Circle((d.centre.x, d.centre.y), d.radius, alpha=0.5)
    subplot[:add_patch](circ)
end

function draw(p::Plane, subplot)  # don't draw anything for plane? This is too naive
end

function draw(table::BilliardTable, subplot)

    # Calculate the minimum and maximum of the points used to define the planes:
    xmax, xmin = -Inf, Inf
    ymax, ymin = -Inf, Inf

    for obstacle in table.obstacles
        if typeof(obstacle) == Plane
            xmin = min(xmin, obstacle.c.x)
            xmax = max(xmax, obstacle.c.x)
            ymin = min(ymin, obstacle.c.y)
            ymax = max(ymax, obstacle.c.y)
        else
            draw(obstacle, subplot)  # draw the non-plane objects
        end
    end

    # Draw a square bounding box according to the points found:
    subplot[:plot]( [xmin, xmax, xmax, xmin, xmin],
                    [ymin, ymin, ymax, ymax, ymin],"-k",lw=2)
end

function draw{T}(xs::Array{Vector2D{T},1}, subplot)

    x = [pt.x for pt in xs]
    y = [pt.y for pt in xs]
    subplot[:plot](x, y, "-", alpha=0.5)
    subplot[:axis]("image")

end

function draw(p::Particle, subplot; draw_velocity=true)
    pos = p.x
    subplot[:plot]([pos.x], [pos.y], "o")

   # if draw_velocity, draw an arr

end
