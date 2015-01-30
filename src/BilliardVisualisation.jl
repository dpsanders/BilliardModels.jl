using PyPlot
using PyCall
using BilliardModels

@pyimport matplotlib.patches as patches

#import PyPlot.draw
export bdraw

@doc """The various `bdraw` functions draw parts of the billiard(hence the `b`
at the start of the name, to avoid using `draw` which conflicts with PyPlot.)""" ->
function bdraw(d::Disc, subplot)  # default is the current axis

    circ = patches.Circle((d.centre.x, d.centre.y), d.radius, alpha=0.5)
    subplot[:add_patch](circ)
end

function bdraw(p::AbstractPlane, subplot)  # don't draw anything for plane? This is too naive
end

function bdraw(table::BilliardTable, subplot)

    # Calculate the minimum and maximum of the points used to define the planes:
    xmax, xmin = -Inf, Inf
    ymax, ymin = -Inf, Inf

    for obstacle in table.obstacles
        if isa(obstacle, AbstractPlane)
            xmin = min(xmin, obstacle.c.x)
            xmax = max(xmax, obstacle.c.x)
            ymin = min(ymin, obstacle.c.y)
            ymax = max(ymax, obstacle.c.y)
        else
            bdraw(obstacle, subplot)  # draw the non-plane objects
        end
    end

    # Draw a square bounding box according to the points found:
    subplot[:plot]( [xmin, xmax, xmax, xmin, xmin],
                    [ymin, ymin, ymax, ymax, ymin],"-k",lw=2)
end

function bdraw(xs::Array, subplot)

    x = [pt.x for pt in xs]
    y = [pt.y for pt in xs]
    subplot[:plot](x, y, "-", alpha=0.5)
    subplot[:axis]("image")

end

function bdraw(p::AbstractParticle, subplot, draw_velocity=true)
    pos = p.x
    vel = p.v

    subplot[:plot]([pos.x], [pos.y], "o")

    dt = 0.05

    if draw_velocity
        subplot[:arrow](pos.x, pos.y, dt*vel.x, dt*vel.y, head_width=0.05, head_length=0.05)
    end

end
