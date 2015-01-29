using PyPlot
using PyCall
using BilliardModels
using Vector2d

@pyimport matplotlib.patches as patches



function visualize!(d::Disc, subplot)

    circ = patches.Circle((d.centre.x, d.centre.y), d.radius)
    subplot[:add_patch](circ)
end

function visualize!(p::Plane, subplot)
end

function visualize!(table::BilliardTable, subplot)

    xmax, xmin = -Inf, Inf
    ymax, ymin = -Inf, Inf

    for obstacle in table.obstacles
        if typeof(obstacle) == Plane
            xmin = min(xmin, obstacle.c.x)
            xmax = max(xmax, obstacle.c.x)
            ymin = min(ymin, obstacle.c.y)
            ymax = max(ymax, obstacle.c.y)
        else
            visualize!(obstacle, subplot)
        end
    end
    
    subplot[:plot]([xmin,xmax,xmax,xmin,xmin],
                 [ymin,ymin,ymax,ymax,ymin],"-k",lw=2)
end

function visualize!{T}(xs::Array{Vector2D{T},1}, subplot)

    x = [pt.x for pt in xs]
    y = [pt.y for pt in xs]
    subplot[:plot](x, y, "-")
    axes[:axis]("image")

end
