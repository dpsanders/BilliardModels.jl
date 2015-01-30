
#module BilliardModels

#include("Vector2d.jl")

#using Vector2d
using Docile


#export run_Lorentz_gas
#export run_many_particles

export Particle, Disc, collision_time, Plane, BilliardTable
export Sinai_billiard
export calculate_next_collision, billiard_dynamics, initial_condition

Vector2D

type Particle
    x::Vector2D
    v::Vector2D
end


abstract Obstacle

immutable Disc <: Obstacle
    centre::Vector2D
    radius::Real
end

@doc """Compute *time of collision* of a particle and a disc,
        assuming the particle starts outside the disc and the
        particle speed is one (using the quadratic formula).

        Returns -1 if no collision.
        """ ->
function collision_time(p::Particle, disc::Disc)

    disp = p.x-disc.centre

    B = dot(p.v, disp)
    C = dot(disp, disp) - disc.radius^2

    discriminant = B*B - C

    if discriminant < 0
        return -1.
    end

    return -B - √discriminant  # NB: this supposes that the particle is *outside* the disc

end

@doc """Compute normal vector to disc boundary, at a point x that must lie *on the boundary*""" ->
function normal(disc::Disc, x)
    return (x - disc.centre) / disc.radius
end

@doc """Check whether a given position is inside a disc.""" ->
isvalid(x, disc::Disc) = norm(x - disc.centre) > disc.radius  # must change if


immutable Plane <: Obstacle
    c::Vector2D  # an aribtrary point on the plane
    normal::Vector2D
end


collision_time(p::Particle, plane::Plane) =
    dot(plane.c - p.x, plane.normal) / dot(p.v, plane.normal)

normal(plane::Plane, x) = plane.normal  # same normal vector for all positions x on plane

@doc """Convention: the normal points to the allowed part of the billiard tabl""" ->
isvalid(x, plane::Plane) = dot(x - plane.c, plane.normal) > 0.0


type BilliardTable
  obstacles::Vector{Obstacle}
end

function isvalid(x, table::BilliardTable)
    for obstacle in table.obstacles
        if !isvalid(x, obstacle)
            return false
        end
    end

    true
end


function calculate_next_collision(p::Particle, billiard_table, previous_obstacle_hit)

    obstacles = billiard_table.obstacles

    first_collision_time = Inf
    which_obstacle_hit = nothing

    for obstacle in obstacles
        obstacle === previous_obstacle_hit && continue

        t = collision_time(p, obstacle)

        if t > 0.0 && t < first_collision_time
            which_obstacle_hit, first_collision_time = obstacle, t
        end
    end

    x_collision = p.x + p.v*first_collision_time

    v_new = post_collision_velocity(x_collision, p.v, which_obstacle_hit)

    return x_collision, v_new, first_collision_time, which_obstacle_hit
end

function post_collision_velocity(x_collision, v, obstacle)
    n = normal(obstacle, x_collision)

    v_new = v - 2.0*dot(n,v)*n  # reflejar solo si es disco; si es plano, pasar a traves

    speed = norm(v_new)
    v_new /= speed
    v_new
end


# design so that the normal vector points into the allowed region!
function initial_condition(table, xmin, xmax, ymin, ymax)
    x, y = xmin, ymin

    while true
        x = xmin + rand()*(xmax-xmin)
        y = ymin + rand()*(ymax-ymin)

        if isvalid(Vector2D(x,y), table)  # valid if outside all discs
            break
        end
    end

    # velocity in 2D:
    theta = rand() * 2*pi

    xx = Vector2D(x, y)
    vv = Vector2D(cos(theta), sin(theta))

    return xx, vv
end




### Main program

function Sinai_billiard(radius)

    obstacles = Obstacle[]

    push!(obstacles, Disc([0., 0.], radius) )

    push!(obstacles, Plane([-0.5, -0.5], [0., 1.]) )
    push!(obstacles, Plane([0.5, -0.5], [-1., 0.]) )
    push!(obstacles, Plane([0.5, 0.5], [0., -1.]) )
    push!(obstacles, Plane([-0.5, 0.5], [1., 0.]) )

    return BilliardTable(obstacles)
end


@doc """Simulate a single particle p on a billiard table for a given of collisions"""->
function billiard_dynamics(p, table, num_collisions)

    xs = [p.x]

    which_obstacle_hit = nothing

    for t in 1:num_collisions
        x_collision, v_new, first_collision_time, which_obstacle_hit = calculate_next_collision(p, table, which_obstacle_hit)

        p.x, p.v = x_collision, v_new

        push!(xs, p.x)

    end

    xs

end


