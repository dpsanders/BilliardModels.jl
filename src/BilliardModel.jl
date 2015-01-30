
#module BilliardModels

#include("Vector2d.jl")

#using Vector2d
using Docile
using BilliardModels

#export run_Lorentz_gas
#export run_many_particles

export Particle, Disc, collision_time, Plane, BilliardTable, AbstractPlane, AbstractParticle
export Sinai_billiard
export calculate_next_collision, billiard_dynamics, initial_condition

# lattice stuff:
export ParticleOnLattice, CellBoundary, billiard_dynamics_on_lattice

abstract AbstractParticle

type Particle <: AbstractParticle
    x::Vector2D
    v::Vector2D
end

type ParticleOnLattice <: AbstractParticle
    x::Vector2D
    v::Vector2D
    lattice_vector::Vector2D{Int}
end


abstract Obstacle

immutable Disc <: Obstacle
    centre::Vector2D
    radius::Real
end

abstract AbstractPlane <: Obstacle

immutable Plane <: AbstractPlane
    c::Vector2D  # an aribtrary point on the plane
    normal::Vector2D
end

@doc """A CellBoundary is a boundary between parts of a lattice.
        It is assumed to be planar.  The constructor creates an incomplete
        object, so the other_side must be added afterward.""" ->
type CellBoundary <: AbstractPlane
    c::Vector2D
    normal::Vector2D
    lattice_increment::Vector2D{Int}
    other_side::CellBoundary

    CellBoundary(c,normal,lattice_increment) = new(c,normal,lattice_increment)
end

get_lattice_increment(o::Obstacle) = [0, 0]
get_lattice_increment(o::CellBoundary) = o.lattice_increment

@doc """Compute *time of collision* of a particle and a disc,
        assuming the particle starts outside the disc and the
        particle speed is one (using the quadratic formula).

        Returns -1 if no collision.
        """ ->
function collision_time(p::AbstractParticle, disc::Disc)

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


collision_time(p::AbstractParticle, plane::AbstractPlane) =
    dot(plane.c - p.x, plane.normal) / dot(p.v, plane.normal)

normal(plane::AbstractPlane, x) = plane.normal  # same normal vector for all positions x on plane

@doc """Convention: the normal points to the allowed part of the billiard tabl""" ->
isvalid(x, plane::AbstractPlane) = dot(x - plane.c, plane.normal) > 0.0


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


function calculate_next_collision_on_lattice(p::ParticleOnLattice, billiard_table, previous_obstacle_hit)

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

    lattice_increment = get_lattice_increment(which_obstacle_hit)

    x_new, v_new, which_obstacle_hit = collide(x_collision, p.v, which_obstacle_hit)

    return x_new, v_new, first_collision_time, which_obstacle_hit, lattice_increment
end

function collide(x_collision, v, obstacle::Obstacle)
    n = normal(obstacle, x_collision)

    v_new = v - 2.0*dot(n,v)*n  # reflejar solo si es disco; si es plano, pasar a traves

    speed = norm(v_new)
    v_new /= speed
    return x_collision, v_new, obstacle
end

@doc """This function assumes that the distance between opposite
        CellBoundary objects is 1 and the normal is a unit normal
        pointing inwards.""" ->
function collide(x_collision, v, boundary::CellBoundary)

    x_new = x_collision + boundary.normal

    return x_new, v, boundary.other_side
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

function Sinai_billiard(radius,periodic_x=false,periodic_y=false)

    obstacles = Obstacle[]

    push!(obstacles, Disc([0., 0.], radius) )

    if periodic_x
        right = CellBoundary([0.5, -0.5], [-1., 0.],Vector2D(1,0))
        left  = CellBoundary([-0.5, 0.5], [1., 0.],Vector2D(-1,0))
        right.other_side = left
        left.other_side = right
        push!(obstacles, right, left)
    else
        push!(obstacles, Plane([0.5, -0.5], [-1., 0.]) )
        push!(obstacles, Plane([-0.5, 0.5], [1., 0.]) )
    end

    if periodic_y
        top = CellBoundary([0.5, 0.5], [0., -1.],Vector2D(0,1))
        bottom  = CellBoundary([-0.5, -0.5], [0., 1.],Vector2D(0,-1))
        top.other_side = bottom
        bottom.other_side = top
        push!(obstacles, top, bottom)
    else
        push!(obstacles, Plane([-0.5, -0.5], [0., 1.]) )
        push!(obstacles, Plane([0.5, 0.5], [0., -1.]) )
    end

    return BilliardTable(obstacles)
end


@doc """Simulate a single particle p on a billiard table for a given of collisions"""->
function billiard_dynamics(p, table, num_collisions)

    xs = [p.x]

    which_obstacle_hit = nothing

    for t in 1:num_collisions
        x_collision, v_new, first_collision_time, which_obstacle_hit =
            calculate_next_collision(p, table, which_obstacle_hit)

        p.x, p.v = x_collision, v_new

        push!(xs, p.x)

    end

    xs

end

@doc """Simulate a single particle p on a billiard table for a given of collisions on
        a periodic billiard table."""->
function billiard_dynamics_on_lattice(p, table, num_collisions)

    xs = [p.x]
    lattice_vectors = [p.lattice_vector]

    which_obstacle_hit = nothing

    for t in 1:num_collisions
        x_new, v_new, first_collision_time, which_obstacle_hit, lattice_increment =
            calculate_next_collision_on_lattice(p, table, which_obstacle_hit)

        # Position where it hit the obstacle: (This is a bit of a hack)
        x_collision = p.x + p.v*first_collision_time
        push!(xs,x_collision)
        push!(lattice_vectors, p.lattice_vector)


        p.x, p.v = x_new, v_new
        p.lattice_vector += lattice_increment
        push!(xs, p.x)
        push!(lattice_vectors, p.lattice_vector)

    end

    return xs, lattice_vectors

end
