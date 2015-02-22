
export Particle, Disc, collision_time, Plane, BilliardTable, AbstractPlane, AbstractParticle
export Sinai_billiard
export calculate_next_collision, billiard_dynamics, initial_condition

# exports of lattice functionality:
export ParticleOnLattice, CellBoundary, billiard_dynamics_on_lattice,
        calculate_next_collision_on_lattice

export continuous_time


# particle types:

abstract AbstractParticle{T}

type Particle{T} <: AbstractParticle{T}
    x::Vector2D{T}
    v::Vector2D{T}
end

Particle(x::Vector, v::Vector) = Particle(Vector2D(x), Vector2D(v))

type ParticleOnLattice{T} <: AbstractParticle{T}
    x::Vector2D{T}
    v::Vector2D{T}
    lattice_vector::Vector2D{Int}
end

# obstacle types:

abstract Obstacle{T}

immutable Disc{T} <: Obstacle{T}
    centre::Vector2D{T}
    radius::T
end

Disc(v::Vector, r) = Disc(Vector2D(v), r)

abstract AbstractPlane{T} <: Obstacle{T}

immutable Plane{T} <: AbstractPlane{T}
    point_on_plane::Vector2D{T}  # an aribtrary point on the plane
    normal::Vector2D{T}
end

Plane(v1::Vector, v2::Vector) = Plane(Vector2D(v1), Vector2D(v2))

@doc """A `CellBoundary` is a boundary between neighbouring cells on a lattice.
        It is assumed to be planar.  The constructor creates an incomplete
        object, so the `other_side` field must be added later
        (once the corresponding `CellBoundary` on the other side of the next cell has been constructed).

        `lattice_increment` is a vector which says which new lattice cell a particle enters when it crosses this
        `CellBoundary`.
        """ ->
type CellBoundary{T} <: AbstractPlane{T}  # immutable not (currently?) allowed due to the partial construction
    point_on_plane::Vector2D{T}
    normal::Vector2D{T}
    lattice_increment::Vector2D{Int}
    other_side::CellBoundary{T}

    CellBoundary(point, normal, lattice_increment) = new(point, normal, lattice_increment)  # partial constructor leaving other_side undefined
end

CellBoundary{T}(point::Vector{T}, normal::Vector{T}, lattice_increment::Vector) =
    CellBoundary{T}(Vector2D(point), Vector2D(normal), Vector2D(lattice_increment))


# BilliardTable type.

@doc """A `BilliardTable` is just a list of obstacles.""" ->
immutable BilliardTable{T}
  obstacles::Vector{Obstacle{T}}
end


@doc """`get_lattice_increment` returns the lattice_increment of an obstacle. This is a zero vector
except for `CellBoundary`s.""" ->

get_lattice_increment(o::Obstacle) = Vector2D(0, 0)  # non-boundary obstacles do not change the lattice cell
get_lattice_increment(o::CellBoundary) = o.lattice_increment



@doc """Compute *time of collision* of a particle and a disc,
        assuming the particle starts outside the disc and the
        particle speed is one (using the quadratic formula).

        Returns -Inf if no collision.
        """ ->
function collision_time{T}(p::AbstractParticle{T}, disc::Disc{T})

    disp = p.x - disc.centre

    B = dot(p.v, disp)
    C = dot(disp, disp) - disc.radius^2

    discriminant = B^2 - C

    if discriminant < zero(T)
        return -convert(T, Inf)
    end

    return convert(T, -B - √discriminant)  # NB: this supposes that the particle is *outside* the disc

end

@doc """Compute normal vector to disc boundary, at a point x that must lie *on the boundary*""" ->
function normal{T}(disc::Disc{T}, x::Vector2D{T})
    return (x - disc.centre) / disc.radius
end

@doc """Check whether a given position is inside a disc.""" ->
isoutside{T}(x::Vector2D{T}, disc::Disc{T}) = norm(x - disc.centre) > disc.radius


@doc """Calculate collison time of particle with plane.""" ->
collision_time{T}(p::AbstractParticle{T}, plane::AbstractPlane{T}) =
    dot(plane.point_on_plane - p.x, plane.normal) / dot(p.v, plane.normal)

normal(plane::AbstractPlane, x) = plane.normal  # same normal vector for all positions x on plane

@doc """Convention: the normal points towards the allowed part of the billiard table""" ->
isoutside(x, plane::AbstractPlane) = dot(x - plane.point_on_plane, plane.normal) > 0.0

@doc """`isoutside` finds if the particle is in the allowed part of the billiard table,
i.e. is "outside" the billiard obstacles.""" ->
function isoutside(x, table::BilliardTable)
    for obstacle in table.obstacles
        if !isoutside(x, obstacle)
            return false
        end
    end

    true
end

@doc """The two versions of `calculate_next_collision` return the information of what
the result is of doing the next collision, *without* updating the particle.""" ->
function calculate_next_collision_on_lattice{T}(p::ParticleOnLattice{T},
                                             billiard_table::BilliardTable{T},
                                             previous_obstacle_hit::Obstacle{T})

    obstacles = billiard_table.obstacles

    first_collision_time = convert(T, Inf)
    artificial_obstacle = Plane{T}([zero(T), zero(T)], [zero(T), zero(T)])
    which_obstacle_hit = artificial_obstacle
    # artificial non-existent obstacle to avoid type instability

    for obstacle in obstacles
        obstacle === previous_obstacle_hit && continue

        t = collision_time(p, obstacle)

        if 0.0 < t < first_collision_time
            which_obstacle_hit, first_collision_time = obstacle, t
        end
    end

    #@assert which_obstacle_hit != artificial_obstacle  # hit a real obstacle

    x_collision = p.x + p.v*first_collision_time
    lattice_increment = get_lattice_increment(which_obstacle_hit)

    x_new, v_new, which_obstacle_hit = collide(x_collision, p.v, which_obstacle_hit)

    return x_new, v_new, first_collision_time, which_obstacle_hit, lattice_increment
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

    x_collision, v_new, which_obstacle_hit = collide(x_collision, p.v, which_obstacle_hit)
    # Actually only need to return new velocity

    return x_collision, v_new, first_collision_time, which_obstacle_hit
end


@doc """`collide` *implements* an elastic collision""" ->
function collide{T}(x_collision::Vector2D{T}, v::Vector2D{T}, obstacle::Obstacle{T})
    n = normal(obstacle, x_collision)

    v_new = v - 2.0*dot(n,v)*n  # reflejar solo si es disco; si es plano, pasar a traves
    v_new /= norm(v_new)
    return x_collision, v_new, obstacle
end


@doc """This version of `collide` implements a "collision" with a `CellBoundary`.
        This just moves the particle to the opposite boundary (updating its position).
        This assumes that the distance between opposite
        CellBoundary objects is 1 and the normal is a unit normal
        pointing inwards (i.e. towards the available space on the billiard table).""" ->
function collide{T}(x_collision::Vector2D{T}, v::Vector2D{T}, boundary::CellBoundary{T})

    x_new = x_collision + boundary.normal

    return x_new, v, boundary.other_side  # velocity stays the same
end

@doc """Generate a random initial condition, uniformly in the allowed region of the billiard table
("outside") the billiard obstacles. Here, "outside" for a plane is taken to mean that
the allowed region lies *in the direction of the plane's normal vector*.

Velocities are generated uniformly on the unit circle S^1, so that in total the initial condition
is generated uniformly with respect to Liouville measure.
""" ->
function initial_condition(table, xmin, xmax, ymin, ymax)
    x, y = xmin, ymin

    while true
        x = xmin + rand()*(xmax-xmin)
        y = ymin + rand()*(ymax-ymin)

        if isoutside(Vector2D(x,y), table)  # valid if outside all discs
            break
        end
    end

    # velocity in 2D:
    theta = rand() * 2*pi

    xx = Vector2D(x, y)
    vv = Vector2D(cos(theta), sin(theta))

    return xx, vv
end


@doc """Create a Sinai billiard table (a square with a disc in the centre).
It may be periodic in the x and/or y directions.""" ->
function Sinai_billiard{T}(radius::T, periodic_x=false, periodic_y=false)

    obstacles = Obstacle{T}[]

    push!(obstacles, Disc{T}([0., 0.], radius) )

    if periodic_x
        # create a pair of opposite `CellBoundary`s:
        right = CellBoundary([0.5, -0.5], [-1., 0.], [1, 0])
        left  = CellBoundary([-0.5, 0.5], [1., 0.], [-1, 0])
        right.other_side = left
        left.other_side = right
        push!(obstacles, right, left)
    else
        push!(obstacles, Plane([0.5, -0.5], [-1., 0.]) )
        push!(obstacles, Plane([-0.5, 0.5], [1., 0.]) )
    end

    if periodic_y
        top = CellBoundary([0.5, 0.5], [0., -1.], [0, 1])
        bottom  = CellBoundary([-0.5, -0.5], [0., 1.], [0, -1])
        top.other_side = bottom
        bottom.other_side = top
        push!(obstacles, top, bottom)
    else
        push!(obstacles, Plane([-0.5, -0.5], [0., 1.]) )
        push!(obstacles, Plane([0.5, 0.5], [0., -1.]) )
    end

    return BilliardTable(obstacles)
end


@doc """Simulate a single particle p on a billiard table for a given number of collisions."""->
function billiard_dynamics{T}(p::Particle{T}, table::BilliardTable{T}, num_collisions::Integer)

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

@doc """Simulate discrete-time (collision) dynamics for a single particle
        p on a billiard table for a given number of obstacle collisions
        (excluding "collisions" on boundaries) on a *periodic* billiard table.\

        p contains the initial conditions.

        Returns the positions in the unit cells (*xs*) and the corresponding lattice cells (*ls*),
        as well as the free paths."""->
function billiard_dynamics_on_lattice{T}(p::ParticleOnLattice{T}, table::BilliardTable{T}, num_collisions::Integer)

    xs = [p.x]
    lattice_vectors = [p.lattice_vector]

    which_obstacle_hit = Plane([-Inf, -Inf], [0., 0.])
    # which_obstacle_hit = nothing

    free_path_length = 0.0
    free_paths = Float64[]
    num_obstacle_collisions = 0

    first_collision_time = Inf

    #for t in 1:num_collisions
    while num_obstacle_collisions < num_collisions
        x_new, v_new, first_collision_time, which_obstacle_hit, lattice_increment =
            calculate_next_collision_on_lattice(p, table, which_obstacle_hit)

        # version to draw periodised orbits:
        # outputs the position twice, once before and once after the collision

        ## Position where it hit the obstacle at the collision: (This is a bit of a hack)
        # x_collision = p.x + p.v*first_collision_time
        # push!(xs,x_collision)
        # push!(lattice_vectors, p.lattice_vector)

        # Periodise:
        p.x, p.v = x_new, v_new
        p.lattice_vector += lattice_increment


        # push!(xs, p.x)
        # push!(lattice_vectors, p.lattice_vector)


        ## version which takes account only of collisions with objects that are not CellBoundary

        free_path_length += first_collision_time

        #@show typeof(which_obstacle_hit)

        if !isa(which_obstacle_hit, CellBoundary)

            push!(xs, p.x)
            push!(lattice_vectors, p.lattice_vector)
            push!(free_paths, free_path_length)
            free_path_length = 0.0

            num_obstacle_collisions += 1

        end
    end

    return xs, lattice_vectors, free_paths

end


@doc """Extract the particle position for given continuous times that are multiples of delta_t
        from a discrete-time trajectory.

        Use delta_t that are representable floating-point numbers
        e.g. 0.125 instead of 0.1""" ->
function continuous_time(xs, lattice_vectors, free_paths, delta_t)
    t = 0.0
    collision_positions = xs + lattice_vectors
    collision_times = cumsum([0; free_paths])

    t_final = collision_times[end]
    push!(collision_times, t_final)  # duplicate -- necessary?


    collision = 1
    time_counter = 1

    positions = [xs[1];]  # initial position at time 0
    times = [0.0;]

    next_output_time = delta_t

    time_left_for_step = delta_t

    # velocity:
    direction = collision_positions[collision+1] - collision_positions[collision]
    direction /= norm(direction)

    while next_output_time < t_final

        if next_output_time < collision_times[collision+1]

            flight_time = next_output_time - collision_times[collision]

            push!(positions, collision_positions[collision] + flight_time * direction)
            push!(times, next_output_time)
            time_counter += 1
            next_output_time = time_counter * delta_t
            # more accurate than adding delta_t if delta_t is not representable

        else

            collision += 1

            direction = collision_positions[collision+1] - collision_positions[collision]
            direction /= norm(direction)


        end
    end

    positions, times



end
