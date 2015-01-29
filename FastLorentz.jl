
module Lorentz

using Vector2d
using Docile


#using ImmutableArrays
#import Base.convert

export run_Lorentz_gas
export run_many_particles

export Particle, Disc, collision_time, Plane


type Particle
    x::Vector2D
    v::Vector2D
end


immutable Disc
    centre::Vector2D
    radius::Float64
end

@doc """Compute *time of collision* of a particle and a disc,
        assuming the particle starts outside the disc and the
        particle speed is one (using the quadratic formula).

        Returns -1 if no collision.
        """ ->
function collision_time(p::Particle, disc::Disc)

    disp = p.x-disc.centre

    B = dot(p.v, disp)
    C = dot(disp, disp) - disc.radius*disc.radius

    discriminant = B*B - C

    if discriminant < 0
        return -1.
    end

    return -B - sqrt(discriminant)  # NB: this supposes that the particle is *outside* the disc

end

@doc """Compute normal vector to disc boundary, at a point x that must lie *on the boundary*""" ->
function normal(disc::Disc, x)
    return (x - disc.centre) / disc.radius
end


immutable Plane
    c::Vector2D  # point on plane
    normal::Vector2D
end

#Plane(c::Vector{Float64}, normal::Vector{Float64}) =
#    Plane(Vector2D(c), Vector2D(normal))


collision_time(p::Particle, plane::Plane) =
    dot(plane.c - p.x, plane.normal) / dot(p.v, plane.normal)

normal(plane::Plane, x) = plane.normal


function collision(p::Particle, current_cell, boundaries, jump_directions, previous_hit)


    min_collision_time = 100.0
    which_hit = -1

    for i in 1:5
        if i == previous_hit
            continue
        end

        @inbounds collision_t = collision_time(p, boundaries[i])


        if collision_t > 0.0 && collision_t < min_collision_time
            which_hit, min_collision_time = i, collision_t
        end
    end


    x_new = p.x + p.v*min_collision_time
    new_cell = current_cell

    @assert which_hit > 0

    if which_hit == 1  # reflect off disc
        @inbounds n = normal(boundaries[1], x_new)

        v_new = p.v - 2.0*n*dot(n, p.v)  # reflejar solo si es disco; si es plano, pasar a traves

        speed = sqrt(dot(v_new, v_new))
        v_new /= speed

    else
        # hit plane, so periodise
        @inbounds jump = jump_directions[which_hit]
        # print "jump = ", jump

        new_cell += jump
        x_new -= jump  # this is supposing that we are in exactly unit cell

        v_new = p.v

        which_hit += 2
        if which_hit > 5
            which_hit -= 4
        end
    end

    return x_new, v_new, new_cell, min_collision_time, which_hit
end


function initial_condition(radius)

    x, y = 0.0, 0.0

    while true
        x, y = rand(2) .- 0.5

        if (x*x + y*y) >= radius^2
            break
        end
    end

    # velocity in 2D:
    theta = rand() * 2*pi

    xx = Vector2D(x, y)
    vv = Vector2D(cos(theta), sin(theta))

    return xx, vv
end


# fast version:
function collision_time_with_disc(x, v, r)
    B = dot(v, x);
    C = dot(x, x) - r*r;

    discriminant = B*B - C;

    if (discriminant < 0)
        return -1.0;

    else
        return -B - sqrt(discriminant);
    end
end

intersection_horiz_boundary(xx, vv, y) =
    @inbounds return (y - xx[2]) / vv[2]

intersection_vert_boundary(xx, vv, x) =
    @inbounds return (x - xx[1]) / vv[1]


function new_collision(xx, vv, r, cell, previous_hit)

    min_collision_time = 1000.0

    collision_times =
        ( collision_time_with_disc(xx, vv, r),
            intersection_vert_boundary(xx, vv, -0.5),
            intersection_horiz_boundary(xx, vv, 0.5),
            intersection_vert_boundary(xx, vv, 0.5),
            intersection_horiz_boundary(xx, vv, -0.5)
        )

    object_hit = -1

    for i in 1:5  # num of collision times

        if i==previous_hit
            continue
        end

        @inbounds coll_time = collision_times[i]

        if coll_time > 0.0 && coll_time < min_collision_time

            min_collision_time = coll_time
            object_hit = i
        end
    end

    #@assert object_hit > 0

    #println(min_collision_time)

    xx += vv*min_collision_time

    @inbounds begin
        x, y = xx[1], xx[2]
        x_cell, y_cell = cell[1], cell[2]
    end

    if object_hit == 1
        n = xx / r  # vector normal

        vv -= 2. * dot(vv, n) * n
        speed = sqrt(dot(vv, vv))
        vv /= speed

    elseif object_hit == 2
        x = 0.5
        x_cell -= 1
        object_hit = 4

    elseif object_hit == 3
        y = -0.5
        y_cell += 1
        object_hit = 5

    elseif object_hit == 4
        x = -0.5
        x_cell += 1
        object_hit = 2

    elseif object_hit == 5
        y = 0.5
        y_cell -= 1
        object_hit = 3
    end

    cell = Vector2D(x_cell, y_cell)
    xx = Vector2D(x, y)

    return xx, vv, cell, min_collision_time, object_hit


end






### Main program

function make_boundaries(radius)

    boundaries = {}

    push!(boundaries, Disc([0., 0.], radius) )

    push!(boundaries, Plane([-0.5, -0.5], [0., -1.]) )
    push!(boundaries, Plane([0.5, -0.5], [1., 0.]) )
    push!(boundaries, Plane([0.5, 0.5], [0., 1.]) )
    push!(boundaries, Plane([-0.5, 0.5], [-1., 0.]) )

    jump_directions = [Vector2D(0.,0.)]  # dummy

    for boundary in boundaries[2:end]
        push!(jump_directions, normal(boundary, 0))
    end

    println("# Boundaries: ", boundaries)
    println("# Jump directions: ", jump_directions)

    return boundaries, jump_directions
end


function Lorentz_gas(radius, boundaries, jump_directions, t_final)

    current_cell = Vector2D(0., 0.)

    #collision_times = -ones(5)


    x0, v0 = initial_condition(radius)
    p = Particle(x0, v0)

    total_time = 0.0

    xs = [x0[1]]
    ys = [x0[2]]

    which_hit = -1
    time_since_disc_collision = 0.0
    disc_collision_times = Float64[]

    output_time = 10.
    next_output_time = output_time
    num_data = int(round(t_final / output_time))

    println("num_data = ", num_data)


    step = 0

    positions = [x0]
    times = [0.0]

    #while total_time <= t_final
    while step <= num_data+1
        x_new, v_new, new_cell, collision_t, which_hit =
            #collision(p, current_cell, boundaries,
            #        jump_directions, which_hit)
        new_collision(p.x, p.v, radius, current_cell, which_hit)

        next_time = total_time + collision_t

        #println("next_time = ", next_time, "; next_output_time = ", next_output_time)

        #while (next_time > next_output_time)
        if (next_time >= next_output_time)
            # not necessary to use a while

            new_position = p.x + (next_output_time - total_time) * p.v
            new_position += new_cell

            #println(next_output_time, "\t", new_position[1], "\t", new_position[2])



            push!(positions, new_position)
            push!(times, next_output_time)

            next_output_time += output_time
            step += 1

            if step == num_data
                break
            end


        end


        total_time = next_time

        if total_time + collision_t < t_final
            total_time += collision_t
            current_cell = new_cell

            p.x, p.v = x_new, v_new


            #println(x_new[1], "\t", x_new[2], "\t", v_new[1], "\t", v_new[2])

            time_since_disc_collision += collision_t

            if which_hit == 1  # disc

                x_new += new_cell

                push!(xs, x_new[1])
                push!(ys, x_new[2])

                push!(disc_collision_times, time_since_disc_collision)
                time_since_disc_collision = 0.0
            end


        else
            collision_t = t_final - total_time

            x_new = p.x + p.v*collision_t
            x_new += current_cell

            println("# Position at time ", t_final, " = ", x_new)
            println("# Cell: ", current_cell)
            println("# Distance from origin = ", (dot(x_new, x_new))^0.5)
#            println("# Distance squared from origin = ", dot(x_new, x_new))

            break
        end

    end

    return xs, ys, disc_collision_times, times, positions
end

function many_particles(N, radius, boundaries, jump_directions, t_final)

    xs, ys, disc_collision_times, times, positions = Lorentz_gas(radius, boundaries, jump_directions, t_final)

 println("times = ", times)
    println(length(positions))

    all_positions = { { positions[i] } for i in 1:length(positions) }


    for n in 1:N-1
        xs, ys, disc_collision_times, times, positions = Lorentz_gas(radius, boundaries, jump_directions, t_final)

         println("times = ", times)
    println("all_positions = ", all_positions)
    println(length(positions))
    println(length(all_positions))

        for i in 1:length(positions)
            push!(all_positions[i], positions[i])
        end
    end

    all_positions
end

# function continuous_time_trajectory()

#     p = Particle(initial_condition()...)

#     t = 0.0;
#     next_output_time = output_time;

#     step = 1;  /// location in data array

#     while (t < max_time) {

#         next_t = t + min_collision_time;

#         while (next_t > next_output_time) {
#             /// use a while in case have very long path

#             new_position = p.position + (next_output_time - t)*p.velocity;

#             new_position += Vec(x_cell, y_cell);

#             displacement = new_position - p.initial_position;
#             disp_norm = dot(displacement, displacement);

#             /// |Delta x|
#             basic_moment = sqrt(sqrt(disp_norm)); /// |Delta x|^{1/2}
#             /// the basic unit is |Delta x|^{1/2}

#             moment = basic_moment;
#             /// moment is multiplied over and over by basic_moment to give the powers

#             for (int m=1; m<=num_moments; m++) {
#                 moments_data[m][step] += moment;
#                 moment *= basic_moment;
#             }

#             vel_vel_correlation_data[step] += dot(p.velocity, p.initial_velocity);
#             vel_pos_correlation_data[step] += dot(displacement, p.initial_velocity);
#             current_vel_pos_correlation_data[step] += dot(displacement, p.velocity);

#             step++;
#             next_output_time += output_time;

#             if (step > num_data) return;

#         }

#         implementCollision();

#         t = next_t;

#     //cout << "After collision:\n";
#     //cout << "Pos: ";

#         findCollision();

#     }


# }


function run_Lorentz_gas(radius, max_time)

    println("# Using radius: ", radius)

    boundaries, jump_directions = make_boundaries(radius)
    @time xs, ys = Lorentz_gas(radius, boundaries, jump_directions, max_time)
end


function run_many_particles(N, radius, max_time)

    println("# Using radius: ", radius)
    println("# Using $N particles")

    boundaries, jump_directions = make_boundaries(radius)
    @time all_positions = many_particles(N, radius, boundaries, jump_directions, max_time)
end

if length(ARGS) > 1
    radius = float(ARGS[1])
    run_Lorentz_gas(radius)
end


# using PyPlot
# plot(xs, ys)

end


function elastic_collision(obstacle, particle)
end

function periodic_boundary(obstacle, particle)  # needs *more* information and returns more information
end

# There's a fundaental difference betweeen "nromal"elastic objects and periodic boundaries
# There could also be other types of boundary, e.g. ones where crossing it randomizes the position of the disc in the new cell or just updates which cell the particle
# lieves in.

# E.g. think of a sequence of cells with pre-randomised (quenched) disc positions.
