
# FOR PERIODISING ON A LATTICE:
#  else
#         # hit plane, so periodise
#         @inbounds jump = jump_directions[which_hit]
#         # print "jump = ", jump

#         new_cell += jump
#         x_new -= jump  # this is supposing that we are in exactly unit cell

#         v_new = p.v

#         which_hit += 2
#         if which_hit > 5
#             which_hit -= 4
#         end
#     end


function Lorentz_gas(init_cond, table, t_final)


    x0, v0 = initial_condition(radius)  # CHANGE TO MAKE MORE GENERAL FOR GENERAL TABLE
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
            #collision(p, current_cell, obstacles,
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


function run_Lorentz_gas(radius, max_time)

    println("# Using radius: ", radius)

    table = Sinai_billiard(radius)
    @time xs, ys = Lorentz_gas(radius, obstacles, jump_directions, max_time)
end


function run_many_particles(N, radius, max_time)

    println("# Using radius: ", radius)
    println("# Using $N particles")

    obstacles, jump_directions = make_obstacles(radius)
    @time all_positions = many_particles(N, radius, obstacles, jump_directions, max_time)
end

# if length(ARGS) > 1
#     radius = float(ARGS[1])
#     run_Lorentz_gas(radius)
# end


# using PyPlot
# plot(xs, ys)

#end
