# Try to use Wang-Landau to estimate "density of states" of "energy"=path length

using BilliardModels

function first_free_path(billiard_table, x, v)
    l = Vector2D(0, 0)
    p = ParticleOnLattice(x, v, l)
    num_collisions = 1
    xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, num_collisions)

    free_paths[1]
end


unif_rand(δ) = δ*(rand()-0.5)  # uniform real ∈ [-δ, δ]



function propose(x, y, θ, radius, δ)  # random walk proposal
    x_new = x + unif_rand(δ)
    y_new = y + unif_rand(δ)

     # check outside disc:
    if x*x + y*y < radius^2
        x_new, y_new = x, y   # return to where it was
    else
        # periodize:
        (x_new > 0.5) && (x_new -= 1)
        (x_new < -0.5) && (x_new += 1)
        (y > 0.5) && (y -= 1)
        (y < -0.5) && (y += 1)
    end

    θ += unif_rand(δ)  # in any case change θ

    x, y, θ
end

acceptance_prob(x, y, θ) = 1.0  # standard random walk

function explore_phase_space(radius=0.1, δ=0.1, num_steps=1000)

    billiard_table = Sinai_billiard(radius, true, true)  # periodic in x and y

    xx, vv = initial_condition(billiard_table, -.5, .5, -.5, .5)

    θ = atan2(vv[2], vv[1])
    # this is how I generated the velocity to start with (from an angle) -- it's a bit ridiculous to undo that

    xs = [xx[1]]
    ys = [xx[2]]
    θs = [θ]
    free_paths = [find_free_path(billiard_table, xx, vv)]

    num_free_path_bins = 1000
    free_path_distribution = zeros(num_free_path_bins)

    for i in 1:num_steps
        x, y = xx[1], xx[2]

        x_new, y_new, θ_new = propose(x, y, θ, radius, δ)

        if rand() < acceptance_prob(x, y, θ)
            x, y, θ = x_new, y_new, θ_new
        end

        xx_new = Vector2D(x, y)
        vv = Vector2D(cos(θ), sin(θ))


        free_path = find_free_path(billiard_table, xx_new, vv)

        push!(xs, x)
        push!(ys, y)
        push!(θs, θ)
        push!(free_paths, free_path)

        free_path_bin = ceil(free_path)
        (free_path_bin > num_free_path_bins) && (free_path_bin = num_free_path_bins)

        free_path_distribution[free_path_bin] += 1

        xx = xx_new

    end

    xs, ys, θs, free_paths, free_path_distribution
end



