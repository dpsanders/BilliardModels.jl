# Accumulate statistics (mean, variance etc.) at different times

# Stores statistics at one time

using StreamStats

# Stores statistics at given times:
type StatisticsObject
    data::Dict{Float64, StreamStats.Moments}
end

StatisticsObject() = StatisticsObject(Dict())


# Adds the data x at time t
function add_data!(o::StatisticsObject, t, x)
    if haskey(o.data, t)
        #add_data!(o.data[t], x)
        update!(o.data[t], x)
    else
        o.data[t] = StreamStats.Moments() #StatisticsContainer(x)
        update!(o.data[t], x)
    end
end

function get_statistics(o::StatisticsObject)
    times = sort(collect(keys(o.data)))

    [(t, state(o.data[t])) for t in times]

end


