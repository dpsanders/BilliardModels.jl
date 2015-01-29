module Vector2d

importall Base
export Vector2D, convert, dot, getindex

immutable Vector2D{T}
	x::T
	y::T
end

Vector2D{T}(v::Array{T,1}) = Vector2D(v[1], v[2])

+(v::Vector2D, w::Vector2D) = Vector2D(v.x+w.x, v.y+w.y)
-(v::Vector2D, w::Vector2D) = Vector2D(v.x-w.x, v.y-w.y)
*(v::Vector2D, lamb::Number) = Vector2D(v.x*lamb, v.y*lamb)
*(lamb::Number, v::Vector2D) = v*lamb
/(v::Vector2D, lamb::Number) = Vector2D(v.x/lamb, v.y/lamb)

dot(v::Vector2D, w::Vector2D) = v.x*w.x + v.y*w.y

convert{T}(::Type{Vector2D}, v::Array{T, 1}) = Vector2D(v[1], v[2])

getindex(v::Vector2D, i) = (i==1) ? v.x : v.y

end
