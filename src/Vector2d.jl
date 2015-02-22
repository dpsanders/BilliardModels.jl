#module Vector2d

importall Base
export Vector2D, convert, dot, getindex

@doc """A simple type for a 2-component vector.
Previous tests seemed to show that this version was faster than
e.g. ImmutableArrays. But these performance tests should be redone.

This should presumably be a subtype of an abstract array type to avoid
so much rewriting of preexisting methods. """ ->
immutable Vector2D{T} # <: AbstractArray{T,1}
	x::T
	y::T
end

Vector2D{T}(v::Vector{T}) = Vector2D{T}(v[1], v[2])
#Vector2D{T}(v::Vector) = Vector2D{T}(v[1], v[2])
#Vector2D(v::Vector) = Vector2D(v[1], v[2])

+(v::Vector2D, w::Vector2D) = Vector2D(v.x+w.x, v.y+w.y)
+{T}(v::Vector2D{T}, w::Vector2D{T}) = Vector2D{T}(v.x+w.x, v.y+w.y)
#+{S,T}(v::Vector2D{S}, w::Vector2D{T}) = Vector2D{promote_rule(T,S)}(v.x+w.x, v.y+w.y)
-{T}(v::Vector2D{T}, w::Vector2D{T}) = Vector2D{T}(v.x-w.x, v.y-w.y)
*{T}(v::Vector2D{T}, lamb::Number) = Vector2D{T}(v.x*lamb, v.y*lamb)
*{T}(lamb::Number, v::Vector2D{T}) = v*lamb
/{T}(v::Vector2D{T}, lamb::Number) = Vector2D{T}(v.x/lamb, v.y/lamb)

=={T}(v::Vector2D{T}, w::Vector2D{T}) = v.x==w.x && v.y==w.y

size{T}(v::Vector2D{T}) = 2
dot{T}(v::Vector2D{T}, w::Vector2D{T}) = v.x*w.x + v.y*w.y
norm{T}(v::Vector2D{T}) = √(v ⋅ v)

getindex{T}(v::Vector2D{T}, i) = (i==1) ? v.x : v.y

#There are problems with this `convert`
convert{T}(::Type{Vector2D{T}}, v::Vector{T}) = Vector2D{T}(v[1], v[2])
convert(::Type{Vector2D}, v::Vector) = Vector2D(v[1], v[2])



#end
