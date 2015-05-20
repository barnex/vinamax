//This file contains some useful vector operators

package vinamax

import "math"

type Vector [3]float64

//Dot product between two vectors
func (x Vector) Dot(y Vector) float64 {
	return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])
}

//cross product between two vectors
func (x *Vector) Cross(y Vector) Vector {
	return Vector{x[1]*y[2] - x[2]*y[1], y[0]*x[2] - y[2]*x[0], x[0]*y[1] - x[1]*y[0]}
}

// Set norm of a vector to one
// TODO: -> normalize
func norm(x Vector) Vector {
	magnitude := math.Sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
	return x.Times(1. / magnitude)
}

//returns the norm of a vector
// TODO: rm
func size(x Vector) float64 {
	return math.Sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
}

func (x Vector) Norm() float64 {
	return math.Sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
}

//multiply each component of a vector by a float
func (x Vector) Times(i float64) Vector {
	return Vector{x[0] * i, x[1] * i, x[2] * i}
}

// point-wise multiplication of components
func (v Vector) Mul3(a Vector) Vector {
	return Vector{a[0] * v[0], a[1] * v[1], a[2] * v[2]}
}

// Returns (1/a)*v.
func (v Vector) Div(a float64) Vector {
	return v.Times(1 / a)
}

//Add two vectors
func (x Vector) Add(i Vector) Vector {
	x[0] += i[0]
	x[1] += i[1]
	x[2] += i[2]
	return x
}

//cubes a number
func cube(x float64) float64 {
	return x * x * x
}

//squares a number
func sqr(x float64) float64 {
	return x * x
}

// absolute value of all components
func (v Vector) Abs() Vector {
	return Vector{math.Abs(v[0]), math.Abs(v[1]), math.Abs(v[2])}
}

// Returns the uniform norm of v
// (maximum of absolute values of components)
func (v Vector) MaxNorm() float64 {
	x := math.Abs(v[0])
	y := math.Abs(v[1])
	z := math.Abs(v[2])
	return math.Max(math.Max(x, y), z)
}

// Returns a+s*b.
func (a Vector) MAdd(s float64, b Vector) Vector {
	return Vector{a[0] + s*b[0], a[1] + s*b[1], a[2] + s*b[2]}
}

// Returns a-b.
func (a Vector) Sub(b Vector) Vector {
	return Vector{a[0] - b[0], a[1] - b[1], a[2] - b[2]}
}

const (
	X = 0
	Y = 1
	Z = 2
)
