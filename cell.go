package vinamax

import "fmt"

// Cell in the FMM tree
type Cell struct {
	child            [8]*Cell    // octree of child cells
	partner          []*Cell     // I receive field taylor expansions from these cells
	near             []*Cell     // I recieve brute-force field contributions form these cells
	center           Vector      // my position
	size             Vector      // my diameter (x, y, z)
	m                Vector      // sum of child+particle magnetizations
	b0               Vector      // Field in cell center
	dbdx, dbdy, dbdz Vector      // Derivatives for Taylor expansion of field around center
	particle         []*Particle // If I'm a leaf cell, particles inside me (nil otherwise)
}

func (c *Cell) UpdateB(parent *Cell) {
	if c == nil {
		return
	}

	// propagete parent field expansion to this cell,
	// (applies shift to Taylor expansion)
	if parent == nil {
		c.b0 = Vector{0, 0, 0}
		c.dbdx = Vector{0, 0, 0}
		c.dbdy = Vector{0, 0, 0}
		c.dbdz = Vector{0, 0, 0}
	} else {
		sh := c.center.Sub(parent.center)
		c.b0 = parent.b0.MAdd(sh[X], parent.dbdx).MAdd(sh[Y], parent.dbdy).MAdd(sh[Z], parent.dbdz)
		c.dbdx = parent.dbdx
		c.dbdy = parent.dbdy
		c.dbdz = parent.dbdz
	}

	// add expansions of fields of partner sources
	for _, p := range c.partner {
		r := c.center.Sub(p.center)
		if r.Dot(r) == 0 {
			panic("self")
		}
		B := DipoleField(p.m, r)
		c.b0 = c.b0.Add(B)
		//c.bx.add(derivative of B)
		// ..
	}

	// if !leaf:
	// propagete field to children
	for _, ch := range c.child {
		ch.UpdateB(c)
	}

	// if leaf cell:
	for _, dst := range c.particle {

		// start with field from cell's Taylor expansion:
		sh := dst.Center().Sub(c.center)
		dst.Bdemag = c.b0.MAdd(sh[X], c.dbdx).MAdd(sh[Y], c.dbdy).MAdd(sh[Z], c.dbdz)

		// then add, in a brute-force way, the near particle's fields
		for _, n := range c.near {
			for _, src := range n.particle {
				r := dst.Center().Sub(src.Center())
				if r.Dot(r) != 0 { // exclude self
					B := DipoleField(src.m, r)
					dst.Bdemag = dst.Bdemag.Add(B)
				}
			}
		}
	}

}

// recursively update this cell's m as the sum
// of its children's m.
func (c *Cell) UpdateM() {
	c.m = Vector{0, 0, 0}

	// leaf node: sum particle m's.
	if c.particle != nil {
		for _, p := range c.particle {
			c.m = c.m.Add(p.m)
		}
		return
	}

	// non-leaf: update children, then add to me.
	for _, ch := range c.child {
		if ch != nil {
			ch.UpdateM()
			c.m = c.m.Add(ch.m)
		}
	}
}

// Recursively find partner cells from a list candidates.
// To be called on the root cell with level[0] as candidates.
// Partners are selected from a cell's own level, so they have the same size
// (otherwise it's just anarchy!)
// TODO: when we have unit tests, it can be optimized not to do so many allocations
func (c *Cell) FindPartners(candidates []*Cell) {

	// select partners based on proximity,
	// the rest goes into "near":
	var near []*Cell
	for _, cand := range candidates {
		if IsFar(c, cand) {
			c.partner = append(c.partner, cand)
			//totalPartners++
		} else {
			near = append(near, cand)
		}
	}

	// leaf cell uses near cells for brute-force,
	// others recursively pass near cells children as new candidates
	if c.IsLeaf() {
		c.near = near
		//totalNear += len(near)
	} else {
		// children of my near cells become parter candidates
		newCand := make([]*Cell, 0, 8*len(near))
		for _, n := range near {
			newCand = append(newCand, n.child[:]...)
		}

		// recursively find partners in new candidates
		for _, c := range c.child {
			c.FindPartners(newCand)
		}
	}
}

// Is this cell a leaf cell?
func (c *Cell) IsLeaf() bool {
	for _, c := range c.child {
		if c != nil {
			return false
		}
	}
	return true
}

// Are the cells considered far separated?
func IsFar(a, b *Cell) bool {
	// TODO: this is more or less a touch criterion: improve!
	dist := a.center.Sub(b.center).Norm()
	return dist > 1.1*a.size.Norm()
}

// Create child cells to reach nLevels of levels and add to global level array.
// nLevels == 1 stops creating children (we always already have at least 1 level),
// but still adds the cell to the global level array.
func (c *Cell) Divide(nLevels int) {

	// add to global level array
	myLevel := len(level) - nLevels
	level[myLevel] = append(level[myLevel], c)
	totalCells++

	if nLevels == 1 {
		return
	}

	// create children
	for i := range c.child {
		newSize := c.size.Div(2)
		newCenter := c.center.Add(direction[i].Mul3(newSize.Div(2)))
		c.child[i] = &Cell{center: newCenter, size: newSize}
	}

	// recursively go further
	for _, c := range c.child {
		c.Divide(nLevels - 1)
	}
}

func (c *Cell) String() string {
	if c == nil {
		return "nil"
	} else {
		typ := "node"
		if c.IsLeaf() {
			typ = "leaf"
		}
		return fmt.Sprint(typ, "@", c.center, len(c.partner), "partners, ", len(c.near), "near")
	}
}

// unit vectors for left-back-bottom, left-back-top, ...
var direction = [8]Vector{
	{-1, -1, -1}, {-1, -1, +1}, {-1, +1, -1}, {-1, +1, +1},
	{+1, -1, -1}, {+1, -1, +1}, {+1, +1, -1}, {+1, +1, +1}}
