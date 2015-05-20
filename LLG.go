package vinamax

//Calculates the torque working on the magnetisation of a particle
//using the Landau Lifshitz equation
func (p *Particle) tau(temp Vector) Vector {
	pm := &p.m
	mxB := pm.Cross(p.b_eff(temp))
	amxmxB := pm.Cross(mxB).Times(Alpha)
	mxB = mxB.Add(amxmxB)
	return mxB.Times(-gammaoveralpha)
}
