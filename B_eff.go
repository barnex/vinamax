package vinamax

import (
	"math"
	"math/rand"
)

var rng = rand.New(rand.NewSource(0))

//Sums the individual fields to the effective field working on the particle
func (p *Particle) b_eff(temp Vector) Vector {
	return p.Bdemag.Add(p.anis().Add(p.zeeman().Add(temp)))
}

//Set the randomseed for the temperature
func Setrandomseed(a int64) {
	randomseedcalled = true
	rng = rand.New(rand.NewSource(a))
}

// factor 4/3pi in "number" because they are spherical
func (p *Particle) calculatetempnumber() {
	p.tempnumber = math.Sqrt((2. * kb * Alpha * Temp) / (gamma0 * p.msat * 4. / 3. * math.Pi * cube(p.r) * Dt))
}

func calculatetempnumbers(lijst []*Particle) {
	for i := range lijst {
		lijst[i].calculatetempnumber()
	}
}

//Calculates the the thermal field B_therm working on a particle
func (p *Particle) temp() Vector {
	B_therm := Vector{0., 0., 0.}
	if Brown {
		if Temp != 0. {
			etax := rng.NormFloat64()
			etay := rng.NormFloat64()
			etaz := rng.NormFloat64()

			B_therm = Vector{etax, etay, etaz}
			B_therm = B_therm.Times(p.tempnumber)
		}
	}
	return B_therm
}

//Calculates the Zeeman field working on a particle
func (p *Particle) zeeman() Vector {
	x, y, z := B_ext(T)
	x2, y2, z2 := B_ext_space(T, p.x, p.y, p.z)
	return Vector{x + x2, y + y2, z + z2}
}

//Calculates the anisotropy field working on a particle
func (p *Particle) anis() Vector {
	//2*Ku1*(m.u)*u/p.msat

	mdotu := p.m.Dot(p.u_anis)
	uniax := p.u_anis.Times(2. * Ku1 * mdotu / p.msat)

	cubic := Vector{0., 0., 0.}
	if Kc1 != 0 {
		c1m := p.m.Dot(p.c1_anis)
		c2m := p.m.Dot(p.c2_anis)
		c3m := p.m.Dot(p.c3_anis)
		firstterm := p.c1_anis.Times(c1m * (c3m*c3m + c2m*c2m))
		secondterm := p.c2_anis.Times(c2m * (c3m*c3m + c1m*c1m))
		thirdterm := p.c3_anis.Times(c3m * (c2m*c2m + c1m*c1m))

		cubic = firstterm.Add(secondterm.Add(thirdterm)).Times(-2. * Kc1 / p.msat)
	}
	return uniax.Add(cubic)
}
