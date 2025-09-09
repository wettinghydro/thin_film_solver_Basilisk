/**
Header file used for the definition of the functions used to define the
disjoining pressure related functions*/

double disp (const double h)
{
	const double hi       = 1.0 / h;
	const double eByH     = precFilm * hi;
	const double fracMpow = pow(eByH, m);
	const double fracMn   = pow(eByH, n - m);

	return c1 * fracMpow * (1.0 - fracMn);
}

/**
Then, we define the disjoining pressure derivative w.r.t. the film height
$h$, i.e.

$$ \Pi^\prime(h) \equiv \left[\Pi(h)\right]_h
= \mathscr{c_1} \left(\frac{n}{h}\left( \frac{\varepsilon}{h}\right)^n
- \frac{m}{h}  \left(\frac{\varepsilon}{h}\right)^m \right) $$

or

$$ \left[\Pi(h)\right]_h = \frac{\mathscr{c_1}}{h}
\left(\frac{\varepsilon}{h}\right)^m
\left(n \left(\frac{\varepsilon}{h}\right)^{n-m}  -m \right) $$
with $\mathscr{c_1}$ = ```c1``` (see [Global Constants](#physical-constants)).
*/
double dpdh (const double hf)
{
	const double hi       = 1.0 / hf;
	const double eByH     = precFilm * hi;
	const double fracMpow = pow(eByH, m);
	const double fracMn   = pow(eByH, n - m);

	return c1 * hi * fracMpow * (n * fracMn - m);
}

/**
Accordingly, the expression of the $2^{nd}$ derivative w.r.t. $h$ is computed,
$$ \Pi^{\prime\prime}(h) \equiv \left[\Pi(h)\right]_{hh} =...=
\frac{\mathscr{c_1}}{h^2} \left(\frac{\varepsilon}{h}\right)^m
\left(m(m+1) - n (n+1) \left(\frac{\varepsilon}{h}\right)^{n-m} \right) $$
*/
double dpdh_der (const double hf)
{
	const double hi       = 1.0 / hf;
	const double eByH     = precFilm * hi;
	const double fracMpow = pow(eByH, m);
	const double fracMn   = pow(eByH, n - m);

	return c1 * sq(hi) * fracMpow * (m * (m+1.0) - n * (n+1.0) * fracMn);
}

/**
The next two functions are used to compute
$\left[h^3\left[\Pi(h)\right]_h\right]_h$
*/
//computes the derivative of a power term
double pow_der (const double x, const double exponent)
{
	return exponent * pow(x, exponent-1.0);
}

//computes the by parts derivative of [h^3 * dpdh]' at using h (here at faces)
double pressureDerTermByParts (const double hf)
{
	return pow_der(hf, 3.0) * dpdh(hf) + pow(hf, 3.0) * dpdh_der(hf);
}

//computes the by parts derivative of [h^3 * ph]' at using h (here at faces)
double pressureByParts(const double hf, const double h_grad)
{
	return pow_der(hf, 3.0) * disp(hf) + pow(hf, 3.0) * dpdh(hf) * h_grad;
}

/**
## Substrate Heterogeneity
The substrate heterogeneity is defined by prescribing the hamaker constant
accros the computational domain.

Then, the face substrate heterogeneities ```hetf``` are used to compute
the $\tilde\Pi(\mathbf x, h)$ terms.
*/
double ham_disp (const double het, const double h)
{
	return het * disp(h);
}

/**
$$ \mathbf \nabla h^3 \tilde{\Pi}(\mathbf{x}, h)
= \mathbf \nabla \left( h^3 \mathcal{H}(\mathbf{x}) \Pi(h) \right)
= \mathbf \nabla \left( \mathcal{H}(\mathbf{x})\right) h^3 \Pi(h)
+ \mathbf \nabla \left( h^3 \Pi(h) \right) \mathcal{H} (\mathbf{x})
$$
*/
double het_disPressureDer(double hetf, double hf,
		double h_grad, double het_grad)
{
	return pressureDerTermByParts(hf) * hetf * h_grad
		    + pressureByParts(hf, h_grad) * het_grad;
}

