/**
# Precursor Film Model (PFM) for the simulation of thin fluid flows under the Long-Wave (LW) Approximation
This work extends work done by [Dallaston et al.](#references) (see
[Singularities in the thin 2D film equation](http://basilisk.fr/sandbox/lopez/layer.c))
and [Lam et al.](#references) on the PFM.

Within the LW limit, thin film flow can be described via a fourth-order
non-linear parabolic differential equation. Considering a Coordinate system
$O(\mathrm{x_1}, \mathrm{x_2})$, the said PFM is written as:

$$ h_t + \mathbf\nabla \cdot \left( \mathbf \nabla
F(h, \mathbf{\nabla} h, \mathbf{\nabla} \Delta h) \right)= 0 $$

$$ F(h, \mathbf{\nabla} h, \mathbf{\nabla} \Delta h) =
h^3 \left(\mathbf\nabla\Delta h
- \Pi^\prime(h) \mathbf\nabla h\right) $$
where $h(\mathbf x, t)$ is the droplet thickness,
$\mathbf{\nabla}$ is the gradient operator,
and $\Delta=\mathbf{\nabla}\cdot\mathbf{\nabla}$ is the Laplace operator.

The disjoining pressure term $\Pi(h)$, introduced to treat the stress
singularity at the moving contact line, is given as
$$ \Pi(h)=\frac{(m-1)(n-1)}{2\varepsilon(n-m)}
\left[ \left(\frac{\varepsilon}{h}\right)^m -
\left(\frac{\varepsilon}{h}\right)^n \right] $$

with $\varepsilon << 1$ being the precursor film height, typically in the
$1-100nm$ range, and $n$, $m$, are the mobility and interaction exponents
$(n>m)$, respectively; here, $n=3$, $m=2$ are used.
The first term ($m-$exponent) describes the attractive forces, while the second
one ($n-$exponent) the liquid-solid repulsive ones.

~~~gnuplot Disjoining Pressure $\Pi(h)$
set grid;
set ylabel '{/Symbol P}(h)'; set xlabel 'h/{/Symbol e}'; set log x;
e=0.001; n=3.0; m=2.0; C=(m-1.0)*(n-1.0) / (2.0 * e * (n-m) );
f(x) = 1/x; p(x)  = C * (f(x)**m - f(x)**n);

p[1:10**2] p(x) w l lw 2 not
~~~
As one can see, the disjoining pressure term is zero at the precursor film
$h = \varepsilon$ and as $h\rightarrow\infty$, while in between the
said forces are introduced.

The expression for the pressure field throughout the droplet is given as the
balance of the surface curvature approximation and the disjoining pressure
term, i.e.
$$
    p = -\Delta h + {\Pi}(h)
$$

## Heterogeneity-augmented PFM
The PFM presented, can be extended to include additional effects, such as
substrate heterogeneities, of both chemical and topographic nature, body forces,
evaporative/condensing effects, and account for the presence of electric fields.
Before continuing, we assume that for inclined substrates, this manifests along
the $\mathrm{x_1}$ axis with an inclination angle $\phi$, whereas substrate
vibrations occures at anangle $\psi$ relative to the substrate vertical.

Based on the above, the extended PFM reads:

$$ h_t + \mathbf\nabla \cdot \left[h^3 \left(\mathbf\nabla\Delta h
- \underbrace{\tilde\Pi(\mathbf x, h)}_
{\text{substrate}\atop\text{heterogeneities}}
- \underbrace{w_1\mathbf\nabla h}_{\text{body forces}} \right)\right]
+ \underbrace{\mathbf\nabla \left(h^3w_2\right)\cdot \mathbf{e}_{\mathrm{x_1}}}
_{\text{body forces}}
+ \underbrace{\mathcal J}_{\text{evaporation}} = 0 $$

Substrate **chemical** heterogeneities are introduced via the Hamaker coefficient
$\mathcal{H}(\mathbf x)$ [Brusch et al.](#references), which gives rise to
$$ \tilde\Pi(\mathbf x, h) = \mathcal{H}(\mathbf x) \Pi(h) $$

**Body forcing effects** are included via the following two terms
[Savva and Kalliadasis](#references)
$$ w_1 = B_o \left(\cos \phi + r \gamma(t) \cos \psi \right) $$
$$ w_2 = \frac{B_o}{\tan \alpha_s}\left(\sin \phi+r \gamma(t)\sin\psi\right) $$
with $B_o=\frac{\varrho gR^2}{\sigma}$ being the Bond number, $r$ the virbration
strength, $\gamma(t)$ its acceleration and $\alpha_s$ the contact angle.

**Evaporative/condensing effects** are introduced via the source term $\mathcal J$
that quantifies the evaporation/condensation flux [Pham and Kumar](#references)

$$ \mathcal J = \frac{\mathcal E + \delta(p - p_v)}{\mathcal{K} + h}
=\frac{\mathcal E -\delta(\Delta h-\tilde\Pi(\mathbf x, h) - w_1 h)}
{\mathcal K+h} $$

with $\mathcal E$ being the evaporation number and $\delta$ the liquid
temperature variation dependency (Kelvin effect) and $\mathcal{K}$ the kinetic
parameter.

## Surface Roughness
Considering a substrate with topographic features (c.f. [Gaskell et al.](#references)),
defined by the distribution $s(\mathbf{x})$, the PFM is re-written as

$$ h_t + \mathbf\nabla \cdot \left[h^3 \left(\mathbf\nabla
\Delta {\color{crimson} z} - \tilde\Pi(\mathbf x, h)
- w_1 \mathbf\nabla {\color{crimson}z} \right) \right]
+ \mathbf\nabla \left(h^3w_2\right)\cdot \mathbf{e}_{\mathrm{x_1}}
+ \frac{\mathcal E -\delta(\Delta {\color{crimson} z} -\tilde\Pi(\mathbf x, h)
		-w_1{\color{crimson} z})} {\mathcal K+{\color{crimson} z}} =0$$

where $z = h + s$.
*/

/**
# Linearization of the PFM
The PFM is marched in time using Newton's iterative method to solve the
linearized PDE via a MultiGrid (MG) method.
<span style="color:red"> Write a bit cleaner</span>

The resulting linearized equation, that need be discretized and solved,
takes the form
$$ h_t
+ \mathbf \nabla\cdot ( \hat{\mathbf C} \hat h )
+ 2\ \nabla_{\mathrm{x_1}} (\hat G \hat h)
+ \hat E (\mathcal K + \hat h + \hat z)
= \mathbf\nabla \cdot \left( \hat A \mathbf\nabla (-\Delta z) \right)
+ \mathbf\nabla \cdot \left( \hat B \mathbf\nabla h \right)
+ \mathbf\nabla \cdot \left( \hat B_2 \mathbf\nabla s \right)
+ \mathbf\nabla \cdot \left( \hat{\mathbf C} h  \right)
+ 3\ \nabla_{\mathrm{x_1}} (\hat G h)
+ \hat E h
$$
or by using $\mathscr n$ to denote physical timesteps and $k$ to denote Newton
iterations
$$
\begin{aligned}
&h_t
+ \mathbf \nabla\cdot (\mathbf C_k^{\mathscr n+1} h_k^{\mathscr n+1} )
+ 2\ \nabla_{\mathrm{x_1}} (G_k^{\mathscr n+1} h_k^{\mathscr n+1})
+ E_k^{\mathscr n+1} (\mathcal K + h_k^{\mathscr n+1} + z_k^{\mathscr n+1} )
\\=&
 \mathbf\nabla\cdot\left(A^{\mathscr n+1}_{k} \mathbf\nabla(-\Delta h_{k+1}^{\mathscr n + 1}) \right)
+\mathbf\nabla\cdot\left(A^{\mathscr n+1}_{k} \mathbf\nabla(-\Delta s) \right)
+\mathbf\nabla\cdot\left(B_k^{\mathscr n+1} \mathbf\nabla h_{k+1}^{\mathscr n+1}\right)
+\mathbf\nabla\cdot\left({B_2}_k^{\mathscr n+1} \mathbf\nabla s\right)
+\mathbf\nabla\cdot \left(\mathbf C_k^{\mathscr n+1}h_{k+1}^{\mathscr n+1}\right)
+3\ \nabla_{\mathrm{x_1}} (G_k^{\mathscr n+1}h_{k+1}^{\mathscr n+1})
+ E_k^{\mathscr n+1} h_{k+1}^{\mathscr n+1}
\end{aligned}
$$

with
$$
\begin{aligned}
A & = h^3 \\
B & = h^3 \left( \tilde\Pi^\prime(\mathbf x, h) + w_1 \right) ,\qquad B_2 = h^3 w_1\\
\mathbf C & = 3h^2\mathbf\nabla (-\Delta z )+
\left\{ \left[h^3\tilde\Pi^\prime(\mathbf x, h)\right]_h + 3h^2 w_1 \right\}
\mathbf\nabla h  + 3h^2 w_1 \mathbf\nabla s  \\
G  &= -h^2 w_2\\
E  &= \frac{\mathcal{E} - \delta(\Delta z - \Pi(\mathbf{x}, h) - w_1 z)}{(\mathcal{K} + z)^2}
\end{aligned}
$$
*/

/**
# Basilisk Implementation
*/
#define BGHOSTS 2

#include "grid/multigrid.h"
#include "run.h"
#include "poisson.h"
#include "fractions.h"
#include "view.h"
#include "maxruntime.h"
#include "outputManager.h"

#if 0
#define LOG (fprintf(stderr, "# %s  %d\n",__func__, __LINE__))
#define LOGSTATS (fields_stats())
#else
#define LOG
#define LOGSTATS
#endif

/**
Utility Macro to define scalars that are inserted in lists automatically
*/
#define createList(sc1, sc2, scl, type) type sc1[], *scl = {sc1}

#define createScalars(sc1, sc2, scl) createList(sc1, sc2, scl, scalar)
/**
## Variables
### Droplet thickness and Interface charge density
We define the ```scalar``` fields $h$ and $h_{old}$ to track the evolution of
the film equation.  If $1^{st}$ order time integration is performed
```hold``` is only used to monitor convergence. If ```BDF2``` is defined,
the solution at last iteration is stored in ```hold```.
*/

createScalars(h, q, fList);
createScalars(hold, qold, fListOld);

/**
### Heterogeneity and Roughness scalars
```scalar hamaker``` is used to define the sustrate heterogeneity. Finally, if
surface roughness is modelled, the scalars ```gs``` and ```lapS``` are
also introduced.
*/
scalar hamaker[];
#if ROUGH
face vector gs[];
scalar s[], lapS[];
#endif
/**
# Problem Setup
### Physical Constants
We define some problem specific input parameters, such as the precursor
film height ```precFilm```, and the disjoining pressure mobility and iteraction
exponents ```n, m```.
*/
//double precFilm = 5.e-3, n = 3.0, m = 2.0;
double precFilm = 1.e-2, n = 3.0, m = 2.0;

/**
For convenience, the constant disjoining pressure non--dimensional prefactor
$\mathscr{c_1}$ is stored as
*/
#define c1 ( (m - 1.0) * (n - 1.0) / (2. * precFilm * (n - m)) )

/**
To account for cases of inclined substrates the Bond number, and inclination
and contact angles are defined (in degrees).
*/
#if INCLINED
double bond = 1.00, inclineA_d = 6., contactA_d = 45.0;
#endif

/**
In addition, if the substrate is vibrated, the vibration strength ```rvib```,
acceleration ```gam(t)``` with period ```Tper``` and the angle w.r.t. the
substrate vertical ```psiA_d``` are also defined.
*/
#if FORCING
double rvib = 25.00, psiA_d = 8., Tper = 60.0;
#define gam(t) ( sin(2.0 * pi * t / Tper) )
#endif

/**
Evaporative/condensing effects depend on three non-dimensionalized
parameters, $\mathcal{K}$, $\delta$, and $\mathcal E$ defined as ```kappa,
delta, Efactor```.
*/
#if EVAPORATION
double kappa = 3.e-1, delta = 1.e-2, Efactor = -1.e-1;
#endif

/**
### Numerics Constants
Following, some constants are defined that are associated with the numerical
solution process, such as the timestep bounds ```dt```
$\in [$```minDT, maxDT```$]$ and the total simulation time, ```Tend```.
*/
double maxDT = 5.e-4, minDT = 1.e-8, Tend = 100.00;

/**
In addition, the tolerances of the MG solver ```tolMG```, the Newton scheme
```tolNewton``` and the maximum iteration of the Newton scheme ```nitermax```
are provided.
*/
const double tolMG = 1.e-3, tolNewton = 1.e-10;
const int nitermax = 12;

/**
Next, the intervals for printing convergence statistics and the number of
visual outputs to be generated (flowfields and png figures) are defined.
Also, the output directories relative paths are defined.
*/
const int printStatsEvery = 100, fieldSolToSave = 150, noOfFigsToSave= 250;
char dirFigures[]    = "png/";
char dirFlowFields[] = "flows/";
char dirGnu[]        = "gnuData/";

/**
Mesh resolution is bounded by a maximum level.
*/
// -fixme: obsolete without Quadtrees (only Max is taken into account)
int LEVEL_IN  = 10; // initial resolution (optional input arg)
int LEVEL_MAX = 10; // maximum resolution

/**
During the solution process, we track the sum of MG iterations within one
Newton iteration.
*/
int mgIterSum = 0;

/**
## Auxilliary Initialization Functions
The initialization functions used to setup the studied problems can be found
in the ```initializationFunctions.h``` header file.
*/

#include "initializationFunctions.h"

/**
## Disjoining Pressure Computation Functions
Here, we define some functions that are used to compute the disjoining
pressure terms that arise for the linearization and discretization of the
equation (see [Linearization](#linearization-of-the-tfm)).
First, a function to compute the disjoining pressure term is defined as
*/

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
double het_disPressureDer(double hetf, double hf, double h_grad, double het_grad)
{
	return pressureDerTermByParts(hf) * hetf * h_grad
		    + pressureByParts(hf, h_grad) * het_grad;
}

/**
## Vibration Forces
The contributions due to body forcing are computed via the following macros.
*/
#if INCLINED
#define w1I ( bond*cos(inclineA_d*pi/180.) )
#define w2I ( bond*sin(inclineA_d*pi/180.)/tan(contactA_d*pi/180.) )
#endif

#if FORCING
#define fprefac ( bond * rvib * gam(t) ) // forcing variables
#define w1F ( fprefac * cos(psiA_d*pi/180.) )
#define w2F ( fprefac * sin(psiA_d*pi/180.)/tan(contactA_d*pi/180.) )
#endif

#if (INCLINED && FORCING)
#define grav1 ( w1I + w1F )
#define grav2 ( w2I + w2F )
#elif INCLINED
#define grav1 ( w1I )
#define grav2 ( w2I )
#elif FORCING
#define grav1 ( w1F )
#define grav2 ( w2F )
#endif

/**
## Evaporation / Condensation

Evaporative/condensing effects are incorporated as an additional source term
to the governing equations, extending relevant literature [V.S. Azaev](#references).

Here, first the macro ```equilPressure``` is used to compute the equilibrium
pressure depending on the additional effects included
($p-p_v=-\Delta z + \tilde{\Pi}(\mathbf{x}, h) + w_1 z$).
The macros ```evaporationSource``` and ```evaporationTerm``` compute $\mathcal J$
and $E$, respectively (see [Linearization](#linearization-of-the-tfm)).
*/

#if EVAPORATION
#if ROUGH
#define z ( h[] + s[] )
#define lapZ ( lapH[] + lapS[] )
#else
#define z ( h[] )
#define lapZ ( lapH[] )
#endif //ROUGH

#if ( INCLINED || FORCING )
#define T ( -grav1 * z )
#else
#define T ( 0 )
#endif //( INCLINED || FORCING )

// Pressure Term based on heterogeneous Effects
#define equilPressure ( -lapZ - ham_disp(hamaker[], h[]) + T )

// J-term
#define evaporationSource ( (Efactor - delta *  equilPressure ) / (kappa + z) )

// E-term
#define evaporationTerm ( evaporationSource / (kappa + z) )
#endif //EVAPORATION

/**
Based on the defined physical constants, three regimes are observed based on the
value of $\alpha(\mathbf{x})= \frac{\varepsilon \mathcal{E}}
{\delta\mathcal H(\mathbf{x})}$. The three regimes describe the evaporation,
weak condensation (dropwise) or strong condensation (filmwise).

~~~gnuplot Evaporation / Condensation Regimes
reset;
set grid;
e = 1.0; p(x) = (e/x)**3-(e/x)**2; f(x, a) = -a + p(x*e)
set xlabel 'h/{/Symbol e}';
set title 'F= -a + {/Symbol P}(h)/{/Symbol e}';
p[0:4][-2:4] f(x, 1    ) w l lw 2 t 'Evaporation - a = 1', \
             f(x, -0.08) w l lw 2 t 'Weak Condensation - a = -0.08', \
             f(x, -1   ) w l lw 2 t 'Strong Condensation - a = -1', \
                  0 w l ls 1 lc 'black' lw 1 not
~~~
*/

/**
## Surface Roughness
Since surface roughness is predefined at each iteration, we opt to define $s(\mathbf{x})$
and compute and store $\mathbf{\nabla} s$ and $-\Delta s$. Their contributions
are then computed and added on-the-fly.
*/

/**
# PFM solver
Here, we define the data structure to store the data needed to be passed
to the MG solver.
 */
struct ThinData
{
	//mandatory
	scalar *f;
	double dt;

	scalar *fold;  //only for 2nd Order Time Integration
	double  dtold; //only for 2nd Order Time Integration

	//optional
	scalar *lambda;    // diagonal term
	double  timeDiag;  // diagonal constant term
	scalar *r;				 // rhs
	double  tolerance; // MG tolerance
/**
The discretized equation coefficients are stored as in the struct to be passed
to the MG solver (see the [linearized equation](#linearization-of-the-tfm))
*/
	(const) face vector  A;   // height Laplacian coefficient
	(const) face vector *Bl;  // height gradient  coefficient
	(const) face vector *Cl;  // height term coefficient (Jacobian)
};

/**
## PFM Multigrid residual and relaxation functions
Definitions of the MG residual and relaxation functions are provided in the
header file ```filmEquationMG.h```.
 */
#include "filmEquationMG.h"

/**
## PFM Newton Scheme
We initiate the Newton scheme and compute the linearized terms in each
iteration. We limit the Newton scheme to a specific tolerance or ```nitermax```
iterations.  Simultaneously, we track the number of total Multigrid
solutions as an indicator to change the timestep size.
*/
int newton_solve (struct ThinData p)
{
	int niter = 0;
	double err = 0.0;

	mgIterSum = 0;

/**
The procedure is as follows. First, recover the *constant* (within the Newton
Iterations) terms from the data structure
*/
	scalar *rl = p.r;
	scalar *fl = p.f;

/**
and define the ```scalar``` fields $\mathcal{RHS}$ and last Newton iteration
solution
*/
	createScalars(b, bq, bl);
	createScalars(fPr, fqPr, fPrl);

	createScalars(lambdah, lambdaQ, lambdal);

	scalar f = fl[0];

	face vector gradHet[]; // could pre-compute
	foreach_face()
		gradHet.x[] = face_gradient_x(hamaker, 0);

	restriction({hamaker, gradHet});

	do
	{

/**
compute $\mathbf \nabla h$ and $-\Delta h$
*/
		face vector gradH[];
		scalar lapH[];
		computeGradAndNegLaplace(f, gradH, lapH);

/**
and ```face vector``` $A$, $B$, $\mathbf C$ and $\mathbf{gC} = \mathbf C h_f$
*/
		face vector A[], gC[];

		createList(B, Bq, Bl, face vector);
		createList(C, Cq, Cl, face vector);
		foreach_face()
		{
			const double hf   = face_value(f, 0);
			const double hetf = face_value(hamaker, 0);

			A.x[]  = cube(hf);
			B.x[]  = A.x[] * hetf * dpdh(hf);

			C.x[]  = het_disPressureDer(hetf, hf, gradH.x[], gradHet.x[]);
			C.x[] += pow_der(hf, 3.0) * face_gradient_x(lapH, 0);

			gC.x[] = C.x[] * hf  - pow(hf, 3.0) * disp(hf) * gradHet.x[];
		}

/**
Here, the $\mathcal{RHS}$ ```b``` is computed with the additional terms arising from
the linearization
*/
		foreach()
		{
			scalar b, r;
			for (b, r in bl, rl) b[] = r[]; // update time history

			scalar fp, f;
			for (fp, f in fPrl, fl) fp[] = f[]; // update last iteration field

			for (scalar lambda in lambdal) lambda[] = p.timeDiag; // revert diag

			b = bl[0]; // h-res
			foreach_dimension()
			{
				b[] += dflux(gC);
			}
		}

/**
The struct holding the coefficients is updated. This also allows storing the
and utilizing them in the following functions.
*/
		p.lambda = lambdal;
		p.A  = A;
		p.Bl = Bl;
		p.Cl = Cl;

/**
Compute inclination terms, i.e. ```face vector``` $Gr, gGr$
*/
#if INCLINED
		newton_Contributions_ITFM(p, gradH, bl);
#endif

#if ROUGH
		newton_Contributions_Rough(p, bl);
#endif

/**
Compute the evaporation source term, i.e. ```scalar``` $E$
*/
#if EVAPORATION
		newton_Contributions_ETFM(p, lapH, bl);
#endif

/**
The coefficients are restricted to all grid levels. We first group them in a
scalar list to avoid multiple calls to the restriction function.
*/
		scalar *coefs = NULL;
		for (scalar sc in p.lambda)
			coefs = list_add(coefs, sc);

		for (scalar sc in (scalar *)p.Bl)
			coefs = list_add(coefs, sc);

		for (scalar sc in (scalar *)p.Cl)
			coefs = list_add(coefs, sc);

		for (scalar sc in (scalar *){p.A})
		coefs = list_add(coefs, sc);

		restriction(coefs);

		free(coefs);
/**
Now, the MG solver is called that updated the ```fl``` values
*/
		mgstats s = mg_solve(fl, bl, residual_thin, relax_thin, &p);
/**
After the MG cycle, the maximum absolute difference between
$h^{\mathscr n + 1}_k$ and $h^{\mathscr n + 1}_{k+1}$ is computed to check
the convergence of the Newton scheme.
*/
		err = 0.0;
		foreach( reduction(max:err) )
		{
			scalar fPr, f;
			for (fPr, f in fPrl, fl)
			{
				err = max(err, fabs(fPr[] - f[]));
			}

		}

		mgIterSum += s.i; // add internal MG iterations
	} while (err > tolNewton && niter++ < nitermax-1);

/**
After the Newton iterations, check if ```nitermax``` is reached to warn
the user.
*/

	if (niter == nitermax)
		fprintf(ferr, "# niter: %d/%d, err: %g \n", niter, nitermax, err);

	return niter;
}
/**
## The PFM solver
Sets up the PFM solver. Computes the constant coefficients
and initiates the Newton iterations to treat the non-linear terms.
In each Newton step the non--linear terms are recomputed and passed to
```ThinData``` struct.
*/
trace
int thinFilmSolver (struct ThinData p)
{
	int newtonIterations = 0;

	if (p.dt == 0)
	{
		return newtonIterations;
	}

	const double dTol = TOLERANCE;

	if (p.tolerance)
	{
		TOLERANCE = p.tolerance;
	}

/**
  Obtain the current solution and timestep coefficient
*/
	const double invDt = -1.0 / p.dt;

/**
and define the ```scalar``` fileds to store the flow solution history
```f_his``` and the linear(constant) diagonal term ``` timeDiag``` that carries
time information, computed prior to the Newton scheme
 */

	createScalars(f_his, fq_his, fl_his);

	double timeDiag = 0.0;

/**
If ```p.fold``` is given, one can perform a second--order time integration
using a variable timestep BDF2 (a-BDF2), i.e.
$$ \frac{1+2w_{\mathscr n}}{1+w_{\mathscr n}} h^{\mathscr n+1}
- \frac{(1+w_{\mathscr n})^2}{1+w_{\mathscr n}} h^{\mathscr n}
+ \frac{w_{\mathscr n}^2}{1+w_{\mathscr n}} h^{\mathscr n-1}
= dt^{\mathscr n} \mathcal{RHS}^{\mathscr n+1} $$
where $w_{\mathscr n} = \frac{dt^{\mathscr n}}{dt^{\mathscr n-1}}$, and
$\mathcal{RHS}$ incorporates the *linearized* spatial terms.

This is re-written as
$$ \mathcal{RHS}^{\mathscr n+1} - \frac{\alpha}{dt^{\mathscr n}}
h^{\mathscr n+1} = -\frac{1}{dt^{\mathscr n}} \mathcal F
(h^{\mathscr n},h^{\mathscr n-1}) $$

where $$ \mathcal F(h^{\mathscr n},h^{\mathscr n-1})= \beta h^{\mathscr n}
- \gamma h^{\mathscr n-1}$$ and the a-BDF2 coefficients are:
$$ \alpha = \frac{1 + 2w_{\mathscr n}}{1+w_{\mathscr n}},\quad
\beta  = 1+w_{\mathscr n},\quad
 \gamma = \frac{w_{\mathscr n}^2}{1+w_{\mathscr n}}$$
*/

	if (p.fold)
	{
		if (!p.dtold)
		{
		 	fprintf(ferr, "Second Order requires dtold as well\n");
			exit(1);
		}

		double w = p.dt / p.dtold;

		double alpha = (2.0 * w + 1.0) / (w + 1.0);
		double beta  = w + 1.0;
		double gamma = sq(w) / (w + 1.0);

		timeDiag = alpha * invDt;
		foreach()
		{
			scalar f_his, f, fold;
			for (f_his, f, fold in fl_his, p.f, p.fold)
			{
				f_his[]  = (beta * f[] - gamma * fold[]) * invDt;
			}
		}
	}
	else /** Otherwise, perform BDF1.  */
	{
		timeDiag = invDt;

		foreach()
		{
			scalar f_his, f;
			for (f_his, f in fl_his, p.f)
			{
				f_his[]  = f[] * invDt;
			}
		}
	}

/**
Store the computed ```scalar``` fields to the ```ThinData struct```
*/
	p.r  = fl_his;
	p.timeDiag = timeDiag;

/**
and pass it to the Newton scheme
*/
	newtonIterations = newton_solve(p);

/**
Finally, ```TOLERANCE``` is reset if it was changed.
*/
	if (p.tolerance)
	{
		TOLERANCE = dTol;
	}

	return newtonIterations;
}

/**
# Simulations Setup
At the beginning, the boundary conditions, domain size, and computational mesh
are initialized.
*/

#include "inputReader.h"

int main (int argc, char *argv[])
{
	//check for maximum run time
	maxruntime(&argc, argv);
	//check for input file with extension "*.dat"
	setupFromFile(&argc, argv);

	//BC
	foreach_dimension()
		periodic(right);

	//domain size
	L0 = 10.0;

	//center computational domain
	origin( -4.0, -6.0); 

	init_grid( 1<<LEVEL_MAX );

/**
We refine and prolongate linearly to ensure conservation of $h$ during
refinement.
*/
#if TREE
	h.gradient = minmod2;
	h.refine = h.prolongation = refine_linear;
	hamaker.gradient = minmod2;
	hamaker.refine = hamaker.prolongation = refine_linear;
#endif //!TREE

/**
![Comparison of linear and bilinear refinement on a coarse mesh](conservation.png "Comparison of linear and bilinear refinement on a coarse mesh")
*/

	run();
}

/**
We initialize either from a previous solution, or using a combination of the
[initialization](#auxilliary-initialization-functions) functions.
*/
event init (t=0)
{
	// flowfields
	createDirectory(dirFlowFields);
	// animations & pngs
	createDirectory(dirFigures);
	// gnp data
	createDirectory(dirGnu);

	if ( !restore(file="restart") )
	{
/**
Initially, ```h``` is set equal to the precursor film and ```hamaker``` to unity
(homogeneous substrate).
*/
		foreach()
		{
			h[]       = precFilm ;
			hamaker[] = 1.0 ;
		}

/**
Next, then the film height distribution by defining the initial droplets radii
and centers, and their heights using the initialization functions. In this step,
additional refinement is introduced at the initial contact lines.
*/
		struct ContactLineInitializer
					cli1 = {.r = 1.0, .a=4.0, .c = {.x = 4.0, .y =-4.0} };

		{
			char b[100]; sprintf(b, "%scircles.dat", dirGnu);
			FILE *fp = fopen(b, "w");

			simpleParapoloidInitialization(fp, cli1);

			fclose(fp);
		}

/**
Then, the Hamaker field values are defined to introduce substrate
heterogeneities.
*/
		//defaultHeterogeneityExpression(&c_Ham);

		heterogeneity(true, 
		 sq(2.0 + 0.1*x - 0.2*(tanh(cos(pi*(x+2.0*y))) + tanh(cos(pi*(2.0*x-y)))))
			 );

		restriction({hamaker});

/**
Similarly, substrate roughness is accordingly defined.
*/

#if ROUGH
		foreach()
		{
			s[] = 0.0;
		}

		roughness( true , cos(10.0 * x) );

		boundary({s});

		computeGradAndNegLaplace(s, gs, lapS);
		restriction({s, lapS, gs});
		gs.x.nodump = gs.y.nodump=true;
#endif

		boundary({h, hamaker});
		hold.nodump=true;

/**
The initialization fields are written based on the selected format, defined in
header file ```outputManager.h```.
*/
		outputData(fList, dirFlowFields, suffix="Init");
		outputData({hamaker}, dirFlowFields, suffix="Init");

/**
and create a directory to store the output files and write the generated
Hamaker field.
*/
#if 1
		view(width=800, height=800);
		{
 			stats x = statsf(hamaker); //needed for homogeneous case
			squares("hamaker", min=x.min, max=x.max, linear=1);
			draw_string("t=0", lc={0,0,0}, lw = 2);

			char b[100]; sprintf(b, "%shamakerInitial.png", dirFigures);
			save(b);
		}
		{
			squares("h", spread = -1, linear=1);
			draw_string("t=0", lc={0,0,0}, lw = 2);

			char b[100]; sprintf(b, "%shinitial.png", dirFigures);
			save(b);
		}
		//clear();
#else
		outputPNGSnapshot(hamaker, dirFigures, "Init");
		outputPNGSnapshot(h, dirFigures, "Init");
#if ROUGH
		outputPNGSnapshot(s, dirFigures, "Init", spread=-1);
#endif

#endif

		fprintf(fout, "# New Simulation with %ld cells\n", grid->tn);
	}
	else
	{
		fprintf(fout, "# Resuming Previous Simulation ith %ld cells\n", grid->tn);
	}

/**
The initial ```dt``` is computed based on ```minDT```.
*/
	dt = 1.e1 * minDT;

/**
We output statistics and simulated effects, along with the headers of the
variables that are printed on screen
*/
	fields_stats();

	reportSimulationEffects();
}


/**
## Variable Timestep Calculation
```dt``` is computed at each timestep based on the number of MG iteration
performed in each Newton iteration. Starting from an initial timestep, this is
progressively varied by a predetermined as indicated by the number of MG
iterations performed.

This takes the following form $$ dt_{MG} = \begin{cases} dt * factor, & mgIterSum < I_L \\
dt / factor, & mgIterSum > I_U\\ dt, &otherwise \end{cases} $$
where $I_L, I_U$ are lower and upper iteration bounds; the selected values were
found to adequately facilitate convergence.

Finally, the minimum of the two is taken, bounded by the ```minDT``` and
```maxDT```, i.e. $dt= \max(\min(dt_{max}, dt_{MG}), dt_{min})$.
*/

struct TimestepCalculator
{
	int lb,ub;
	double factor;
};

double timestep (struct TimestepCalculator p)
{
	const int lb = p.lb ? p.lb : 10;
	const int ub = p.ub ? p.ub : lb + 10;

  const double factor = p.factor ? p.factor : 1.04;

	double dtmg = dt;

	     if ( mgIterSum < lb ) dtmg *= factor;
	else if ( mgIterSum > ub ) dtmg /= factor;

	return clamp(dtmg, minDT, maxDT);
}

/**
The typical growth of $dt$ is shown in the next figure for an isolated droplet
reaching the equilibrium condition.

![BDF(1) vs BDF(2)](bdfDiff.png)
*/

/**
Here, convergence statistics are overviewed.
*/
void monitorFieldStats (scalar* fl, scalar* flold)
{
	scalar f, fold;
	for (f, fold in fl, flold)
	{
 		stats sf = statsf(f);

/**
A check is performed in case of divergence
*/
		if ( isnan(sf.sum) )
		{
			fprintf(ferr, "# Simulation crashed at i=%d %g. field: %s\n",
					iter, dt, f.name); exit(1);
		}

		const double df = change(f, fold) / precFilm;
/**
Here, the relative change of each field is outputted.
*/
		fprintf(fout, "\t% 05.2e", log10(df+1.e-30));

		fprintf(fout, "\t% 07.6f\t% 07.6f\t% 07.4e", sf.min, sf.max, sf.sum);

	}
}

// do every printStatsEvery
bool printEveryPrintStats()
{
 	return !(iter%printStatsEvery);
}

void computeConvergenceStatistics ()
{
/**
Every ```printStatsEvery``` iterations the convergence statistics are outputted
to the log screen for each field.
*/
	if ( !printEveryPrintStats() ) return;

	fprintf(fout, "%-7d", iter);

	monitorFieldStats(fList, fListOld);

	fprintf(fout, "\t%-10.7f\t%-6.3e\t%-d\n", t, dt, mgIterSum);
	fflush(fout);
}

/**
Time integration is performed up to ```Tend```
*/


event integration (t <= Tend; t += dt)
{
/**
If a-BDF2 is active, we need to store the previous ```dt``` to compute the
[differencing weights](#the-tfm-solver)
*/

#if BDF2
	double dtold = dt;
	dt = dtnext(timestep());

/**
For the first iteration, we perform an BDF1 since no ```hold``` field
is available. Prior to each solver call the old values are updated (see below).
*/
	if (!mgIterSum)  // if (1st iteration or continuation)
	{
		thinFilmSolver(fList, dt, tolerance=tolMG);
	}
	else
	{
		thinFilmSolver(fList, dt, fListOld, dtold, tolerance=tolMG);
	}
#else /* BDF1 */
	dt = dtnext(timestep());

	thinFilmSolver(fList, dt, tolerance=tolMG);
#endif

	computeConvergenceStatistics();
}

event update(i++, last)
{
#if BDF2
	foreach()
	{
		scalar cur, old;
		for (cur, old in fList, fListOld)
		{
			old[] = cur[];
		}
	}

#endif
}

/**
# Outputs
The contact lines at each time are written in separate gnp files.
*/
event outputContactLines(t+=0.5, last)
{
	char fname[80]; sprintf(fname, "%spostlines%06.02f.gnp", dirGnu, t);

	writeContactLinesFile(h, fname, 1.2*precFilm);
}

/**
A movie of the simulation is also generated
*/
#if 0
event buDumps(t+=0.1, last)
{
	char b[100];  sprintf(b, "%sdump%06.02f", dirFlowFields, t);
	dump(file = b);
}
#endif

event movie(t+=1.0, last)
{
	char timestring[100]; sprintf(timestring, "t=%06.02f", t);

	//create vertex scalar hv
	vertex scalar hv[];
	foreach_vertex()
		hv[] = (h[] + h[-1] + h[0,-1] + h[-1,-1])/4.;

	view(width=800, height=800);
	squares("h", min=precFilm, max=1.5, linear=1, map=blue_white_red);
	isoline("hv", 1.2 * precFilm, linear=1, lc={1,0,0}, lw=2.0);
//	isoline("hv", 1.4 * precFilm, linear=1, lc={0,1,0}, lw=2.0);
//  isoline("hv", 2.0 * precFilm, linear=1, lc={1,1,1}, lw=2.0);

	draw_string(timestring, lc={0,0,0}, lw = 2.0);

#if 0
	mirror(n={1,0}, alpha=0.5*L0) {
		squares("q", spread=-1, linear = 1, map = jet);
		isoline("q", 0, lc={0,0,0}, lw=2.0);
		cells(lw=0.5);
	}
#endif
	char b[80]; sprintf(b, "%s/animation.mp4", dirFigures);
	save(b);
}

/**
As well as slices
*/
#if 0
event slicer(t+=1.0, last)
{
	{
		char b[100]; sprintf(b, "%sslice%06.02f.gnp", dirGnu, t);
		outputSlice(h, b, linear = 1, n={0,1,0}, alpha=0.0);
	}

#if ROUGH
	{
		scalar eta[];
		foreach()
		{
			eta[] = z;
		}

		char b[100]; sprintf(b, "%szslice%06.02f.gnp", dirGnu, t);
		outputSlice(eta, b, linear = 1, n={0,1,0}, alpha=0.0);
	}
#endif
}
#endif

/**
Binary files of the current solution can be written
*/

event outputSolution (t+=1., last)
{
	outputData({h}, dirFlowFields, t);
}

/**
This event prepares the simulation termination by creating a dump file and
computing the *final* contact lines.
*/

event finalize (t=end)
{

	fprintf(fout, "#------------------------------------------\n");

	fields_stats();

	dump(file = "last_dump");

	fprintf(fout, "# Terminated after i=%d iterations", i);
	fprintf(fout, " and a simulation time of t=%g\n", t);
	fprintf(fout, "#------------------------------------------\n");
}

/**
# Adaptive Mesh Refinement
Mesh adaptation is triggered based on the wavelet error estimates on $\frac1{h}$.
This sensor was chosen after a series of test regarding different sensors, shown
below.

Several adaptation criteria that were tried out

| Criterion               | Remarks                                                                                               |
| -----                   | ------------                                                                                          |
| $h$                     | doesn't follow the contact line                                                                       |
| sigmoid($h$)            | very localized adaptation changes. Seems too weak irrespective of tolerance                           |
| $\Delta h$              | very strong adaptation sensor. tolerance should be scaled based on $\Delta$ to avoid over--adaptation |
| $h\Delta h$             | very strong adaptation sensor. h is not enough to attenuate                                           |
| $\frac1{h}$             | most potent sensor. Fails slightly to focus only on contact line.Should increase the tolerance                                      |
| sigmoid($\frac1{h}$)    | very localized and stong adaptation changes. Difficult to find tolerances                             |
| $\frac1{h}$, $\Delta h$ | most potent combination. tolerances need to be better adjusted

Table: Adaptation Techniques

<img src="adaptation_png/adaptation_h.png"           width="200" height="200" >
<img src="adaptation_png/adaptation_sigma(h).png"    width="200" height="200" >
<img src="adaptation_png/adaptation_lap(h).png"      width="200" height="200" >
<img src="adaptation_png/adaptation_hlap(h).png"     width="200" height="200" >
<img src="adaptation_png/adaptation_hInv.png"        width="200" height="200" >
<img src="adaptation_png/adaptation_sigma(hInv).png" width="200" height="200" >
<img src="adaptation_png/adaptation_hInvlaph.png"    width="200" height="200" >

**UPDATE:** In most cases adaptation seems to be slower that running on a
uniform mesh. Especially in cases with initial pertrubation, which result in
large meshes, the uniform meshes seems to be the preferred choice. In any case,
the setup can be changed just by including the quadtree.h/multigrid.h file.
*/

#if TREE
event adapt (i++)
{
	scalar invH[];
	foreach ()
	{
		invH[] = 1.0 / h[];
	}

#if ROUGH
	astats s = adapt_wavelet({invH, lapS},(double[]){0.3, 0.5}, LEVEL_MAX, LEVEL_IN-1);
#else
	astats s = adapt_wavelet({invH},(double[]){0.3}, LEVEL_MAX, LEVEL_IN-1);
#endif

/**
If mesh adaptation occurs, we output statistics to keep track of the mesh size
*/
	if (s.nf || s.nc )
	{
		fprintf(ferr , "#Grid adaptation %10d"             , i);
		fprintf(ferr , "\tRefined %6d\tCoarsened: %6d"     , s.nf    , s.nc);
		fprintf(ferr , "\tNumber of (L/G) Cells:%ld/%ld\n" , grid->n , grid->tn);
	}
}
#endif

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif
/**
# References
~~~bib
@article{Dallaston2021,
	title     = {Regular and complex singularities of the generalized thin film
	equation in two dimensions},
	volume    = {917},
	DOI       = {10.1017/jfm.2021.286},
	journal   = {Journal of Fluid Mechanics},
	publisher = {Cambridge University Press},
	author    = {Dallaston, M.C. and Fontelos, M.A. and Herrada, M.A.
	and Lopez-Herrera, J.M. and Eggers, J.},
	year      = {2021},
	pages     = {A20},
	pdf       = {https://people.maths.bris.ac.uk/~majge/regular_and_complex_singularities_of_the_generalized_thin_film_equation_in_two_dimensions.pdf}
}
@article{Lam2019,
	author    = {Michael Angelo Y.H. Lam and Linda J. Cummings and Lou Kondic},
	doi       = {10.1016/j.jcpx.2018.100001},
	issn      = {25900552},
	journal   = {Journal of Computational Physics: X},
	keywords  = {Film instabilities,Finite difference simulations,GPU computing,
	Long-wave approximation,Thin films},
	month     = {3},
	publisher = {Academic Press Inc.},
	title     = {Computing dynamics of thin films via large scale GPU-based simulations},
	volume    = {2},
	year      = {2019},
}
@article{Ajaev2005,
	title     = {Spreading of thin volatile liquid droplets on uniformly heated surfaces},
	volume    = {528},
	DOI       = {10.1017/S0022112005003320},
	journal   = {Journal of Fluid Mechanics},
	publisher = {Cambridge University Press},
	author    = {Ajaev, Vladimir S.},
	year      = {2005},
	pages     = {279–296},
	pdf       = {https://www.researchgate.net/publication/231787844_Spreading_of_thin_liquid_droplets_on_uniformly_heated_substrates}
}
@article{Gaskell2004,
  title={Efficient and accurate time adaptive multigrid simulations of droplet spreading},
  author={Gaskell, PH and Jimack, PK and Sellier, M and Thompson, HM},
  journal={International Journal for Numerical Methods in Fluids},
  volume={45},
  number={11},
  pages={1161--1186},
  year={2004},
  publisher={Wiley Online Library}
}
@article{Pillai2021,
   author = {Dipin S. Pillai and Kirti Chandra Sahu and Ranga Narayanan},
   doi = {10.1103/PhysRevFluids.6.073701},
   issn = {2469-990X},
   issue = {7},
   journal = {Physical Review Fluids},
   month = {7},
   pages = {073701},
   title = {Electrowetting of a leaky dielectric droplet under a time-periodic electric field},
   volume = {6},
   url = {https://link.aps.org/doi/10.1103/PhysRevFluids.6.073701},
   year = {2021},
}
@article{Kainikkara2021,
  title={Equivalence of sessile droplet dynamics under periodic and steady electric fields},
  author={Kainikkara, Muhamed Ashfak and Pillai, Dipin S and Sahu, Kirti Chandra},
  journal={npj Microgravity},
  volume={7},
  number={1},
  pages={47},
  year={2021},
  publisher={Nature Publishing Group UK London}
}
@article{Brusch2002,
  title={Dewetting of thin films on heterogeneous substrates: Pinning versus coarsening},
  author={Brusch, Lutz and K{\"u}hne, Heiko and Thiele, Uwe and B{\"a}r, Markus},
  journal={Physical Review E},
  volume={66},
  number={1},
  pages={011602},
  year={2002},
  publisher={APS}
}
@article{Savva2014,
  title={Low-frequency vibrations of two-dimensional droplets on heterogeneous substrates},
  author={Savva, Nikos and Kalliadasis, Serafim},
  journal={Journal of fluid mechanics},
  volume={754},
  pages={515--549},
  year={2014},
  publisher={Cambridge University Press}
}
@article{Pham2017,
  title={Drying of droplets of colloidal suspensions on rough substrates},
  author={Pham, Truong and Kumar, Satish},
  journal={Langmuir},
  volume={33},
  number={38},
  pages={10061--10076},
  year={2017},
  publisher={ACS Publications}
}
~~~
*/
