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
- \underbrace{\mathbf\nabla\tilde\Pi(\mathbf x, h)}_
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
with $B_o=\frac{\varrho gL^2}{\sigma}$ being the Bond number, $r$ the virbration
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
## Electro--Wetting
When electric fields are introduced, ...
we consider a leaky dielectric droplet
<span style="color:red"> TODO</span>
[Pillai et al.](#references)
[Kainikkara et al.](#references)
*/

#include "poisson.h"
/**
First some utility macros are defined.
A macro is used to define scalars that are inserted in lists automatically
to be able to switch between solving a system of equations or a single
equation.
*/
#if ELECTRO
#define createList(sc1, sc2, scl, type) type sc1[], sc2[], *scl = {sc1, sc2}
#else
#define createList(sc1, sc2, scl, type) type sc1[], *scl = {sc1}
#endif

#define createScalars(sc1, sc2, scl) createList(sc1, sc2, scl, scalar)

#define dflux(s) ( (s.x[1] - s.x[0])/Delta )

/**
# Variables
## Droplet thickness and Interface charge density
We define the ```scalar``` fields $h$ and $h_{old}$ to track the evolution of
the film equation.  If $1^{st}$ order time integration is performed
```hold``` is only used to monitor convergence. If ```-DBDF2``` is defined,
the solution at last iteration is stored in ```hold```.
*/

createScalars(h, q, fList);
createScalars(hold, qold, fListOld);

/**
## Heterogeneity and Roughness scalars
```scalar hamaker``` is used to define the sustrate heterogeneity. Finally, if
surface roughness is modelled, the scalars ```gs``` and ```lapS``` are
also introduced.
*/
scalar hamaker[];
face vector gradHet[];

#if ROUGH
face vector gs[];
scalar s[], lapS[];
#endif

/**
# Problem Setup
## Numerics Constants - Solver attributes
Some constants are defined that are associated with the numerical
solution process, such as the timestep bounds ```dt```
$\in [$```minDT, maxDT```$]$ and the total simulation time, ```Tend```.
*/
double maxDT = 1.e-2, minDT = 1.e-7, Tend = .025;

/**
The next $3$ constants are used to control the adaptive time-stepping process,
i.e. ```timestepFactor``` defines the factor to alter the timestep, while
```lMGcycles```, ```uMGcycles``` are the lower and upper bounds of total number
of MG iterations that are used to control the timestep
([See](#variable-timestep-calculation)).
*/
double timestepFactor = 1.04;
int lMGcycles = 10, uMGcycles = 20;

/**
In addition, the tolerances of the MG solver ```tolMG```, the Newton scheme
```tolNewton``` and the maximum iteration of the Newton scheme ```nitermax```
are provided.
*/
double tolMG = 1.e-6, tolNewton = 1.e-10;
int nitermax = 10;

/**
Simulation overview is controlled by ```printStatsEvery``` that defines the
iteration intervals between each statistics output.
*/
int	printStatsEvery = 100;

/**
## Physical Constants
We define some problem specific input parameters, such as the precursor
film height ```precFilm```, and the disjoining pressure mobility and iteraction
exponents ```n, m```.
*/
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
double bond = 1.00, inclineA_d = 40., contactA_d = 45.0;
#endif

/**
In addition, if the substrate is vibrated, the vibration strength ```rvib```,
acceleration ```gam(t)``` with period ```Tper``` and the angle w.r.t. the
substrate vertical ```psiA_d``` are also defined.
*/
#if FORCING
double rvib = 4.25, psiA_d = 40., Tper = 8.0;
#define gam(t) ( sin(2.0 * pi * t / Tper) )
#endif

/**
Evaporative/condensing effects depend on three non-dimensionalized
parameters, $\mathcal{K}$, $\delta$, and $\mathcal E$ defined as ```kappa,
delta, Efactor```.
*/
#if EVAPORATION
double kappa = 1.e-1, delta = 1.e-2, Efactor = 1.e-1;
#endif

/**
<span style="color:red">
TODO W-I-P
</span>
*/

#if ELECTRO
double hRatio = 2.0; // beta
double eRatio = 9.45, sRatio = 4.0, psiC = 13.0; // epsilon, sigma, psiC
double ePotential = 4.9, eAmp = 12.0, Tele = 005.00; // E, EforcingPeriod
double eAmp2 = 1.0, Tele2=0.1;
#endif

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

/*
## Disjoining Pressure Functions
Here, we define some functions that are used to compute the disjoining
pressure terms that arise for the linearization and discretization of the
equation (see [Linearization](#linearization-of-the-tfm)).
First, a function to compute the disjoining pressure term is defined as
*/

#include "disjoiningPressure.h"

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
## Electro--Wetting
[Pillai et al.](#references)
[Kainikkara et al.](#references)
Separate header file
<span style="color:red">
WORK-IN-PROGRESS
</span>
*/

/**
# PFM solver Interface
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
#if ELECTRO
	(const) face vector *JQ;  // Q term coefficient (Jacobian) // maybe merge with Cl (?)
	(const) scalar JSource;  // EW-Source Jacobian -fixme: pass better
#endif
};

/**
Given a ```scalar``` field $f$, the ```face vector``` field $\mathbf\nabla f$
and ```scalar``` field $-\Delta f$ are computed
*/
static void computeGradAndNegLaplace(scalar f, face vector grad, scalar lap)
{
	LOG;
	foreach_face()
	{
		grad.x[] = face_gradient_x(f, 0);
	}
	boundary({grad});

	foreach()
	{
		lap[] = 0.0;

		foreach_dimension()
		{
			lap[] -= dflux(grad);
		}
	}
	boundary({lap});
	LOG;
}


/**
## PFM Multigrid residual and relaxation functions
The definitions of the MG residual and relaxation functions that include
contributions depending on the heterogeneities defined are provided next.
*/

#if ELECTRO
#include "electrowettingUtilities.h"
#endif

/**
Then this function is used to compute the residual of the linearized equations.
*/
static double residual_thin(scalar *al, scalar *bl, scalar *resl, void *data)
{
	LOG;
	scalar f=al[0];
#if ELECTRO
	scalar fq=al[1];
#endif

	struct ThinData *thData = data;

	(const) face vector A    = thData->A;
	(const) face vector Bh   = thData->Bl[0];
	// h- eq Jacobians(C)
	(const) face vector Cdh  = thData->Cl[0];
#if ELECTRO
	(const) face vector Bq   = thData->Bl[1];
	(const) face vector Cdq  = thData->Cl[1];
	// Q- eq Jacobians(J)
	(const) face vector Jqdh = thData->JQ[0];
	(const) face vector Jqdq = thData->JQ[1];
#endif

/**
Compute $\mathbf \nabla f$ and $\Delta h = f$
*/
	face vector gradf[]; scalar p[];
	computeGradAndNegLaplace(f, gradf, p);

/**
and $A \mathbf \nabla p$, $B \mathbf \nabla h$, and $\mathbf C h$
*/

	createList(gh, gq, g, face vector);

/**
the inclination term $ 3 \nabla_x ( G h )$ is included in C,
*/
	foreach_face()
	{
		gh.x[] = A.x[]   * face_gradient_x(p, 0)
					 + Bh.x[]  * gradf.x[]
					 + Cdh.x[] * face_value(f, 0);
#if ELECTRO
		Variables pv = {face_value(f, 0), face_value(q, 0), {x, y}};

		gh.x[] += Bq.x[] * face_gradient_x(fq, 0)
					 + Cdq.x[] * face_value(fq, 0)
					 + psiPotentialSpaceTerm(pv) * gradPsiPotential_x(pv);

		gq.x[] = Jqdh.x[] * face_value(f , 0)
			     + Jqdq.x[] * face_value(fq, 0);
#endif
	}

/**
the substrate roughness term
*/
#if ROUGH
	foreach_face()
	{
		gh.x[] += A.x[] * face_gradient_x(lapS, 0);
#if ( INCLINED || FORCING )
		gh.x[] += A.x[] * grav1 * gs.x[];
#endif
	}
#endif

/**
and the evaporation term $ ( E )$ is already included in the diagterm
*/

/**
Here, the Residual is computed
*/
	double maxres = 0.0;
	foreach(reduction(max:maxres))
	{
		// f: var, b: RHS of f, res: residual, l:diagonal of var
		scalar f, b, res, l;

		// geq: fluxes of var equation
	 	face vector geq;
		for ( f, b, res, l, geq in al, bl, resl, thData->lambda, g)
		{
			res[] = b[] - l[] * f[];
/**
The following terms are also added to the residual, <span style="color:red">
depending on the considered equation -FIX TEXT </span>,
$\mathbf \nabla \cdot (A \mathbf \nabla p)$,
$\mathbf \nabla \cdot (B \mathbf \nabla h)$,
and $\mathbf \nabla \cdot (\mathbf C h)$
*/
			foreach_dimension()
			{
				res[] -= dflux(geq);
			}

			maxres = max(maxres, fabs(res[]));
		}
	}
	boundary(resl);

	LOG;
	return maxres;
}

/**
The relaxation function is obtained by taking the ```residual_thin```
function and inverting for $h$
*/
static void relax_thin(scalar *al, scalar *bl, int l, void *data)
{
	LOG;
	scalar f = al[0], b = bl[0];
	struct ThinData *thData = data;

	(const) scalar lambdah   = thData->lambda[0];
	(const) face vector A    = thData->A;
	(const) face vector B    = thData->Bl[0];
	(const) face vector Cdh  = thData->Cl[0];


#if ELECTRO
	(const) face vector Bq   = thData->Bl[1];
	// q- field
	scalar fq = al[1];
	// q- Contributions
	(const) face vector Cdq = thData->Cl[1];
#endif

#if JACOBI
	scalar c[];
#else
	scalar c = f;
#endif

	foreach_level_or_leaf (l) {
		double n = sq(sq(Delta)) * b[], d = lambdah[] * sq(sq(Delta));
		foreach_dimension() {
			// DIV(GRAD(LAPLACIAN)) TERM
			n += A.x[1] * (f[ 2] - 5.0 * f[ 1] - f[-1] - f[0,-1] - f[0,1] + f[ 1, 1] + f[ 1,-1]);
			n += A.x[ ] * (f[-2] - 5.0 * f[-1] - f[ 1] - f[0,-1] - f[0,1] + f[-1, 1] + f[-1,-1]);
			d -= 5.0 * (A.x[1] + A.x[]);
			// DIV(GRADIENT) TERM
			n -= (B.x[1] * f[1] + B.x[] * f[-1]) * sq(Delta);
			d -= (B.x[1]        + B.x[]        ) * sq(Delta);
			// DIV(FACEVAL) TERM
			n -= 0.5 * (Cdh.x[1] * f[1] - Cdh.x[] * f[-1]) * cube(Delta);
			d += 0.5 * (Cdh.x[1]        - Cdh.x[]        ) * cube(Delta);
#if ELECTRO
			// DIV(GRADIENT) TERM
			n -= (Bq.x[1] * fq[1] + Bq.x[] * fq[-1]) * sq(Delta);
			d -= (Bq.x[1]         + Bq.x[]         ) * sq(Delta);
			// DIV(FACEVAL) TERM - Q
			n -= 0.5 * (Cdq.x[1] * fq[1] - Cdq.x[] * fq[-1]) * cube(Delta);
			d += 0.5 * (Cdq.x[1]         - Cdq.x[]         ) * cube(Delta);
#endif
		}

		c[] = n / d;
	}

	LOG;
#if ELECTRO
	relax_electroThin(al, bl, l, data);
#endif

#if JACOBI
  foreach_level_or_leaf (l)
    f[] = (f[] + 2.*c[])/3.;
#endif

#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = f[];
  trash ({f});
  foreach_level_or_leaf (l)
    f[] = a1[];
#endif
}

#if INCLINED
static void newton_Contributions_ITFM (struct ThinData p, face vector gradH, scalar *bl)
{
	scalar f = p.f[0];
/**
Add Inclination Contributions
*/
	face vector A = p.A;
	face vector B = p.Bl[0];
	face vector C = p.Cl[0];

	face vector gCi[];
	foreach_face()
	{
		const double hf = face_value(f, 0);

		B.x[] += A.x[] * grav1;

		double iContribution = pow_der(hf, 3.0) * grav1 * gradH.x[];

		C.x[]   += iContribution;
		gCi.x[]  = iContribution * hf;
	}

/**
Compute inclination terms, i.e. ```face vector``` $Gr, gGr$
*/
	face vector gGr[];

	foreach_face(x)
	{
		const double hf   = face_value(f, 0);
		const double whf2 = sq(hf) * grav2;

		// Gr is immediately included in C to avoid extra face vector definitions
		C.x[]   += -3.0 * whf2;
		gGr.x[] = -2.0 * whf2 * hf;
	}
/**
Here, the $\mathcal{RHS}$ ```b``` is computed with the additional terms arising from
the linearization
*/

	scalar b = bl[0];

	foreach()
	{
		b[] += dflux(gGr);

		foreach_dimension()
		{
			b[] += dflux(gCi);
		}
	}
}
#endif

#if ROUGH
static void newton_Contributions_Rough(struct ThinData p, scalar *bl)
{
	scalar f = p.f[0];
/**
Add Roughness Contributions
*/
	face vector C = p.Cl[0];
	face vector gCr[];

	foreach_face()
	{
		const double hf = face_value(f, 0);

		double rContribution = pow_der(hf, 3.0) * face_gradient_x(lapS, 0);
#if INCLINED
		rContribution += pow_der(hf, 3.0) * grav1 * gs.x[];
#endif

		C.x[]   += rContribution;
		gCr.x[]  = rContribution * hf;
	}

	scalar b = bl[0];
	foreach()
	{
		foreach_dimension()
		{
			b[] += dflux(gCr);
		}
	}
}
#endif

#if EVAPORATION
static void newton_Contributions_ETFM(struct ThinData p, scalar lapH, scalar *bl)
{
/**
Add Evaporation Contributions
*/
	//Source Term
	scalar b = bl[0];
	scalar f = p.f[0];
	scalar l = p.lambda[0];

	foreach()
	{
		const double eTerm = evaporationTerm;

		// Evap Source - LHS
		l[] += eTerm;

		//evaporation linearization term - RHS
		b[] += ( kappa + z + f[] ) * eTerm;
	}
}
#endif

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
Compute the electro terms,
*/
#if ELECTRO
		createList(Jh, Jq, JQ, face vector); // store Jacobians of Q
	 	p.JQ = JQ;

		scalar JSource[];
		p.JSource = JSource;
		newton_Contributions_ElectroTFM(p, gradH, lapH, bl);
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

#if ELECTRO
		for (scalar sc in (scalar *)p.JQ)
			coefs = list_add(coefs, sc);

			coefs = list_add(coefs, p.JSource);
#endif

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


event init (t=0);

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
	dt = dtnext(timestep(lb = lMGcycles, ub = uMGcycles, factor=timestepFactor));

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
	dt = dtnext(timestep(lb = lMGcycles, ub = uMGcycles, factor=timestepFactor));

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

