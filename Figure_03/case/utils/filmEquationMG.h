/**
Header file used for the definition of the [TFE](#thinFilm.h) Multigrid (MG)
solver residual and relaxation functions.

First a utility macro is defined:
*/
#define dflux(s) ( (s.x[1] - s.x[0])/Delta )

/**
Given a ```scalar``` field $f$, the ```face vector``` field $\mathbf\nabla f$
and ```scalar``` field $-\Delta f$ are computed
*/
static void computeGradAndNegLaplace(scalar f, face vector grad, scalar lap)
{
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
}

/**
This function is used during the MG iterations and
computes the residual of the linearized equations.
*/
static double residual_thin(scalar *al, scalar *bl, scalar *resl, void *data)
{
	scalar f=al[0];

	struct ThinData *thData = data;

	(const) face vector A    = thData->A;
	(const) face vector Bh   = thData->Bl[0];
	// h- eq Jacobians(C)
	(const) face vector Cdh  = thData->Cl[0];

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

	return maxres;
}

/**
The relaxation function is obtained by taking the ```residual_thin```
function and inverting for $h$
*/
static void relax_thin(scalar *al, scalar *bl, int l, void *data)
{
	scalar f = al[0], b = bl[0];
	struct ThinData *thData = data;

	(const) scalar lambdah   = thData->lambda[0];
	(const) face vector A    = thData->A;
	(const) face vector B    = thData->Bl[0];
	(const) face vector Cdh  = thData->Cl[0];



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
		}

		c[] = n / d;
	}

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

