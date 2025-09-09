/**
# Header file that deals with the electric terms arising in the PFM.
*/
#if !ELECTRO // include guard
#error Include only when ELECTRO is specified
#endif

#if ROUGH
#error Rough Electro not yet implemented
#endif

/**
Utility definitions
*/
typedef struct
{
	double hf, qf;
	coord  c;
} Variables;

#define eRatioMinus ( 1.0 - eRatio )

double denomInv(double hf)
{
	return 1.0 / (hf * eRatioMinus + eRatio * hRatio);
}

/**
## Potential Excitation
*/
double psiPotential(Variables p)
{
	return eAmp * sin(2.0 * pi * t / Tele) + eAmp2 * cos( (p.c.x + p.c.y) / Tele2 );
}

double gradPsiPotential_x(Variables p)
{
	return -eAmp2 / Tele2 * sin( (p.c.x + p.c.y) / Tele2 );
}

double gradPsiPotential_y(Variables p)
{
	return -eAmp2 / Tele2 * sin( (p.c.x + p.c.y ) / Tele2);
}

double df1_dpsi(double hf)
{
	return 1.0 - hf * denomInv(hf);
}

double df2dpsi(Variables p)
{
	return (hRatio - p.hf ) * eRatio * denomInv(p.hf);
}

/**
## Potential 1 - $\Phi_1$

### $2^{nd}$ order derivatives
*/
double ddf1_dhz(Variables p)
{
	return ( eRatioMinus * psiPotential(p) - hRatio * p.qf ) * sq(denomInv(p.hf));
}

double ddf1_dhh(Variables p)
{
	return -2.0 * p.hf * eRatioMinus * ddf1_dhz(p) * denomInv(p.hf);
}

double ddf1_dqh(Variables p)
{
	return -p.hf * hRatio * sq(denomInv(p.hf));
}

double ddf1_dqz(Variables p)
{
	return (hRatio - p.hf) * denomInv(p.hf);
}

/**
### $1^{st}$ order derivatives
*/
double df1_dh(Variables p)
{
	return ( p.hf * ddf1_dhz(p) );
}

double df1_dz(Variables p)
{
	return ( (hRatio - p.hf) * p.qf - psiPotential(p) ) * denomInv(p.hf);
}

double df1_dq(double hf)
{
	return hf * ( hf - hRatio ) * denomInv(hf);
}

/**
## Potential 2 - $\Phi_2$
*/
double df2_dz(Variables p)
{
	return -( eRatio * psiPotential(p) + p.hf * p.qf ) * denomInv(p.hf);
}

double ddf2_dhz(Variables p)
{
	return eRatio * ddf1_dhz(p);
}

double ddf2_dqz(double hf)
{
	return -hf * denomInv(hf);
}

/**
## Maxwell Stress Components computations
The maxwell tangential contributions $\mathcal C^{(h, Q)}$ in
$\mathcal C^{h} \mathbf \nabla h$ and $\mathcal C^{Q} \mathbf \nabla Q$ are
computed using the following functions.

### Utility Definitions
*/
double tangentialMaxwellMultiplier(double hf)
{
	return 0.5 * sq(hf) * ePotential;
}

/**
Term $\frac{\mathcal C^{(h)}}{Q}$
*/
double cHTermToQ(Variables p)
{
	//introduced to avoid division by zero when Qf = 0.
	return -tangentialMaxwellMultiplier(p.hf) * (df1_dh(p) + df1_dz(p));
}

/**
Term $\mathcal C^{(h)}$
*/
double cHTerm(Variables p)
{
	return cHTermToQ(p) * p.qf;
}

/**
Term $\mathcal C^{(Q)}$
*/
double cQTerm(Variables p)
{
	return -tangentialMaxwellMultiplier(p.hf) * p.qf * df1_dq(p.hf);
}


/**
Here, the functions that compute the Jacobians i.e. $\mathcal C^{(h)}$ and $\mathcal C^{(Q)}$
derivatives w.r.t.  $\mathbf U$ are computed. First the function that computes
$\partial_h \mathcal C^{(h)}$ is defined as
*/
double cHTerm_dh(Variables p)
{
	// derivatives w.r.t. h
	const double df1_ders = ddf1_dhh(p) + ddf1_dhz(p);

	const double dcdh_1 = 2.0 * cHTerm(p) / p.hf;
	const double dcdh_2 = tangentialMaxwellMultiplier(p.hf) * p.qf * df1_ders;

	return dcdh_1 - dcdh_2;
}

/**
And then $\partial_Q \mathcal C^{(h)}$
*/

double cHTerm_dq(Variables p)
{
	// derivatives w.r.t. q
	const double df1_ders = ddf1_dqh(p) + ddf1_dqz(p);

	const double dcdq_1 = cHTermToQ(p);
	const double dcdq_2 = tangentialMaxwellMultiplier(p.hf) * p.qf * df1_ders;

	return dcdq_1 - dcdq_2;
}

/**
$\partial_h \mathcal C^{(Q)}$
*/

double cQTerm_dh(Variables p)
{
	// derivatives w.r.t. h
	const double dcdh_1 = 2.0 * cQTerm(p) / p.hf;
	const double dcdh_2 = tangentialMaxwellMultiplier(p.hf) * p.qf * ddf1_dqh(p);

	return dcdh_1 - dcdh_2;
}

/**
Finally, $\partial_Q \mathcal C^{(Q)} = 0$.

The same procedure is followed regarding the maxwell normal component i.e.
$\mathcal D^{(h,Q)}$
*/
double normalMaxwellMultiplier(double hf)
{
	return eRatio * ePotential * cube(hf) / 3.0;
}

/**
Term $\mathcal D^{(h)}$
*/
double dHTerm(Variables p)
{
	const double df1df2 = df1_dz(p) - df2_dz(p);
	const double ders   = ddf1_dhz(p) * df1df2;

	return -normalMaxwellMultiplier(p.hf) * ders;
}

/**
$\mathcal D^{(Q)}$
*/
double dQTerm(Variables p)
{
	const double df1 = eRatio * df1_dz(p) * ddf1_dqz(p);
	const double df2 = df2_dz(p) * ddf2_dqz(p.hf);

	return -cube(p.hf) * ePotential * ( df1 - df2 ) / 3.0;
}

/**
Then the Jacobians of $\mathcal D^{(h,Q)}$ are defined.

```ders``` are the derivaties of the potential when computing
$\frac{\partial \mathcal D}{\partial h}$, i.e.
$$
\frac{\partial}{\partial h}
\left[
\left(\frac{\partial}{\partial h} \frac{\partial \Phi_1}{\partial z}\right)
\left(\frac{\partial\Phi_1}{\partial z}-\frac{\partial\Phi_2}{\partial z}\right)
\right]
= 3 (1-\epsilon)
\left(
\frac{\partial^2\Phi_1}{\partial z\partial h}
\right) ^2
$$

First, $\partial_h \mathcal{D}^{(h)}$
*/
double dHTerm_dh(Variables p)
{
	const double ders = eRatioMinus * sq(ddf1_dhz(p));

	const double dDdh_1 = dHTerm(p) / p.hf;
	const double dDdh_2 = normalMaxwellMultiplier(p.hf) * ders;

	return 3.0 * ( dDdh_1 - dDdh_2 );
}


/**
$\partial_Q \mathcal{D}^{(h)}$
*/
double dHTerm_dq(Variables p)
{
// derivatives w.r.t. q
  const double preFac = -2.0 * normalMaxwellMultiplier(p.hf) * hRatio;
	return preFac * ddf1_dhz(p) * denomInv(p.hf);
}

/**
$\partial_h \mathcal{D}^{(Q)}$
*/
double dQTerm_dh(Variables p)
{
	const double ders = 2.0 * hRatio * ddf1_dhz(p) * denomInv(p.hf);
	// f1 to give eRatio to normalMaxwellMultiplier

	const double dDdh_1 = 3.0 * dQTerm(p) / p.hf;
	const double dDdh_2 = normalMaxwellMultiplier(p.hf) * ders;

	return dDdh_1 - dDdh_2;
}

/**
And, $\partial_Q \mathcal{D}^{(Q)}$
*/
double dQTerm_dq(Variables p)
{
	const double prefactor = -ePotential * cube(p.hf) * sq(denomInv(p.hf)) / 3.0;
	const double der_dq = eRatio * hRatio * (hRatio - 2.0 * p.hf)
		                  - sq(p.hf) * eRatioMinus;

	return prefactor * der_dq;
}

/**
### Contributions Computations

$\mathcal C^{(h)} + \mathcal D^{(h)}$
*/

double electroFluxH(Variables p)
{
	return -( cHTerm(p) + dHTerm(p) );
}

/**
$B^Q= - (\mathcal C^{(Q)} + \mathcal D^{(Q)} )$
*/
double electroFluxQ(Variables p)
{
	return -( cQTerm(p) + dQTerm(p) );
}

/**
Components for $\mathbf \Gamma ^{(\delta h)}$. (need to be multiplied by
$\nabla \mathbf U$)
*/
double electroJacobiandhGradH(Variables p)
{
	return -( cHTerm_dh(p) + dHTerm_dh(p) );
}

double electroJacobiandhGradQ(Variables p)
{
	return -( cQTerm_dh(p) + dQTerm_dh(p) );
}

double electroJacobianH(Variables p, double h_grad, double q_grad)
{
	return electroJacobiandhGradH(p) * h_grad
		   + electroJacobiandhGradQ(p) * q_grad;
}

/**
Components for $\mathbf \Gamma ^{(\delta q)}$. (need to be multiplied by
$\nabla \mathbf U$)
*/

double electroJacobiandQGradH(Variables p)
{
	return -( cHTerm_dq(p) + dHTerm_dq(p) );
}

double electroJacobiandQGradQ(Variables p)
{
	return -dQTerm_dq(p); // cTerm == 0
}

double electroJacobianQ(Variables p, double h_grad, double q_grad)
{
	return electroJacobiandQGradH(p) * h_grad
		   + electroJacobiandQGradQ(p) * q_grad;
}

/**
## Charge influx rate - Source term $S(\mathbf{U})$
*/
double chargeInflux(Variables p)
{
	// computed at center
	return psiC * ( df2_dz(p) - sRatio * df1_dz(p) );
}

double chargeInfluxDerivativeH(Variables p)
{
	// computed at center
	return psiC * ( eRatio - sRatio ) * ddf1_dhz(p);
}

double chargeInfluxDerivativeQ(Variables p)
{
	// computed at center
	return psiC * ( sRatio * ( p.hf - hRatio ) - p.hf ) * denomInv(p.hf);
}

/**
## Potential Excitation Contribution
*/
//fixme: -move frome here
double psiPotentialSpaceTerm(Variables p)
{
	const double Cpsi = -tangentialMaxwellMultiplier(p.hf) * p.qf * df1_dpsi(p.hf);

	const double df1df2 = df1_dz(p) - df2_dz(p);
	const double Dpsi = -normalMaxwellMultiplier(p.hf) * denomInv(p.hf) * df1df2;

	return Cpsi + Dpsi;
}


/**
## Newton Contributions
*/

static void newton_Contributions_ElectroTFM
	(struct ThinData p, face vector gradH, scalar lapH, scalar *bl)
{
	scalar fh = p.f[0];
	scalar fq = p.f[1];

	face vector A    = p.A;
	face vector Bh   = p.Bl[0];
	face vector Bq   = p.Bl[1];
	face vector Cdh  = p.Cl[0], Cdq = p.Cl[1];
	face vector Jqdh = p.JQ[0], Jqdq = p.JQ[1];

	scalar Js = p.JSource;

	face vector gC[], gG[];
	face vector gradQ[];

	//RHS - Contributions
	foreach_face()
	{
		const double hf = face_value(fh, 0);
		const double qf = face_value(fq, 0);
		const double ihf = 1.0 / hf;

		Variables p = {hf, qf, {x, y}};

		gradQ.x[] = face_gradient_x(fq, 0);

		// Maxwell Stresses
		Bh.x[] += electroFluxH(p);
		Bq.x[]  = electroFluxQ(p);

		// Maxwell Stresses Jacobian -> dh
		const double jacobianDh = electroJacobianH(p, gradH.x[], gradQ.x[]);

		// Add Contribution to Gamma_dh
		Cdh.x[] += jacobianDh;

		// Maxwell Stresses Jacobian -> dQ
		Cdq.x[] = electroJacobianQ(p, gradH.x[], gradQ.x[]);

		// h-residual
		const double G =  A.x[] * face_gradient_x(lapH, 0)
                   + Bh.x[] * gradH.x[]
                   + Bq.x[] * gradQ.x[];

		//store Jacobians for MG method
		Jqdh.x[] = 0.5 * qf * sq(ihf) * ( Cdh.x[] * hf - G );
		Jqdq.x[] = 0.5 * ihf          * ( Cdq.x[]      + G );

		// Maxwell contributions to h-eq - We only add ew contributions
		// since the others are already included
		gC.x[] = jacobianDh * hf + Cdq.x[] * qf;

		// Contributions to Q-eq
		gG.x[]  = 0.5 * qf * ihf * ( Cdh.x[] * hf + Cdq.x[] - G );
	}

	//Source Term for Q
	scalar bh = bl[0], bq = bl[1];

	scalar lambdaq = p.lambda[1];

	foreach()
	{
		Variables p = {fh[], fq[], {x, y}};

		const double Su    = chargeInflux(p);
		const double dSudh = chargeInfluxDerivativeH(p);
		const double dSudq = chargeInfluxDerivativeQ(p);

		lambdaq[] += dSudq;
		//Source Term - Diag
		Js[] = dSudh;

		//Source Term - RHS
		bq[] = dSudh * fh[] + dSudq * fq[] - Su;

		foreach_dimension()
		{
			bq[] += dflux(gG);

			bh[] += dflux(gC);
		}
	}
}

/**
Relaxation function for Q-Equation
*/

static void relax_electroThin(scalar *al, scalar *bl, int l, void *data)
{
	scalar fh = al[0];
	scalar f = al[1], b = bl[1];

	struct ThinData *thData = data;

	(const) scalar lambda = thData->lambda[1];
	(const) scalar Js     = thData->JSource;
	// Q- eq Jacobians(J)
	(const) face vector Jqdh = thData->JQ[0];
	(const) face vector Jqdq = thData->JQ[1];

#if JACOBI
	scalar c[];
#else
	scalar c = f;
#endif

	foreach_level_or_leaf (l) {
		double n = Delta * b[], d = lambda[] * Delta;

		n -= Js[] * fh[];

		foreach_dimension() {
			// DIV(FACEVAL) TERM - Q
			n -= 0.5 * (Jqdq.x[1] * f[1] - Jqdq.x[] * f[-1]);
			d += 0.5 * (Jqdq.x[1]        - Jqdq.x[]        );
			// DIV(FACEVAL) TERM - h
			n -= 0.5 * (Jqdh.x[1] * fh[1] - Jqdh.x[] * fh[-1]);
			d += 0.5 * (Jqdh.x[1]         - Jqdh.x[]         );
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
