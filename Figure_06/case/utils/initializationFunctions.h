/**
Header file used for the definition of the functions used to define the
initial conditions for the studied problems.

# Interface
First, an interface ```struct``` is introduced to help generate the initial
conditions for $h$.
The struct that holds relevant information is defined as
*/


typedef struct
{
	double a; // amplitude (max height). Default = r
	double r; // radius
	coord c;  // center of circular initialization
} ContactLineInitializer;

ContactLineInitializer * dropsList;

// defined in "inputReader.h" - fixme
void printDropletParameters(ContactLineInitializer * cl);

/**
Merges Input Droplets and those given from file.
*/
ContactLineInitializer *mergeInputDrops(ContactLineInitializer * clIn, int sizeIn)
{
	Array *dropAr = array_new();

	for (int i = 0; i < sizeIn; i++)
	{
		array_append(dropAr, &(clIn[i]), sizeof(ContactLineInitializer));
	}


	int i = 0;
	if (dropsList)
	{
		while (dropsList[i].a != nodata)
		{
			array_append(dropAr, &(dropsList[i]), sizeof(ContactLineInitializer));
			i++ ;
		}
		free(dropsList);
	}


	ContactLineInitializer cli = {.a=nodata, .r=nodata, .c.x=nodata, .c.y=nodata};
	array_append(dropAr, &cli, sizeof(ContactLineInitializer));

	return array_shrink(dropAr);
}
/**
A helper function is defined that pre--refines the mesh near the defined
contact lines to ensure adequate mesh resolution at early iterations.
*/

void prerefine (const double r, const coord c)
{
#if TREE
	refine( sq( 0.95 * r ) < sq(x-c.x) + sq(y-c.y)
			 && sq( 1.05 * r ) > sq(x-c.x) + sq(y-c.y)
			 && level < LEVEL_MAX-1 );
#endif
}

/**
# Free-surface Initialization
Simple parapoloid functions or smoothened parapoloid mappings are used to
define the $h$ field values at predefined locations.

## Parapoloid Function
The $h$ values are defined for a radius ```r``` around center
$\mathbf c=[c_x, c_y]^T$ based on the following function
$$ h = \max\left( \varepsilon, A (1.0 -
\frac{||\mathbf{x} -\mathbf c||_2}{r^2} \right) $$
*/

void simpleParapoloidInitialization (FILE *fp, ContactLineInitializer cl)
{
	if (!cl.r) return;
	if (!fp) fp = fopen("circles.dat", "w");
	if (!cl.a) cl.a = 0.5 * cl.r;

	prerefine(cl.r, cl.c); // refine near initial contact line

	const double r2 = cl.r * cl.r;
	const double r2Inv = 1.0 / (r2);
	foreach()
	{
		const double r2Cell = sq(x-cl.c.x) + sq(y-cl.c.y);
		if (r2Cell < r2)
			h[] = max(precFilm, cl.a * (1.0 - (r2Cell)*r2Inv)); // add height at collocated initialization
	}

	if (pid() == 0) fprintf(fp, "%g %g %g \n", cl.c.x, cl.c.y, cl.r);
}

/**
## Parapoloid Mapping
In this case, the $h$ values are defined via the following mapping
$$ h(\mathbf x, t=0) \leftarrow G\left(H\left(F\left(\mathbf x, a, r\right)
\right)\right)$$
with $a$ and $r$ being
the smoothness parameter and radius, respectively.
*/
double funH (const double f, const double r)
{
	return -f * f + r * r;
}

double funF (FILE *fp, const double x, const double y, const double a, const double r)
{
	const double posNorm = sqrt( sq(x) + sq(y) );
	const double sgxy = sign(x) * sign(y);

	const double nom = 1.0 + exp(-2.0 * a * sgxy * ( posNorm + r) );
	const double den = 1.0 + exp(-2.0 * a * sgxy * ( posNorm - r) );

	if (pid() == 0) fprintf(fp, "%g %g %g \n", x, y, r);

	return sgxy * r + 0.5 * log( nom/den ) / a;
}

//Initialize h with a parapoloid function that connects smoothly with the
//precursor film
void smoothParapoloidInitialization (FILE *fp, ContactLineInitializer cl)
{
	if (!cl.r) return;
	if (!fp) fp = fopen("circles.dat", "w");

	if (!cl.a) cl.a = 50.0;

	prerefine(cl.r, cl.c); // refine near initial contact line

	foreach()
		h[] = (1.0 - precFilm) *
					funH( funF(fp,x-cl.c.x, y-cl.c.y, cl.a, cl.r), cl.r );
}

/**
## Usage
*/

void initializeInputDrops(ContactLineInitializer * clIn, int sizeIn)
{
	dropsList = mergeInputDrops(clIn, sizeIn);

	char b[100]; sprintf(b, "circles.dat");
	FILE *fp = fopen(b, "w");
	if (fp == NULL)
	{
		perror(b);
		exit(1);
	}

	int i = 0;
	while (dropsList[i].a != nodata)
	{
		fprintf(fout, "# %g %g %g %g\n",
				dropsList[i].a, dropsList[i].r, dropsList[i].c.x, dropsList[i].c.y);
		simpleParapoloidInitialization(fp, dropsList[i]);
		i++;
	}

	fclose(fp);

	printDropletParameters(dropsList);

	free(dropsList);
}

/**
# Gnuplotting
~~~gnuplot Initialization Functions (r=1, A=1, c=[1.5, (0.0)] )
reset;
set grid; set size ratio -1;
e=0.01; A=1; r=1; r2=r*r; c=1.5;
max(x,y)=x>y?x:y;
f(x) = max(e, A/r2 * (1.0 - (x-c)**2/r2))

H(f, r) = -f*f+r*r;
FN(x,a,r) = 1.0 + exp(-2.0 * a * sgn(x) * (abs(x) + r));
FD(x,a,r) = 1.0 + exp(-2.0 * a * sgn(x) * (abs(x) - r));
F(x,a,r) = sgn(x) * r + 0.5 * log( FN(x,a,r) / FD(x,a,r) ) / a;
g(x, a) = e + (1.0 - e) * H(F(x-c,a,r),r)

p[0:3] f(x)    w l lw 2 t 'Simple Parapoloid', \
       g(x,5)   w l lw 2 t 'Smooth - a=5', \
       g(x,10)  w l lw 2 t 'Smooth - a=10', \
       g(x,100) w l lw 2 t 'Smooth - a=100'
~~~

The initial configuration can be visualized via gnuplot by calling

~~~gnuplot Initial Droplet Position
stats 'circles.dat' nooutput
stats 'circles.dat' u 3 name "R" nooutput
rmax = 1.5 * R_max;
xp=STATS_max_x+rmax;
xm=STATS_min_x-rmax;
yp=STATS_max_y+rmax;
ym=STATS_min_y-rmax;
set style fill transparent solid 0.2
set grid; set size ratio -1;
set ylabel 'y'; set xlabel 'x';
p[xm:xp][ym:yp] 'circles.dat' with circles not
~~~
*/

/**
# Heterogeneity and Roughness Initialization
A utility function is also defined to facilitate the initialization of
prescribed fields such as substrate heterogeneity and roughness.
*/
#define initializeField(field, cond, expr) do { \
	foreach() {         \
		if(cond) {        \
			field[] = expr; \
		}                 \
	}                   \
}	while (0)

/**
The Hamaker profile is by default defined as
$$\Phi(\mathbf x) = 1 + p_0\tilde{x} + p_1\tilde{y}
+ p_2\tanh\bigl[p_3\cos\bigl(p_4\pi\left(\sin\left(p_6\right)\tilde{x} +
\cos\left(p_6\right)\tilde{y}\right)\bigr)\cos\bigl(p_5\pi\tilde{x}\bigr)\bigr]
$$
with
$$\tilde{x} = \cos(p_7)(x/L0) - \sin(p_7)(y/L0)$$
$$\tilde{y} = \sin(p_7)(x/L0) + \cos(p_7)(y/L0)$$

The coefficients lie in the ranges of
$$p_0, p_1 \in\left[ 0,1    \right],
p_2      \in\left[ 0,0.2  \right],
p_3      \in\left[-5,5    \right],
p_4, p_5  \in\left[ 0,20		 \right],
p_6      \in\left[ -\tfrac{\pi}{2}, \tfrac{\pi}{2},     \right],
p_7      \in\left[ -\pi,\pi    \right]$$
*/
union HamakerCoefficients
{
	// ISO C99
	struct HC_struct {
		double p0, p1, p2, p3, p4, p5, p6, p7;
	} hcs_Data;

	double hamC[8];
} c_Ham;

double hamakerExpression (double x, double y, const struct HC_struct *coefs_ptr )
{
	const struct HC_struct coefs = *coefs_ptr;

	const double s6 = sin(coefs.p6), c6 = cos(coefs.p6);
	const double s7 = sin(coefs.p7), c7 = cos(coefs.p7);

	const double xTilde = ( c6 * x - s6 * y ) / L0;
	const double yTilde = ( s7 * x + c7 * y ) / L0;

	const double c4 = cos(coefs.p4 * pi * (s6 * xTilde + c6 * yTilde));
	const double c5 = cos(coefs.p5 * pi * xTilde);

	const double tanhTerm = tanh(coefs.p3 * c4 * c5);

	return 1.0 + coefs.p0 * xTilde + coefs.p1 * yTilde + coefs.p2 * tanhTerm;
}

void defaultHeterogeneityExpression(const union HamakerCoefficients *hc)
{
	const struct HC_struct *coefs_ptr = &(hc->hcs_Data);
	foreach()
	{
			hamaker[] = hamakerExpression(x, y, coefs_ptr);
	}
}

/**
However, it can also be initialized using
any expression provided during compile-time.
*/
#define heterogeneity(cond, expr) initializeField(hamaker, cond, expr)

#if ROUGH
#define roughness(cond, expr) initializeField(s, cond, expr)
#endif
