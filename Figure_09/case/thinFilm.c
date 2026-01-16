/**
# Example for the Precursor Film Model solver
*/
#if 0
#define LOG (fprintf(stderr, "# %s  %d\n",__func__, __LINE__))
#define LOGSTATS (fields_stats())
#else
#define LOG
#define LOGSTATS
#endif

#define BGHOSTS 2

#include "grid/multigrid.h"
#include "run.h"
#include "precursorFilmModel.h"
#include "fractions.h"
#include "view.h"
#include "maxruntime.h"
#include "outputManager.h"

/**
# Problem Setup
## Numerics Constants
Following, some constants are defined that are associated with the numerical
solution process, such as the timestep bounds ```dt```
$\in [$```minDT, maxDT```$]$ and the total simulation time, ```Tend```.
*/

/**
In addition, the tolerances of the MG solver ```tolMG```, the Newton scheme
```tolNewton``` and the maximum iteration of the Newton scheme ```nitermax```
are provided.
*/

/**
Here we define, the simulation output directories relative paths.
*/
char dirFigures[]    = "png/";
char dirFlowFields[] = "flows/";
char dirGnu[]        = "gnuData/";

/**
Mesh resolution is bounded by a maximum level.
*/
int LEVEL_IN  = 11; // initial resolution (optional input arg)
int LEVEL_MAX = 11; // maximum resolution

/**
# Simulations Setup
At the beginning, the boundary conditions, domain size, and computational mesh
are initialized.
*/

#include "inputReader.h"

/**
We initialize either from a previous solution, or using a combination of the
[initialization](#auxilliary-initialization-functions) functions.
*/
event init (t=0)
{
	LOG;
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
			h[]       = precFilm;// + precFilm * (noise() + 1.0);
			hamaker[] = 1.0;
		}

/**
Next, then the film height distribution by defining the initial droplets radii
and centers, and their heights using the initialization functions. In this step,
additional refinement is introduced at the initial contact lines.
*/
		//ContactLineInitializer cls[0];
		ContactLineInitializer cls[] = {
			{.r = 0.35, .c = {.x =-1.5, .y =-1.5} },
			{.r = 0.75, .c = {.x =-1.25,.y = 0.5} },
			{.r = 0.15, .c = {.x = 0.25,.y =0.25} },
			{.r = 0.50, .c = {.x = 1.3, .y = 1.3} },
			{.r = 0.85, .c = {.x = 1.4, .y =-1.4} },
			//{.r = 0.2, .c = {.x = 1.0, .y = 2.0} },
		};

		int sizeCL = sizeof(cls) / sizeof(cls[0]);
		initializeInputDrops(cls, sizeCL);
/**
Then, the Hamaker field values are defined to introduce substrate
heterogeneities.
*/
		defaultHeterogeneityExpression(&c_Ham);
		//heterogeneity(true, 4. * (noise() + 1.0) + 1.0 );

		{
			scalar tempH[], hdist[];
			srand(300*(mpi_rank+1));
			initializeField(tempH, true, 4. * (noise() + 1.0) + 1.0);


			/**
				Preform 30x laplacian smoothing
			*/
			for (int i = 0; i < 30; i++)
			{
				foreach()
				{
					hdist[] = (
							tempH[1] +
							tempH[-1] +
							tempH[ 0, 1] +
							tempH[ 0,-1] +
							tempH[ 1, 1] +
							tempH[-1,-1] +
							tempH[-1, 1] +
							tempH[ 1,-1] ) / 8.0;
				}

				foreach()
				{
					tempH[] = hdist[];
				}
			}

			stats sh = statsf(hdist);
			// range min, max
			double rmax = 4.0 * precFilm;
			double rmin = 0.0;
			double invDiff = (rmax - rmin) / (sh.max - sh.min);

			foreach()
			{
				//normalize hdist between [rmax, rmin]
				hdist[] = (hdist[] - sh.min) * invDiff + rmin;
				// add the disturbance to h field
				h[] += hdist[]; 
			}
		}


		{
			scalar tempHam[];
			srand(100*(mpi_rank+1));
			initializeField(tempHam, true,4. * (noise() + 1.0) + 1.0);


			/**
				Preform 2x laplacian smoothing
			*/
			for (int i = 0; i < 22; i++)
			{
				foreach()
				{
					hamaker[] = (
							tempHam[1] +
							tempHam[-1] +
							tempHam[ 0, 1] +
							tempHam[ 0,-1] +
							tempHam[ 1, 1] +
							tempHam[-1,-1] +
							tempHam[-1, 1] +
							tempHam[ 1,-1] ) / 8.0;
				}

				foreach()
				{
					tempHam[] = hamaker[];
				}
			}

			stats sh = statsf(hamaker);
			// range min, max
			double rmax = 9.0;
			double rmin = 1.0;
			double invDiff = (rmax - rmin) / (sh.max - sh.min);

			foreach()
			{
				hamaker[] = (hamaker[] - sh.min) * invDiff + rmin;
			}
		}

		foreach_face()
		{
			gradHet.x[] = face_gradient_x(hamaker, 0);
		}

		restriction({hamaker, gradHet});

/**
Similarly, substrate roughness is accordingly defined.
*/

#if ROUGH
		foreach()
		{
			s[] = 0.0;
		}

		roughness( true , 2.0 * precFilm * (noise() + 1.0) );

		boundary({s});

		computeGradAndNegLaplace(s, gs, lapS);
		restriction({s, lapS, gs});
		gs.x.nodump = gs.y.nodump=true;
#endif

		boundary({h, hamaker});
		hold.nodump=true;
		gradHet.x.nodump = gradHet.y.nodump=true;

/**
The initialization fields are written based on the selected format, defined in
header file ```outputManager.h```.
*/
		outputData(fList, dirFlowFields, suffix="Init");
		//outputData({hamaker}, dirFlowFields, suffix="Init");
		{
		    output_matrix_MPI(hamaker, "hamakerInitial.bin", dirFlowFields, linear=true);
		}

		{
		    char b[100]; sprintf(b, "%shamakerInitial.dat", dirFlowFields);
		    FILE * fp = fopen (b, "w");
		    output_field({hamaker}, fp, linear=true);
		    fclose(fp);
		}
/**
and create a directory to store the output files and write the generated
Hamaker field.
*/
		outputInitalPNGs(dirFigures);

		fprintf(fout, "# New Simulation with %ld cells\n", grid->tn);
	}
	else
	{
		backUpContinuationMovie(dirFigures);

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
	LOG;
}

int main (int argc, char *argv[])
{
 	Tend = 60.025;
	maxDT = 1.e-4, minDT = 1.e-9;

	// defined default values
	precFilm = 3.e-3;
	inclineA_d = 30.0;
	contactA_d = 30.0;
	bond = 2.2;

	// timestepFactor = 1.04;
	lMGcycles = 12, uMGcycles = 22;
	tolMG = 1.e-3, tolNewton = 1.e-10;
	nitermax = 25;
	printStatsEvery = 10;

	//check for maximum run time
	maxruntime(&argc, argv);
	//check for input file with extension "*.dat"
	setupFromFile(&argc, argv);

	//BC
	foreach_dimension()
		periodic(right);
 	//h[right ] = dirichlet(precFilm);
 	//h[top   ] = dirichlet(precFilm);
 	//h[bottom] = dirichlet(precFilm);
 	//h[left  ] = dirichlet(precFilm);
	//domain size

	L0 = 6.0;

	//center computational domain
	origin( -0.5 * L0 , -0.5 * L0 );

	N = 1<<LEVEL_MAX;

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
	LOG;
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

event movie(t+=0.5, last)
{
#if 1
	char timestring[100]; sprintf(timestring, "t=%06.02f", t);

	//create vertex scalar hv
	vertex scalar hv[];
	foreach_vertex()
		hv[] = (h[] + h[-1] + h[0,-1] + h[-1,-1])/4.;

	view(width=800, height=800);
	squares("h", spread=-1, linear=1, map=gray);
	isoline("hv", 1.2 * precFilm, linear=1, lc={1,0,0}, lw=2.0);
//	isoline("hv", 1.4 * precFilm, linear=1, lc={0,1,0}, lw=2.0);
//	isoline("hv", 2.0 * precFilm, linear=1, lc={1,1,1}, lw=2.0);

	draw_string(timestring, lc={0,0,0}, lw = 2.0);

	char b[80]; sprintf(b, "%s/animation.mp4", dirFigures);
	save(b);
#endif
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
	LOG;
}
#endif

/**
Binary files of the current solution can be written
*/

event outputSolution (t+=0.5, last)
{
	outputData({h}, dirFlowFields, t);
}

/**
This event terminates the simulation and creates a dump file.
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
	astats s = adapt_wavelet({invH},(double[]){0.25}, LEVEL_MAX, LEVEL_IN-1);
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

