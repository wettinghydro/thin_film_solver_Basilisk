/**
# Output Manager
Interface to output simulation data.
*/
#include <sys/stat.h>

#ifdef TEC
	#error -DTEC is obsolete
	//#include "tecplot.h" // obsolete
#elif PRINT
	#include "myoutput.h"
#endif /* TEC || PRINT */

/**
Store the Directory name. Give a default value for safety
*/
struct DirCreator
{
	char *dirName;
};

void createDirectory(struct DirCreator p)
{
	char *file = p.dirName ? p.dirName : "default";

/**
Append directory suffix if not there
*/
	if (file[strlen(file) - 1] != '/')
		strcat(file, "/");

	struct stat st = {0};

/**
Create if does not exists
*/
	if (stat(file, &st) == -1)
		mkdir(file, 0700);
}

void addOutputFileHeader()
{
	int n = 0;
  fprintf(fout, "#%d:i", ++n);
  fprintf(fout, "\t%d:log10(dh/e)", ++n);
  fprintf(fout, "\t%d:min(h)\t%d:max(h)\t%d:sum(h*Omega)", n+1, n+2, n+3);
  n+=3;
#if ELECTRO
  fprintf(fout, "\t%d:log10(dq/e)", ++n);
  fprintf(fout, "\t%d:min(q)\t%d:max(q)\t%d:sum(q*Omega)", n+1, n+2, n+3);
  n+=3;
#endif
  fprintf(fout, "\t%d:t\t%d:dt\t%d:mgIterSum\n", n+1, n+2, n+3);

	fflush(fout);
}

void reportSimulationEffects()
{
	fprintf(fout, "# Simulation with \n");
#if BDF2
	fprintf(fout, "# \t 2nd Order Time accuracy \n");
#endif
#if ROUGH
	fprintf(fout, "# \t Surface Roughness \n");
#endif
#if INCLINED
	fprintf(fout, "# \t Inclination effects\n");
#endif
#if FORCING
	fprintf(fout, "# \t Vibration effects\n");
#endif
#if EVAPORATION
	fprintf(fout, "# \t Evaporation effects\n");
#endif
#if ELECTRO
	fprintf(fout, "# \t ElectroWetting effects\n");
#endif
#if PRINT
	fprintf(fout, "# \t Binary field Outputs \n");
#endif

	addOutputFileHeader();
}

/**
# User Interface
These could go in myoutput.h
*/

struct OutputInterface
{
	const scalar *fields; //data to output
	const char   *dir;		//optional for filename.
	const double  time;   //time
	const char   *suffix; //optional for filename.
};

void outputData(struct OutputInterface p)
{
	if (p.dir[strlen(p.dir) - 1] != '/')
	{
		fprintf(ferr, "dir member in function %s should end in \'/\'\n", __func__);
		exit(1);
	}

#ifdef TEC
	if (p.suffix) output_tec (p.fields, p.dir, suffix=p.suffix);
	else 					output_tec (p.fields, p.dir, p.time);
#elif PRINT
	if (p.suffix) output_data(p.fields, p.dir, suffix=p.suffix);
	else 					output_data(p.fields, p.dir, p.time);
#endif /* TEC || PRINT */
}

struct Slicer
{
	scalar f;
  char * fileName;
	int    count;  //number of point
	bool linear;
	coord  n;	    //normal plane
	double alpha; // n \cdot x = a
};

void outputSlice(struct Slicer p)
{
	if (p.count == 0 ) p.count = N;

	char def[] = "defSlice";
	if (!p.fileName) p.fileName = def;
	if (!p.alpha) p.alpha = 0.0;

	FILE *fp = fopen(p.fileName, "w");

	if (fp == NULL)
	{
		perror(p.fileName); exit(1);
	}

	normalize(&p.n);
	double fn = p.count;
	double Delta = L0/fn;

  if (p.linear) {
    scalar f = p.f;
    boundary ({f});
  }

	double *field = qmalloc(2 * p.count, double); //distance along plane + value

	//compute slice values
	for (int i = 0; i < p.count; i++)
	{
		double xC, yC;
		if (p.n.y > p.n.x) // fixme: general and optimize
		{
			xC = Delta * i + X0 + 0.5 * Delta;
			yC = (p.alpha - p.n.x * xC) / p.n.y;

			field[2*i] = sign(xC) * sqrt( sq(xC) + sq(yC) );
		} else {
			yC = Delta * i + Y0 + 0.5 * Delta;
			xC = (p.alpha - p.n.y * yC) / p.n.x;
			field[2*i] = sign(yC) * sqrt( sq(xC) + sq(yC) );
		}

		if (0) // interpolate crashes profiling
		{
			field[2*i +1] = interpolate(p.f, xC, yC);
		}
		else
		{
			Point point = locate (xC, yC);
			field[2*i+1] = point.level >= 0 ? val(p.f) : nodata;
		}
	}
	// writing ASCI file
	if (pid() == 0) //master
	{
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field, 2*p.count, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif
		fprintf(fp, "#n:(%g %g), alpha=%g\n", p.n.x, p.n.y, p.alpha);
		fprintf(fp, "#1:t, 2:val:\n");
		//t-coords
	  for (int j = 0; j < p.count; j++)
		{
			fprintf(fp, "%g %g\n", field[2*j], field[2*j+1]);
	  }
		fflush(fp);
	}
@if _MPI
	else //slave
	{
    MPI_Reduce (field, NULL, 2*p.count, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
	}
@endif

  free(field);
	fclose(fp);
}

struct Snapper
{
	scalar f;
	char  *dir;
  char  *suffix;
	int spread;
};

void outputPNGSnapshot(struct Snapper p)
{
	char b[100];
	sprintf(b, "%s%s%s.png", p.dir, p.f.name, p.suffix);

	if (p.spread < 0)
	{
		output_ppm(p.f, spread=-1, file=b);
	}
	else
	{
		stats x = statsf(p.f);

		double defmin = (x.min == x.max) ? x.max -1 : x.min;
		output_ppm(p.f, min=defmin, max=x.max, file=b);
	}
}

/**
At different time intervals, the current solution is processed to compute and
output the contact lines
at $h=\varepsilon+\epsilon \approx 1.2 \varepsilon$ using ```fractions```. The
extra $0.2$ is used to avoid numerical oscilations on the contact lines due to
the existence of the precursor film.
Differences in the *exact* contact lines are negligible.
*/

struct PostProcess
{
	scalar f;     // scalar field used to compute contact lines
	char  *fname; // file
	double value; // intersection limit
};

void writeContactLinesFile (struct PostProcess pos)
{
	scalar f = pos.f;

	vertex scalar cl[];
	scalar tc[];

	foreach_vertex()
	{
		cl[] = (f[] + f[-1] + f[0,-1] + f[-1,-1])/4.;
	}

	face vector tf[];
	fractions(cl, tc, tf, pos.value);

/**
If ```MPI``` is used, each processor writes its own part of the contact line
and are subsequently merged.
*/
#if _MPI
	char tempName[25]; sprintf(tempName, "tempOut%d", pid());

	FILE *fp = fopen(tempName, "w");

	if (fp == NULL)
	{
		perror(tempName);
		exit(1);
	}

	output_facets(tc, fp, tf);

	fclose(fp);

	if (pid() == 0) // master
	{
		MPI_Barrier(MPI_COMM_WORLD);

		char com[100]; sprintf(com, "cat tempOut* > %s && rm tempOut* ", pos.fname);

		system(com); // merge parallel output and cleanup.
	}
	else // slave
	{
		MPI_Barrier(MPI_COMM_WORLD);
	}
#else // !_MPI

	FILE *fp = fopen(pos.fname, "w");
	if (fp == NULL)
	{
		perror(pos.fname);
		exit(1);
	}

	output_facets(tc, fp, tf);

	fclose(fp);
#endif
}
