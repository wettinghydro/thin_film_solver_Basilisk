/**
# Purpose
Small header file that defined specific output functions
1. ```output_data``` : output binary adapted meshes \
(lvl, x, y, z, $\Delta x$, $\phi$) using MPI
1. ```outputSlice``` : output slice of $2D$ scalar
1. ```outputPNGSnapshot``` : output_ppm wrapper
1. ```writeContactLinesFile``` : outputs current contact lines
*/
struct OutputData
{
	const scalar *f_list; //scalars to print
	const char   *dir;    //directory Name
	const double	t;		  //time
	const char	 *suffix; //optional for filename.
};

@if !_MPI
trace
void output_data (struct OutputData p)
{
/**
outputs are written in a separate directory (```p.dir```) to avoid pollution
*/
	for (scalar f in p.f_list)
	{
		char filename[60];
		sprintf(filename, "%sdata_%s_%06.02f.bin", p.dir, f.name, p.t);

		FILE * fh = fopen(filename,"w");

		boundary({f});

		foreach()
		{
			const double lvl = (double) level;
			fwrite(&lvl,sizeof(double),1,fh);
			fwrite(&x  ,sizeof(double),1,fh);
			fwrite(&y  ,sizeof(double),1,fh);
#if dimension==3
			fwrite(&z,sizeof(double),1,fh);
#endif
			fwrite(&Delta,sizeof(double),1,fh);
			fwrite(&f[]  ,sizeof(double),1,fh);
		}
		fclose(fh);
	}
}
@else // _MPI
trace
void output_data(struct OutputData p)
{
	for (scalar f in p.f_list)
	{
		char filename[60];
		sprintf(filename, "%sdata_%s_%06.02f.bin", p.dir, f.name, p.t);

		FILE * fh = fopen(filename,"w");

		boundary({f});

		scalar index = {-1};
		index = new scalar;
		z_indexing (index,true);

#if dimension==2
		const int len = 5*sizeof(double);
#elif dimension==3
		const int len = 6*sizeof(double);
#endif

		foreach()
		{
		  if (is_local(cell))
			{
		  	const int offset = index[]*len;
		  	fseek(fh,offset,SEEK_SET);
		  	const double lvl = (double) level;

		  	fwrite(&lvl,sizeof(double),1,fh);
		  	fwrite(&x,sizeof(double),1,fh);
		  	fwrite(&y,sizeof(double),1,fh);
#if dimension==3
				fwrite(&z,sizeof(double),1,fh);
#endif
	  		fwrite(&Delta,sizeof(double),1,fh);
	  		fwrite(&f[],sizeof(double),1,fh);
			}
		}

		delete({index});
		fclose(fh);
	}
}
@endif // _MPI


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
/**
# Output of initial setup.
*/

void outputInitalPNGs(char *dirFigures)
{
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
		cells(lw=0.5);
		draw_string("t=0", lc={0,0,0}, lw = 2);

		char b[100]; sprintf(b, "%shinitial.png", dirFigures);
		save(b);
	}
#else
		outputPNGSnapshot(hamaker, dirFigures, "Init");
		outputPNGSnapshot(h, dirFigures, "Init");
#if ROUGH
		outputPNGSnapshot(s, dirFigures, "Init", spread=-1);
#endif
#endif
}
