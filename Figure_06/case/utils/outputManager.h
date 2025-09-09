/**
# Output Manager
Interface to output simulation data.
*/
#include <sys/stat.h>
#include "myoutput.h"

#if TEC
	#error -DTEC is obsolete
	//#include "tecplot.h" // obsolete
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

#if TEC
	if (p.suffix) output_tec (p.fields, p.dir, suffix=p.suffix);
	else 					output_tec (p.fields, p.dir, p.time);
#elif PRINT
	if (p.suffix) output_data(p.fields, p.dir, suffix=p.suffix);
	else 					output_data(p.fields, p.dir, p.time);
#endif /* TEC || PRINT */
}
