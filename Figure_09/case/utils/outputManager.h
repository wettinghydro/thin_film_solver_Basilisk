/**
# Output Manager
Interface to output simulation data.
*/
#include <sys/stat.h>
#include "myoutput.h"

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

#if PRINT
	if (p.suffix) output_data(p.fields, p.dir, suffix=p.suffix);
	else 					output_data(p.fields, p.dir, p.time);
#endif /* PRINT */
}

void backUpContinuationMovie(const char* dir)
{
		/* Backup Saved animation */
		if (pid() == 0)
		{
			FILE * fp = NULL;
			int i = 0;
			bool done = false;

			while(i < 10 && !done)
			{
				char b[80]; sprintf(b, "%s/animation%d.mp4", dir, i);
				if ((fp = fopen(b, "r")))
				{
					fclose(fp);
					i++;
				}
				else
				{
					char com[100]; sprintf(com, "mv %s/animation.mp4 %s", dir, b);
					system(com); 
					done = true;
				}
			}

		}
}
