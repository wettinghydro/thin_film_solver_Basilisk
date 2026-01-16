union CaseParameters {

	// ISO C99
	struct ParameterData {
		// general constants
		double precFilm, n, m;
		// inclination constants
		double bond, theta, contact;
		// forcing constants
		double famp, psiA_d, Tper;
		// evaporation constants
		double kappa, delta, Efactor;
		// simulation parameters
		double Tend;
		double lmin, lmax; // store as double to simplify
	} d;

	double v[15];
};

void readCaseParameters(char *inName, union CaseParameters *params)
{
	FILE* fp = fopen(inName, "r");

	if (fp == NULL) {
		perror(inName);
		exit(1);
	}

	size_t bsize = 256;
	char * b = (char *) malloc(bsize * sizeof(char));

	int i = 0;
	//read case parameters
	while (i < 15 && getline(&b, &bsize, fp) != -1)
	{
		if (b[0] == '#') continue; //skip comments

		params->v[i] = atof(b);
		i ++;
	}

	free(b);
	fclose(fp);
}

void setAndPrintCaseParameters(const struct ParameterData *p)
{
	fprintf(fout, "# Simulation Parameters:\n");
	if (p->precFilm || p->n || p->m)
	{
		precFilm = p->precFilm;
		n        = p->n;
		m        = p->m;
	}
	fprintf(fout, "#Phys:\n# \t\t (ε, n, m) \n#\t\t %g %g %g\n", precFilm, n, m);
#if INCLINED
	if (p->bond  || p->theta || p->contact)
	{
		bond       = p->bond;
		inclineA_d = p->theta;
		contactA_d = p->contact;
	}
	fprintf(fout, "#Incl:\n# \t\t (B, θ, a) \n#\t\t %g %g %g\n", bond, inclineA_d, contactA_d);
#endif
#if FORCING
	if (p->famp || p->psiA_d  || p->Tper)
	{
		rvib   = p->famp;
		psiA_d = p->psiA_d;
		Tper   = p->Tper;
	}
	fprintf(fout, "#Vibration:\n# \t\t (r, ψ, T) \n#\t\t %g %g %g\n", rvib, psiA_d, Tper);
#endif
#if EVAPORATION
	if (p->kappa || p->delta || p->Efactor)
	{
		kappa   = p->kappa;
		delta   = p->delta;
		Efactor = p->Efactor;
	}
	fprintf(fout, "#Evap:\n# \t\t (k, δ, E) \n#\t\t %g %g %g\n" , kappa, delta, Efactor);
#endif
	if (p->Tend || p->lmin || p->lmax)
	{
		Tend      =       p->Tend;
		LEVEL_IN  = (int) p->lmin;
		LEVEL_MAX = (int) p->lmax;
	}
	fprintf(fout, "#Mesh:\n# \t\t (Lmin, Lmax, Tend) \n#\t\t %d %d %g\n", LEVEL_IN, LEVEL_MAX, Tend);
}

void readHeterogeneityParameters(char *inName, struct HC_struct *p)
{
	FILE* fp = fopen(inName, "r");

	if (fp == NULL) {
		perror(inName);
		exit(1);
	}

	size_t bsize = 256;
	char * b = (char *) malloc(bsize * sizeof(char));

	int toSkip = 1;
	while (getline(&b, &bsize, fp) != -1)
	{
		if( strstr(b, "Heterogeneity") != NULL ) toSkip = 0;

		if (toSkip || b[0] == '#') continue; //skip comments

		sscanf(b, "%lf %lf %lf %lf %lf %lf %lf %lf",
				&p->p0, &p->p1, &p->p2, &p->p3, &p->p4, &p->p5,
				&p->p6, &p->p7);
	}

	free(b);
	fclose(fp);
}

void printHetParameters()
{
	fprintf(fout, "#Heterogeneity:\n# ");
	for (int i = 0; i < 8; i ++)
	{
		fprintf(fout, "%g\t", c_Ham.hamC[i]);
	}
	fprintf(fout, "\n");
}


void setupFromFile(int *argc, char *argv[])
{
	for (int i = 0; i < *argc; i++)
	{
		int len = strlen(argv[i]);
		if (!strcmp (&argv[i][len-4], ".dat") )
		{
			fprintf(fout, "#Reading Input file: %s\n", argv[i]);

			union CaseParameters params;

			readCaseParameters(argv[i], &params);

			setAndPrintCaseParameters(&(params.d));

			readHeterogeneityParameters(argv[1], &(c_Ham.hcs_Data));

			printHetParameters();

			// reset command line arguments
			*argc -= 1;
			for (int j = i; j < *argc; j++)
				argv[i] = argv[i+1];

		}
	}
}
