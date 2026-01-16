
#include <stdio.h>
#include <stdlib.h>


int main()
{
	fprintf(stdout, "Creating Input File\n");


	char fN[] = "parameters.dat";
	FILE *fp = fopen(fN, "w");

	fprintf(fp, "#Fill in the following variables\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#-------General Constants-------\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#1. Precursor film height (precFilm) \n");
	fprintf(fp, "\n");
	fprintf(fp, "#2. Mobility exponent n\n");
	fprintf(fp, "\n");
	fprintf(fp, "#3. Interaction exponent m\n");
	fprintf(fp, "\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#-----Inclination Constants-----\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#1. Bond number \n");
	fprintf(fp, "\n");
	fprintf(fp, "#2. Substrate Angle (Degrees)\n");
	fprintf(fp, "\n");
	fprintf(fp, "#3. Contact Angle (Degrees)\n");
	fprintf(fp, "\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#-----Vibration Constants-----\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#1. Vibration amplitude \n");
	fprintf(fp, "\n");
	fprintf(fp, "#2. Force Angle ψ (Degrees)\n");
	fprintf(fp, "\n");
	fprintf(fp, "#3. Force Period \n");
	fprintf(fp, "\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#-----Evaporation Constants-----\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#1. kappa \n");
	fprintf(fp, "\n");
	fprintf(fp, "#2. delta \n");
	fprintf(fp, "\n");
	fprintf(fp, "#3. Evaporation factor (E) \n");
	fprintf(fp, "\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#-----Simulation Parameters-----\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#1. Tend \n");
	fprintf(fp, "\n");
	fprintf(fp, "#2. Level Min \n");
	fprintf(fp, "\n");
	fprintf(fp, "#3. Level Max \n");
	fprintf(fp, "\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#----Heterogeneity Parameters---\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#0-7. Parameters (same line)\n");
	fprintf(fp, "#1.p0 2.p1 3.p2 4.p3 5.p4 6.p5 7.p6 8.p7\n");
	fprintf(fp, "\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#----Droplet Initialization---\n");
	fprintf(fp, "#-------------------------------\n");
	fprintf(fp, "#1.a 2.r 3.c_x 4.c_y\n");
	fprintf(fp, "\n");
	fclose(fp);

	fprintf(stdout, "****Done\n");
	fprintf(stdout, "\"%s\" has been created!\n", fN);

	return EXIT_SUCCESS;
}
