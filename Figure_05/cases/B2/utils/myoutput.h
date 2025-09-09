/**
# Purpose
Small header file that outputs the *unstructured* ```scalar``` fields in binary
form on adapted meshes.
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
