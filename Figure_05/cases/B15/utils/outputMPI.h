struct OutputMatrixMPI
{
	scalar f;
	FILE *fp;
	int n;
	bool linear;
};

trace
void output_matrix_MPI (struct OutputMatrixMPI p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;

  float fn = p.n;
  float Delta = (float) L0/fn;

  if (p.linear) {
    scalar f = p.f;
    boundary ({f});
  }
  float ** field = (float **) matrix_new (p.n, p.n, sizeof(float));

	//computing field
  for (int i = 0; i < p.n; i++)
	{
    float x = Delta*i + X0 + Delta / 2.0;
    for (int j = 0; j < p.n; j++)
		{
      float y = Delta*j + Y0 + Delta / 2.0;

			if (p.linear)
			{
				field[i][j] = interpolate(p.f, x, y);
			}
			else
			{
				Point point = locate (x, y);
	  		field[i][j] = point.level >= 0 ? val(p.f) : nodata;
			}
		}
	}

	// writing binary file
	if (pid() == 0) //master
	{
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif
		// size
  	fwrite (&fn, sizeof(float), 1, p.fp);

		//y-coords
	  for (int j = 0; j < p.n; j++)
		{
	    float yp = (float) (Delta*j + Y0 + Delta/2.);
	    fwrite (&yp, sizeof(float), 1, p.fp);
	  }
		//x-coord + all vals
  	for (int i = 0; i < p.n; i++)
		{
  	  float xp = (float) (Delta*i + X0 + Delta/2.);
  	  fwrite (&xp, sizeof(float), 1, p.fp);
  	  for (int j = 0; j < p.n; j++)
			{
  	    fwrite (&(field[i][j]), sizeof(float), 1, p.fp);
			}
		}
		fflush(p.fp);
	}
@if _MPI
	else //slave
	{
    MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,
		MPI_COMM_WORLD);
	}
@endif
  matrix_free (field);
}
