/**
# Intro
Header file that can be used to output [Tecplot](https://www.tecplot.com)
compatible binary files of multiple ```scalar``` and ```vertex scalar``` fields.
Still under--development with some missing features.

For example:

	* only for 2D cases.

	* difficulty with periodic boundary conditions that have adaptation along the boundaries.

	* Testing was mainly focused on scalar fields. Might have some bugs for other fields.

	* not MPI combatible.

	* Hanging nodes along the bigger neighbor are not included in its connectivity. That would require additional connectivity for each face and a FEpolyhedron logic.

Performance considerations:

	* Currently based on the generation and storage of a dual connectivity mesh - memory intensive and can be slow.
*/
#include <TECIO.h>

/**
The main idea is to find the dominant node at each vertex. When adaptation
occurs, vertices between refinement levels are duplicated and need be
identified and then *unified*.
*/
#define uniqueNode(cell, point) (\
			(!is_leaf(aparent(0)) && is_leaf(cell)) \
		  || (is_prolongation(cell) && (point.i + point.j)%2) \
		  || is_boundary(cell))

/**
A ```struct``` that imitates the structure of the quad/octree and stores the
uniquely enumerated vertices.
*/
struct Ar_wrapper
{
	int** _p;
	int _lay;
};

/**
Multiple utility functions that are used to traverse the connectivity tree.
We define 3 *counters*. One using ```i,j,k```, one using a single ```ID``` and
a ```_Point``` based one.
*/
//utility functions
int  lineN(int i, int j, int l)
{
	return _size(l) * j + i;
}
//getters
int  getID(struct Ar_wrapper *ar, int ID, int l)
{
	return ar->_p[l][ID];
}
int  getIJ(struct Ar_wrapper *ar, int i, int j, int l)
{
	return getID(ar, lineN(i, j, l), l);
}
int  getp (struct Ar_wrapper *ar, struct _Point *p)
{
	return getID(ar, lineN(p->i, p->j, p->level), p->level);
}
//setters
void setID(struct Ar_wrapper *ar, int ID, int l, int val)
{
	ar->_p[l][ID] = val;
}
void setIJ(struct Ar_wrapper *ar, int i, int j, int l, int val)
{
	setID(ar, lineN(i, j, l), l, val);
}
void setp (struct Ar_wrapper *ar, struct _Point *p, int val)
{
	setID(ar, lineN(p->i, p->j, p->level), p->level, val);
}

/**
Initialize the tree structure. *Could be done more memory efficient*
*/
struct Ar_wrapper *ar_init()
{
	struct Ar_wrapper *p = qmalloc(1, struct Ar_wrapper);
	if (!p) { fprintf(ferr, "AR Allocation failed\n"); exit(1);}

	p->_lay = grid->maxdepth+1;
	p->_p = qmalloc(p->_lay, int*);

	if (!p->_p) { fprintf(ferr, "AR Allocation failed\n"); exit(1);}

	for (int i = 0; i < p->_lay; i++)
	{
		p->_p[i] = qmalloc(poolsize(i, 1), int);
		if (!p->_p[i]) {fprintf(ferr, "AR[%d] Allocation failed\n", i); exit(1);}
	}

	for (int l = 0; l < p->_lay; l++)
	{
			for (int i = 0; i<poolsize(l, 1); i++)
			{
				setID(p, i, l, -1);
			}
		}
	return p;
}

/**
delete the tree structure pointers.
*/
void ar_del(struct Ar_wrapper *ar)
{

	for (int i = 0; i < ar->_lay; i++)
	{
		free(ar->_p[i]);
		ar->_p[i] = NULL;
	}

	free(ar->_p); ar->_p = NULL;

	free(ar); ar = NULL;
}



/**
# User Interface
A ```struct``` that facilitates the exchange of data. tec_cc is the
```scalar``` list to be outputted.
*/
struct OutputTec
{
  const scalar *cc; 		//scalar list
	const char   *dir;    //directory Name
	const double	t;			//time.
	const char	 *suffix; //optional for filename.
};

trace
void output_tec(struct OutputTec tec)
{
	INTEGER4 nNodes,nCells,
	         idebug,itec,dIsDouble,vIsDouble,zoneType,
	         strandID,parentZn,isBlock,
	         nFConns, fNMode, shrConn, fileType, nFaces,
	         connectivityCount, iCellMax, jCellMax, kCellMax;
	INTEGER4 fileFormat;

#if _MPI
	assert(1==0);
#endif

	char varNames[256] = "x y";
	for (scalar s in tec.cc)
	{
		strcat(varNames, " ");
		strcat(varNames, s.name);
	}

	struct Ar_wrapper* treeStructure = ar_init();

	nCells = tree->leaves.n;

/**
# Connectivity creation
enumerate all vertices and store their number in our tree.
Both unique and non--unique.
*/
	int vertSize = 0;
	foreach_vertex(serial)
	{
	 		setp(treeStructure, &point, vertSize);
	 		vertSize++;
	}

	int pos = vertSize;
/**
Deal with Periodic boundaries. When ```periodic``` is defined,
[basilisk](Front Page) eliminates a row/column of vertices and connects the
interior nodes. Here we need to re--create those rows.
*/

	foreach_vertex(serial)
	{
		if ( Period.x && (point.i == GHOSTS) )
			vertSize ++;

		if ( Period.y && (point.j == GHOSTS) )
			vertSize ++;

		if ( (point.i == GHOSTS && point.j == GHOSTS)
				&& Period.y && Period.x) vertSize ++ ;

#if dimension >2
		assert(1==0);
		//TODO
		if ( Period.z && (point.k == GHOSTS) )
			vertSize ++;
#endif
	}
	int * vID = qmalloc(vertSize, int);

/**
We need to find non--dominant vertices and *associate* them to the dominant one.
*/
	Array *swaps = array_new();

	int counter = 0;
	foreach_vertex(serial) {
		if ( uniqueNode(cell,point) )
/**
If the node is unique, give it an ID(```vpos```).
*/
		{
			const int vpos = getp(treeStructure, &point);
			vID[vpos] = counter;
			counter ++;
		} else if (is_leaf(aparent(0))) {
/**
If the node is not unique, check its parent. If unique, store both in our array
to associate them.
*/
			array_append(swaps, &point,  sizeof(struct _Point));
			array_append(swaps, &parent, sizeof(struct _Point));
		} else {
/**
If neither the node or its parent is unique, then it is the child.
*/
			struct _Point c = {0};
			c.i = 2*point.i-GHOSTS;
			c.j = 2*point.j-GHOSTS;
#if dimension>2
			c.k = 2*point.k-GHOSTS;
#endif
			c.level = point.level + 1;

			array_append(swaps, &point, sizeof(struct _Point));
			array_append(swaps, &c,     sizeof(struct _Point));
		}
	}

	const int size = swaps->len / sizeof(struct _Point);
	struct _Point * ps = (struct _Point *)array_shrink(swaps);

/**
Unify numbering. Only keep dominant vertex ID by overwritting non--dominant one
enumeration in the tree structure.
*/
	for (int _i = 0; _i < size; _i +=2)
	{
		struct _Point duplicate = ps[_i];
		struct _Point unique    = ps[_i+1];

		const int vpos = getp(treeStructure, &unique);
		setp(treeStructure, &duplicate, vpos);
	}
	free(ps); ps = NULL;

/**
Again, periodicity complicates the process.
*/
	foreach_vertex(serial)
	{
		const int idxMax = _size(point.level) - GHOSTS;
		if (Period.x)
		{
			if (point.i == GHOSTS)
		 	{
				const int vpos = getIJ(treeStructure, idxMax, point.j, point.level);
				if (vpos < 0 )
				{ // unspecified
					setIJ(treeStructure, idxMax, point.j, point.level, pos);
					vID[pos] = counter;
					counter ++; pos ++;
				}
				//corner node
				if (point.j == GHOSTS && Period.y)
				{
					const int vpos = getIJ(treeStructure, idxMax, idxMax, point.level);
					if (vpos < 0 )
					{ // unspecified
						setIJ(treeStructure, idxMax, idxMax, point.level, pos);
						vID[pos] = counter;
						counter ++; pos ++;
					}
				}
			}
		}
		if (Period.y)
		{
			if (point.j == GHOSTS )
			{
				const int vpos = getIJ(treeStructure, point.i, idxMax, point.level);
				if (vpos < 0 )
				{ // unspecified
					setIJ(treeStructure, point.i, idxMax, point.level, pos);
					vID[pos] = counter;
					counter ++; pos ++;
				}
			}
		}
#if dimension>2
		assert(1==0);
		//TODO
		fprintf(ferr,"periodicity not done\n");
		exit(1);
#endif
	}
	nNodes = counter;

	const int numcc = list_len((scalar *)tec.cc);

/**
"dimension" corresponds to the number of the coordinate axes,
numcc is the number of variables at cell center.
*/
	int valueLocation[dimension + numcc];

	for(int k=0;k<dimension;k++)
		valueLocation[k]=1; //x,y (,z) are given at vertexes.

	for(int k=dimension;k<dimension+numcc;k++)
		valueLocation[k]=0; //Scalars are given at cell centers.

	int *connectivity;
	double *xn,*yn,*sc;
	#if dimension ==3
	double *zn;
	#endif

	fileFormat = 0;
	fileType   = 0;
	idebug     = 0;
	vIsDouble  = 1;

	#if dimension==2
	zoneType  = 3;      /* FEQuadrilateral */
	#elif dimension==3
	zoneType  = 5;      /* FEBrick */
	#endif
	nFaces    = 0;
	iCellMax  = 0;
	jCellMax  = 0;
	kCellMax  = 0;
	strandID  = 1;     /* StaticZone */
	parentZn  = 0;      /* No Parent */
	isBlock   = 1;      /* Block */
	nFConns   = 0;
	fNMode    = 0;
	dIsDouble = 1;
	shrConn   = 0;

/**
# Output File
Create the filename.
*/
	char filename[60];

	if (tec.suffix)
		sprintf(filename,"%sdata_%s.plt", tec.dir, tec.suffix);
	else
		sprintf(filename,"%sdata_step%06.02f_zone%d.plt", tec.dir, tec.t, pid());

	itec = TECINI142((char*)"BASILISK DATASET",
	                 varNames,
	                 filename,
	                 (char*)".",
	                 &fileFormat,
	                 &fileType,
	                 &idebug,
	                 &vIsDouble);

	if (itec)
	{
		fprintf(ferr, "Tecplot Header Failed\n");
		exit(1);
	}

	itec = TECZNE142((char*)"Primary Zone",
	                 &zoneType,
	                 &nNodes,
	                 &nCells,
	                 &nFaces,
	                 &iCellMax,
	                 &jCellMax,
	                 &kCellMax,
	                 &t,
	                 &strandID,
	                 &parentZn,
	                 &isBlock,
	                 &nFConns,
	                 &fNMode,
	                 0,              /* TotalNumFaceNodes */
	                 0,              /* NumConnectedBoundaryFaces */
	                 0,              /* TotalNumBoundaryConnections */
	                 0,//NULL,           /* PassiveVarList */
	                 valueLocation,  /* ValueLocation = Nodal */
	                 0,//NULL,           /* SharVarFromZone */
	                 &shrConn);

	if (itec)
 	{
		fprintf(ferr, "Tecplot Zone Failed\n");
		exit(1);
	}

/**
# Vertices Position
Store the unique vertex positions to temporary arrays.
 */
	xn  = (double*)malloc(nNodes * sizeof(double));
	yn  = (double*)malloc(nNodes * sizeof(double));
	#if dimension>2
	zn  = (double*)malloc(nNodes * sizeof(double));
	#endif
	int ID = 0;
	foreach_vertex(serial)
	{
		if ( uniqueNode(cell, point) )
		{
			xn[ID]=x;
			yn[ID]=y;
#if dimension>2
			zn[ID]=z;
#endif
			ID++;
		}
	}

/**
If periodic, build the missing nodes.
 */
	foreach_vertex(serial)
	{
		if (Period.x)
	 	{ // add Coord
			if (point.i == GHOSTS)
			{
				const int mirrorID = getp(treeStructure, &point);
				xn[ID] = X0 + L0;
				yn[ID] = yn[vID[mirrorID]];
				ID ++;
				if (point.j == GHOSTS && Period.y)
				{
					xn[ID] = X0 + L0;
					yn[ID] = Y0 + L0;
					ID ++;
				}
			}
		}
		if (Period.y)
		{ // add Coord
			if (point.j == GHOSTS )
			{
				const int mirrorID = getp(treeStructure, &point);
				xn[ID] = xn[vID[mirrorID]];
				yn[ID] = Y0 + L0;
				ID ++;
			}
		}
	}

/**
Use Tecplot functions to print the data and ```free``` the arrays.
 */
	if (TECDAT142(&nNodes, xn, &dIsDouble))
	{
		fprintf(ferr, "X-nodes Failed\n"); exit(1);
	}
	if (TECDAT142(&nNodes, yn, &dIsDouble))
	{
		fprintf(ferr, "Y-nodes Failed\n");
		exit(1);
	}
	#if dimension>2
	if (TECDAT142(&nNodes, zn, &dIsDouble))
	{
		fprintf(ferr, "Z-nodes Failed\n");
		exit(1);
	}
	#endif
	free(xn); xn = NULL;
	free(yn); yn = NULL;
	#if dimension>2
	free(zn); zn = NULL;
	#endif

/**
# Scalars
Write the ```scalar``` fields
 */
	sc  = (double*)malloc(nCells * sizeof(double));
	for (scalar s in tec.cc)
	{
		int IDc = 0;
	  foreach(serial)
		{
	     sc[IDc]=s[];
			 IDc ++ ;
	  }
	  if(TECDAT142(&nCells, sc, &dIsDouble))
		{
			fprintf(ferr, "Writing Scalar Field Failed\n");
			exit(1);
		}
	}
	free(sc); sc = NULL;

	#if dimension==2
	connectivityCount = 4 * nCells;
	#elif dimension==3
	connectivityCount = 8 * nCells;
	#endif
	connectivity = qmalloc(connectivityCount, INTEGER4);

/**
# Connectivity Usage
Now we use the connectivity that was computed.
Specify the connectivity between the IDs of cell center and vertex
counter-clockwise orientation
*/
	int _ID=0;
	foreach(serial)
	{
		INTEGER4 itecc=(INTEGER4)_ID;
		_ID ++;
		const int ii = point.i;
		const int jj = point.j;
		const int ll = point.level;
		#if dimension==2
		//add one since node numbering starts from 1
		const int v1 = vID[getIJ(treeStructure, ii  ,jj  ,ll)] + 1;
		const int v2 = vID[getIJ(treeStructure, ii+1,jj  ,ll)] + 1;
		const int v3 = vID[getIJ(treeStructure, ii+1,jj+1,ll)] + 1;
		const int v4 = vID[getIJ(treeStructure, ii  ,jj+1,ll)] + 1;

		if (sign(v1) * sign(v2) * sign(v3) * sign(v4) < 0)
		{
			fprintf(ferr,
					"Unassigned node %g %g \n %d %d %d\n %d %d %d %d\n %d %d %d %d\n",
				x,y, ii, jj, ll,
				getIJ(treeStructure, ii  ,jj  ,ll),
				getIJ(treeStructure, ii+1,jj  ,ll),
				getIJ(treeStructure, ii+1,jj+1,ll),
				getIJ(treeStructure, ii  ,jj+1,ll),
				v1, v2, v3, v4);
			exit(1);
		}

		//cast to INTEGER4
		connectivity[itecc*4]   = (INTEGER4)v1;
		connectivity[itecc*4+1] = (INTEGER4)v2;
		connectivity[itecc*4+2] = (INTEGER4)v3;
		connectivity[itecc*4+3] = (INTEGER4)v4;
		#elif dimension==3
		assert(1==0);
		// TODO
		//connectivity[itecc*8]   =(INTEGER4)idv[0,0,0];
		//connectivity[itecc*8+1] =(INTEGER4)idv[1,0,0];
		//connectivity[itecc*8+2] =(INTEGER4)idv[1,1,0];
		//connectivity[itecc*8+3] =(INTEGER4)idv[0,1,0];
		//connectivity[itecc*8+4] =(INTEGER4)idv[0,0,1];
		//connectivity[itecc*8+5] =(INTEGER4)idv[1,0,1];
		//connectivity[itecc*8+6] =(INTEGER4)idv[1,1,1];
		//connectivity[itecc*8+7] =(INTEGER4)idv[0,1,1];
		#endif
	}

/**
Memory Cleanup
*/
	if (TECNODE142(&connectivityCount, connectivity))
	{
		fprintf(ferr, "Writing Connectivity Failed\n");
		exit(1);
	}
	free(connectivity); connectivity = NULL;

	if (TECEND142())
	{
		fprintf(ferr, "Closing plt file Failed\n");
		exit(1);
	}

	free(vID); vID = NULL;
	ar_del(treeStructure);
}
