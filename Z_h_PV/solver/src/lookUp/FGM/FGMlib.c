# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include <unistd.h>
# include <sys/stat.h>
# include "FGMlib.h"

# define MAX_LINE_LENGTH 256
# define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a > _b ? _a : _b; })
# define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; })

FGM * readFGM (const char filename[])
{
    char line[MAX_LINE_LENGTH];
    FILE *fid;
  
    int Ncv = 0;
    int *Ngrid = NULL;
    double *gridpower = NULL;
    int Nvar = 0;
    char (*varname)[VAR_NAME_LENGTH];
    double *data = NULL;

    FGM *fgm = NULL;
  
    int i,j;
    double f;

    /* Open the file */
    fid = fopen(filename, "r");
    if (fid == NULL) {
        fprintf(stderr, "Error: Failed to open FGM database file\n");
        exit(EXIT_FAILURE);
    }
  
    /* Check the identifier [FGM] */
    fgets(line,MAX_LINE_LENGTH,fid);
    if (strncmp(line,"[FGM]",5) != 0) {
        fprintf(stderr, "Error: Incorrect FGM database format\n");
        exit(EXIT_FAILURE);
    }
  
    /* Read the number of cv's */
    if (locateKeyWord(fid,"[DIMENSION]")) {
        fprintf(stderr, "Error: Couldn't find keyword DIMENSION\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fid,"%i",&Ncv);
    printf("Ncv = %i\n",Ncv);
  
    /* Read space for the data size */
    Ngrid = (int *)malloc(sizeof(int) * Ncv);
    if (Ngrid == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for data size\n");
        exit(EXIT_FAILURE);
    }
  
    /* Read the data size */
    if (locateKeyWord(fid,"[DATASIZE]")) {
        fprintf(stderr, "Error: Couldn't find keyword DATASIZE\n");
        exit(EXIT_FAILURE);
    }
    int Ntotal = 1;
    for (i = 0; i < Ncv; i++) {
        fscanf(fid,"%i",&j);
        Ngrid[i] = j;
        printf("Ngrid[%i] = %i\n",i,Ngrid[i]);
        Ntotal = Ntotal * Ngrid[i];
    }
  
    /* Read the number of variables */
    if (locateKeyWord(fid,"[VARIABLES]")) {
        fprintf(stderr, "Error: Couldn't find keyword VARIABLES\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fid,"%i",&Nvar);
    printf("Nvar = %i\n",Nvar);

    /* Allocate space for the variable names */
    varname = malloc(sizeof *varname * Nvar);
    if (varname == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for variable names\n");
        exit(EXIT_FAILURE);
    }
  
    /* Read the variable names */
    for (j = 0; j < Nvar; j++) {
        fscanf(fid,"%s",varname[j]);
        printf("%s\n",varname[j]);
    }
  
    /* Allocate space for the data */
    data = malloc(sizeof(double) * Ntotal * Nvar);
    if (data == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for data\n");
        exit(EXIT_FAILURE);
    }
  
    /* Read the data */  
    if (locateKeyWord(fid,"[DATA]")) {
        fprintf(stderr, "Error: Couldn't find keyword DATA\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < Ntotal; i++) {
        for (j = 0; j < Nvar; j++) {
            fscanf(fid,"%lf",&f);
            data[i*Nvar + j] = f;
        }
    }
  
    /* Allocate space for grid powers */
    gridpower = (double *)malloc(sizeof(double) * Ncv);
    if (gridpower == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for gridpower\n");
        exit(EXIT_FAILURE);
    }
  
    /* Read grid powers for non-linear distribution of mesh points */
    if (locateKeyWord(fid,"[GRIDPOWER]")) {
        fprintf(stderr, "Warning: Couldn't find keyword GRIDPOWER\n");
        for (i = 0; i < Ncv; i++) {
            gridpower[i] = -1.0; // Set to negative value
        }
    }
    else {
        for (i = 0; i < Ncv; i++) {
            fscanf(fid,"%lf",&f);
            // If the power is equal to 1 then set it to a negative value
            if ( fabs(f-1.0) < 1.e-6 ) {
                gridpower[i] = -1.0;
            }
            else {
                gridpower[i] = 1.0 / f; // Store the reciprocal values
            }
            printf("Gridpower[%i] = %f, %f\n", i, gridpower[i], f);
        }
    }
  
    /* Close file */
    fclose(fid);
  
    /* Allocate memory for FGM */
    fgm = (FGM *)malloc(sizeof(FGM));
    if (fgm == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for fgm\n");
        exit(EXIT_FAILURE);
    }
  
    /* Assign values */
    fgm->Ncv = Ncv;
    fgm->Ngrid = Ngrid;
    fgm->gridpower = gridpower;
    fgm->Nvar = Nvar;
    fgm->varname = varname;
    fgm->data = data;
  
    return fgm;
};

int freeFGM(FGM *fgm)
{
    if (fgm != NULL) {
        if (fgm->Ngrid != NULL) { free(fgm->Ngrid); }
        if (fgm->gridpower != NULL) { free(fgm->gridpower); }
        if (fgm->data != NULL) { free(fgm->data); }
        free(fgm);
    } else {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
};

int locateKeyWord(FILE *fid, char keyword[])
{
    char line[MAX_LINE_LENGTH];
    int keywordlength;
  
    keywordlength = strlen(keyword);
  
    rewind(fid);
    do {
        if (fgets(line,MAX_LINE_LENGTH,fid) == NULL) {
            return EXIT_FAILURE;
        }
    } while (strncmp(line,keyword,keywordlength) != 0);
  
    return EXIT_SUCCESS;
};

int lookupFGM_1D (FGM *fgm, double *x, double *f)
{
  
    int N1 = fgm->Ngrid[0]-1;
    int Nvar = fgm->Nvar;
    int ivar, i1m, i1p;
    double eta1, w1m, w1p;
    
    // Compute extremes of cv1
    double xmin = fgm->data[0];
    double xmax = fgm->data[N1*Nvar];
    
    // Compute eta1 by normalising cv1
    eta1 = (*x - xmin) /  (xmax - xmin);
    if (fgm->gridpower[0] >= 0.0) {
      eta1 = pow(eta1, fgm->gridpower[0]);
    }

    // Determine mesh point
    i1m = (int) (eta1 * N1);
    // Limit to [0,N-1]
    i1m = max(0,i1m);
    i1m = min(N1-1,i1m);
    i1p = i1m + 1;
    
    // Determine weight factors
    w1p = eta1 * N1 - i1m;
    w1m = 1 - w1p;
  
    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = w1m * fgm->data[i1m*Nvar+ivar];
    }
    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = *(f + ivar) + w1p * fgm->data[(i1p)*Nvar+ivar];
    }
  
    return EXIT_SUCCESS;
};

int lookupFGM_2D (FGM *fgm, double *x, double *f)
{    
    int N1 = fgm->Ngrid[0]-1;
    int N1p = N1 + 1;
    int N2 = fgm->Ngrid[1]-1;
    int Nvar = fgm->Nvar;
    int ivar, i1m, i1p, i2m, i2p;
    double eta1, eta2, w1m, w1p, w2m, w2p;
    
    // Compute extremes of cv1
    double xmin = fgm->data[0];
    double xmax = fgm->data[N1*Nvar]; 

    // Compute eta1 by normalising cv1
    eta1 = (*x - xmin) /  (xmax - xmin + 1e-8);
    if (fgm->gridpower[0] >= 0.0) {
      eta1 = pow(eta1, fgm->gridpower[0]);
    }
    
    // Determine mesh point
    i1m = (int) (eta1 * N1);
    // Limit to [0,N-1]
    i1m = max(0,i1m);
    i1m = min(N1-1,i1m);
    i1p = i1m + 1;
    
    // Determine weight factors
    w1p = eta1 * N1 - i1m;
    w1m = 1 - w1p;
    
    // Compute extremes of cv2
    xmin = w1m * fgm->data[i1m*Nvar+1]
         + w1p * fgm->data[i1p*Nvar+1];
    xmax = w1m * fgm->data[(N2*N1p+i1m)*Nvar+1]
         + w1p * fgm->data[(N2*N1p+i1p)*Nvar+1];

    // Compute eta2 from cv2
    eta2 = (*(x+1) - xmin) /  (xmax - xmin + 1e-8);
    if (fgm->gridpower[1] >= 0.0) {
      eta2 = pow(eta2, fgm->gridpower[1]);
    }
    
    // Determine mesh point
    i2m = (int) (eta2 * N2);
    i2m = max(0,i2m);
    i2m = min(N2-1,i2m);
    i2p = i2m + 1;

    // Determine weight factors
    w2p = eta2 * N2 - i2m;
    w2m = 1 - w2p;

    // Do the actual interpolation in eta space
    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = w1m * w2m * fgm->data[(i2m*N1p+i1m)*Nvar+ivar];
    }

    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = *(f + ivar) + w1p * w2m * fgm->data[(i2m*N1p+i1p)*Nvar+ivar];
    }

    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = *(f + ivar) + w1m * w2p * fgm->data[(i2p*N1p+i1m)*Nvar+ivar];
    }

    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = *(f + ivar) + w1p * w2p * fgm->data[(i2p*N1p+i1p)*Nvar+ivar];
    }
    
    return EXIT_SUCCESS;
};

int lookupFGM_2Dr (FGM *fgm, double *x, double *f)
{
    
    int N1 = fgm->Ngrid[0]-1;
    int N1p = N1 + 1;
    int N2 = fgm->Ngrid[1]-1;
    int Nvar = fgm->Nvar;
    int ivar, i1m, i1p, i2m, i2p;
    double eta1, eta2, w1m, w1p, w2m, w2p;
    
    /* Compute extremes of cv2 */
    double xmin = fgm->data[1];
    double xmax = fgm->data[N2*N1p*Nvar+1];
    
    /* Compute eta2 by normalising cv2 */
    eta2 = (*(x+1) - xmin) /  (xmax - xmin);
    if (fgm->gridpower[1] >= 0.0) {
      eta2 = pow(eta2, fgm->gridpower[1]);
    }
    /* printf("CV2: %e %e %e %e \n",*(x+1),xmin,xmax,eta2); */
    
    /* Determine mesh point*/
    i2m = (int) (eta2 * N2);
    /* Limit to [0,N-1]*/
    i2m = max(0,i2m);
    i2m = min(N2-1,i2m);
    i2p = i2m + 1;
    
    /* Determine weight factors*/
    w2p = eta2 * N2 - i2m;
    w2m = 1 - w2p;
    
    /* Compute extremes of cv1 */
    xmin = w2m * fgm->data[i2m*N1p*Nvar]
         + w2p * fgm->data[i2p*N1p*Nvar];
    xmax = w2m * fgm->data[(i2m*N1p+N1)*Nvar]
         + w2p * fgm->data[(i2p*N1p+N1)*Nvar];

    /* Compute eta1 from cv1 */
    eta1 = (*x - xmin) /  (xmax - xmin);
    if (fgm->gridpower[0] >= 0.0) {
      eta1 = pow(eta1, fgm->gridpower[0]);
    }
    /* printf("CV1: %e %e %e %e \n",*x,xmin,xmax,eta1); */

    /* Determine mesh point*/
    i1m = (int) (eta1 * N1);
    i1m = max(0,i1m);
    i1m = min(N1-1,i1m);
    i1p = i1m + 1;

    /* Determine weight factors*/
    w1p = eta1 * N1 - i1m;
    w1m = 1 - w1p;

    /* Do the actual interpolation in eta space*/
    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = w1m * w2m * fgm->data[(i2m*N1p+i1m)*Nvar+ivar];
    }

    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = *(f + ivar) + w1p * w2m * fgm->data[(i2m*N1p+i1p)*Nvar+ivar];
    }

    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = *(f + ivar) + w1m * w2p * fgm->data[(i2p*N1p+i1m)*Nvar+ivar];
    }

    for (ivar = 0; ivar < Nvar; ivar++) {
        *(f + ivar) = *(f + ivar) + w1p * w2p * fgm->data[(i2p*N1p+i1p)*Nvar+ivar];
    }
    
    return EXIT_SUCCESS;
};

int lookupFGM_3Dr(FGM *fgm, double *x, double *f)
{

	int N1 = fgm->Ngrid[0] - 1;
	int N1p = N1 + 1;
	int N2 = fgm->Ngrid[1] - 1;
	int N2p = N2 + 1;
	int N3 = fgm->Ngrid[2] - 1;
	int Nvar = fgm->Nvar;
	int ivar, i1m, i1p, i2m, i2p, i3m, i3p;
	double eta1, eta2, eta3, w1m, w1p, w2m, w2p, w3m, w3p;


        /* i1 = 0:N1, i2 = 0:N2, i3 = 0:N3, ivar = 0,Nvar-1 */
        /* Convert ijk index to linear index */
        /* (i1,i2,i3,ivar) => ((i3*N2p + i2)*N1p +i1)*Nvar + ivar */

	/* Compute extremes of cv3 (ivar = 2) */
        /* i1 = 0, i2 = 0, i3 = 0 */
	double xmin = fgm->data[2];
        /* i1 = 0, i2 = 0, i3 = N3 */
	double xmax = fgm->data[N3*N2p*N1p*Nvar + 2];

        
        //printf("x 3 %f\n",*(x + 2));
        //printf("xmin 3 %f\n",xmin);
        //printf("xmax 3 %f\n",xmax);

	/* Compute eta3 by normalising cv3 */
	eta3 = (*(x + 2) - xmin) / (xmax - xmin);
	if (fgm->gridpower[2] >= 0.0) {
		eta3 = pow(eta3, fgm->gridpower[2]);
	}

        //printf("eta3 %f\n",eta3);
	
        /* Determine mesh point*/
	i3m = (int)(eta3 * N3);
	/* Limit to [0,N-1]*/
	i3m = max(0, i3m);
	i3m = min(N3 - 1, i3m);
	i3p = i3m + 1;

	/* Determine weight factors*/
	w3p = eta3 * N3 - i3m;
	w3p = max(w3p, 0.0); /* No extrapolation */
        w3p = min(w3p, 1.0);
	w3m = 1 - w3p;
	
	/* printf("CV3: %e %e %e %e %d %d %e %e\n",*(x+2),xmin,xmax,eta3, i3p, i3m, w3p, w3m); */

	/* Compute extremes of cv2 (ivar = 1)*/
        /* i1 = 0, i2 = 0, i3 = i3m,i3p */
	xmin = w3m * fgm->data[i3m*N2p*N1p*Nvar + 1]
             + w3p * fgm->data[i3p*N2p*N1p*Nvar + 1];
        /* i1 = 0, i2 = N2, i3 = i3m,i3p */
	xmax = w3m * fgm->data[(i3m*N2p+N2)*N1p*Nvar + 1] 
             + w3p * fgm->data[(i3p*N2p+N2)*N1p*Nvar + 1];


        //printf("x 2 %f\n",*(x + 1));
        //printf("xmin 2 %f\n",xmin);
        //printf("xmax 2 %f\n",xmax);

	/* Compute eta2 from cv2 */
	eta2 = (*(x + 1) - xmin) / (xmax - xmin);
	if (fgm->gridpower[1] >= 0.0) {
		eta2 = pow(eta2, fgm->gridpower[1]);
	}

        //printf("eta2 %f\n",eta2);

	/* Determine mesh point*/
	i2m = (int)(eta2 * N2);
	i2m = max(0, i2m);
	i2m = min(N2 - 1, i2m);
	i2p = i2m + 1;

	/* Determine weight factors*/
	w2p = eta2 * N2 - i2m;
	w2p = max(w2p, 0.0); /* No extrapolation */
        w2p = min(w2p, 1.0);
	w2m = 1 - w2p;
	
	/* printf("CV2: %e %e %e %e %d %d %e %e\n",*(x+1),xmin,xmax,eta2, i2p, i2m, w2p, w2m); */

	/* Compute extremes of cv1 (ivar = 0) */
	/* i1 = 0, i2 = i2m,i2p, i3 = i3m,i3p */
	xmin = w3m * w2m * fgm->data[((i3m*N2p + i2m)*N1p)*Nvar]
             + w3p * w2m * fgm->data[((i3p*N2p + i2m)*N1p)*Nvar]
             + w3m * w2p * fgm->data[((i3m*N2p + i2p)*N1p)*Nvar]
             + w3p * w2p * fgm->data[((i3p*N2p + i2p)*N1p)*Nvar];
	/* i1 = N1, i2 = i2m,i2p, i3 = i3m,i3p */		 
	xmax = w3m * w2m * fgm->data[((i3m*N2p + i2m)*N1p + N1)*Nvar]
             + w3p * w2m * fgm->data[((i3p*N2p + i2m)*N1p + N1)*Nvar]
             + w3m * w2p * fgm->data[((i3m*N2p + i2p)*N1p + N1)*Nvar]
             + w3p * w2p * fgm->data[((i3p*N2p + i2p)*N1p + N1)*Nvar];	  

        //printf("x 1 %f\n",*x);
        //printf("xmin 1 %f\n",xmin);
        //printf("xmax 1 %f\n",xmax);
		 
	/* Compute eta1 from cv1 */
	eta1 = (*x - xmin) / (xmax - xmin);
	if (fgm->gridpower[0] >= 0.0) {
		eta1 = pow(eta1, fgm->gridpower[0]);
	}

        //printf("eta1 %f\n",eta1);

	/* Determine mesh point*/
	i1m = (int)(eta1 * N1);
	i1m = max(0, i1m);
	i1m = min(N1 - 1, i1m);
	i1p = i1m + 1;

	/* Determine weight factors*/
	w1p = eta1 * N1 - i1m;
	w1p = max(w1p, 0.0); /* No extrapolation */
        w1p = min(w1p, 1.0);
        w1m = 1 - w1p;
	
	/* printf("CV1: %e %e %e %e %d %d %e %e\n",*(x),xmin,xmax,eta1, i1p, i1m, w1p, w1m); */

	/* Do the actual interpolation in eta space */
        /* Sum contributions of all 8 corner points */

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = w1m * w2m * w3m * fgm->data[((i3m*N2p + i2m)*N1p + i1m)*Nvar + ivar];
	}

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = *(f + ivar) 
                            + w1p * w2m * w3m * fgm->data[((i3m*N2p + i2m)*N1p + i1p)*Nvar + ivar];
	}

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = *(f + ivar) 
                            + w1m * w2p * w3m * fgm->data[((i3m*N2p + i2p)*N1p + i1m)*Nvar + ivar];
	}

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = *(f + ivar) 
                            + w1p * w2p * w3m * fgm->data[((i3m*N2p + i2p)*N1p + i1p)*Nvar + ivar];
	}

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = *(f + ivar) 
                            + w1m * w2m * w3p * fgm->data[((i3p*N2p + i2m)*N1p + i1m)*Nvar + ivar];
	}

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = *(f + ivar) 
                            + w1p * w2m * w3p * fgm->data[((i3p*N2p + i2m)*N1p + i1p)*Nvar + ivar];
	}

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = *(f + ivar) 
                            + w1m * w2p * w3p * fgm->data[((i3p*N2p + i2p)*N1p + i1m)*Nvar + ivar];
	}

	for (ivar = 0; ivar < Nvar; ivar++) {
		*(f + ivar) = *(f + ivar) 
                            + w1p * w2p * w3p * fgm->data[((i3p*N2p + i2p)*N1p + i1p)*Nvar + ivar];
	}

	return EXIT_SUCCESS;
};

int lookupFGM_ND (FGM *fgm, double *cv, double *f)
{
    /* void free (void* f); */
    int dd, j, ivar;
    int nk[2]; /* iData: index for the FGM data */
    int *S = NULL, *jd = NULL; /* auxiliary vectors for look-up */
    int *jmin = NULL, *jmax = NULL; /* min & max indices */
    int dummy0=0;  /* dummies for array computations */
    int icv=0, index_cv=0;
    int jleft, jright, mid;
    double etaleft, etaright, etamid;
    double cvmin, cvmax, eta, wp;
    double cvhlp[2]; /*,cvminn[fgm->Ncv],cvmaxx[fgm->Ncv]; */
    /* [2] is the # of points required for the interpolation */

    /*  Gridpoint index and weight */
    int l0 = pow(2,fgm->Ncv); /* # of points */
    int *jj = NULL;
    double *wj = NULL;

    S = malloc(sizeof(int) * fgm->Ncv);      // int 3
    jd = malloc(sizeof(int) * fgm->Ncv+1);   // int 4
    jmin = malloc(sizeof(int) * fgm->Ncv);   // int 3
    jmax = malloc(sizeof(int) * fgm->Ncv);   // int 3
    jj = malloc(sizeof(int) * l0);           // int 8
    wj = malloc(sizeof(double) * l0);        // int 8

    int l1=0;
    S[0] = 1;
    for( dummy0=1; dummy0<fgm->Ncv; dummy0++ )
    {
        l1 = fgm->Ngrid[dummy0-1];
        //printf("l1 %i\n",l1); 
        S[dummy0] = S[dummy0-1]*l1;
    }
    //printf("S[0] %i\n",S[0]);
    //printf("S[1] %i\n",S[1]);
    //printf("S[2] %i\n",S[2]);

    jd[0] = 1;
    for( dummy0=0; dummy0<fgm->Ncv; dummy0++ )
    {
        jd[dummy0+1] = pow(2, dummy0+1);
    }
    //printf("jd[0] %i\n",jd[0]);
    //printf("jd[1] %i\n",jd[1]);
    //printf("jd[2] %i\n",jd[2]);
    //printf("jd[3] %i\n",jd[3]);


     /* Min and max indices for each dimension */
    for( dummy0=0; dummy0<fgm->Ncv; dummy0++ )
    {
        jmin[dummy0] = 1;
        jmax[dummy0] = fgm->Ngrid[dummy0]-1;
    }
    //printf("jmin[0] %i\n",jmin[0]);
    //printf("jmax[0] %i\n",jmax[0]);
    //printf("jmin[1] %i\n",jmin[1]);
    //printf("jmax[1] %i\n",jmax[1]);
    //printf("jmin[2] %i\n",jmin[2]);
    //printf("jmax[2] %i\n",jmax[2]);


    /*  Evaluate the weights */
    for ( dummy0=0; dummy0<l0; dummy0++ ) /* setting to initial values for each look-up */
    {
        jj[dummy0] = 1;
        wj[dummy0] = 1.0;
    }
    //printf("Initial jj %i\n",jj[0]);
    //printf("Initial wj %f\n",wj[0]);

    for ( dd=0; dd<fgm->Ncv; dd++ ) /* loop for each CV */
    {
        //printf("In loop with dd = %i\n",dd);
        /* Set the initial values to zero */
        cvmin=0.0, cvmax=0.0;
        eta = 0.0, wp = 0.0;
        cvhlp[0] = 0.0, cvhlp[1] = 0.0;

        for ( icv=0; icv<jd[dd]; icv++ )
        {   
            //printf("In innter loop where icv = %i\n",icv);
            //printf("wj[icv] = %f\n",wj[icv]);
            //printf("jj[icv] = %i\n",jj[icv]);
            //printf("(jj[icv] - 1)*fgm->Nvar + dd = %i\n",(jj[icv]-1)*fgm->Nvar+dd);
            //printf("fgm->data([jj[icv]-1) = %f\n",fgm->data[(jj[icv]-1)*fgm->Nvar+dd]);
            cvmin = cvmin + wj[icv]*fgm->data[(jj[icv]-1)*fgm->Nvar+dd];
            cvmax = cvmax + wj[icv]*fgm->data[(jj[icv]+jmax[dd]*S[dd]-1)*fgm->Nvar+dd];
            //printf("cvmin = %f\n",cvmin);
            //printf("cvmax = %f\n",cvmax);

        }
        eta = (*(cv+dd) - cvmin)/(cvmax - cvmin + 1e-8); /* scaling the CVs */
        //printf("eta = %f\n",eta);
        //eta = (*(cv+dd) - cvmin)/(cvmax - cvmin); /* scaling the CVs */
        eta = max(eta, 0.0);
        eta = min(eta, 1.0);
        //printf("eta clipped = %f\n",eta);

        /* Find grid point index in the current dimension */
        if (dd > -1) /* for now grid inversion is used for all variables */
        {
            /* Analitycal inversion of grid distribution */
            if ( fgm->gridpower[dd] >= 0.0 )
            {
                eta = pow(eta, fgm->gridpower[dd]);
            }

            j = floor(eta * jmax[dd] + 1);
            
            j = max(j, jmin[dd]);
            j = min(j, jmax[dd]);
        }
        else
        {
            /* Binary search */
            if (fabs(cvmin-cvmax)<1e-10)
            {
                j = jmin[dd];
            }
            else
            {
                jleft = 0;
                jright = jmax[dd];

                etaleft = 0.0;
                etaright = 1.0;
                while (jleft<jright-1)
                {
                    if ((jright - jleft) == 2)
                    {
                        mid = jleft +1;
                    }
                    else
                    {
                        mid = ceil( ((double) (jleft + jright)) / 2.0);
                    }

                    etamid = 0.0;
                    for ( icv=0; icv<jd[dd]; icv++ )
                    {
                        etamid = etamid + wj[icv]*fgm->data[(jj[icv]+mid*S[dd]-1)*fgm->Nvar+dd];
                    }
                    etamid = (etamid - cvmin ) / (cvmax - cvmin);

                    if (eta < etamid)
                    {
                        jright = mid;
                        etaright = etamid;
                    }
                    else
                    {
                        jleft = mid;
                        etaleft = etamid;
                    }
                }
                j = jright;
            }

        }



        /* ****************************
        Linear interpolation without extrapolation.
        Even if custom interpolation is used, the following needs to be
        computed for the cv values calculation in the subsequent dimensions.

        Find the indexes of the gridpoints in the table */
        for ( icv=0; icv<jd[dd]; icv++ )
        {
            jj[jd[dd]+icv] = jj[icv] + j*S[dd];

            jj[icv] = jj[icv] + (j-1)*S[dd];

        }

        /* Interpolate cv values at j and j+1 */
        for ( icv=0; icv<jd[dd]; icv++ ) /*for ( icv=dd; icv<fgm->Ncv; icv++ ) */
        {
            cvhlp[0] = cvhlp[0] + wj[icv]*fgm->data[(jj[icv]-1)*fgm->Nvar+dd];
            cvhlp[1] = cvhlp[1] + wj[icv]*fgm->data[(jj[jd[dd]+icv]-1)*fgm->Nvar+dd];

        }

        /* Compute weights for the search in the next dimension */
        wp = cvhlp[1] - cvhlp[0];
        
        /* these weights are for the purpose of search in the next
        dimension, so clipping them to the boundaries of this cell */
        wp = min(max((*(cv+dd) - cvhlp[0])/(wp + 1e-8), 0.0), 1.0);
        //wp = min(max((*(cv+dd) - cvhlp[0])/wp, 0.0), 1.0);

        for ( icv=0; icv<jd[dd]; icv++ )
        {
            wj[jd[dd]+icv] = wj[icv]*wp;
            wj[icv] = wj[icv]*(1.0-wp);
        }
    }

    /*  Actual interpolation */
    for ( ivar = 0; ivar<fgm->Nvar; ivar++ )
    {
        *(f + ivar) = 0.0;
        for ( dummy0=0; dummy0<l0; dummy0++ )
        {
            *(f + ivar) = *(f + ivar) + wj[dummy0]*fgm->data[(jj[dummy0]-1)*fgm->Nvar+ivar];
        }
    }

    free(S);
    free(jd);
    free(jmin);
    free(jmax);
    free(jj);
    free(wj);

    return EXIT_SUCCESS;
};
