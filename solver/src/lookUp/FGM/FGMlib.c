/*---------------------------------------------------------------------------*\

███████╗ ██████╗ ███╗   ███╗███████╗ ██████╗  █████╗ ███╗   ███╗
██╔════╝██╔════╝ ████╗ ████║██╔════╝██╔═══██╗██╔══██╗████╗ ████║
█████╗  ██║  ███╗██╔████╔██║█████╗  ██║   ██║███████║██╔████╔██║
██╔══╝  ██║   ██║██║╚██╔╝██║██╔══╝  ██║   ██║██╔══██║██║╚██╔╝██║
██║     ╚██████╔╝██║ ╚═╝ ██║██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║
╚═╝      ╚═════╝ ╚═╝     ╚═╝╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝

FGMFoam - Combustion solver using Flamelet Generated Manifolds 
Developed by Stijn Schepers - Eindhoven University of Technology

FGMFoam is developed for research and experimental purposes.
Use, modification, and redistribution are subject to the
applicable OpenFOAM license terms.

\*---------------------------------------------------------------------------*/

# include "FGMlib.h"

FGM * readFGM (const char filename[])
{
    char line[MAX_LINE_LENGTH];
    FILE *fid;
  
    int Ncv    = 0;
    int *Ngrid = NULL;
    int Nvar   = 0;

    char (*varname)[VAR_NAME_LENGTH];
    float *data = NULL;

    int pressure = 101325; // Default pressure in Pa, can be read from file if needed

    FGM *fgm = NULL;
  
    int i,j;
    float f;

    printf("Reading FGM database file %s\n",filename);

    /* Open the file */
    fid = fopen(filename, "r");
    if (fid == NULL) {
        fprintf(stderr, "Error: Failed to open FGM database file\n");
        exit(EXIT_FAILURE);
    }
  
    /* Check the identifier [FGM] */
    if (!fgets(line, MAX_LINE_LENGTH, fid)) {
        fprintf(stderr, "Error reading line from file\n");
        exit(EXIT_FAILURE);
    }
    if (strncmp(line,"[FGM]",5) != 0) {
        fprintf(stderr, "Error: Incorrect FGM database format\n");
        exit(EXIT_FAILURE);
    }
  
    if (!locateKeyWord(fid, "[PRESSURE]"))  // found
    {
        if (fscanf(fid, "%d", &pressure) != 1) {
            fprintf(stderr, "Error: Failed to read pressure\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        fprintf(stderr, "Warning: Couldn't find keyword PRESSURE, using default ambient value\n");
    }

    /* Read the number of cv's */
    if (locateKeyWord(fid,"[DIMENSION]")) {
        fprintf(stderr, "Error: Couldn't find keyword DIMENSION\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(fid,"%i", &Ncv) != 1) {
        fprintf(stderr, "Error: Failed to read Ncv\n");
        exit(EXIT_FAILURE);
    }
    
    // fscanf(fid,"%i",&Ncv);
    // printf("Ncv = %i\n",Ncv);
  
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
        if (fscanf(fid, "%i", &j) != 1) {
            fprintf(stderr, "Error: Failed to read Ngrid[%d] from FGM file\n", i);
            exit(EXIT_FAILURE);
        }
        Ngrid[i] = j;
        Ntotal *= j;
    }
  
    /* Read the number of variables */
    if (locateKeyWord(fid,"[VARIABLES]")) {
        fprintf(stderr, "Error: Couldn't find keyword VARIABLES\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(fid,"%i",&Nvar) != 1) {
        fprintf(stderr, "Error: Failed to read Nvar\n");
        exit(EXIT_FAILURE);
    }

    /* Allocate space for the variable names */
    varname = malloc(sizeof *varname * Nvar);
    if (varname == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for variable names\n");
        exit(EXIT_FAILURE);
    }
  
    /* Read the variable names */
    for (j = 0; j < Nvar; j++) {
        if (fscanf(fid,"%s",varname[j]) != 1) {
            fprintf(stderr, "Error: Failed to read varname[%d] from FGM file\n", j);
            exit(EXIT_FAILURE);
        }
        // fscanf(fid,"%s",varname[j]);
        // printf("%s\n",varname[j]);
    }
  
    /* Allocate space for the data */
    data = malloc(sizeof(float) * Ntotal * Nvar);
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
            // fscanf(fid,"%f",&f);
            if (fscanf(fid,"%f",&f) != 1) {
                fprintf(stderr, "Error: Failed to read data[%d][%d] from FGM file\n", i, j);
                exit(EXIT_FAILURE);
            }
            data[i*Nvar + j] = f;
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
    fgm->pressure = pressure;
    fgm->Ncv = Ncv;
    fgm->Ngrid = Ngrid;
    fgm->Ntotal = Ntotal;
    fgm->Nvar = Nvar;
    fgm->varname = varname;
    fgm->data = data;

    printf("FGM loaded successfully.\n");
  
    return fgm;
};

int freeFGM(FGM *fgm)
{
    if (fgm != NULL) {
        if (fgm->Ngrid != NULL) { free(fgm->Ngrid); }
        if (fgm->data != NULL) { free(fgm->data); }
        if (fgm->varname != NULL) { free(fgm->varname); }
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

int lookupFGM_2D (FGM *fgm, float *x, float *f)
{
    int N1  = fgm->Ngrid[0]-1;
    int N1p = N1 + 1;
    int N2  = fgm->Ngrid[1]-1;
    int ivar, i1m, i1p, i2m, i2p;
    float eta1, eta2;
    float w1m, w1p, w2m, w2p;
    float min_cv1, max_cv1;
    float min_cv2, max_cv2;

    // CV2 bounds are constant over CV1 (ivar=1, i1=0, at i2=0 and i2=N2)
    min_cv2 = fgm->data[1];
    max_cv2 = fgm->data[(N2*N1p) * fgm->Nvar + 1];

    // CV2 is the second input
    eta2 = (*(x+1) - min_cv2) / (max_cv2 - min_cv2);

    i2m = (int) (eta2 * N2);
    i2m = fgm_max(0, i2m);
    i2m = fgm_min(N2-1, i2m);
    i2p = i2m + 1;

    w2p = eta2 * N2 - i2m;
    w2m = 1 - w2p;

    // CV1 bounds depend on CV2 row (ivar=0, i1=0 for min, i1=N1 for max)
    min_cv1 = w2m * fgm->data[(i2m*N1p) * fgm->Nvar] + w2p * fgm->data[(i2p*N1p) * fgm->Nvar];
    max_cv1 = w2m * fgm->data[(i2m*N1p + N1) * fgm->Nvar] + w2p * fgm->data[(i2p*N1p + N1) * fgm->Nvar];

    // CV1 is the first input
    eta1 = (*x - min_cv1) / (max_cv1 - min_cv1);

    i1m = (int) (eta1 * N1);
    i1m = fgm_max(0, i1m);
    i1m = fgm_min(N1-1, i1m);
    i1p = i1m + 1;

    w1p = eta1 * N1 - i1m;
    w1m = 1 - w1p;

    for (ivar = 0; ivar < fgm->Nvar; ivar++) {
        f[ivar] =
            w1m * w2m * fgm->data[(i1m + i2m*N1p) * fgm->Nvar + ivar] +
            w1p * w2m * fgm->data[(i1p + i2m*N1p) * fgm->Nvar + ivar] +
            w1m * w2p * fgm->data[(i1m + i2p*N1p) * fgm->Nvar + ivar] +
            w1p * w2p * fgm->data[(i1p + i2p*N1p) * fgm->Nvar + ivar];
    }

    return EXIT_SUCCESS;
}