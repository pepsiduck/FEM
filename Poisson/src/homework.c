#include "fem.h"

# ifndef NOPOISSONCREATE

int isInArray(int *array, int size, int arg)
{
    for(int i = 0; i < size; ++i)
    {
        if(array[i] == arg)
            return 1;
    }
    return 0;
}

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges; 
    int nBoundary = theEdges->nElem;
    
    //  A completer :-)

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    for(int i = 0; i < nBoundary; ++i)
        theBoundary->elem[i] = -1;
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");
 
    int j = 0;
    for(int i = 0; i < nBoudary * 2; ++i)
    {
        if(!isInArray(theBoundary->elem,nBoundary,theEdges->elem[i]))
        {
            theBoundary->elem[j] = theEdges->elem[i];
            ++j;
        }
    }

    // A completer :-)

}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    geoMeshFree(theProblem->geo);
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->rule);
    femFullSystemFree(theProblem->system);
    free(theProblem);
    // A completer :-)
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    
    for(int i = 0; i < theProblem->space.n; ++i)
    {
        map[i] = theMesh->elem[theProblem->space.n*iElem + i];
        x[i] = theMesh->nodes->X[map[i]];
        y[i] = theMesh->nodes->Y[map[i]];
    }

    //  A completer :-)

}

# endif
# ifndef NOPOISSONSOLVE

double jacobian(double *x, double *y)
{
    return abs((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]));
}

void matrixSolve(double **A, double *B, int size)
{ 
    int i, j, k;
    /* Gauss elimination */
    for (k=0; k < size; k++)
    {
        if (A[k][k] == 0) 
            Error("zero pivot");
        for (i = k+1 ; i < size; i++) 
        {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++)
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor;
        }
    }
    /* Back-substitution */
    for (i = (size)-1; i >= 0 ; i--) 
    {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i];
    }
}

void femPoissonSolve(femPoissonProblem *theProblem)
{
    //Mesh contenant tout les points du problème + une liste reprenant les points contituants les triangle [tr11, tr12, tr13, tr21, ...]
    femMesh *theMesh = theProblem->geo->theElements;

    //Ensemble contenant tous les neouds à la frontière
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");

    //Matrice A et vecteur B
    femFullSystem *theSystem = theProblem->system;

    //règle d'intégration
    femIntegration *theRule = theProblem->rule;

    //3 pour triangles, 4 pour quads
    femDiscrete *theSpace = theProblem->space;
 
    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
            
    //triangles
    if (theSpace->n == 3)
    {
        for(int i = 0; i < theSystem->size; ++i)
        {
            for(int j = 0; j < theSystem->size; ++j)
                theSystem->A[i][j] = 0;
            theSystem->B[i][j] = 0;
        }

        double **Ae = (double *) malloc(3*sizeof(int*));
        for(int i = 0; i < 3; ++i)
            Ae[i] = (double *) malloc(3*sizeof(int*));
        double *Be = malloc(3*sizeof(int));
        int *map = (int *) malloc(3*sizeof(int));
        double *x = (double *) malloc(3*sizeof(double));
        double *y = (double *) malloc(3*sizeof(double));

        double Je = 0.0;

        for(int n = 0; n < theMesh->nElem; ++n)
        {
            femPoissonLocal(theProblem,n,map,x,y);
            
            Je = jacobian(x,y);
        }
    
        for(int i = 0; i < 3; ++i)
            free(Ae[i]);
        free(Ae);
        free(Be);
        free(map);
        free(x);
        free(y);
    }

    // A completer :-)
}

# endif



