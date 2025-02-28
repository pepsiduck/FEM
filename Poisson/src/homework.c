#include "fem.h"
#include "stdio.h"

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
    for(int i = 0; i < nBoundary * 2; ++i)
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
    
    for(int i = 0; i < theProblem->space->n; ++i)
    {
        map[i] = theMesh->elem[theProblem->space->n*iElem + i];
        x[i] = theMesh->nodes->X[map[i]];
        y[i] = theMesh->nodes->Y[map[i]];
    }

    //  A completer :-)

}

# endif
# ifndef NOPOISSONSOLVE

double absolute_value(double arg)
{
    if(arg > 0)
        return arg;
    return -1 * arg;
}

double jacobian(double *x, double *y)
{
    return absolute_value((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]));
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


    for(int i = 0; i < theSystem->size; ++i)
    {
        for(int j = 0; j < theSystem->size; ++j)
            theSystem->A[i][j] = 0;
        theSystem->B[i] = 0;
    }

    double Ae[4][4] = {{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0}};
    double Be[4] = {0.0,0.0,0.0,0.0};
    int map[4] = {-1,-1,-1,-1};
    double x[4] = {0.0,0.0,0.0,0.0};
    double y[4] = {0.0,0.0,0.0,0.0};

    double phi[4] = {0.0,0.0,0.0,0.0};
    double dphideta[4] = {0.0,0.0,0.0,0.0};
    double dphidxsi[4] = {0.0,0.0,0.0,0.0};

    double Je = 0.0;


    for(int n = 0; n < theMesh->nElem; ++n)
    {
        femPoissonLocal(theProblem,n,map,x,y);
        Je = jacobian(x,y);

        //A
        for(int i = 0; i < theSpace->n; ++i)
        {
            for(int j = 0; j < theSpace->n; ++j)
            {
                Ae[i][j] = 0;
                for(int k = 0; k < theRule->n; ++k)
                {
                    theSpace->dphi2dx(theRule->xsi[k],theRule->eta[k],dphidxsi,dphideta);
                    Ae[i][j] += theRule->weight[k] * ((dphidxsi[i]*(y[2]-y[0]) + dphideta[i]*(y[0]-y[1]))*(dphidxsi[j]*(y[2]-y[0]) + dphideta[j]*(y[0]-y[1])) + (dphidxsi[i]*(x[0]-x[2]) + dphideta[i]*(x[1]-x[0]))*(dphidxsi[j]*(x[0]-x[2]) + dphideta[j]*(x[1]-x[0])));
                }
                Ae[i][j] /= Je;
            }
        }
            
        //B
        for(int i = 0; i < theSpace->n; ++i)
        {
            Be[i] = 0;
            for(int k = 0; k < theRule->n; ++k)
            {
                theSpace->phi2(theRule->xsi[k],theRule->eta[k],phi);
                Be[i] += theRule->weight[k] * phi[i];
            }
            Be[i] *= Je; //f * jacobian
        }

        //Assemblage
        for(int i = 0; i < theSpace->n; ++i)
        {
            for(int j = 0; j < theSpace->n; ++j)
                theSystem->A[map[i]][map[j]] += Ae[i][j];
        }
        for(int i = 0; i < theSpace->n; ++i)
            theSystem->B[map[i]] += Be[i];
            
    }

    //A et B sont complets

    //Contrainte
    for(int u = 0; u < theBoundary->nElem; ++u)
        femFullSystemConstrain(theSystem,theBoundary->elem[u],0.0); 

    /*
    for(int i = 0; i < theSystem->size; ++i)
    {
            printf("%f\n",theSystem->B[i]);
    }*/
            
    //Résolution
    femFullSystemEliminate(theSystem);
    
    // A completer :-)
}

# endif



