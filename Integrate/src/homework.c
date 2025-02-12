#include <stdio.h>
#include <math.h>
#include "glfem.h"

double jacobien(double x[3], double y[3])
{
    return (x[1] - x[0])*(y[2] - y[0]) - (x[2] - x[0])*(y[1] - y[0]);
}

double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    //Xloc et Yloc points utilisés pour l'intégration

    xLoc[0] = x[0] * 0.66667 + (x[1] + x[2]) * 0.16667; //  = X0 * 2/3 + X1 * 1/6 + X2 * 1/6 
    xLoc[1] = (x[0] + x[1]) * 0.16667 + x[2] * 0.66667; //  = X0 * 1/6 + X1 * 1/6 + X2 * 2/3
    xLoc[2] = (x[0] + x[2]) * 0.16667 + x[1] * 0.66667; //  = X0 * 1/6 + X1 * 2/3 + X2 * 1/6
    yLoc[0] = y[0] * 0.66667 + (y[1] + y[2]) * 0.16667; //  = Y0 * 2/3 + Y1 * 1/6 + Y2 * 1/6
    yLoc[1] = (y[0] + y[1]) * 0.16667 + y[2] * 0.66667; //  = Y0 * 1/6 + Y1 * 1/6 + Y2 * 2/3
    yLoc[2] = (y[0] + y[2]) * 0.16667 + y[1] * 0.66667; //  = Y0 * 1/6 + Y1 * 2/3 + Y2 * 1/6

    for(unsigned short int i = 0; i < 3; ++i)
        I += f(xLoc[i],yLoc[i]);

    I *= 0.16667 * jacobien(x,y); //poids tous égaux à 1/6

    
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    

    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

    double I = 0.0;

    //Si pas de sous-triangles, on calcule directement l'intégrale
    //Sinon, on calcules les 4 sous-triangles et on recommence sur chacun d'entre eux

    if(n == 0)
        I = integrate(x,y,f);
    else
    {
        double x2[3];
        double y2[3];
        x2[0] = (x[0] + x[1]) / 2;
        x2[1] = (x[1] + x[2]) / 2;
        x2[2] = (x[2] + x[0]) / 2;

        y2[0] = (y[0] + y[1]) / 2;
        y2[1] = (y[1] + y[2]) / 2;
        y2[2] = (y[2] + y[0]) / 2;

        //points des milieux des segments


        double trx[3] = {x[0],x2[0],x2[2]}; 
        double try[3] = {y[0],y2[0],y2[2]}; 
        I += integrateRecursive(tr1x,tr1y,f,n - 1);
        trx[3] = {x2[0],x[1],x2[1]}; 
        try[3] = {y2[0],y[1],y2[1]}; 
        I += integrateRecursive(tr2x,tr2y,f,n - 1);
        trx[3] = {x2[2],x2[1],x[2]}; 
        try[3] = {y2[2],y2[1],y[2]}; 
        I += integrateRecursive(tr3x,tr3y,f,n - 1);
        trx[3] = {x2[0],x2[1],x2[2]}; 
        try[3] = {y2[0],y2[1],y2[2]}; 
        I += integrateRecursive(tr4x,tr4y,f,n - 1);

        //somme des intégrales sur les 4 sous-triangles
    }       
     
    return I;
}
