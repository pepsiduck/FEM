#include <stdio.h>
#include <math.h>
#include "glfem.h"

double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    //Xloc et Yloc points utilisés pour l'intégration

    xLoc[0] = x[0] * 0.66667 + (x[1] + x[2]) * 0.16667;
    xLoc[1] = (x[0] + x[1]) * 0.16667 + x[2] * 0.66667;
    xLoc[2] = (x[0] + x[2]) * 0.16667 + x[1] * 0.66667;
    yLoc[0] = y[0] * 0.66667 + (y[1] + y[2]) * 0.16667;
    yLoc[1] = (y[0] + y[1]) * 0.16667 + y[2] * 0.66667;
    yLoc[2] = (y[0] + y[2]) * 0.16667 + y[1] * 0.66667;

    for(unsigned short int i = 0; i < 3; ++i)
        I += f(xLoc[i],yLoc[i]);
    I *= 1.6667;
    
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    

    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

    double I = 0.0;

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


        double tr1x[3] = {x[0],x2[0],x2[2]}; 
        double tr1y[3] = {y[0],y2[0],y2[2]}; 
        I += integrateRecursive(tr1x,tr1y,f,n - 1);
        double tr2x[3] = {x2[0],x[1],x2[1]}; 
        double tr2y[3] = {y2[0],y[1],y2[1]}; 
        I += integrateRecursive(tr2x,tr2y,f,n - 1);
        double tr3x[3] = {x2[2],x2[1],x[2]}; 
        double tr3y[3] = {y2[2],y2[1],y[2]}; 
        I += integrateRecursive(tr3x,tr3y,f,n - 1);
        double tr4x[3] = {x2[0],x2[1],x2[2]}; 
        double tr4y[3] = {y2[0],y2[1],y2[2]}; 
        I += integrateRecursive(tr4x,tr4y,f,n - 1);
    }       
     
    return I;
}
