// Copyright (c) 2014  
// University of Illinois at Urbana-Champaign, All rights reserved. 
//
// This file is a CGAL version to calculate the smallest enclosing circle 
// using CGAL library (https://www.cgal.org/). The function prototype is:
//     diameter=CGAL_minR(c);     
//     Input:   
//       c - the binary matrix of the single particle image
//     Output:  
//       diameter - the diameter of circle
//
// Author(s)     : Wei Wu   weiwu3@illinois.edu
//
// Compile command: mex CGAL_minR.cpp
// 
// Code History:
//   * First version created by Wei Wu, July 24, 2014
//   * Added support to take input of tas ratio and output circle center, Wei Wu, June 3, 2019

#include "mex.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include "Miniball.hpp"
#include <math.h>
#include "matrix.h"

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{   
    typedef double mytype;            // coordinate type
  
    #define B_OUT plhs[0]
    #define A_IN prhs[0]
    #define A2_IN prhs[1] // Added support for TAS ratio between probe & aircraft ~ Joe Finlon 5/12/19

    //double *B[3];
    double *B;
    double *tasRatio; // Added support for TAS ratio ~ Joe Finlon 5/12/19

    double *B1;
    double *B2;
    double *B3;

    int m,n;
    mxChar *ICEPIC;
    
    B_OUT=mxCreateDoubleMatrix(1,3,mxREAL);
    //B_OUT=mxCreateDoubleScalar(0.0);
    B=mxGetPr(B_OUT);
    //B2=mxGetPr(B_OUT[1]);
    //B3=mxGetPr(B_OUT[2]);
    
    m=mxGetM(A_IN);
    n=mxGetN(A_IN);    
    tasRatio=mxGetPr(A2_IN); // Added support for TAS ratio ~ Joe Finlon 5/12/19  
    ICEPIC=mxGetChars(A_IN);
    
    if ( 0==m*n ) //"No Illuminated Doide, return 0
    {
        cout << "No Illuminated Doide!\n" << endl;
        *(++B) = 0.0;
        *(++B) = 0.0;
        *B = 0.0;
        return;
    }

    int i,j,k=0;
    std::list<std::vector<mytype> > lp;
   
    for( i = 0; i < n; i++ )
    {
        for (j=0; j<m; j++)
        {
            //cout<< ICEPIC[i*m+j] << ", ";
            if(48==ICEPIC[i*m+j])
            {
                std::vector<mytype> p(2);
                p[0]=(double)j/ *tasRatio;
                p[1]=(double)i ; // Added support for TAS ratio ~ Joe Finlon 5/12/19
                lp.push_back(p);
            }
        }        
    }

    //cout<< k << endl;
    
      // define the types of iterators through the points and their coordinates
    // ----------------------------------------------------------------------
    typedef std::list<std::vector<mytype> >::const_iterator PointIterator; 
    typedef std::vector<mytype>::const_iterator CoordIterator;

    // create an instance of Miniball
    // ------------------------------
    typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > 
    MB;
    MB mb (2, lp.begin(), lp.end());

    double radius;
    radius=sqrt(mb.squared_radius());

    //float radius;
    //minEnclosingCircle(Mat(points), center, radius);
    
    //cout << "Center: " << center << " Radius: " << radius <<endl;

    
    // center
    const double* center = mb.center(); 
    *B = 1.+radius*2; //Assume round doide, add 1 for the boundary
    B++;
    *B = *center;
    center++;
    B++;
    *B = *center;
            
    return;
   
}