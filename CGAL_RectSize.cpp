// Copyright (c) 2014  
// University of Illinois at Urbana-Champaign, All rights reserved. 
//
// This file is a CGAL version to calculate the smallest enclosing rectangle 
// using CGAL library (https://www.cgal.org/). The function prototype is:
//     [length,width,angle]=CGAL_RectSize(c);     
//     Input:   
//       c - the binary matrix of the single particle image
//     Output:  
//       length - the lenght of the rectangle
//       width - the width of the rectangle
//       angle - the angle of the rectangle in radian
//
// Author(s)     : Wei Wu   weiwu3@illinois.edu
//
// Compile command: mex CGAL_RectSize.cpp -I/opt/local/include -L/opt/local/lib -lCGAL -lgmp
// 
// Code History:
//   * First version created by Wei Wu, July 24, 2014
//
#include "mex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/convex_hull_2.h>
typedef CGAL::Simple_cartesian<double>            Kernel;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef Kernel::Point_2                           Point_2;
typedef Polygon_2::Vertex_iterator VertexIterator;
typedef Polygon_2::Edge_const_iterator EdgeIterator;

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{   
      
    #define B_OUT1 plhs[0]
    #define B_OUT2 plhs[1]
    #define B_OUT3 plhs[2]

    #define A_IN prhs[0]

    double *B0,*B1,*B2;
    int m,n;
    mxChar *ICEPIC;
    

    //B_OUT=mxCreateDoubleMatrix(1,3,mxREAL);
    B_OUT1=mxCreateDoubleScalar(0.0);
    B_OUT2=mxCreateDoubleScalar(0.0);
    B_OUT3=mxCreateDoubleScalar(0.0);

    B0=mxGetPr(B_OUT1); 
    B1=mxGetPr(B_OUT2); 
    B2=mxGetPr(B_OUT3); 
        
    m=mxGetM(A_IN);
    n=mxGetN(A_IN);    
    ICEPIC=mxGetChars(A_IN);

    int i,j,k=0,mm=0,nn=0;
    double L,W,anglex1,angley1,anglex2,angley2;
    std::vector<Point_2> points;
    for( i = 0; i < n; i++ )
    {
        for (j=0; j<m; j++)
        {
            //cout<< ICEPIC[i*m+j] << ", ";
            if(ICEPIC[i*m+j]==48) //48
            {
                points.push_back(Point_2(i,j));
            }
        }
        //cout<<"\n"<<endl;
        
    }
    //cout<< k << endl;

    if (points.size()==0)
    {
        //cout << "No Illuminated Doide!\n" << endl;
        *B0 = 0;
        *B1 = 0;
        *B2 = 0;

        return;
    }
          

    Polygon_2 ch;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch));

    Polygon_2 p_m1;
    CGAL::min_rectangle_2(ch.vertices_begin(), ch.vertices_end(),std::back_inserter(p_m1));
    //std::cout << p_m1 << std::endl;
    
    // traverse the vertices and the edges
    //     n=0;
    //     for (VertexIterator vi = p_m1.vertices_begin(); vi != p_m1.vertices_end(); ++vi)
    //       std::cout << "vertex " << n++ << " = " << *vi << std::endl;
    //     std::cout << std::endl;
    n=0;
    for (EdgeIterator ei = p_m1.edges_begin(); ei != p_m1.edges_end(); ++ei)
    {
      //std::cout << "edge " << n << " = " << *ei << ", ";
      //std::cout << "edge length " << n << " = " << sqrt(ei->squared_length())  << ", " ;
      //std::cout << "edge orientation " << n << " = " << ei->direction()  << ", ";
      //std::cout << "edge orientation x " << n << " = " << ei->direction().dx() << ", ";
      //std::cout << "edge orientation y " << n << " = " << ei->direction().dy()  << std::endl;

      if (0==n)
      {
          L=1.+sqrt(ei->squared_length()); 
          anglex1=ei->direction().dx();
          angley1=ei->direction().dy();
          //std::cout << "L: " << L << std::endl;

      }
      else if (1==n)
      {
          W=1.+sqrt(ei->squared_length()); 
          anglex2=ei->direction().dx();
          angley2=ei->direction().dy();
          //std::cout << "W: " << W << std::endl;
      }
      
      n++;
    }

    *B0=L>W? L:W;
    *B1=L>W? W:L;
    *B2=L>W? (0==anglex1? 3.1415926/2:atan(angley1/anglex1)):(0==anglex2? 3.1415926/2:atan(angley2/anglex2));
    
    return;

   
}
