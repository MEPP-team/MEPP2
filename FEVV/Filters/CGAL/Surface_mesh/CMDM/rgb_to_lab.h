// Copyright (c) 2012-2020 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#pragma once 

inline double F(double input) // function f(...), which is used for defining L, a and b
                       // changes within [4/29,1]
{
  if (input > 0.008856)
    return std::cbrt(input); // maximum 1 --- prefer cbrt to pow for cubic root
  else
    return ((double(841.0) / 108.0) * input +
            double(4.0) / 29.0); // 841/108 = 29*29/36*16
}

inline void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z)
{
    // RGB Working Space: sRGB 
    // R 0 -- 1, G 0 -- 1, B 0 -- 1

    //Gamma correction
	R = ((R > 0.0404482362771076) ? std::pow((R + 0.055) / 1.055, 2.4) : (R / 12.92)) ;
    G = ((G > 0.0404482362771076) ? std::pow((G + 0.055) / 1.055, 2.4) : (G / 12.92));
    B = ((B > 0.0404482362771076) ? std::pow((B + 0.055) / 1.055, 2.4) : (B / 12.92));

    //Transform
    /*X = 0.4124 *R + 0.3576*G + 0.1805*B; 
    Y = 0.2126 *R + 0.7152*G + 0.0722*B; 
    Z = 0.0193 *R + 0.1192*G + 0.9505*B;*/

    // Coefficients more precise 
   X = 0.412396 *R + 0.357583*G + 0.180493*B;
   Y = 0.212586*R + 0.71517*G + 0.0722005*B;
   Z = 0.0192972*R + 0.119184 *G + 0.950497*B;

}

// XYZ to CIELab
inline void XYZtoLab(double X, double Y, double Z, double *lab)
{

    //D65 white
    //const double Xo = 0.95047; 
    //const double Yo = 1.00000;
    //const double Zo = 1.08883;

	// matlab White point
    const double Xo = 0.950456;
    const double Yo = 1.000000;
    const double Zo = 1.088754;

    lab[0] = 116.0 * F(Y / Yo) - 16.0; // maximum L = 100
    lab[1] = 500.0 * (F(X / Xo) - F(Y / Yo)); // maximum 
    lab[2] = 200.0 * (F(Y / Yo) - F(Z / Zo));

}

// RGB to CIELab
inline void rgb_to_lab(double R, double G, double B, double *lab)
{
    double X, Y, Z;
    RGBtoXYZ(R, G, B, X, Y, Z);
    XYZtoLab(X, Y, Z, lab);
}
