// +-------------------------------------------------------------------------
// | color.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2012
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the TopTop library.
// |
// |    TopTop is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    TopTop is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with TopTop.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#include "vec.h"

Vec4f blendBins(float index, const Vec4f *bins, int NBins);
inline Vec4f rgb(uint r, uint g, uint b) {
    return Vec4f(r/255.0, g/255.0, b/255.0, 1.0);
}

// Ramps
//  please give values in range 0.0 to 1.0
Vec4f brewerBlueRamp(float value)
{
    static const Vec4f ramp[] = {
        rgb(3,78,123),
        rgb(5,112,176),
        rgb(54,144,192),
        rgb(116,169,207),
        rgb(166,189,219),
        rgb(208,209,230),
        rgb(241,238,246),
    };
    
    return blendBins(clamp(value, 0.0f, 1.0f), ramp,7);
}
Vec4f brewerRedRamp(float value)
{
    static const Vec4f ramp[] = {
        rgb(145,0,63),
        rgb(206,18,86),
        rgb(231,41,138),
        rgb(223,101,176),
        rgb(201,148,199),
        rgb(212,185,218),
        rgb(241,238,246),
    };
    
    return blendBins(clamp(value, 0.0f, 1.0f), ramp,7);
}
// Diverging Ramps
//  please give values in range -1.0 to 1.0
Vec4f brewerRedBlueDiverging(float value)
{
    static const Vec4f ramp[] = {
        rgb(215,25,28),
        rgb(253,174,97),
        rgb(255,255,191),
        rgb(171,217,233),
        rgb(44,123,182),
    };
    
    return blendBins(clamp((value + 1.0f)/2.0f, 0.0f, 1.0f), ramp, 5);
}




// index is in range 0.0 to 1.0
inline
Vec4f blendBins(float index, const Vec4f *bins, int NBins)
{
    float amplified = index * float(NBins-1);
    int bin = int(amplified);
    if(bin > (NBins-2)) bin = NBins-2; // push top value down...
    float fraction = amplified - float(bin);
    return (1.0f-fraction) * bins[bin] + fraction * bins[bin+1];
}





