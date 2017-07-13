// Copyright (C) Geom Software e.U, Bernhard Kornberger, Graz/Austria
//
// This file is part of the Fade2D library. The student license is free
// of charge and covers personal non-commercial research. Licensees
// holding a commercial license may use this file in accordance with
// the Commercial License Agreement.
//
// This software is provided AS IS with NO WARRANTY OF ANY KIND,
// INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE.
//
// Please contact the author if any conditions of this licensing are
// not clear to you.
//
// Author: Bernhard Kornberger, bkorn (at) geom.at
// http://www.geom.at



#pragma once

#include "common.h"
#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif

/** \brief Colorname
*
* A few colors are predefined for convenience
*/
enum Colorname
{
	CRED,
	CGREEN,
	CBLUE,
	CBLACK,
	CYELLOW,
	CPINK,
	CGRAY
};

/**  \brief %Color
*
* \see Visualizer2 where Color objects are used for visualizations
*
*/
class CLASS_DECLSPEC Color
{
public:
/**
*
* Constructor for a color to be used together with the Visualizer2 class
* @param r_ red
* @param g_ green
* @param b_ blue
* @param width_ linewidth
* @param bFill_ fill-color (default value: \e false)
*
*/

	Color(double r_,double g_,double b_,double width_,bool bFill_=false);

/**
 *
 * Constructor for a color to be used together with the Visualizer2 class
 *
 * @param c is a predefined color name
 * @param width_ linewidth
 * @param bFill_ fill-color (default value: \e false)
 */
	Color(Colorname c,float width_=0.001,bool bFill_=false);
	Color();


	bool operator!=(Color& other) const;
	bool operator==(Color& other) const;
	float r,g,b,width;
	bool bFill;
	friend std::ostream &operator<<(std::ostream &stream, const Color& c);
};



inline std::ostream &operator<<(std::ostream &stream, const Color& c)
{
	stream<<"Color (r,g,b): "<<c.r<<","<<c.g<<","<<c.b<<", linewidth="<<c.width<<", fill="<<c.bFill;
	return stream;
}


} // (namespace)
