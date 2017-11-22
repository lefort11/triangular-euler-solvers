
// (c) 2010 Geom e.U. Bernhard Kornberger, Graz/Austria. All rights reserved.
//
// This file is part of the Fade2D library. The licensee mentioned below may use
// this license file for evaluation purposes during the agreed period
// Commercial use requires a valid commercial license, this applies also to
// inhouse usage. This license file is personalized, DO NOT SHARE IT.

// This software is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING 
// THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Please contact the author if any conditions of this licensing are not clear 
// to you.
// 
// Author: Bernhard Kornberger, bkorn (at) geom.at 
//         C++ Freelancer
// http://www.geom.at/products/fade2d/


#pragma once 
#include "Fade_2D.h"

namespace{
#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_LIC=GEOM_FADE25D;
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_LIC=GEOM_FADE2D;
#endif

struct License
{
License()
{
	GEOM_LIC::setLic(
		"Lomonosov Moscow State University;Faculty of Computational Mathematics and Cybernetics;Russia;119991, Moscow, GSP-1, 1-52;Leninskiye Gory;;Student: Gleb Plaksin;Supervisor: Mikhail Abakumov;Research license, valid until 11/2019;;Not valid for use in commercial software;Do not share this License File;",
		"[LicType,eval],[2D,1e6],[25D,5e5],[MeshGen,5e5],[SegCheck,5e5],[CutFill,1e4]",
		"[LF:F/C]",
		"940ebf7d",
		"7637e127");
	}
};
License lic;
}
