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
#include "Segment2.h"
#include "ConstraintSegment2.h"
#include "Edge2.h"
#include <map>


#include "common.h"
#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif

// FWD
class Dt2;
class ConstraintSegment2;
class GeomTest;
class Visualizer2;
class Color;

/** \brief ConstraintGraph2 is a set of Constraint Edges to be enforced. These are ConstraintSegment2 objects.
*
* \see \ref createConstraint in the Fade2D class
*
* \image html crop_ignoreBike.jpg "Constraints in a triangulation"
* \image latex crop_ignoreBike.eps "Constraints in a triangulation" width=12cm
*/
class ConstraintGraph2
{
public:
/// @private
	CLASS_DECLSPEC
	ConstraintGraph2(
		Dt2* pDt2_,
		std::vector<ConstraintSegment2*>& vCSegments,
		ConstraintInsertionStrategy
		);

/// @private
	CLASS_DECLSPEC
	ConstraintGraph2(
		Dt2* pDt2_,
		std::vector<ConstraintSegment2*>& vCSegments_,
		std::map<std::pair<Point2*,Point2*>,bool > mPPReverse,
		ConstraintInsertionStrategy cis_
	);

/// @private
	void init(std::vector<ConstraintSegment2*>& vCSegments_);


/** \brief Does the constraint graph form a closed polygon?
*
* \return true when the present ConstraintGraph forms a closed polygon.
*
* @note This method won't check if it is a simple polygon (one without
* self-intersections).
*/
	CLASS_DECLSPEC
	bool isPolygon() const;

/** \brief Are the segments of the constraint graph oriented?
*
* @return true if the constraint graph has been created with
* bOrientedSegments=true or if automatic reorientation was possible
* which is the case for simple polygons.
*/
	CLASS_DECLSPEC
	bool isOriented() const;

/** \brief Get the vertices of the constraint segments
*
* Use this method to retrieve the segments of the present object
* in form of a vector of vertices. If *this is a closed polygon, then
* the points are ordered and oriented in counterclockwise direction,
* e.g. (a,b,b,c,c,d,d,a). If the *this is not a polygon, then the
* segments are returned in the original direction.
*
* @note If it was necessary to split the constraint segments, then the
* splitted segments are returned. If, in the above example, the constraint
* segment (a,b) crosses some previously inserted point x, then the
* result is (a,x,x,b,b,c,c,d,d,a).
*
*/
	CLASS_DECLSPEC
	void getPolygonVertices(std::vector<Point2*>& vTriangulationPoints_) ;


/** \brief Get the constraint insertion strategy
*
* \if SECTION_FADE25D
* @return CIS_CONFORMING_DELAUNAY, CIS_CONFORMING_DELAUNAY_SEGMENT_LEVEL or @n
* CIS_CONSTRAINED_DELAUNAY
* \else
* @return CIS_CONFORMING_DELAUNAY or CIS_CONSTRAINED_DELAUNAY
* \endif
*
*/
	CLASS_DECLSPEC
	ConstraintInsertionStrategy getInsertionStrategy() const;

/** \brief Check if an edge is a constraint
*
* Checks if the edge (p0,p1) is a constraint of the present
* ConstraintGraph2 object.
*/
	CLASS_DECLSPEC
	bool isConstraint(Point2* p0,Point2* p1) const;

/** \brief Check if a ConstraintSegment2 is a member
*
* The present ConstraintGraph2 has been created using a set of edges
* and this method checks if the ConstraintSegment2 \p pCSeg is one of
* them. Original edges that have been splitted are not alive anymore
* and are no members. But their child segments are members.
*
*/
	CLASS_DECLSPEC
	bool isConstraint(ConstraintSegment2* pCSeg) const;


/** \brief Visualization
*
*/
	CLASS_DECLSPEC
	void show(const std::string& name);

/** \brief Visualization
*
*/
	CLASS_DECLSPEC
	void show(Visualizer2* pVis,const Color& color);


/** \brief Get the original ConstraintSegment2 objects
*
* Get the original, ConstraintSegment2 objects. They are not subdivided
* but may be dead and have child segments (which may also be dead and
* have child segments...)
*
*/
	CLASS_DECLSPEC
	void getOriginalConstraintSegments(std::vector<ConstraintSegment2*>& vConstraintSegments_) const;

/** \brief Get child ConstraintSegment2 objects
*
* Returns the current constraint segments, i.e., the original ones or,
* if splitted, their child segments.
*
*/
	CLASS_DECLSPEC
	void getChildConstraintSegments(std::vector<ConstraintSegment2*>& vConstraintSegments_) const;

/// @private
bool dbg_hasDirection(ConstraintSegment2* pCSeg) const;


///// @private
//void updateSplittedConstraintSegment(
		//ConstraintSegment2* pCSeg,
		//bool bDirChange0,
		//bool bDirChange1,
		//ConstraintSegment2* pChild0,
		//ConstraintSegment2* pChild1,
		//bool bUpdateCMGR);
/// @private
void updateSplittedConstraintSegment(
		ConstraintSegment2* pCSeg,
		bool bUpdateCMGR);
/**
* \return the Delaunay class it belongs to
*/
	Dt2* getDt2();
/** \brief Get direct children
*
* \param [in] pParent is a ConstraintSegment that may have been splitted
* \param [out] pChild0, pChild1 are the direct child segments of
* \p pParent. They can be alive or dead (splitted).
*
* The children are returned in the correct order of the present
* ConstraintGraph2.
*/
	CLASS_DECLSPEC
	void getDirectChildren(ConstraintSegment2* pParent,ConstraintSegment2*& pChild0,ConstraintSegment2*& pChild1);
/// @private
	void getAliveConstraintChain(std::vector<ConstraintSegment2*>& vAliveCSeg) ; // For debugging
	/// @private
	void setDirectionsRecursive(std::vector<ConstraintSegment2*>& vCSegments);
/** Get the orientation of a ConstraintSegment2
*
* A ConstraintSegment2 \p pCSeg is unoriented because it may participate
* (with different orientations) in more than just one ConstraintGraph2
* and thus the vertices returned by pCSeg->getSrc() and pCSeg->getTrg()
* do not carry any orientation information. However, the orientation of
* \p pCSeg is stored in the ConstraintGraph2 objects where \p pCSeg is
* a member and this method returns if the source and target vertex must
* be exchanged to match the present graph's direction.
*
*
*/
	CLASS_DECLSPEC
	bool isReverse(ConstraintSegment2* pCSeg) const;
	//void getAliveConstraintChain_old(std::vector<ConstraintSegment2*>& vAliveCSeg);
protected:

/// @private
	bool checkAndSortPolygon(std::vector<ConstraintSegment2*>& vCSegments_);
/// @private
	bool checkAndSortPolygonSub(std::vector<ConstraintSegment2*>& vCSegments);
/// @private
	void makeSelfOwner(std::vector<ConstraintSegment2*>& vCSeg);

	// Data
	Dt2* pDt2;
	GeomTest* pGeomPredicates;
	ConstraintInsertionStrategy cis;
	std::vector<ConstraintSegment2*> vCSegParents;
	bool bIsPolygon;
	std::map<ConstraintSegment2*,bool,func_ltDerefPtr<ConstraintSegment2*> > mCSegReverse;
	std::map<Point2*,size_t> mSplitPointNum;
	bool bIsOriented;
};

} // (namespace)
