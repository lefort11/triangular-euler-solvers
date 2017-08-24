#include "Area.h"

using namespace euler;


std::vector<GEOM_FADE2D::Point2> ConstraintFunction::Discretize(int const pointNumber) const
{
	auto const h = 1.0 / pointNumber;
	std::vector<GEOM_FADE2D::Point2> vPoints;

	for(int i = 0; i < pointNumber; ++i)
	{
		auto point = (*this)(i*h);
		vPoints.push_back(point);
	}

	return vPoints;

}




TriangularMesh Area::Triangulate(int const discrPointNumber,
								 std::array<double, 3> const& triangleProperties,
								 std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const &initStateFunc)
{

//	GEOM_FADE2D::Fade_2D globalArea;

	//making bounding rectangle
	GEOM_FADE2D::Point2 p1(0.0, 0.0), p2(0.0, 1.0), p3(1.0, 0.0), p4(1.0, 1.0);
	globalArea.insert(p1);
	globalArea.insert(p2);
	globalArea.insert(p3);
	globalArea.insert(p4);


	std::vector<std::vector<GEOM_FADE2D::Segment2>> vvSegment(m_zones.size());
	std::vector<GEOM_FADE2D::ConstraintGraph2 *> vCSG;
	std::vector<GEOM_FADE2D::Zone2 *> vpZonesDealunay;

	for (int i = 0; i < m_zones.size(); ++i)
	{
		//discretizing zone's constraints
		auto vPoints = m_zones[i].Discretize(discrPointNumber);

		//creating segments
		for (int j = 0; j < vPoints.size(); ++j)
		{
			vvSegment[i].emplace_back(GEOM_FADE2D::Segment2(vPoints[j], vPoints[(j + 1) % vPoints.size()]));
		}

		//creating constraint graphs
		vCSG.push_back(globalArea.createConstraint(vvSegment[i], GEOM_FADE2D::CIS_CONSTRAINED_DELAUNAY));


		//creating fade2D delaunay zones
		m_zones[i].IsInside() ? vpZonesDealunay.push_back(globalArea.createZone(vCSG[i], GEOM_FADE2D::ZL_INSIDE)) :
		vpZonesDealunay.push_back(globalArea.createZone(vCSG[i], GEOM_FADE2D::ZL_OUTSIDE));


	}


	auto pGrowZone = globalArea.createZone(vCSG, GEOM_FADE2D::ZL_GROW, p1);

	//calculating final zone
	auto const size = static_cast<int>(vpZonesDealunay.size() - 1);
	for (int i = 0; i < size; ++i)
	{
		pGrowZone = GEOM_FADE2D::zoneIntersection(vpZonesDealunay[i], vpZonesDealunay[i + 1]);
	}


	globalArea.applyConstraintsAndZones();



//	pGrowZone->show("kekas.ps", false, true);


	//refining final zone
	auto pBoundedZone(pGrowZone->convertToBoundedZone());

	GEOM_FADE2D::MeshGenParams params(pBoundedZone);
	params.minAngleDegree = triangleProperties[0];
	params.minEdgeLength = triangleProperties[1];
	params.maxEdgeLength = triangleProperties[2];

//	globalArea.refine(pBoundedZone, triangleProperties[0], triangleProperties[1], triangleProperties[2], true);
	globalArea.refineAdvanced(&params);

	pBoundedZone->show("lul.ps", false, true);


//	globalArea.show("kek.ps");


	std::vector<GEOM_FADE2D::Triangle2*> vTriangles2;

	pBoundedZone->getTriangles(vTriangles2);


	std::vector<Triangle*> vTriangle;

	CopyGraph(vTriangle, vTriangles2);



	//setting init state for each triangle
	for (auto &triangle : vTriangle)
	{
		auto const point = triangle->getBarycenter();
		triangle->SetState(initStateFunc(point));
	}


	return TriangularMesh(vTriangle);
//	return TriangularMesh(triangles);


}

