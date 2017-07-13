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
								 std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const &initStateFunc) const
{

	GEOM_FADE2D::Fade_2D globalArea;

	//making bounding rectangle
	GEOM_FADE2D::Point2 p1(-5, -5), p2(-5, 5), p3(5, -5), p4(5, 5);
	globalArea.insert(p1);
	globalArea.insert(p2);
	globalArea.insert(p3);
	globalArea.insert(p4);


	std::vector<std::vector<GEOM_FADE2D::Segment2>> vvSegment(m_Zones.size());
	std::vector<GEOM_FADE2D::ConstraintGraph2*> vCSG;
	std::vector<GEOM_FADE2D::Zone2*> vpZonesDealunay;

	for(int i = 0; i < m_Zones.size(); ++i)
	{
		//discretizing zone's constraints
		auto vPoints = m_Zones[i].Discretize(discrPointNumber);

		//creating segments
		for(int j = 0; j < vPoints.size(); ++j)
		{
			vvSegment[i].push_back(GEOM_FADE2D::Segment2(vPoints[j], vPoints[(j+1) % vPoints.size()]));
		}

		//creating constraint graphs
		vCSG.push_back(globalArea.createConstraint(vvSegment[i], GEOM_FADE2D::CIS_CONSTRAINED_DELAUNAY));


		//creating fade2D delaunay zones
		m_Zones[i].IsInside() ? vpZonesDealunay.push_back(globalArea.createZone(vCSG[i], GEOM_FADE2D::ZL_INSIDE)):
								vpZonesDealunay.push_back(globalArea.createZone(vCSG[i], GEOM_FADE2D::ZL_OUTSIDE));



	}



	auto pGrowZone = globalArea.createZone(vCSG, GEOM_FADE2D::ZL_GROW, p1);

	//calculating final zone
	for(int i = 0; i < vpZonesDealunay.size() - 1; ++i)
	{
		pGrowZone = GEOM_FADE2D::zoneIntersection(vpZonesDealunay[i], vpZonesDealunay[i+1]);
	}


	globalArea.applyConstraintsAndZones();

	std::vector<GEOM_FADE2D::Triangle2*> triangles;

	pGrowZone->getTriangles(triangles);

	pGrowZone->show("kekas.ps", false, true);


	//refining final zone
	auto pBoundedZone(pGrowZone->convertToBoundedZone());

	globalArea.refine(pBoundedZone, triangleProperties[0], triangleProperties[1], triangleProperties[2], true);

	pBoundedZone->show("lul.ps", false, true);


	globalArea.show("kek.ps");


	std::vector<GEOM_FADE2D::Triangle2*> vTriangles;

	pBoundedZone->getTriangles(vTriangles);

	//setting init state for each triangle
	std::vector<std::array<double, 4>> vInitStates;

	for(int i = 0; i < vTriangles.size(); ++i)
	{
		auto const point = vTriangles[i]->getBarycenter();
		vInitStates.push_back(initStateFunc(point));
	}


	//i-th element of vInitStates corresponing to i-th triangle

	return TriangularMesh(vTriangles, vInitStates);


}

