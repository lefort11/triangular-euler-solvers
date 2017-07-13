#ifndef TRIANGULAR_SOLVERS_AREA_H
#define TRIANGULAR_SOLVERS_AREA_H

#include <vector>
#include <functional>

#include "TriangularMesh.h"

namespace euler
{

	class ConstraintFunction: public std::function<GEOM_FADE2D::Point2(double)> // function of argument t from 0 to 1
	{

	public:

		ConstraintFunction(std::function<GEOM_FADE2D::Point2(double)> func): std::function<GEOM_FADE2D::Point2(double)>(func)
		{}


		std::vector<GEOM_FADE2D::Point2> Discretize(int const pointNumber) const;


	};


	//!@brief Zone class describes a zone bounded by some functions.
	class Zone
	{
		const ConstraintFunction m_Function;
		const bool m_IsInside;
		//const std::array<double, 4> m_InitialState;

	public:
		Zone(ConstraintFunction const& func, bool isInside, std::array<double, 4> const& initState): m_Function(func),
																									 m_IsInside(isInside)
		{}


		std::vector<GEOM_FADE2D::Point2> Discretize(int const pointNumber) const
		{
			return m_Function.Discretize(pointNumber);
		}

		bool IsInside() const
		{
			return m_IsInside;
		}


	};

	//!@brief Area class describes an area inside the bounding rectangle
	class Area
	{
		std::vector<Zone> m_Zones;

	public:
		Area(std::vector<Zone> const& constraints): m_Zones(constraints)
		{}

		//!@param triangleProperties = {minAngleDegree, minEdgeLength, maxEdgeLength}
		TriangularMesh Triangulate(int const discrPointNumber,
								   std::array<double, 3> const& triangleProperties,
								   std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunc) const;

	private:


	};

}

/* сначала дискретизируем границы, получая вектор точек для каждой границы. потом эти точки пихаем куда нужно,
 * сегменты, зоны все дела, потом ZL_INSIDE или OUTSIDE и пересекаем */


#endif //TRIANGULAR_SOLVERS_AREA_H
