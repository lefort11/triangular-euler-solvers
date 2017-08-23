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

		explicit ConstraintFunction(std::function<GEOM_FADE2D::Point2(double)> const& func):
				std::function<GEOM_FADE2D::Point2(double)>(func)
		{}


		std::vector<GEOM_FADE2D::Point2> Discretize(int pointNumber) const;


	};


	//!@brief Zone class describes a zone bounded by some functions.
	class Zone
	{
		const ConstraintFunction m_function;
		const bool m_isInside;


		const std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> m_boundaryConditions; //returns {ro, u, v, P}

	public:
		Zone(ConstraintFunction const& func, bool isInside, std::array<double, 4> const& initState): m_function(func),
																									 m_isInside(isInside)
		{}


		std::vector<GEOM_FADE2D::Point2> Discretize(int const pointNumber) const
		{
			return m_function.Discretize(pointNumber);
		}

		bool IsInside() const
		{
			return m_isInside;
		}


	};

	//!@brief Area class describes an area inside the bounding rectangle
	class Area
	{
		std::vector<Zone> m_zones;

	public:

		Area(std::vector<Zone> const& constraints): m_zones(constraints)
		{}


		/**@brief trianglulate area consisted of zones
		 * @param triangleProperties = {minAngleDegree, minEdgeLength, maxEdgeLength}
		 * @param discrPointNumber - point number of constraint function discretezation
		 * @param initStateFunc - initial function of the problem
		**/
		TriangularMesh Triangulate(int discrPointNumber,
								   std::array<double, 3> const& triangleProperties,
								   std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunc) const;

	private:


	};

}


#endif //TRIANGULAR_SOLVERS_AREA_H
