#ifndef TRIANGULAR_SOLVERS_SOLVER_H
#define TRIANGULAR_SOLVERS_SOLVER_H

#include "../SpaceMesh/Area.h"


namespace euler
{

	class Solver
	{
		Area m_area;

		int m_discrPointNumber;

		std::array<double, 3> m_triangleProperties;

		TriangularMesh m_mesh;

		const std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> initStateFunc;

	public:
		Solver(std::vector<Zone> const& constraints, int const discrPointNumber = 0,
			   std::array<double, 3> const& triangleProp = {0.0, 0.0, 0.0}): m_area(constraints),
																			 m_discrPointNumber(discrPointNumber),
																			 m_triangleProperties(triangleProp)
		{}

		void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& function)
		{
			m_mesh = m_area.Triangulate(m_discrPointNumber, m_triangleProperties, function);
		}

		void Calculate();


	protected:

		/**@brief calculates flux between triangleNumber-th triangle and its edgeNumber-th opposite triangle
		 *
		 * @param triangleNumber - number of triangle
		 *
		 * @param edgeNumber - number of edge - 'cell's boundary'
		**/
		virtual void CalculateFlux(int const triangleNumber, int const edgeNumber) const = 0;

	};

}
#endif //TRIANGULAR_SOLVERS_SOLVER_H
