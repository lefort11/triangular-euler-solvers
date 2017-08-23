#ifndef TRIANGULAR_SOLVERS_SOLVER_H
#define TRIANGULAR_SOLVERS_SOLVER_H

#include "../SpaceMesh/Area.h"

#include "../Maths/Algebra.h"


namespace euler
{

	class Solver
	{
		Area m_area;

		int m_discrPointNumber;

		std::array<double, 3> m_triangularizationProperties;

		TriangularMesh m_triangles;

		double const m_gamma;

		//const std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> m_initStateFunc;

	public:



		explicit Solver(std::vector<Zone> const& constraints, int const discrPointNumber = 0,
			   std::array<double, 3> const& triangleProp = {0.0, 0.0, 0.0},
						double gamma = 5.0/3.0): m_area(constraints),
												 m_discrPointNumber(discrPointNumber),
												 m_triangularizationProperties(triangleProp),
												 m_gamma(gamma)
		{}

		/**@brief method initializes the mesh with initial state function
		 *
		 * @param initStateFunction - initial state funtion returning state vector corresponding to some point
		 *
		**/

		void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunction)
		{
			m_triangles = m_area.Triangulate(m_discrPointNumber, m_triangularizationProperties, initStateFunction);
		}

		/**@brief method does all calculations
		 *
		 * @param time - time of the result
		 *
		**/

		void Calculate(double time);


	protected:

		/**@brief method calculates flux between triangleNumber-th triangle and its edgeNumber-th opposite triangle
		 *
		 * @param triangleNumber - number of triangle
		 *
		 * @param edgeNumber - number of triangle's edge - 'cell's boundary'
		**/
		virtual Vec4 CalculateFlux(Vec4 const& qVec, int triangleNumber, int edgeNumber) const = 0;


		/**@brief method makes an iteration of Runge Kutta method modernized for TVD schemes
		 * 		Gottlieb, S., Shu, C.-W.: Total Variation Diminishing Runge-Kutta Schemes.
		 * 		Math of Computation, 73–85 (1998)
		 *
		 * @param current_q - q_0
		 * @param f	- right part of the diff equation
		 * @param delta_t - time step
		 * @return q-vector at next time layer
		**/
		Vec4 RungeKuttaTVDStep(Vec4 const& current_q, double delta_t, std::function<Vec4(Vec4)> const& f) const;

		/**
		 * @return triangleNumber-th triangle's area
		**/

		double CalculateTriangleArea(int triangleNumber)
		{
			return m_triangles[triangleNumber]->getArea2D();
		}

		/**
		 * @return length of the edgeNumber-th edge of triangleNumber-th triangle
		**/
		double CalculateTriangleEdgeLength(int triangleNumber, int edgeNumber)
		{
			return std::sqrt(m_triangles[triangleNumber]->getSquaredEdgeLength(edgeNumber));
		}


		/**
		 * @return outer normal's coordinates
		**/
		std::array<double, 2> CalculateNormal(int triangleNumber, int edgeNumber)
		{
			auto const vertex1 = m_triangles[triangleNumber]->getCorner((edgeNumber + 1) % 3);
			auto const vertex2 = m_triangles[triangleNumber]->getCorner((edgeNumber + 2) % 3);

			auto const p_x = vertex2->x() - vertex1->x();
			auto const p_y = vertex2->y() - vertex1->y();

			return {p_y, -p_x};
		};



	};




}
#endif //TRIANGULAR_SOLVERS_SOLVER_H
