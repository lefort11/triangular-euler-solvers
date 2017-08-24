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


		//const std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> m_initStateFunc;

	protected:

		TriangularMesh m_triangles;

		double const m_gamma;

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
		 * @param initStateFunction initial state funtion returning state vector corresponding to some point
		 *
		**/

		void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunction)
		{
			m_triangles = m_area.Triangulate(m_discrPointNumber, m_triangularizationProperties, initStateFunction);

//			std::cout<<m_triangles[0]->getCorner(0)->x() << " " << m_triangles[0]->getCorner(0)->y() << std::endl;
		}

		/**@brief method does all calculations
		 *
		 * @param time time of the result
		 *
		**/

		void Calculate(double time) const;


	protected:

		/**@brief method calculates flux between triangleNumber-th triangle and its edgeNumber-th opposite triangle
		 *
		 * @param qVec q Vector corresponing triangleNumber-th triangle. It is needed for Runge-Kutta method
		 * implementation
		 *
		 * @param triangleNumber number of triangle
		 *
		 * @param edgeNumber number of triangle's edge - 'cell's boundary'
		**/
		virtual Vec4 CalculateFlux(Vec4 const& qVec, int triangleNumber, int edgeNumber) const = 0;

//		virtual Vec4 GaussianIntegration(Vec4 const& qVec) const = 0;


		/**@brief method makes an iteration of Runge Kutta method modernized for TVD schemes
		 * 		Gottlieb, S., Shu, C.-W.: Total Variation Diminishing Runge-Kutta Schemes.
		 * 		Math of Computation, 73–85 (1998)
		 *
		 * @param current_q q_0
		 * @param f	right part of the diff equation
		 * @param delta_t time step
		 * @return q-vector at next time layer
		**/
		Vec4 RungeKuttaTVDStep(Vec4 const& current_q, double delta_t, std::function<Vec4(Vec4)> const& f) const;

		/**
		 * @return triangleNumber-th triangle's area
		**/

		double CalculateTriangleArea(int triangleNumber) const
		{
			return m_triangles[triangleNumber]->getArea2D();
		}

		/**
		 * @return length of the edgeNumber-th edge of triangleNumber-th triangle
		**/
		double CalculateTriangleEdgeLength(int triangleNumber, int edgeNumber) const
		{
			return std::sqrt(m_triangles[triangleNumber]->getSquaredEdgeLength(edgeNumber));
		}


		/**
		 * @return outer normal's coordinates
		**/
		std::array<double, 2> CalculateNormal(int triangleNumber, int edgeNumber) const
		{
			auto const vertex1 = m_triangles[triangleNumber]->getCorner((edgeNumber + 1) % 3);
			auto const vertex2 = m_triangles[triangleNumber]->getCorner((edgeNumber + 2) % 3);

			auto const p_x = vertex2->x() - vertex1->x();
			auto const p_y = vertex2->y() - vertex1->y();

			auto const p_abs = std::sqrt(p_x * p_x + p_y * p_y);

			return {p_y/p_abs, -p_x/p_abs};
		};

		/**
		 *
		 * @return Time step calculated according to Courant criterion
		 */

		virtual double CalculateTimeStep() const
		{

			double const sigma = 0.1;

			auto min_area = m_triangles[0]->getArea2D();
			auto denominator = 0.0;
			for(int i = 0; i < m_triangles.size(); ++i)
			{
				//lookin' for smallest area
				auto const area = m_triangles[i]->getArea2D();
				if(min_area > area)
					min_area = area;

				auto const longestEdge = sqrt(m_triangles[i]->getSquaredEdgeLength(m_triangles[i]->getMaxIndex()));
				auto const c0 = sqrt(m_gamma * m_triangles[i]->pressure / m_triangles[i]->density);
				auto const velocity_abs = sqrt(sqr(m_triangles[i]->velocityX) + sqr(m_triangles[i]->velocityY));

				auto const lambda = fabs(velocity_abs + c0) > fabs(velocity_abs - c0) ?
									fabs(velocity_abs + c0) : fabs(velocity_abs - c0); //eigenvalue for each cell

				denominator += longestEdge * lambda;

//				if(area == 0)
//					double d = 0;


			}

			return sigma * min_area / denominator;

		}


		/**
		 * @brief Method forms a q-vector using triangle's data members
		 *
		 * @param qVec Q vector wanted to being formed
		 * @param triangle
		 */
		void FormQVector(Vec4 & qVec, Triangle* const triangle) const
		{
			auto const density = triangle->density;
			auto const velocityX = triangle->velocityX;
			auto const velocityY = triangle->velocityY;
			auto const pressure = triangle->pressure;
			auto const eps = pressure / (density * (m_gamma - 1.0));
			auto const E = eps + 0.5 * (sqr(velocityX) + sqr(velocityY)); // E = eps + 1/2 * |v|^2

			qVec[0] = density;
			qVec[1] = density * velocityX;
			qVec[2] = density * velocityY;
			qVec[3] = density * E;
		}


		/**
		 * @brief Method forms a q-vector using gas parameters
		 *
		 * @param qVec Q vector wanted to being formed
		 * @param density
		 * @param velocityX
		 * @param velocityY
		 * @param pressure
		 */
		void FormQVector(Vec4 & qVec, double density, double velocityX, double velocityY, double pressure) const
		{
			auto const eps = pressure / (density * (m_gamma - 1.0));
			auto const E = eps + 0.5 * (sqr(velocityX) + sqr(velocityY)); // E = eps + 1/2 * |v|^2

			qVec[0] = density;
			qVec[1] = density * velocityX;
			qVec[2] = density * velocityY;
			qVec[3] = density * E;
		}
		/**@brief Gets gas parameters from a q-vector
		 *
		 * @param qVec
		 * @param density
		 * @param velocityX
		 * @param velocityY
		 * @param pressure
		 */

		void GetGasParamsFromQ(Vec4 const& qVec, double& density, double& velocityX,
							   double& velocityY, double& pressure) const
		{
			density = qVec[0];
			velocityX = qVec[1] / density;
			velocityY = qVec[2] / density;
			auto const E = qVec[3] / density;
			auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);
			auto const eps = E - 0.5 * velocity_sqr_abs;
			pressure = (m_gamma - 1.0) * eps * density;

		}

	public:

		void DebugOutput(std::string const& filename) const
		{
			std::ofstream densityResultsFile;

			densityResultsFile.open(filename);

			double min_y = 10.0;
			for(int i = 0; i < m_triangles.size(); ++i)
			{
				auto const barycenter = m_triangles[i]->getBarycenter();
				if(barycenter.y() < min_y)
					min_y = barycenter.y();

			}

			for(int i = 0; i < m_triangles.size(); ++i)
			{
				auto const barycenter = m_triangles[i]->getBarycenter();
				if(std::fabs(barycenter.y() - min_y) < 0.3 )
				{
					densityResultsFile << barycenter.x() << " " << m_triangles[i]->density << std::endl;
				}
			}

			densityResultsFile.close();

		}



	};




}
#endif //TRIANGULAR_SOLVERS_SOLVER_H