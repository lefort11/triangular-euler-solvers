#ifndef TRIANGULAR_SOLVERS_SOLVER_H
#define TRIANGULAR_SOLVERS_SOLVER_H

#include "../SpaceMesh/Area.h"

#include "../Maths/Algebra.h"


namespace euler
{

	class Solver
	{
		Area m_area;


		std::array<double, 3> m_triangularizationProperties;


		//const std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> m_initStateFunc;

	protected:

		TriangularMesh m_triangles;

		TriangularMesh m_boundingTriangles;

		double const m_gamma;

		std::function<void(TriangularMesh const& boundaryMesh, TriangularMesh const& mainMesh)>
																					m_boundaryConditionFunction;

	public:



		explicit Solver(std::vector<Zone> const& constraints,
						std::function<void(TriangularMesh const&, TriangularMesh const&)>  const& bcFunc,
						std::array<double, 3> const& triangleProp = {0.0, 0.0, 0.0},
						double gamma = 5.0/3.0): m_area(constraints),
												 m_triangularizationProperties(triangleProp),
												 m_boundaryConditionFunction(bcFunc),
												 m_gamma(gamma)
		{}

		/**@brief method initializes the mesh with initial state function
		 *
		 * @param initStateFunction initial state funtion returning state vector corresponding to some point
		 *
		**/

		virtual void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunction)
		{
			m_triangles = m_area.Triangulate(m_triangularizationProperties, initStateFunction);

			CreateBoundingMesh();
		}

		/**@brief method does all calculations
		 *
		 * @param time time of the result
		 *
		**/

		void Calculate(double time) const;


	private:

		/**@brief Gaussian quadrature of order 3 for the standard quadrilateral element R = [1;1]x[1;1]
		 *
		 * @param f function to be integrated
		 * @return
		 */
		double GaussianQuadrilateralIntegration(std::function<double(double, double)> const& g) const;


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

		/**@brief Method reconstructs the q vector at the gaussian point (x_g, y_g)
		 *
		 * @param qVec
		 * @param triangleNumber
		 * @param x_g
		 * @param y_g
		 * @return
		 */
		virtual Vec4
		Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, Point2 const &gaussianPoint, int edgeNumber) const = 0;


		/**@brief Method generates bounding mesh
		 *
		 */
		virtual void CreateBoundingMesh() = 0;

		/**@brief Method updates bounding triangles according to boundary conditions
		 *
		 */
		void UpdateBoundingMesh() const
		{
			m_boundaryConditionFunction(m_boundingTriangles, m_triangles);
		}


		/**@brief Third-order calculation of an integral of a scalar function
		 * over the triangle pTriangle according to Gaussian formulae
		 *
		 * "Appropriate Gaussian quadrature formulae for triangles" // Farzana Hussain, M. S. Karima, Razwan Ahamada
		 *
		 *
		 * @param func
		 * @param pTriangle
		 * @return
		 */
		double GaussianIntegration(std::function<double(double, double)> const& func, Triangle const* pTriangle) const;


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
		 * @return outward normal's coordinates
		**/
		std::array<double, 2> CalculateNormal(Triangle const* pTriangle, int edgeNumber) const
		{
			return pTriangle->CalculateNormal(edgeNumber);
		};

		/**
		 *
		 * @return Time step calculated according to Courant criterion
		 */

		virtual double CalculateTimeStep() const
		{

			double const sigma = 0.01;
			auto min_area = m_triangles[0]->getArea2D();

			auto lambda_max = 0.0;

			for(int i =0; i < m_triangles.size(); ++i)
			{

				if(min_area > m_triangles[i]->getArea2D())
					min_area = m_triangles[i]->getArea2D();

				auto const velocity_sqr_abs = std::sqrt(sqr(m_triangles[i]->velocityX) + sqr(m_triangles[i]->velocityY));
				auto const sound_speed = std::sqrt(m_gamma * m_triangles[i]->pressure / m_triangles[i]->density);
				if(lambda_max < std::fabs(velocity_sqr_abs + sound_speed))
					lambda_max = std::fabs(velocity_sqr_abs + sound_speed);


			}

			return sigma * std::sqrt(min_area) / (2 * lambda_max);

		}


		/**
		 * @brief Method forms a q-vector using triangle's data members
		 *
		 * @param qVec Q vector wanted to being formed
		 * @param triangle
		 */
		void FormQVector(Vec4 & qVec, Triangle const*  triangle) const
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

/*			double min_y = 10.0;
			for(int i = 0; i < m_triangles.size(); ++i)
			{
				auto const barycenter = m_triangles[i]->getBarycenter();
				if(barycenter.y() < min_y)
					min_y = barycenter.y();

			} */

			for(int i = 0; i < m_triangles.size(); ++i)
			{
				auto const barycenter = m_triangles[i]->getBarycenter();
				if(std::fabs(barycenter.y() - 0.0) < 0.1 )
				{
					densityResultsFile << barycenter.x() << " " << m_triangles[i]->density << std::endl;
				}
			}

			densityResultsFile.close();

		}

		void Output(std::string const& densityFilename, std::string const& velocityFilename,
					std::string const& pressureFilename) const
		{
			std::ofstream densityResultsFile, velocityResultsFile, pressureResultsFile;

			densityResultsFile.open(densityFilename);
			velocityResultsFile.open(velocityFilename);
			pressureResultsFile.open(pressureFilename);

			double maxVelocityAbs = 0.0;

			for(int i = 0; i < m_triangles.size(); ++i)
			{
				auto const velocityAbs = std::sqrt(sqr(m_triangles[i]->velocityX) + sqr(m_triangles[i]->velocityY));
				if(maxVelocityAbs < velocityAbs)
					maxVelocityAbs = velocityAbs;
			}

			for(int i = 0; i < m_triangles.size(); ++i)
			{
				auto const barycenter = m_triangles[i]->getBarycenter();
				densityResultsFile << barycenter.x() << " " << barycenter.y() << " "
								   << m_triangles[i]->density << std::endl;

				velocityResultsFile << barycenter.x() << " " << barycenter.y() << " "
								   << m_triangles[i]->velocityX / (10 * maxVelocityAbs)
									<< " " << m_triangles[i]->velocityY / (10 * maxVelocityAbs) << std::endl;
				pressureResultsFile << barycenter.x() << " " << barycenter.y() << " "
								   << m_triangles[i]->pressure << std::endl;
			}

			densityResultsFile.close();
			velocityResultsFile.close();
		}



	};




}
#endif //TRIANGULAR_SOLVERS_SOLVER_H
