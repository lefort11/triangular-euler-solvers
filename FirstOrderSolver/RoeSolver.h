#ifndef TRIANGULAR_SOLVERS_ROESOLVER_H
#define TRIANGULAR_SOLVERS_ROESOLVER_H

#include "Solver.h"


namespace euler
{


	class RoeSolver: public Solver
	{

	public:


		explicit RoeSolver(std::vector<Zone> const& constraints,
						   std::function<void(TriangularMesh const&, TriangularMesh const&, double)> const& bcFunc,
						   MeshParams const& triangleProp,
						   double gamma = 5.0/3.0): Solver(constraints, bcFunc, triangleProp, gamma)
		{}

	protected:

		/**@brief Calculates a LaxFriedrich numerical flux
		 *
		 * @param qVec current qvector corresponding to this triangle
		 * @param triangleNumber
		 * @param edgeNumber
		 * @return Numerical flux vector
		 */
		virtual Vec4 CalculateFlux(Vec4 const& qVec, int triangleNumber, int edgeNumber) const override;

		Vec4 Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, int edgeNumber, int gPointNumber) const override
		{
			return qVec;
		}

		void CreateBoundingMesh() override;




	};

}

#endif //TRIANGULAR_SOLVERS_ROESOLVER_H
