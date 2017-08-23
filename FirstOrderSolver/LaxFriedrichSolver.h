#ifndef TRIANGULAR_SOLVERS_LAXFRIEDRICHSOLVER_H
#define TRIANGULAR_SOLVERS_LAXFRIEDRICHSOLVER_H

#include "Solver.h"


namespace euler
{

	class LaxFriedrichSolver: public Solver
	{

	public:

		explicit LaxFriedrichSolver(std::vector<Zone> const& constraints, int const discrPointNumber = 0,
				std::array<double, 3> const& triangleProp = {0.0, 0.0, 0.0},
		double gamma = 5.0/3.0): Solver(constraints, discrPointNumber, triangleProp, gamma)
		{}

	protected:

		virtual Vec4 CalculateFlux(Vec4 const& qVec, int triangleNumber, int edgeNumber) const;


	};

}


#endif //TRIANGULAR_SOLVERS_LAXFRIEDRICHSOLVER_H
