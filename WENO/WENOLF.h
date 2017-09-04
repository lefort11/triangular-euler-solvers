#ifndef TRIANGULAR_SOLVERS_WENOLF_H
#define TRIANGULAR_SOLVERS_WENOLF_H

#include "../FirstOrderSolver/LaxFriedrichSolver.h"



namespace euler
{
	class WENOLF : public LaxFriedrichSolver
	{
	private:

		class FOReconstructionPolynomial
		{
			std::array<double, 3> coeff;
			std::array<int, 3> stencil;

		public:
			FOReconstructionPolynomial(int tr_0, int tr_1, int tr_2): coeff{{0,0,0}}, stencil{{tr_0, tr_1, tr_2}}
			{}
		};

	public:

		explicit WENOLF(std::vector<Zone> const &constraints, int const discrPointNumber = 0,
						std::array<double, 3> const &triangleProp = {0.0, 0.0, 0.0},
						double gamma = 5.0 / 3.0) : LaxFriedrichSolver(constraints, discrPointNumber, triangleProp,
																	   gamma)
		{}


	protected:

		Vec4 Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, double x_g, double y_g) const override;

//	Vec4 GetRiemannInvariantVector(Vec4 const& qVec) const;

		Vec4 GetFirstOrderPolynomialReconstruction(Matrix4x4 const &mat,
												   Triangle const *pTriangle, double x_g, double y_g) const;

		Vec4 GetQFromInvariant(Vec4 const &wVec) const;

	};

}


#endif //TRIANGULAR_SOLVERS_WENOLF_H
