#ifndef TRIANGULAR_SOLVERS_WENOLF_H
#define TRIANGULAR_SOLVERS_WENOLF_H

#include "../FirstOrderSolver/LaxFriedrichSolver.h"



namespace euler
{
	class WENOLF : public LaxFriedrichSolver
	{
	private:

		double const m_eps = 1e-4;

		static int const gaussian_points_number = 6;

		struct FOPolynomialCoeffs
		{
			std::array<double, 3> c = {0.0, 0.0, 0.0};
		};

		struct SOPolynomialCoeffs
		{
			std::array<double, 9> gammas = {0, 0, 0, 0, 0, 0, 0, 0, 0};
		};

		struct FOReconstructionPolynomial
		{
			std::array<FOPolynomialCoeffs, gaussian_points_number> coeffsAtPoints;
			std::array<int, 3> stencil = {0, 0, 0};
		};

		struct SOReconstructionPolynomial
		{
			std::array<SOPolynomialCoeffs, gaussian_points_number> coeffsAtPoints;

		};

		struct SmoothIndicatorReconstructionData
		{
			std::array<double, 3> alpha = {0.0, 0.0, 0.0};
			std::array<double, 3> beta = {0.0, 0.0, 0.0};
		};

		struct TriangleReconstructionData
		{
			std::array<FOReconstructionPolynomial, 9> fo_polynomial;
			std::array<Point2, gaussian_points_number> gaussian_points;
			SOReconstructionPolynomial so_polynomial;
			std::array<SmoothIndicatorReconstructionData, 9> smoothIndicatorData;
		};

		std::vector<TriangleReconstructionData> m_vReconstructionData;






	public:

		explicit WENOLF(std::vector<Zone> const &constraints,
						std::array<double, 3> const &triangleProp = {0.0, 0.0, 0.0},
						double gamma = 5.0 / 3.0) : LaxFriedrichSolver(constraints, triangleProp,
																	   gamma)
		{}

		void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunction) override;


	protected:

		Vec4 Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, Point2 const& gaussianPoint) const override;


		void FormL_A(Vec4 const& qVec, arma::mat44& L_A) const;
		void FormL_B(Vec4 const& qVec, arma::mat44& L_B) const;

		void FormR_A(Vec4 const& qVec, arma::mat44& R_A) const;
		void FormR_B(Vec4 const& qVec, arma::mat44& R_B) const;

		void GetStencil(Triangle const* pTriangle, std::array<Triangle const*, 10> &stencil) const;

		void CreateBoundingMesh() override;

		void GetSmoothIndicatorData(TriangleReconstructionData& trRecData, Triangle const* pTriangle) const;

		void GetTriangleReconstructionData(TriangleReconstructionData &trRecData, Triangle const* triangle);

	};

}


#endif //TRIANGULAR_SOLVERS_WENOLF_H
