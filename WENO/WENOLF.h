#ifndef TRIANGULAR_SOLVERS_WENOLF_H
#define TRIANGULAR_SOLVERS_WENOLF_H

#include "../FirstOrderSolver/LaxFriedrichSolver.h"



namespace euler
{
	class WENOLF : public LaxFriedrichSolver
	{
	private:

		static int const gaussian_points_number = 6;

		struct FOReconstructionPolynomial
		{
			std::array<double, 3> coeff;
			std::array<int, 3> stencil;

		public:
			FOReconstructionPolynomial(int tr_0 = 0, int tr_1 = 0, int tr_2 = 0):
					coeff{{0,0,0}}, stencil{{tr_0, tr_1, tr_2}}
			{}
		};

		struct OnePointReconstructionData
		{
			Point2 gaussian_point;
			std::array<FOReconstructionPolynomial, 9> polynomial;

			std::array<double, 9> third_order_coeff; // gammas
		};

		struct TriangleReconstructionData: public std::array<OnePointReconstructionData, gaussian_points_number>
		{


		};


		std::vector<TriangleReconstructionData> m_vReconstructionData;



	public:

		explicit WENOLF(std::vector<Zone> const &constraints, int const discrPointNumber = 0,
						std::array<double, 3> const &triangleProp = {0.0, 0.0, 0.0},
						double gamma = 5.0 / 3.0) : LaxFriedrichSolver(constraints, discrPointNumber, triangleProp,
																	   gamma)
		{}

		void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunction) override;


	protected:

		Vec4 Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, Point2 const& gaussianPoint) const override;


		void FormL_A(Vec4 const& qVec, arma::mat44& L_A) const;
		void FormL_B(Vec4 const& qVec, arma::mat44& L_B) const;

		void FormR_A(Vec4 const& qVec, arma::mat44& R_A) const;
		void FormR_B(Vec4 const& qVec, arma::mat44& R_B) const;

		void GetPointReconstructionData(OnePointReconstructionData &data,
										std::array<Triangle const *, 10> const &stencil, Point2 const &gaussian_point) const;

		void GetStencil(Triangle const* pTriangle, std::array<Triangle const*, 10> &stencil) const;

		void CreateBoundingMesh() override;


		double CalculateSmoothIndicator(int triangle_number, int pol_number,
										std::array<Triangle const*, 10> const& stencil) const;

	};

}


#endif //TRIANGULAR_SOLVERS_WENOLF_H
