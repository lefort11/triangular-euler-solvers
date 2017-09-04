#include "WENOLF.h"


using namespace euler;

Vec4 WENOLF::Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, double x_g, double y_g) const
{

	double ksi_average[10];
	double eta_average[10];

	auto const h = std::sqrt(pTriangle->getArea2D());
	auto const x_0 = pTriangle->getBarycenter().x();
	auto const y_0 = pTriangle->getBarycenter().y();


	double ksi_square_average[10];
	double eta_square_average[10];
	double ksi_eta_average[10];

	Vec4 q[10];



	// ************* Getting the stencil ************** //
	// @todo check if triangle has no neighbour // boundary conditions kek
	Triangle const* triangle[10];
	triangle[0] = pTriangle;
	q[0] = qVec;
	triangle[1] = triangle[0]->GetOppTriangle(0) == nullptr? triangle[0] : triangle[0]->GetOppTriangle(0);
	FormQVector(q[1], triangle[1]);
	triangle[2] = triangle[0]->GetOppTriangle(1) == nullptr? triangle[0] : triangle[0]->GetOppTriangle(1);
	FormQVector(q[2], triangle[2]);
	triangle[3] = triangle[0]->GetOppTriangle(2) == nullptr? triangle[0] : triangle[0]->GetOppTriangle(2);
	FormQVector(q[3], triangle[3]);

	triangle[4] = triangle[1]->GetOppTriangle(triangle[1]->getIntraTriangleIndex(triangle[0]->getCorner(2)));
	if(triangle[4] == nullptr)
		triangle[4] = triangle[1];
	FormQVector(q[4], triangle[4]);
	triangle[5] = triangle[1]->GetOppTriangle(triangle[1]->getIntraTriangleIndex(triangle[0]->getCorner(1)));
	if(triangle[5] == nullptr)
		triangle[5] = triangle[1];
	FormQVector(q[5], triangle[5]);

	triangle[6] = triangle[3]->GetOppTriangle(triangle[3]->getIntraTriangleIndex(triangle[0]->getCorner(1)));
	if(triangle[6] == nullptr)
		triangle[6] = triangle[3];
	FormQVector(q[6], triangle[6]);
	triangle[7] = triangle[3]->GetOppTriangle(triangle[3]->getIntraTriangleIndex(triangle[0]->getCorner(0)));
	if(triangle[7] == nullptr)
		triangle[7] = triangle[3];
	FormQVector(q[7], triangle[7]);

	triangle[8] = triangle[2]->GetOppTriangle(triangle[2]->getIntraTriangleIndex(triangle[0]->getCorner(0)));
	if(triangle[8] == nullptr)
		triangle[8] = triangle[2];
	FormQVector(q[8], triangle[8]);
	triangle[9] = triangle[2]->GetOppTriangle(triangle[2]->getIntraTriangleIndex(triangle[0]->getCorner(2)));
	if(triangle[9] == nullptr)
		triangle[9] = triangle[2];
	FormQVector(q[9], triangle[9]);
	// *************************************************** //


	for(int i = 0; i < 10; ++i)
	{
		auto const area = triangle[i]->getArea2D();
		ksi_average[i] = 1 / area * GaussianIntegration([h, x_0](double x, double y)
														{
															return (x - x_0) / h;
														}, triangle[i]);

		eta_average[i] = 1 / area * GaussianIntegration([h, y_0](double x, double y)
														{
															return (y - y_0) / h;
														}, triangle[i]);

		ksi_square_average[i] = 1 / area * GaussianIntegration([h, x_0](double x, double y)
															   {
																   return sqr(x-x_0) / sqr(h);
															   }, triangle[i]);
		eta_square_average[i] = 1 / area * GaussianIntegration([h, y_0](double x, double y)
															   {
																   return sqr(y-y_0) / sqr(h);
															   }, triangle[i]);
		ksi_eta_average[i] = 1 / area * GaussianIntegration([h, x_0, y_0](double x, double y)
															   {
																   return (x-x_0) * (y - y_0)  / sqr(h);
															   }, triangle[i]);

	}

	arma::vec3 b;
	b << 1 << (x_g - x_0) / h << (y_g - y_0) / h;


	arma::vec3 polynomial_coeff[9];

	// ******* p_0 coeff calculation **** //
	arma::mat33 A;

	A << 1 << 1 << 1 << arma::endr
	  << ksi_average[0] << ksi_average[1] << ksi_average[2] << arma::endr
	  << eta_average[0] << eta_average[2] << eta_average[3] << arma::endr;

	polynomial_coeff[0] = arma::solve(A, b);
	// ********************************** //

	// ******* p_1 coeff calculation ***** //
	A(1,1) = ksi_average[2];
	A(1,2) = ksi_average[3];

	A(2,1) = eta_average[2];
	A(2,2) = eta_average[3];

	polynomial_coeff[1] = arma::solve(A, b);

	// ************************************//

	// ******* p_2 coeff calculation ***** //
	A(1,1) = ksi_average[3];
	A(1,2) = ksi_average[1];

	A(2,1) = eta_average[3];
	A(2,2) = eta_average[1];

	polynomial_coeff[2] = arma::solve(A, b);

	// ************************************//

	// ******* p_3 coeff calculation ***** //
	A(1,1) = ksi_average[3];
	A(1,2) = ksi_average[6];

	A(2,1) = eta_average[3];
	A(2,2) = eta_average[6];

	polynomial_coeff[3] = arma::solve(A, b);

	// ************************************//

	// ******* p_4 coeff calculation ***** //
	A(1,1) = ksi_average[3];
	A(1,2) = ksi_average[7];

	A(2,1) = eta_average[3];
	A(2,2) = eta_average[7];

	polynomial_coeff[4] = arma::solve(A, b);

	// ************************************//

	// ******* p_5 coeff calculation ***** //
	A(1,1) = ksi_average[1];
	A(1,2) = ksi_average[4];

	A(2,1) = eta_average[1];
	A(2,2) = eta_average[4];

	polynomial_coeff[5] = arma::solve(A, b);

	// ************************************//

	// ******* p_6 coeff calculation ***** //
	A(1,1) = ksi_average[1];
	A(1,2) = ksi_average[5];

	A(2,1) = eta_average[1];
	A(2,2) = eta_average[5];

	polynomial_coeff[6] = arma::solve(A, b);

	// ************************************//


	// ******* p_7 coeff calculation ***** //
	A(1,1) = ksi_average[2];
	A(1,2) = ksi_average[8];

	A(2,1) = eta_average[2];
	A(2,2) = eta_average[8];

	polynomial_coeff[7] = arma::solve(A, b);

	// ************************************//


	// ******* p_8 coeff calculation ***** //
	A(1,1) = ksi_average[2];
	A(1,2) = ksi_average[9];

	A(2,1) = eta_average[2];
	A(2,2) = eta_average[9];

	polynomial_coeff[8] = arma::solve(A, b);

	// ************************************//

	// ALL POLYNOMIALS CALCULATED !!!!

	// TIME TO CALCULATE GAMMAS
	/* gleb_ebuchaya_zadrotina_ne_neset_mne_sok_dyadya_misha_primite_meru_protiv_etogo_pidorasa */

	arma::mat::fixed<4, 9> M;
	double p_ksi_square[9];
	double p_eta_square[9];
	double p_ksi_eta[9];

	//forming second row
	p_ksi_square[0] = polynomial_coeff[0][0] * ksi_square_average[0] + polynomial_coeff[0][1] * ksi_square_average[1] +
			polynomial_coeff[0][2] * ksi_square_average[2];
	p_ksi_square[1] = polynomial_coeff[1][0] * ksi_square_average[0] + polynomial_coeff[1][1] * ksi_square_average[2] +
					  polynomial_coeff[1][2] * ksi_square_average[3];
	p_ksi_square[2] = polynomial_coeff[2][0] * ksi_square_average[0] + polynomial_coeff[2][1] * ksi_square_average[3] +
					  polynomial_coeff[2][2] * ksi_square_average[1];
	p_ksi_square[3] = polynomial_coeff[3][0] * ksi_square_average[0] + polynomial_coeff[3][1] * ksi_square_average[3] +
					  polynomial_coeff[3][2] * ksi_square_average[6];
	p_ksi_square[4] = polynomial_coeff[4][0] * ksi_square_average[0] + polynomial_coeff[4][1] * ksi_square_average[3] +
					  polynomial_coeff[4][2] * ksi_square_average[7];
	p_ksi_square[5] = polynomial_coeff[5][0] * ksi_square_average[0] + polynomial_coeff[5][1] * ksi_square_average[1] +
					  polynomial_coeff[5][2] * ksi_square_average[4];
	p_ksi_square[6] = polynomial_coeff[6][0] * ksi_square_average[0] + polynomial_coeff[6][1] * ksi_square_average[1] +
					  polynomial_coeff[6][2] * ksi_square_average[5];
	p_ksi_square[7] = polynomial_coeff[7][0] * ksi_square_average[0] + polynomial_coeff[7][1] * ksi_square_average[2] +
					  polynomial_coeff[7][2] * ksi_square_average[8];
	p_ksi_square[8] = polynomial_coeff[8][0] * ksi_square_average[0] + polynomial_coeff[8][1] * ksi_square_average[2] +
					  polynomial_coeff[8][2] * ksi_square_average[9];

	//forming third row
	p_eta_square[0] = polynomial_coeff[0][0] * eta_square_average[0] + polynomial_coeff[0][1] * eta_square_average[1] +
					  polynomial_coeff[0][2] * eta_square_average[2];
	p_eta_square[1] = polynomial_coeff[1][0] * eta_square_average[0] + polynomial_coeff[1][1] * eta_square_average[2] +
					  polynomial_coeff[1][2] * eta_square_average[3];
	p_eta_square[2] = polynomial_coeff[2][0] * eta_square_average[0] + polynomial_coeff[2][1] * eta_square_average[3] +
					  polynomial_coeff[2][2] * eta_square_average[1];
	p_eta_square[3] = polynomial_coeff[3][0] * eta_square_average[0] + polynomial_coeff[3][1] * eta_square_average[3] +
					  polynomial_coeff[3][2] * eta_square_average[6];
	p_eta_square[4] = polynomial_coeff[4][0] * eta_square_average[0] + polynomial_coeff[4][1] * eta_square_average[3] +
					  polynomial_coeff[4][2] * eta_square_average[7];
	p_eta_square[5] = polynomial_coeff[5][0] * eta_square_average[0] + polynomial_coeff[5][1] * eta_square_average[1] +
					  polynomial_coeff[5][2] * eta_square_average[4];
	p_eta_square[6] = polynomial_coeff[6][0] * eta_square_average[0] + polynomial_coeff[6][1] * eta_square_average[1] +
					  polynomial_coeff[6][2] * eta_square_average[5];
	p_eta_square[7] = polynomial_coeff[7][0] * eta_square_average[0] + polynomial_coeff[7][1] * eta_square_average[2] +
					  polynomial_coeff[7][2] * eta_square_average[8];
	p_eta_square[8] = polynomial_coeff[8][0] * eta_square_average[0] + polynomial_coeff[8][1] * eta_square_average[2] +
					  polynomial_coeff[8][2] * eta_square_average[9];

	//forming fourth row
	p_ksi_eta[0] = polynomial_coeff[0][0] * ksi_eta_average[0] + polynomial_coeff[0][1] * ksi_eta_average[1] +
					  polynomial_coeff[0][2] * ksi_eta_average[2];
	p_ksi_eta[1] = polynomial_coeff[1][0] * ksi_eta_average[0] + polynomial_coeff[1][1] * ksi_eta_average[2] +
					  polynomial_coeff[1][2] * ksi_eta_average[3];
	p_ksi_eta[2] = polynomial_coeff[2][0] * ksi_eta_average[0] + polynomial_coeff[2][1] * ksi_eta_average[3] +
					  polynomial_coeff[2][2] * ksi_eta_average[1];
	p_ksi_eta[3] = polynomial_coeff[3][0] * ksi_eta_average[0] + polynomial_coeff[3][1] * ksi_eta_average[3] +
					  polynomial_coeff[3][2] * ksi_eta_average[6];
	p_ksi_eta[4] = polynomial_coeff[4][0] * ksi_eta_average[0] + polynomial_coeff[4][1] * ksi_eta_average[3] +
					  polynomial_coeff[4][2] * ksi_eta_average[7];
	p_ksi_eta[5] = polynomial_coeff[5][0] * ksi_eta_average[0] + polynomial_coeff[5][1] * ksi_eta_average[1] +
					  polynomial_coeff[5][2] * ksi_eta_average[4];
	p_ksi_eta[6] = polynomial_coeff[6][0] * ksi_eta_average[0] + polynomial_coeff[6][1] * ksi_eta_average[1] +
					  polynomial_coeff[6][2] * ksi_eta_average[5];
	p_ksi_eta[7] = polynomial_coeff[7][0] * ksi_eta_average[0] + polynomial_coeff[7][1] * ksi_eta_average[2] +
					  polynomial_coeff[7][2] * ksi_eta_average[8];
	p_ksi_eta[8] = polynomial_coeff[8][0] * ksi_eta_average[0] + polynomial_coeff[8][1] * ksi_eta_average[2] +
					  polynomial_coeff[8][2] * ksi_eta_average[9];


	for(int i = 0; i < 9; ++i)
		M(0, i) = 1;
	for(int i = 0; i < 9; ++i)
		M(1, i) = p_ksi_square[i];
	for(int i = 0; i < 9; ++i)
		M(2, i) = p_eta_square[i];
	for(int i = 0; i < 9; ++i)
		M(3, i) = p_ksi_eta[i];

	arma::vec4 c;
	c << 1 << sqr(x_g - x_0) / sqr(h) << sqr(y_g - y_0) / sqr(h) << (x_g - x_0) * (y_g - y_0) / sqr(h);

	arma::vec9 gamma = arma::solve(M, c);

	// GAMMAS CALCULATED!

	Vec4 q_reconstructed(0,0,0,0);

	Vec4 p[9];
	p[0] = polynomial_coeff[0][0] * q[0] + polynomial_coeff[0][1] * q[1] +
					  polynomial_coeff[0][2] * q[2];
	p[1] = polynomial_coeff[1][0] * q[0] + polynomial_coeff[1][1] * q[2] +
					  polynomial_coeff[1][2] * q[3];
	p[2] = polynomial_coeff[2][0] * q[0] + polynomial_coeff[2][1] * q[3] +
					  polynomial_coeff[2][2] * q[1];
	p[3] = polynomial_coeff[3][0] * q[0] + polynomial_coeff[3][1] * q[3] +
					  polynomial_coeff[3][2] * q[6];
	p[4] = polynomial_coeff[4][0] * q[0] + polynomial_coeff[4][1] * q[3] +
					  polynomial_coeff[4][2] * q[7];
	p[5] = polynomial_coeff[5][0] * q[0] + polynomial_coeff[5][1] * q[1] +
					  polynomial_coeff[5][2] * q[4];
	p[6] = polynomial_coeff[6][0] * q[0] + polynomial_coeff[6][1] * q[1] +
					  polynomial_coeff[6][2] * q[5];
	p[7] = polynomial_coeff[7][0] * q[0] + polynomial_coeff[7][1] * q[2] +
					  polynomial_coeff[7][2] * q[8];
	p[8] = polynomial_coeff[8][0] * q[0] + polynomial_coeff[8][1] * q[2] +
					  polynomial_coeff[8][2] * q[9];

	for(int i = 0; i < 9; ++i)
		q_reconstructed += gamma[i] * p[i];


	return q_reconstructed;



}