#include "WENOLF.h"


using namespace euler;



void WENOLF::CreateBoundingMesh()
{

	for(int triangle_counter = 0; triangle_counter < m_triangles.size(); ++triangle_counter)
	{

		for(int edge_number = 0; edge_number < 3; ++edge_number)
		{
			if(m_triangles[triangle_counter]->GetOppTriangle(edge_number) == nullptr)
			{
				auto const reflectedTriangle0 = m_triangles[triangle_counter]->ReflectTriangle(edge_number);
//				m_triangles[triangle_counter]->SetOppTriangle(edge_number, reflectedTriangle0);
				m_boundingTriangles.push_back(reflectedTriangle0);


				auto const p1 = m_triangles[triangle_counter]->getCorner((edge_number + 1) % 3);
				auto const p2 = m_triangles[triangle_counter]->getCorner((edge_number + 2) % 3);
				auto const ind_1 = reflectedTriangle0->getIntraTriangleIndex(p1);
				auto const ind_2 = reflectedTriangle0->getIntraTriangleIndex(p2);


				auto const reflectedTriangle1 = reflectedTriangle0->ReflectTriangle(ind_1);
//				reflectedTriangle0->SetOppTriangle(ind_1, reflectedTriangle1);
				auto const reflectedTriangle2 = reflectedTriangle0->ReflectTriangle(ind_2);
//				reflectedTriangle0->SetOppTriangle(ind_2, reflectedTriangle2);
				m_boundingTriangles.push_back(reflectedTriangle1);
				m_boundingTriangles.push_back(reflectedTriangle2);


/*				auto const p1_1 = reflectedTriangle0->getCorner((ind_1 + 1) % 3);
				auto const p2_1 = reflectedTriangle0->getCorner((ind_1 + 2) % 3);
				auto const ind_1_1 = reflectedTriangle1->getIntraTriangleIndex(p1_1);
				auto const ind_2_1 = reflectedTriangle1->getIntraTriangleIndex(p2_1);
				auto const reflectedTriangle3 = reflectedTriangle1->ReflectTriangle(ind_1_1);
				auto const reflectedTriangle4 = reflectedTriangle1->ReflectTriangle(ind_2_1);
				m_boundingTriangles.push_back(reflectedTriangle3);
				m_boundingTriangles.push_back(reflectedTriangle4);

				auto const p1_2 = reflectedTriangle0->getCorner((ind_2 + 1) % 3);
				auto const p2_2 = reflectedTriangle0->getCorner((ind_2 + 2) % 3);
				auto const ind_1_2 = reflectedTriangle2->getIntraTriangleIndex(p1_2);
				auto const ind_2_2 = reflectedTriangle2->getIntraTriangleIndex(p2_2);
				auto const reflectedTriangle5 = reflectedTriangle2->ReflectTriangle(ind_1_2);
				auto const reflectedTriangle6 = reflectedTriangle2->ReflectTriangle(ind_2_2);
				m_boundingTriangles.push_back(reflectedTriangle5);
				m_boundingTriangles.push_back(reflectedTriangle6);


*/

			}
		}


	}

}


void WENOLF::GetStencil(Triangle const* pTriangle, std::array<Triangle const*, 10> & stencil) const
{

	//@todo boundary !!!
	stencil[0] = pTriangle;

	stencil[1] = stencil[0]->GetOppTriangle(0);
	stencil[2] = stencil[0]->GetOppTriangle(1);
	stencil[3] = stencil[0]->GetOppTriangle(2);


	stencil[4] = stencil[1]->GetOppTriangle(stencil[1]->getIntraTriangleIndex(stencil[0]->getCorner(2)));

	stencil[5] = stencil[1]->GetOppTriangle(stencil[1]->getIntraTriangleIndex(stencil[0]->getCorner(1)));


	stencil[6] = stencil[3]->GetOppTriangle(stencil[3]->getIntraTriangleIndex(stencil[0]->getCorner(1)));

	stencil[7] = stencil[3]->GetOppTriangle(stencil[3]->getIntraTriangleIndex(stencil[0]->getCorner(0)));

//	stencil[7] = stencil[4]->GetOppTriangle(stencil[4]->getIntraTriangleIndex(stencil[0]->getCorner(1)));

	stencil[8] = stencil[2]->GetOppTriangle(stencil[2]->getIntraTriangleIndex(stencil[0]->getCorner(0)));

	stencil[9] = stencil[2]->GetOppTriangle(stencil[2]->getIntraTriangleIndex(stencil[0]->getCorner(2)));


}


void WENOLF::GetPointReconstructionData(OnePointReconstructionData &data,
										std::array<Triangle const *, 10> const &stencil,
										Point2 const &gaussian_point) const
{
	data.gaussian_point.x = gaussian_point.x;
	data.gaussian_point.y = gaussian_point.y;



	auto const h = std::sqrt(stencil[0]->getArea2D());
	auto const x_0 = stencil[0]->getBarycenter().x();
	auto const y_0 = stencil[0]->getBarycenter().y();

	double ksi_average[10];
	double eta_average[10];

	double ksi_square_average[10];
	double eta_square_average[10];
	double ksi_eta_average[10];




	for(int i = 0; i < 10; ++i)
	{
		auto const area = stencil[i]->getArea2D();
		ksi_average[i] = 1 / area * GaussianIntegration([h, x_0](double x, double y)
														{
															return (x - x_0) / h;
														}, stencil[i]);

		eta_average[i] = 1 / area * GaussianIntegration([h, y_0](double x, double y)
														{
															return (y - y_0) / h;
														}, stencil[i]);

		ksi_square_average[i] = 1 / area * GaussianIntegration([h, x_0](double x, double y)
															   {
																   return sqr(x-x_0) / sqr(h);
															   }, stencil[i]);
		eta_square_average[i] = 1 / area * GaussianIntegration([h, y_0](double x, double y)
															   {
																   return sqr(y-y_0) / sqr(h);
															   }, stencil[i]);
		ksi_eta_average[i] = 1 / area * GaussianIntegration([h, x_0, y_0](double x, double y)
															{
																return (x - x_0) * (y - y_0)  / sqr(h);
															}, stencil[i]);

	}

	//FOReconstructionPolynomial polynomials[9];
	data.polynomial[0].stencil = {0, 1, 2};
	data.polynomial[1].stencil = {0, 2, 3};
	data.polynomial[2].stencil = {0, 3, 1};
	data.polynomial[3].stencil = {0, 3, 6};
	data.polynomial[4].stencil = {0, 3, 7};
	data.polynomial[5].stencil = {0, 1, 4};
	data.polynomial[6].stencil = {0, 1, 5};
	data.polynomial[7].stencil = {0, 2, 8};
	data.polynomial[8].stencil = {0, 2, 9};

	//Calculating first order polynomials
	arma::vec3 b;
	b << 1 << (gaussian_point.x - x_0) / h << (gaussian_point.y - y_0) / h;

	for(int i = 0; i < 9; ++i)
	{
		arma::mat33 A;

		A << 1 << 1 << 1 << arma::endr
		  << ksi_average[data.polynomial[i].stencil[0]]
		  << ksi_average[data.polynomial[i].stencil[1]] << ksi_average[data.polynomial[i].stencil[2]] << arma::endr
		  << eta_average[data.polynomial[i].stencil[0]]
		  << eta_average[data.polynomial[i].stencil[1]] << eta_average[data.polynomial[i].stencil[2]] << arma::endr;

		arma::vec3 coeffs = arma::solve(A, b);

		data.polynomial[i].coeff[0] = coeffs[0];
		data.polynomial[i].coeff[1] = coeffs[1];
		data.polynomial[i].coeff[2] = coeffs[2];
	}


	double p_ksi_square[9];
	double p_eta_square[9];
	double p_ksi_eta[9];

	arma::mat M(4, 9);
	M.fill(1.0);

	for(int i = 0; i < 9; ++i)
	{
		auto const c_0 = data.polynomial[i].coeff[0];
		auto const c_1 = data.polynomial[i].coeff[1];
		auto const c_2 = data.polynomial[i].coeff[2];

		auto const ind_0 = data.polynomial[i].stencil[0];
		auto const ind_1 = data.polynomial[i].stencil[1];
		auto const ind_2 = data.polynomial[i].stencil[2];

		p_ksi_square[i] = c_0 * ksi_square_average[ind_0] + c_1 * ksi_square_average[ind_1]
						  + c_2 * ksi_square_average[ind_2];
		p_eta_square[i] = c_0 * eta_square_average[ind_0] + c_1 * eta_square_average[ind_1]
						  + c_2 * eta_square_average[ind_2];
		p_ksi_eta[i] = c_0 * ksi_eta_average[ind_0] + c_1 * ksi_eta_average[ind_1]
						  + c_2 * ksi_eta_average[ind_2];

		M(1, i) = p_ksi_square[i];
		M(2, i) = p_eta_square[i];
		M(3, i) = p_ksi_eta[i];

	}

	arma::vec4 c;
	c << 1 << sqr(gaussian_point.x - x_0) / sqr(h) << sqr(gaussian_point.y - y_0) / sqr(h)
	  << (gaussian_point.x - x_0) * (gaussian_point.y - y_0) / sqr(h);

	arma::vec gammas = solve(M, c);

	for(int i = 0; i < 9; ++i)
	{
		data.third_order_coeff[i] = gammas[i];
	}

}

void WENOLF::Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const &initStateFunction)
{
	Solver::Init(initStateFunction);

	m_vReconstructionData.resize(m_triangles.size());


	for(int triangle_number = 0; triangle_number < m_triangles.size(); ++triangle_number)
	{

		static auto const c = 1 / 2 + sqrt(3) / 6;
		//static auto const c = 1/2;

		auto const first_vertex = m_triangles[triangle_number]->getCorner(0);
		auto const second_vertex = m_triangles[triangle_number]->getCorner(1);
		auto const third_vertex = m_triangles[triangle_number]->getCorner(2);

		std::vector<Point2> gaussian_points(gaussian_points_number);

		gaussian_points[0].x = c * first_vertex->x() + (1 - c) * second_vertex->x();
		gaussian_points[0].y = c * first_vertex->y() + (1 - c) * second_vertex->y();

		gaussian_points[1].x = c * second_vertex->x() + (1 - c) * first_vertex->x();
		gaussian_points[1].y = c * second_vertex->y() + (1 - c) * first_vertex->y();

		gaussian_points[2].x = c * second_vertex->x() + (1 - c) * third_vertex->x();
		gaussian_points[2].y = c * second_vertex->y() + (1 - c) * third_vertex->y();

		gaussian_points[3].x = c * third_vertex->x() + (1 - c) * second_vertex->x();
		gaussian_points[3].y = c * third_vertex->y() + (1 - c) * second_vertex->y();

		gaussian_points[4].x = c * third_vertex->x() + (1 - c) * first_vertex->x();
		gaussian_points[4].y = c * third_vertex->y() + (1 - c) * first_vertex->y();

		gaussian_points[5].x = c * first_vertex->x() + (1 - c) * third_vertex->x();
		gaussian_points[5].y = c * first_vertex->y() + (1 - c) * third_vertex->y();



		std::array<Triangle const*, 10> stencil;
		GetStencil(m_triangles[triangle_number], stencil);

		for(int i = 0; i < gaussian_points_number; ++i)
		{
			GetPointReconstructionData(m_vReconstructionData[triangle_number][i],
									   stencil, gaussian_points[i]);
		}


	}


}



Vec4 WENOLF::Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, Point2 const& gaussianPoint) const
{

	std::array<Triangle const*, 10> stencil;
	GetStencil(pTriangle, stencil);

	for(int i = 0; i < 10; ++i)
	{
		if(stencil[i] == nullptr)
			return qVec;
	}

	std::array<Vec4, 10> q;
	q[0] = qVec;
	for(int i = 1; i < 10; ++i)
	{
		FormQVector(q[i], stencil[i]);
	}

	//Reconstruction!
	Vec4 q_reconstructed(0.0, 0.0, 0.0, 0.0);

	auto const triangleReconstructionData = m_vReconstructionData[pTriangle->Index()];

	OnePointReconstructionData pointReconstructionData;

	//searching for corresponding reconstruction data
	for(int j = 0; j < gaussian_points_number; ++j)
	{
		if (triangleReconstructionData[j].gaussian_point == gaussianPoint)
		{
			pointReconstructionData = triangleReconstructionData[j];
			break;
		}
	}

	std::array<Vec4, 9> omega, omega_waved;
	//@todo calculate omega waved here

	auto const h = std::sqrt(pTriangle->getArea2D());
	auto const x_0 = pTriangle->getBarycenter().x();
	auto const y_0 = pTriangle->getBarycenter().y();

	//static double const eps = 1e-4;

	static Vec4 const eps(1e-3, 1e-3, 1e-3, 1e-3);
	for(int i = 0; i < 9; ++i)
	{

		Vec4 smoothIndicator(0.0, 0.0, 0.0, 0.0);

		auto const firstPointData = triangleReconstructionData[0];
		auto const secondPointData = triangleReconstructionData[1];
		auto const thirdPointData = triangleReconstructionData[2];


		auto const ksi_1 = (firstPointData.gaussian_point.x - x_0) / h;
		auto const ksi_2 = (secondPointData.gaussian_point.x - x_0) / h;
		auto const ksi_3 = (thirdPointData.gaussian_point.x - x_0) / h;

		auto const eta_1 = (firstPointData.gaussian_point.y - y_0) / h;
		auto const eta_2 = (secondPointData.gaussian_point.y - y_0) / h;
		auto const eta_3 = (thirdPointData.gaussian_point.y - y_0) / h;


		Vec4 const kappa_1 = firstPointData.polynomial[i].coeff[0] * q[firstPointData.polynomial[i].stencil[0]]
							 + firstPointData.polynomial[i].coeff[1] * q[firstPointData.polynomial[i].stencil[1]]
							 + firstPointData.polynomial[i].coeff[2] * q[firstPointData.polynomial[i].stencil[2]];

		Vec4 const kappa_2 = secondPointData.polynomial[i].coeff[0] * q[secondPointData.polynomial[i].stencil[0]]
							 + secondPointData.polynomial[i].coeff[1] * q[secondPointData.polynomial[i].stencil[1]]
							 + secondPointData.polynomial[i].coeff[2] * q[secondPointData.polynomial[i].stencil[2]];

		Vec4 const kappa_3 = thirdPointData.polynomial[i].coeff[0] * q[thirdPointData.polynomial[i].stencil[0]]
							 + thirdPointData.polynomial[i].coeff[1] * q[thirdPointData.polynomial[i].stencil[1]]
							 + thirdPointData.polynomial[i].coeff[2] * q[thirdPointData.polynomial[i].stencil[2]];


		Vec4 a, b;

		for(int m = 0; m < 4; ++m)
		{
			arma::mat22 A;
			/*A << ksi_1 - ksi_3 << eta_1 - eta_3 << arma::endr
			  << ksi_2 - ksi_3 << eta_2 - eta_3 << arma::endr; */
			A(0, 0) = ksi_1 - ksi_3;
			A(0, 1) = eta_1 - eta_3;
			A(1, 0) = ksi_2 - ksi_3;
			A(1, 1) = eta_2 - eta_3;

			arma::vec2 c;
			c << kappa_1[m] - kappa_3[m] << kappa_2[m] - kappa_3[m];
			arma::vec solution = arma::solve(A, c);
			a[m] = solution[0];
			b[m] = solution[1];
		}

/*		Vec4 const a = 1 / ((eta_2 - eta_3) * (ksi_1 - ksi_3) - (eta_1 - eta_3) * (ksi_2 - ksi_3)) *
				((kappa_1 - kappa_3) * (eta_2 - eta_3) - (kappa_2 - kappa_3) * (eta_1 - eta_3));


		Vec4 const b = 1 / (eta_2 - eta_3) * ( (kappa_2 - kappa_3) - (ksi_2 - ksi_3) * a) ; */


/*
		for(int k = 0; k < 10; ++k)
		{
			smoothIndicator += 1 / sqr(h) * ( sqr(a) + sqr(b) ) * stencil[k]->getArea2D();
		}*/

		smoothIndicator = sqr(a) + sqr(b);


		omega_waved[i] = Vec4(pointReconstructionData.third_order_coeff[i], pointReconstructionData.third_order_coeff[i],
							  pointReconstructionData.third_order_coeff[i], pointReconstructionData.third_order_coeff[i]);

		omega_waved[i] = omega_waved[i] / (eps + smoothIndicator);



	}



	Vec4 o_wave_sum(0.0, 0.0, 0.0, 0.0);
	for(int i = 0; i < 9; ++i)
		o_wave_sum += omega_waved[i];

	for(int i = 0; i < 9; ++i)
	{
		auto const c_0 = pointReconstructionData.polynomial[i].coeff[0];
		auto const c_1= pointReconstructionData.polynomial[i].coeff[1];
		auto const c_2 = pointReconstructionData.polynomial[i].coeff[2];

		auto const ind_0 = pointReconstructionData.polynomial[i].stencil[0];
		auto const ind_1 = pointReconstructionData.polynomial[i].stencil[1];
		auto const ind_2 = pointReconstructionData.polynomial[i].stencil[2];

		omega[i] = omega_waved[i] / o_wave_sum;



		q_reconstructed += omega[i] *
				(c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2]);

	}


	return q_reconstructed;


}


void WENOLF::FormL_A(Vec4 const &qVec, arma::mat44 &L_A) const
{
	double density, velocityX, velocityY, pressure;
	GetGasParamsFromQ(qVec, density, velocityX, velocityY, pressure);

	auto const c = std::sqrt(m_gamma * pressure / density);
	auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);

	L_A(0, 0) = velocity_sqr_abs / 2 + velocityX * c / (m_gamma - 1);
	L_A(0, 1) = -velocityX - c / (m_gamma - 1);
	L_A(0, 2) = -velocityY;
	L_A(0, 3) = 1;

	L_A(1, 0) = c * c / (m_gamma - 1) - velocity_sqr_abs / 2 - velocityY * c /(m_gamma - 1);
	L_A(1, 1) = velocityX;
	L_A(1, 2) = velocityY + c / (m_gamma - 1);
	L_A(1, 3) = -1;

	L_A(2, 0) = c * c / (m_gamma - 1) - velocity_sqr_abs / 2 + velocityY * c /(m_gamma - 1);
	L_A(2, 1) = velocityX;
	L_A(2, 2) = velocityY - c / (m_gamma - 1);
	L_A(2, 3) = -1;

	L_A(3, 0) = velocity_sqr_abs / 2 - velocityX * c / (m_gamma - 1);
	L_A(3, 1) = -velocityX + c / (m_gamma - 1);
	L_A(3, 2) = -velocityY;
	L_A(3, 3) = 1;

	L_A *= (m_gamma - 1) / (2 * c * c);

}

void WENOLF::FormL_B(Vec4 const &qVec, arma::mat44 &L_B) const
{
	double density, velocityX, velocityY, pressure;
	GetGasParamsFromQ(qVec, density, velocityX, velocityY, pressure);

	auto const c = std::sqrt(m_gamma * pressure / density);
	auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);

	L_B(0, 0) = velocity_sqr_abs / 2 + velocityY * c / (m_gamma - 1);
	L_B(0, 1) = -velocityX;
	L_B(0, 2) = -velocityY - c / (m_gamma - 1);
	L_B(0, 3) = 1;

	L_B(1, 0) = c * c / (m_gamma - 1) - velocity_sqr_abs / 2 - velocityX * c / (m_gamma - 1);
	L_B(1, 1) = velocityX + c / (m_gamma - 1);
	L_B(1, 2) = velocityY;
	L_B(1, 3) = -1;

	L_B(2, 0) = c * c / (m_gamma - 1) - velocity_sqr_abs / 2 + velocityX * c / (m_gamma - 1);
	L_B(2, 1) = velocityX - c / (m_gamma - 1);
	L_B(2, 2) = velocityY;
	L_B(2, 3) = -1;

	L_B(3, 0) = velocity_sqr_abs / 2 - velocityY * c / (m_gamma - 1);
	L_B(3, 1) = -velocityX;
	L_B(3, 2) = -velocityY + c / (m_gamma - 1);
	L_B(3, 3) = 1;

	L_B *= (m_gamma - 1) / (2 * c * c);

}

void WENOLF::FormR_A(Vec4 const &qVec, arma::mat44 &R_A) const
{
	double density, velocityX, velocityY, pressure;
	GetGasParamsFromQ(qVec, density, velocityX, velocityY, pressure);

	auto const c = std::sqrt(m_gamma * pressure / density);
	auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);

	auto const eps = pressure / (density * (m_gamma - 1));

	auto const H = eps + pressure / density + velocity_sqr_abs / 2;

	R_A(0, 0) = 1;
	R_A(0, 1) = 1;
	R_A(0, 2) = 1;
	R_A(0, 3) = 1;

	R_A(1, 0) = velocityX - c;
	R_A(1, 1) = velocityX;
	R_A(1, 2) = velocityX;
	R_A(1, 3) = velocityX + c;

	R_A(2, 0) = velocityY;
	R_A(2, 1) = velocityY + c;
	R_A(2, 2) = velocityY - c;
	R_A(2, 3) = velocityY;

	R_A(3, 0) = H - velocityX * c;
	R_A(3, 1) = velocity_sqr_abs / 2 + velocityY * c;
	R_A(3, 2) = velocity_sqr_abs / 2 - velocityY * c;
	R_A(3, 3) = H + velocityX * c;


}

void WENOLF::FormR_B(Vec4 const &qVec, arma::mat44 &R_B) const
{
	double density, velocityX, velocityY, pressure;
	GetGasParamsFromQ(qVec, density, velocityX, velocityY, pressure);

	auto const c = std::sqrt(m_gamma * pressure / density);
	auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);

	auto const eps = pressure / (density * (m_gamma - 1));

	auto const H = eps + pressure / density + velocity_sqr_abs / 2;

	R_B(0, 0) = 1;
	R_B(0, 1) = 1;
	R_B(0, 2) = 1;
	R_B(0, 3) = 1;

	R_B(1, 0) = velocityX;
	R_B(1, 1) = velocityX + c;
	R_B(1, 2) = velocityX - c;
	R_B(1, 3) = velocityX;

	R_B(2, 0) = velocityY - c;
	R_B(2, 1) = velocityY;
	R_B(2, 2) = velocityY;
	R_B(2, 3) = velocityY + c;

	R_B(3, 0) = H - velocityY * c;
	R_B(3, 1) = velocity_sqr_abs / 2 + velocityX * c;
	R_B(3, 2) = velocity_sqr_abs / 2 - velocityX * c;
	R_B(3, 3) = H + velocityY * c;

}

