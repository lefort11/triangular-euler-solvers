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



void WENOLF::Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const &initStateFunction)
{
	Solver::Init(initStateFunction);

	m_vReconstructionData.resize(m_triangles.size());

	for(int triangle_number = 0; triangle_number < m_triangles.size(); ++triangle_number)
	{
		GetTriangleReconstructionData(m_vReconstructionData[triangle_number], m_triangles[triangle_number]);

	}
}

void WENOLF::GetTriangleReconstructionData(TriangleReconstructionData &trRecData, Triangle const *pTriangle)
{

	trRecData.fo_polynomial[0].stencil = {0, 1, 2};
	trRecData.fo_polynomial[1].stencil = {0, 2, 3};
	trRecData.fo_polynomial[2].stencil = {0, 3, 1};
	trRecData.fo_polynomial[3].stencil = {0, 3, 6};
	trRecData.fo_polynomial[4].stencil = {0, 3, 7};
	trRecData.fo_polynomial[5].stencil = {0, 1, 4};
	trRecData.fo_polynomial[6].stencil = {0, 1, 5};
	trRecData.fo_polynomial[7].stencil = {0, 2, 8};
	trRecData.fo_polynomial[8].stencil = {0, 2, 9};

	static auto const c = 1 / 2 + sqrt(3) / 6;

	auto const first_vertex = pTriangle->getCorner(0);
	auto const second_vertex = pTriangle->getCorner(1);
	auto const third_vertex = pTriangle->getCorner(2);

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

	GetStencil(pTriangle, stencil);

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

	for(int g_point_number = 0; g_point_number < gaussian_points_number; ++g_point_number)
	{
		auto& currGPoint = gaussian_points[g_point_number];
		trRecData.gaussian_points[g_point_number] = currGPoint;

		arma::vec3 b;
		b << 1.0 << (currGPoint.x - x_0) / h << (currGPoint.y - y_0) / h;

		arma::vec4 d;
		d << 1.0 << sqr( (currGPoint.x - x_0) / h ) << sqr( (currGPoint.y - y_0) / h )
		  << (currGPoint.x - x_0) * (currGPoint.y - y_0) / sqr(h);

		arma::mat M(4, 9);
		M.fill(1.0);



		for(int polynomial_number = 0; polynomial_number < 9; ++polynomial_number)
		{
			auto& currPolynomial = trRecData.fo_polynomial[polynomial_number];

			auto const ind_0 = currPolynomial.stencil[0];
			auto const ind_1 = currPolynomial.stencil[1];
			auto const ind_2 = currPolynomial.stencil[2];


			arma::mat33 A;
			A << 1.0 << 1.0 << 1.0 << arma::endr
			  << ksi_average[ind_0] << ksi_average[ind_1]
			 			 << ksi_average[ind_2] << arma::endr
			  << eta_average[ind_0] << eta_average[ind_1]
			 			 << eta_average[ind_2] << arma::endr;

			arma::vec3 coeffs = arma::solve(A, b);

			currPolynomial.coeffsAtPoints[g_point_number].c[0] = coeffs[0];
			currPolynomial.coeffsAtPoints[g_point_number].c[1] = coeffs[1];
			currPolynomial.coeffsAtPoints[g_point_number].c[2] = coeffs[2];

			M(1, polynomial_number) = coeffs[0] * ksi_square_average[ind_0] + coeffs[1] * ksi_square_average[ind_1]
									  + coeffs[2] * ksi_square_average[ind_2];

			M(2, polynomial_number) = coeffs[0] * eta_square_average[ind_0] + coeffs[1] * eta_square_average[ind_1]
									  + coeffs[2] * eta_square_average[ind_2];
			M(3, polynomial_number) = coeffs[0] * ksi_eta_average[ind_0] + coeffs[1] * ksi_eta_average[ind_1]
									  + coeffs[2] * ksi_eta_average[ind_2];


		}

		arma::vec9 gammas = arma::solve(M, d);

		for(int i = 0; i < 9; ++i)
			trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas[i] = gammas[i];



	}


#if 0
	for(int g_point_number = 0; g_point_number < gaussian_points_number; ++g_point_number)
	{
		auto& currGPoint = gaussian_points[g_point_number];
		trRecData.gaussian_points[g_point_number] = currGPoint;

		arma::vec3 b;
		b << 1 << (currGPoint.x - x_0) / h << (currGPoint.y - y_0) / h;

		for(int polyn_number = 0; polyn_number < 9; ++polyn_number)
		{
			auto& currPolynomial = trRecData.fo_polynomial[polyn_number];


			arma::mat33 A;

			A << 1 << 1 << 1 << arma::endr
			  << ksi_average[currPolynomial.stencil[0]] << ksi_average[currPolynomial.stencil[1]]
			  								<< ksi_average[currPolynomial.stencil[2]] << arma::endr
			  << eta_average[currPolynomial.stencil[0]] << eta_average[currPolynomial.stencil[1]]
			  								<< eta_average[currPolynomial.stencil[2]] << arma::endr;


			arma::vec3 coeffs = arma::solve(A, b);

			currPolynomial.coeffsAtPoints[g_point_number].c[0] = coeffs[0];
			currPolynomial.coeffsAtPoints[g_point_number].c[1] = coeffs[1];
			currPolynomial.coeffsAtPoints[g_point_number].c[2] = coeffs[2];


		}

		double p_ksi_square[9];
		double p_eta_square[9];
		double p_ksi_eta[9];

		arma::mat M(4, 9);
		M.fill(1.0);

		for(int i = 0; i < 9; ++i)
		{
			auto const c_0 = trRecData.fo_polynomial[i].coeffsAtPoints[g_point_number].c[0];
			auto const c_1 = trRecData.fo_polynomial[i].coeffsAtPoints[g_point_number].c[1];;
			auto const c_2 = trRecData.fo_polynomial[i].coeffsAtPoints[g_point_number].c[2];;

			auto const ind_0 = trRecData.fo_polynomial[i].stencil[0];
			auto const ind_1 = trRecData.fo_polynomial[i].stencil[1];
			auto const ind_2 = trRecData.fo_polynomial[i].stencil[2];

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
		c << 1 << sqr(currGPoint.x - x_0) / sqr(h) << sqr(currGPoint.y - y_0) / sqr(h)
		  << (currGPoint.x - x_0) * (currGPoint.y - y_0) / sqr(h);

		arma::vec gammas = solve(M, c);

		for(int i = 0; i < 9; ++i)
		{
			trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas[i] = gammas[i];
		}
	}

#endif


	GetSmoothIndicatorData(trRecData, pTriangle);



}

void WENOLF::GetSmoothIndicatorData(TriangleReconstructionData &trRecData, Triangle const *pTriangle) const
{
	std::array<Triangle const*, 10> stencil;
	GetStencil(pTriangle, stencil);

	auto const h = pTriangle->getArea2D();
	auto const x_0 = pTriangle->getBarycenter().x();
	auto const y_0 = pTriangle->getBarycenter().y();

	double ksi_average[10];
	double eta_average[10];


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

	}


	for(int polynomial_number = 0; polynomial_number < 9; ++polynomial_number)
	{
		arma::mat99 A;
		arma::vec9 c;

		auto const ind_0 = trRecData.fo_polynomial[polynomial_number].stencil[0];
		auto const ind_1 = trRecData.fo_polynomial[polynomial_number].stencil[1];
		auto const ind_2 = trRecData.fo_polynomial[polynomial_number].stencil[2];

		int i = 0;
		for (int g_point_number = 0; g_point_number < 3; ++g_point_number)
		{
			auto const c_0 = trRecData.fo_polynomial[polynomial_number].coeffsAtPoints[g_point_number].c[0];
			auto const c_1 = trRecData.fo_polynomial[polynomial_number].coeffsAtPoints[g_point_number].c[1];
			auto const c_2 = trRecData.fo_polynomial[polynomial_number].coeffsAtPoints[g_point_number].c[2];

			auto const ksi = (trRecData.gaussian_points[g_point_number].x - x_0) / h;
			auto const eta = (trRecData.gaussian_points[g_point_number].y - y_0) / h;


			//first row, u = 1;
			A(i, 0) = ksi;
			A(i, 1) = ksi;
			A(i, 2) = ksi;

			A(i, 3) = eta;
			A(i, 4) = eta;
			A(i, 5) = eta;

			A(i, 6) = 1;
			A(i, 7) = 1;
			A(i, 8) = 1;

			c(i) = c_0 + c_1 + c_2;

			//second row, u = ksi
			A(i + 1, 0) = ksi_average[ind_0] * ksi;
			A(i + 1, 1) = ksi_average[ind_1] * ksi;
			A(i + 1, 2) = ksi_average[ind_2] * ksi;

			A(i + 1, 3) = ksi_average[ind_0] * eta;
			A(i + 1, 4) = ksi_average[ind_1] * eta;
			A(i + 1, 5) = ksi_average[ind_2] * eta;

			A(i + 1, 6) = ksi_average[ind_0];
			A(i + 1, 7) = ksi_average[ind_1];
			A(i + 1, 8) = ksi_average[ind_2];
			c(i + 1) = c_0 * ksi_average[ind_0] + c_1 * ksi_average[ind_1] + c_2 * ksi_average[ind_2];

			//third row, u = eta;
			A(i + 2, 0) = eta_average[ind_0] * ksi;
			A(i + 2, 1) = eta_average[ind_1] * ksi;
			A(i + 2, 2) = eta_average[ind_2] * ksi;

			A(i + 2, 3) = eta_average[ind_0] * eta;
			A(i + 2, 4) = eta_average[ind_1] * eta;
			A(i + 2, 5) = eta_average[ind_2] * eta;

			A(i + 2, 6) = eta_average[ind_0];
			A(i + 2, 7) = eta_average[ind_1];
			A(i + 2, 8) = eta_average[ind_2];
			c(i + 2) = c_0 * eta_average[ind_0] + c_1 * eta_average[ind_1] + c_2 * eta_average[ind_2];

			i+=3;
		}

		arma::vec9 solution = arma::solve(A, c);

		trRecData.smoothIndicatorData[polynomial_number].alpha[0] = solution[0];
		trRecData.smoothIndicatorData[polynomial_number].alpha[1] = solution[1];
		trRecData.smoothIndicatorData[polynomial_number].alpha[2] = solution[2];

		trRecData.smoothIndicatorData[polynomial_number].beta[0] = solution[3];
		trRecData.smoothIndicatorData[polynomial_number].beta[1] = solution[4];
		trRecData.smoothIndicatorData[polynomial_number].beta[2] = solution[5];



	}
}


Vec4 WENOLF::Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, Point2 const &gaussianPoint) const
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

/*	std::array<Vec4, 10> w_x;
	std::array<Vec4, 10> w_y;

	arma::mat44 L_A, R_A, L_B, R_B;
	Vec4 q_avrg = 0.1 * (q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7] + q[8] + q[9]);
	FormL_A(q_avrg, L_A);
	FormL_B(q_avrg, L_B);
	FormR_A(q_avrg, R_A);
	FormR_B(q_avrg, R_B);
	auto const n = CalculateNormal(pTriangle, edgeNumber);
	arma::mat44 L = std::fabs(n[0]) * L_A + std::fabs(n[1]) * L_B;
	arma::mat44 R = std::fabs(n[0]) * R_A + std::fabs(n[1]) * R_B;

	for(int i = 0; i < 10; ++i)
	{
		w_x[i] = L_A * q[i];
		w_y[i] = L_B * q[i];
	}


*/
	//Reconstruction!
	Vec4 q_reconstructed(0.0, 0.0, 0.0, 0.0);
	//Vec4 w_reconstructed(0.0, 0.0, 0.0, 0.0);

	auto const triangleReconstructionData = m_vReconstructionData[pTriangle->Index()];

	//Searching for coressponding gaussian point
	int curr_g_point_n = 0;
	for(curr_g_point_n = 0; curr_g_point_n < gaussian_points_number; ++curr_g_point_n)
	{
		if(triangleReconstructionData.gaussian_points[curr_g_point_n] == gaussianPoint)
			break;
	}


	std::array<Vec4, 9> omega, omega_waved;
	//@todo calculate omega waved here

	auto const h = std::sqrt(pTriangle->getArea2D());
	auto const x_0 = pTriangle->getBarycenter().x();
	auto const y_0 = pTriangle->getBarycenter().y();


	Vec4 smoothIndicator(0.0, 0.0, 0.0, 0.0);

	static Vec4 const eps(m_eps, m_eps, m_eps, m_eps);
	for(int i = 0; i < 9; ++i)
	{




		auto const ind_0 = triangleReconstructionData.fo_polynomial[i].stencil[0];
		auto const ind_1 = triangleReconstructionData.fo_polynomial[i].stencil[1];
		auto const ind_2= triangleReconstructionData.fo_polynomial[i].stencil[2];

		smoothIndicator = (sqr(triangleReconstructionData.smoothIndicatorData[i].alpha[0] * q[ind_0]
								+ triangleReconstructionData.smoothIndicatorData[i].alpha[1] * q[ind_1]
								+ triangleReconstructionData.smoothIndicatorData[i].alpha[2] * q[ind_2])
							+ sqr(triangleReconstructionData.smoothIndicatorData[i].beta[0] * q[ind_0]
								  + triangleReconstructionData.smoothIndicatorData[i].beta[1] * q[ind_1]
								  + triangleReconstructionData.smoothIndicatorData[i].beta[2] * q[ind_2]));


/*		for(int k = 0; k < 10; ++k)
		{
			smoothIndicator += (sqr(triangleReconstructionData.smoothIndicatorData[i].alpha[0] * q[ind_0]
								  + triangleReconstructionData.smoothIndicatorData[i].alpha[1] * q[ind_1]
								  + triangleReconstructionData.smoothIndicatorData[i].alpha[2] * q[ind_2])
							  + sqr(triangleReconstructionData.smoothIndicatorData[i].beta[0] * q[ind_0]
									+ triangleReconstructionData.smoothIndicatorData[i].beta[1] * q[ind_1]
									+ triangleReconstructionData.smoothIndicatorData[i].beta[2] * q[ind_2]))
							  * stencil[k]->getArea2D();

			smoothIndicator = 1 / sqr(h) * smoothIndicator;
		}
*/

		omega_waved[i] = Vec4(triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i],
							  triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i],
							  triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i],
							  triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i]);

		omega_waved[i] = omega_waved[i] / (eps + smoothIndicator);



	}



	Vec4 o_wave_sum(0.0, 0.0, 0.0, 0.0);
	for(int i = 0; i < 9; ++i)
		o_wave_sum += omega_waved[i];

	for(int i = 0; i < 9; ++i)
	{
		auto const c_0 = triangleReconstructionData.fo_polynomial[i].coeffsAtPoints[curr_g_point_n].c[0];
		auto const c_1= triangleReconstructionData.fo_polynomial[i].coeffsAtPoints[curr_g_point_n].c[1];
		auto const c_2 = triangleReconstructionData.fo_polynomial[i].coeffsAtPoints[curr_g_point_n].c[2];

		auto const ind_0 = triangleReconstructionData.fo_polynomial[i].stencil[0];
		auto const ind_1 = triangleReconstructionData.fo_polynomial[i].stencil[1];
		auto const ind_2 = triangleReconstructionData.fo_polynomial[i].stencil[2];

		omega[i] = omega_waved[i] / o_wave_sum;



		q_reconstructed += omega[i] *
						   (c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2]);

	}

	//Vec4 q_reconstructed(0.0, 0.0, 0.0, 0.0);
	//q_reconstructed = R * w_reconstructed;

	if(q_reconstructed(3) <= 0)
		throw 1;

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

