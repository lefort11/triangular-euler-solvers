#include "RoeSolver.h"


using namespace euler;

void RoeSolver::CreateBoundingMesh()
{

	for(int triangle_counter = 0; triangle_counter < m_triangles.size(); ++triangle_counter)
	{

		for(int edge_number = 0; edge_number < 3; ++edge_number)
		{
			if(m_triangles[triangle_counter]->GetOppTriangle(edge_number) == nullptr)
			{
				auto const reflectedTriangle = m_triangles[triangle_counter]->ReflectTriangle(edge_number);
				m_triangles[triangle_counter]->SetOppTriangle(edge_number, reflectedTriangle);
				m_boundingTriangles.push_back(reflectedTriangle);
			}
		}

	}

}


Vec4 RoeSolver::CalculateFlux(Vec4 const &qVec, int triangleNumber, int edgeNumber) const
{


//preparing q_minus and q_plus vectors

	//****** Calculating gaussian points ********//

	static auto const gaussian_weight = 1.0 / 2.0 + sqrt(3.0) / 6.0;
	//static auto const gaussian_weight = 1 / 2;

	Point2  gaussian_p_1, gaussian_p_2;
	auto const firstVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 1) % 3);
	auto const secondVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 2) % 3);

	gaussian_p_1.x = gaussian_weight * firstVertex->x() + (1 - gaussian_weight) * secondVertex->x();
	gaussian_p_1.y = gaussian_weight * firstVertex->y() + (1 - gaussian_weight) * secondVertex->y();

	gaussian_p_2.x = gaussian_weight * secondVertex->x() + (1 - gaussian_weight) * firstVertex->x();
	gaussian_p_2.y = gaussian_weight * secondVertex->y() + (1 - gaussian_weight) * firstVertex->y();

	Vec4 flux(0.0, 0.0, 0.0, 0.0);

	auto const normal = CalculateNormal(m_triangles[triangleNumber], edgeNumber);

	for(Point2 g_point: {gaussian_p_1, gaussian_p_2})
	{



		//*****************************************//

		//********** Forming q_minus vector and _minus parameters *******//
		auto const q_minus = Reconstruct(qVec, m_triangles[triangleNumber], g_point, edgeNumber);

//	double density_minus, velocityX_minus, velocityY_minus, pressure_minus;

//	GetGasParamsFromQ(q_minus, density_minus, velocityX_minus, velocityY_minus, pressure_minus);

//	auto const eps_minus = pressure_minus / (density_minus * (m_gamma - 1.0));

		auto const density_minus = q_minus[0];
		auto const velocityX_minus = q_minus[1] / density_minus;
		auto const velocityY_minus = q_minus[2] / density_minus;
		auto const E_minus = q_minus[3] / density_minus;
		auto const velocity_sqr_abs_minus = sqr(velocityX_minus) + sqr(velocityY_minus);
		auto const eps_minus = E_minus - 0.5 * velocity_sqr_abs_minus;
		auto const pressure_minus = (m_gamma - 1.0) * eps_minus * density_minus;
		auto const H_minus = eps_minus + pressure_minus / density_minus + 0.5 * velocity_sqr_abs_minus;
		auto const c_minus = std::sqrt(m_gamma * pressure_minus / density_minus);


		//****************************************************************//


		//********* Forming q_plus vector and _plus parameteres **********//


		Vec4 q_plus;
		Triangle *neighbour_triangle = m_triangles[triangleNumber]->GetOppTriangle(edgeNumber);

		auto density_plus = neighbour_triangle->density;
		auto velocityX_plus = neighbour_triangle->velocityX;
		auto velocityY_plus = neighbour_triangle->velocityY;
		auto pressure_plus = neighbour_triangle->pressure;

		FormQVector(q_plus, density_plus, velocityX_plus, velocityY_plus, pressure_plus);

		//edge_number calculation
		auto const v0_ind =
				neighbour_triangle->getIntraTriangleIndex(m_triangles[triangleNumber]->getCorner((edgeNumber + 1)%3));
		auto const v1_ind =
				neighbour_triangle->getIntraTriangleIndex(m_triangles[triangleNumber]->getCorner((edgeNumber + 2)%3));
		int neighbourEdgeNumber = (v0_ind + 1) % 3;
		if(neighbourEdgeNumber == v1_ind)
			neighbourEdgeNumber = (v1_ind + 1) % 3;


		q_plus = Reconstruct(q_plus, neighbour_triangle, g_point, neighbourEdgeNumber);


		//Updating plus values
		density_plus = q_plus[0];
		velocityX_plus = q_plus[1] / density_plus;
		velocityY_plus = q_plus[2] / density_plus;
		auto const E_plus = q_plus[3] / density_plus;
		auto const velocity_sqr_abs_plus = sqr(velocityX_plus) + sqr(velocityY_plus);
		auto const eps_plus = E_plus - 0.5 * velocity_sqr_abs_plus;
		pressure_plus = (m_gamma - 1.0) * eps_plus * density_plus;
		auto const H_plus = eps_plus + pressure_plus / density_plus + 0.5 * velocity_sqr_abs_plus;
		auto const c_plus = std::sqrt(m_gamma * pressure_plus / density_plus);

		//************************************************************************//



		//****************** Forming flux vectors ********************************//

		//Calculating F_minus
		Vec4 F_minus = {
				density_minus * velocityX_minus,
				density_minus * sqr(velocityX_minus) + pressure_minus,
				density_minus * velocityX_minus * velocityY_minus,
				density_minus * velocityX_minus * H_minus
		};

		//Calculating G_Minus
		Vec4 G_minus = {
				density_minus * velocityY_minus,
				density_minus * velocityY_minus * velocityX_minus,
				density_minus * sqr(velocityY_minus) + pressure_minus,
				density_minus * velocityY_minus * H_minus

		};

		//Calculating F_plus
		Vec4 F_plus = {
				density_plus * velocityX_plus,
				density_plus * sqr(velocityX_plus) + pressure_plus,
				density_plus * velocityX_plus * velocityY_plus,
				density_plus * velocityX_plus * H_plus
		};

		//Calculating G_plus
		Vec4 G_plus = {
				density_plus * velocityY_plus,
				density_plus * velocityY_plus * velocityX_plus,
				density_plus * sqr(velocityY_plus) + pressure_plus,
				density_plus * velocityY_plus * H_plus

		};



		auto const density_star = std::sqrt(density_minus * density_plus);

		auto const velocityX_star = (std::sqrt(density_minus) * velocityX_minus
									 + std::sqrt(density_plus) * velocityX_plus)
									/ (std::sqrt(density_minus) + std::sqrt(density_plus));

		auto const velocityY_star = (std::sqrt(density_minus) * velocityY_minus
									 + std::sqrt(density_plus) * velocityY_plus)
									/ (std::sqrt(density_minus) + std::sqrt(density_plus));

		auto const H_star = (std::sqrt(density_minus) * H_minus
							 + std::sqrt(density_plus) * H_plus)
							/ (std::sqrt(density_minus) + std::sqrt(density_plus));

		auto const velocity_sqr_abs_star = sqr(velocityX_star) + sqr(velocityY_star);

		auto const c_star = std::sqrt( (m_gamma - 1) * (H_star -  0.5 * velocity_sqr_abs_star) );

		std::array<double, 4> lambda_A = {
				velocityX_star - c_star,
				velocityX_star,
				velocityX_star,
				velocityX_star + c_star
		};

		std::array<double, 4> lambda_A_waved = {
				std::min(lambda_A[0], velocityX_minus - c_minus),
				lambda_A[1],
				lambda_A[2],
				std::max(lambda_A[3], velocityX_plus + c_plus)
		};

		std::array<double, 4> lambda_B = {
				velocityY_star - c_star,
				velocityY_star,
				velocityY_star,
				velocityY_star + c_star
		};

		std::array<double, 4> lambda_B_waved = {
				std::min(lambda_B[0], velocityY_minus - c_minus),
				lambda_B[1],
				lambda_B[2],
				std::max(lambda_B[3], velocityY_plus + c_plus)
		};

		std::array<double, 4> lambda = {
				sqrt(velocity_sqr_abs_star) - c_star,
				sqrt(velocity_sqr_abs_star),
				sqrt(velocity_sqr_abs_star),
				sqrt(velocity_sqr_abs_star) + c_star
		};

		std::array<double, 4> lambda_waved = {
				std::min(lambda[0], sqrt(velocity_sqr_abs_minus) - c_minus),
				lambda[1],
				lambda[2],
				std::max(lambda[3], sqrt(velocity_sqr_abs_plus) + c_plus)
		};


		std::array<double, 4> delta_s_A = {
				((pressure_plus - pressure_minus) - density_star * c_star * (velocityX_plus - velocityX_minus))
				/ (2 * sqr(c_star)),

				(sqr(c_star) * (density_plus - density_minus) - (pressure_plus - pressure_minus)
				 + density_star * c_star * (velocityY_plus - velocityY_minus)) / (2 * sqr(c_star)),

				(sqr(c_star) * (density_plus - density_minus) - (pressure_plus - pressure_minus)
				 - density_star * c_star * (velocityY_plus - velocityY_minus)) / (2 * sqr(c_star)),

				((pressure_plus - pressure_minus) + density_star * c_star * (velocityX_plus - velocityX_minus))
				/ (2 * sqr(c_star))


		};

		std::array<double, 4> delta_s_B = {

				((pressure_plus - pressure_minus) - density_star * c_star * (velocityY_plus - velocityY_minus))
				/ (2 * sqr(c_star)),

				(sqr(c_star) * (density_plus - density_minus) - (pressure_plus - pressure_minus)
				 + density_star * c_star * (velocityX_plus - velocityX_minus)) / (2 * sqr(c_star)),

				(sqr(c_star) * (density_plus - density_minus) - (pressure_plus - pressure_minus)
				 - density_star * c_star * (velocityX_plus - velocityX_minus)) / (2 * sqr(c_star)),

				((pressure_plus - pressure_minus) + density_star * c_star * (velocityY_plus - velocityY_minus))
				/ (2 * sqr(c_star))

		};

		std::array<Vec4, 4> r_A = {
				Vec4(1.0, velocityX_star - c_star, velocityY_star, H_star - velocityX_star * c_star),
				Vec4(1.0, velocityX_star, velocityY_star + c_star, 0.5 * velocity_sqr_abs_star + velocityY_star * c_star),
				Vec4(1.0, velocityX_star, velocityY_star - c_star, 0.5 * velocity_sqr_abs_star - velocityY_star * c_star),
				Vec4(1.0, velocityX_star + c_star, velocityY_star, H_star + velocityX_star * c_star)

		};

		std::array<Vec4, 4> r_B = {
				Vec4(1.0, velocityX_star, velocityY_star - c_star, H_star - velocityY_star * c_star),
				Vec4(1.0, velocityX_star + c_star, velocityY_star, 0.5 * velocity_sqr_abs_star + velocityX_star * c_star),
				Vec4(1.0, velocityX_star - c_star, velocityY_star, 0.5 * velocity_sqr_abs_star - velocityX_star * c_star),
				Vec4(1.0, velocityX_star, velocityY_star + c_star, H_star + velocityY_star * c_star)

		};


		Vec4 F(0.0, 0.0, 0.0, 0.0), G(0.0, 0.0, 0.0, 0.0);

		F = 0.5 * (F_minus + F_plus);
		G = 0.5 * (G_minus + G_plus);

		flux += normal[0] * F + normal[1] * G;

		std::array<Vec4, 4> l_A = {
				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(velocity_sqr_abs_star / 2 + velocityX_star * c_star / (m_gamma - 1.0),
							 -velocityX_star - c_star / (m_gamma - 1), -velocityY_star, 1.0),
				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(sqr(c_star) / (m_gamma - 1) -
									 velocity_sqr_abs_star / 2 - velocityY_star * c_star /(m_gamma - 1),
							 velocityX_star, velocityY_star + c_star / (m_gamma - 1), -1.0),
				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(sqr(c_star) / (m_gamma - 1) -
									 velocity_sqr_abs_star / 2 + velocityY_star * c_star /(m_gamma - 1),
							 velocityX_star, velocityY_star - c_star / (m_gamma - 1), -1.0),
				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(velocity_sqr_abs_star / 2 - velocityX_star * c_star / (m_gamma - 1.0),
							 -velocityX_star + c_star / (m_gamma - 1), -velocityY_star, 1.0)

		};

		std::array<Vec4, 4> l_B = {
				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(velocity_sqr_abs_star / 2 + velocityY_star * c_star / (m_gamma - 1),
							 -velocityX_star, -velocityY_star - c_star / (m_gamma - 1), 1.0),

				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(sqr(c_star) / (m_gamma - 1) -
									 velocity_sqr_abs_star / 2 - velocityX_star * c_star / (m_gamma - 1),
							 velocityX_star + c_star / (m_gamma - 1), velocityY_star, -1.0),
				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(sqr(c_star) / (m_gamma - 1) -
							 velocity_sqr_abs_star / 2 + velocityX_star * c_star / (m_gamma - 1),
							 velocityX_star - c_star / (m_gamma - 1), velocityY_star, -1.0),

				(m_gamma - 1.0) / (2.0 * sqr(c_star)) *
						Vec4(velocity_sqr_abs_star / 2 - velocityY_star * c_star / (m_gamma - 1),
							 -velocityX_star, -velocityY_star + c_star / (m_gamma - 1), 1.0),

		};

/*		auto const velocityX_local_minus = normal[0] * velocityX_minus + normal[1] * velocityY_minus;
		auto const velocityX_local_plus = normal[0] * velocityX_plus + normal[1] * velocityY_plus;

		auto const velocityY_local_minus = -normal[1] * velocityX_minus + normal[0] * velocityY_minus;
		auto const velocityY_local_plus = -normal[1] * velocityX_plus + normal[0] * velocityY_plus;

		auto const velocityX_local_star = (std::sqrt(density_minus) * velocityX_local_minus
									 + std::sqrt(density_plus) * velocityX_local_plus)
									/ (std::sqrt(density_minus) + std::sqrt(density_plus));

		auto const velocityY_local_star = (std::sqrt(density_minus) * velocityY_local_minus
									 + std::sqrt(density_plus) * velocityY_local_plus)
									/ (std::sqrt(density_minus) + std::sqrt(density_plus));



		std::array<Vec4, 4> r = {
				Vec4(1.0, velocityX_local_star - c_star, velocityY_local_star, H_star - velocityX_local_star * c_star),
				Vec4(1.0, velocityX_local_star, velocityY_local_star + c_star,
					 0.5 * velocity_sqr_abs_star + velocityY_local_star * c_star),

				Vec4(1.0, velocityX_local_star, velocityY_local_star - c_star, 0.5 * velocity_sqr_abs_star -
						velocityY_local_star * c_star),
				Vec4(1.0, velocityX_local_star + c_star, velocityY_local_star, H_star + velocityX_local_star * c_star)

		};

		std::array<double, 4> delta_s= {
				((pressure_plus - pressure_minus) - density_star * c_star * (velocityX_local_plus - velocityX_local_minus))
				/ (2 * sqr(c_star)),

				(sqr(c_star) * (density_plus - density_minus) - (pressure_plus - pressure_minus)
				 + density_star * c_star * (velocityY_local_plus - velocityY_local_minus)) / (2 * sqr(c_star)),

				(sqr(c_star) * (density_plus - density_minus) - (pressure_plus - pressure_minus)
				 - density_star * c_star * (velocityY_local_plus - velocityY_local_minus)) / (2 * sqr(c_star)),

				((pressure_plus - pressure_minus) + density_star * c_star * (velocityX_local_plus - velocityX_local_minus))
				/ (2 * sqr(c_star))


		};
*/


		for(int k = 0; k < 4; ++k)
			flux -= (std::fabs(normal[0] * lambda_A_waved[k]) + std::fabs(normal[1] * lambda_B_waved[k])) *
					(normal[0] * delta_s_A[k] + normal[1] * delta_s_B[k]) *
					(normal[0] * r_A[k] + normal[1] * r_B[k]);



	}

	return 0.5 * flux;
}