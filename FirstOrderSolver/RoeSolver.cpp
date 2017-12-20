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
				m_triangles[triangle_counter]->SetBoundary(true);

				auto const reflectedTriangle = m_triangles[triangle_counter]->ReflectTriangle(edge_number);

				reflectedTriangle->SetVirtual(true);
				reflectedTriangle->SetParentIndex(m_triangles[triangle_counter]->Index());
				reflectedTriangle->SetIndex(m_boundingTriangles.size());

				m_boundingTriangles.push_back(reflectedTriangle);
			}
		}

	}

}



Vec4 RoeSolver::CalculateFlux(Vec4 const &qVec, int triangleNumber, int edgeNumber) const
{


	//preparing q_minus and q_plus vectors


	Vec4 flux{0.0, 0.0, 0.0, 0.0};

	auto const normal = CalculateNormal(m_triangles[triangleNumber], edgeNumber);


	for (int gPointNum = 0; gPointNum < 2; ++gPointNum)
	{
		//*****************************************//

		//********** Forming q_minus vector and _minus parameters *******//
		Vec4 const q_minus = Reconstruct(qVec, m_triangles[triangleNumber], edgeNumber, gPointNum);

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
		//		if (neighbour_triangle == nullptr)
		//			neighbour_triangle = m_triangles[triangleNumber];

		auto density_plus = neighbour_triangle->density;
		auto velocityX_plus = neighbour_triangle->velocityX;
		auto velocityY_plus = neighbour_triangle->velocityY;
		auto pressure_plus = neighbour_triangle->pressure;

		FormQVector(q_plus, density_plus, velocityX_plus, velocityY_plus, pressure_plus);

		//edge_number calculation
		auto const v0_ind =
				neighbour_triangle->getIntraTriangleIndex(m_triangles[triangleNumber]->getCorner((edgeNumber + 1) % 3));

		int neighbourEdgeNumber = (v0_ind + 1) % 3;

		q_plus = Reconstruct(q_plus, neighbour_triangle, neighbourEdgeNumber, (gPointNum + 1) % 2);


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
		Vec4 const F_minus{
				density_minus * velocityX_minus,
				density_minus * sqr(velocityX_minus) + pressure_minus,
				density_minus * velocityX_minus * velocityY_minus,
				density_minus * velocityX_minus * H_minus
		};

		//Calculating G_Minus
		Vec4 const G_minus{
				density_minus * velocityY_minus,
				density_minus * velocityY_minus * velocityX_minus,
				density_minus * sqr(velocityY_minus) + pressure_minus,
				density_minus * velocityY_minus * H_minus

		};

		//Calculating F_plus
		Vec4 const F_plus{
				density_plus * velocityX_plus,
				density_plus * sqr(velocityX_plus) + pressure_plus,
				density_plus * velocityX_plus * velocityY_plus,
				density_plus * velocityX_plus * H_plus
		};

		//Calculating G_plus
		Vec4 const G_plus{
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


		//auto const c_star = std::sqrt( (m_gamma - 1) * (H_star -  0.5 * velocity_sqr_abs_star) );

		auto const delta_v_abs_sqr = sqr(velocityX_plus - velocityX_minus) + sqr(velocityY_plus - velocityY_minus);
		auto const c_star = std::sqrt((std::sqrt(density_minus) * sqr(c_minus) + std::sqrt(density_plus) * sqr(c_plus))
									  / (std::sqrt(density_minus) + std::sqrt(density_plus))
									  + (m_gamma - 1.0) / 2.0
										* std::sqrt(density_minus * density_plus) /
										sqr(std::sqrt(density_minus) + std::sqrt(density_plus))
										* delta_v_abs_sqr);


		auto const vel_n_x_minus = normal[0] * velocityX_minus + normal[1] * velocityY_minus;
//		auto const vel_n_y_minus = -normal[1] * velocityX_minus + normal[0] * velocityY_minus;

		auto const vel_n_x_plus = normal[0] * velocityX_plus + normal[1] * velocityY_plus;
//		auto const vel_n_y_plus = -normal[1] * velocityX_plus + normal[0] * velocityY_plus;

		auto const vel_n_x_star = (std::sqrt(density_minus) * vel_n_x_minus
								   + std::sqrt(density_plus) * vel_n_x_plus)
								  / (std::sqrt(density_minus) + std::sqrt(density_plus));

//		auto const vel_n_y_star = (std::sqrt(density_minus) * vel_n_y_minus
//								   + std::sqrt(density_plus) * vel_n_y_plus)
//								  / (std::sqrt(density_minus) + std::sqrt(density_plus));
		auto const velocity_sqr_abs_star = sqr(velocityX_star) + sqr(velocityY_star);


		std::array<double, 4> const lambdas_n =
				{
						vel_n_x_star - c_star,
						vel_n_x_star,
						vel_n_x_star,
						vel_n_x_star + c_star
				};

//		auto const eta_vl_0 = std::max((vel_n_x_plus - c_plus) - (vel_n_x_minus - c_minus), 0.0);
//		auto const eta_vl_1 = std::max(vel_n_x_plus - vel_n_x_minus, 0.0);

//		auto const eta_vl_3 = std::max((vel_n_x_plus + c_plus) - (vel_n_x_minus + c_minus), 0.0);

		auto const eta_vl_0 = std::max({lambdas_n[0] - (vel_n_x_minus - c_minus),
										(vel_n_x_plus - c_plus) - lambdas_n[0], 0.0});

		auto const eta_vl_1 = std::max(vel_n_x_plus - vel_n_x_minus, 0.0);
		auto const eta_vl_3 = std::max({lambdas_n[3] - (vel_n_x_minus + c_minus),
										(vel_n_x_plus + c_plus) - lambdas_n[3], 0.0});

 /*       auto const lambda_neigh_0 = normal[0] * m_triangles[triangleNumber]->GetOppTriangle((edgeNumber + 1) % 3)->velocityX +
                                    normal[1] * m_triangles[triangleNumber]->GetOppTriangle((edgeNumber + 1) % 3)->velocityY;
        auto const lambda_neigh_1 = normal[0] * m_triangles[triangleNumber]->GetOppTriangle((edgeNumber + 2) % 3)->velocityX +
                                    normal[1] * m_triangles[triangleNumber]->GetOppTriangle((edgeNumber + 2) % 3)->velocityY;

        auto const lambda_neigh_2 = normal[0] * neighbour_triangle->GetOppTriangle((neighbourEdgeNumber + 1) % 3)->velocityX +
                                    normal[1] * neighbour_triangle->GetOppTriangle((neighbourEdgeNumber + 1) % 3)->velocityY;

        auto const lambda_neigh_3 = normal[0] * neighbour_triangle->GetOppTriangle((neighbourEdgeNumber + 2) % 3)->velocityX +
                                    normal[1] * neighbour_triangle->GetOppTriangle((neighbourEdgeNumber + 2) % 3)->velocityY;

        auto const eta_pa = 0.5 *
                std::max({std::fabs(lambdas_n[1] - lambda_neigh_0), std::fabs(lambdas_n[1] - lambda_neigh_1),
                          std::fabs(vel_n_x_plus - lambda_neigh_2), std::fabs(vel_n_x_plus - lambda_neigh_3)}); */


		std::array<double, 4> const lambdas_n_waved =
				{
						std::fabs(lambdas_n[0]) >= 2 * eta_vl_0 ? std::fabs(lambdas_n[0])
																: sqr(lambdas_n[0]) / (4 * eta_vl_0) + eta_vl_0,

                        //std::max(std::fabs(lambdas_n[1]), eta_pa),
                       // std::max(std::fabs(lambdas_n[2]), eta_pa),

                        std::fabs(lambdas_n[1]) >= 2 * eta_vl_1 ? std::fabs(lambdas_n[1])
                                                               	: sqr(lambdas_n[1]) / (4 * eta_vl_1) + eta_vl_1,
                        std::fabs(lambdas_n[1]) >= 2 * eta_vl_1 ? std::fabs(lambdas_n[1])
                                                                : sqr(lambdas_n[1]) / (4 * eta_vl_1) + eta_vl_1,

                        std::fabs(lambdas_n[3]) >= 2 * eta_vl_3 ? std::fabs(lambdas_n[3])
																: sqr(lambdas_n[3]) / (4 * eta_vl_3) + eta_vl_3
				};

		std::array<Vec4, 4> r_n =
				{
						Vec4{1.0,
							 velocityX_star - c_star * normal[0],
							 velocityY_star - c_star * normal[1],
							 H_star - c_star * vel_n_x_star
						},
						Vec4{1.0,
							 velocityX_star,
							 velocityY_star,
							 0.5 * velocity_sqr_abs_star
						},
						Vec4{0.0,
							 normal[1],
							 -normal[0],
							 velocityX_star * normal[1] - velocityY_star * normal[0]
						},
						Vec4{1.0,
							 velocityX_star + c_star * normal[0],
							 velocityY_star + c_star * normal[1],
							 H_star + c_star * vel_n_x_star
						}
				};

		std::array<Vec4, 4> l_n =
				{
						Vec4{ ( (m_gamma - 1) * 0.5 * velocity_sqr_abs_star + c_star * vel_n_x_star ) / (2 * sqr(c_star)),
							  ( (1 - m_gamma) * velocityX_star - c_star * normal[0] ) / (2 * sqr(c_star)),
							  ( (1 - m_gamma) * velocityY_star - c_star * normal[1] ) / (2 * sqr(c_star)),
							  (m_gamma - 1) / (2 * sqr(c_star))
						},

						Vec4{ ( sqr(c_star) - (m_gamma - 1) * 0.5 * velocity_sqr_abs_star ) / sqr(c_star),
							  (m_gamma - 1) * velocityX_star / sqr(c_star),
							  (m_gamma - 1) * velocityY_star / sqr(c_star),
							  (1 - m_gamma) / sqr(c_star)
						},

						Vec4{ (velocityY_star * normal[0] - velocityX_star * normal[1]),
							 normal[1],
							 -normal[0],
							 0.0
						},

						Vec4{ ( (m_gamma - 1) * 0.5 * velocity_sqr_abs_star - c_star * vel_n_x_star ) / (2 * sqr(c_star)),
							  ( (1 - m_gamma) * velocityX_star + c_star * normal[0] ) / (2 * sqr(c_star)),
							  ( (1 - m_gamma) * velocityY_star + c_star * normal[1] ) / (2 * sqr(c_star)),
							  (m_gamma - 1) / (2 * sqr(c_star))
						}

				};


        flux += 0.5 * (normal[0] * (F_minus + F_plus) + normal[1] * (G_minus + G_plus));
        for (int k = 0; k < 4; ++k)
        {
            flux -= 0.5 * std::fabs(lambdas_n_waved[k]) * arma::dot(l_n[k], q_plus - q_minus) * r_n[k];
        }


	}

	return 0.5 * flux;
}