#include "LaxFriedrichSolver.h"

using namespace euler;


void LaxFriedrichSolver::CreateBoundingMesh()
{

	for(int triangle_counter = 0; triangle_counter < m_triangles.size(); ++triangle_counter)
	{

		for(int edge_number = 0; edge_number < 3; ++edge_number)
		{
			if(m_triangles[triangle_counter]->GetOppTriangle(edge_number) == nullptr)
			{
				auto const reflectedTriangle = m_triangles[triangle_counter]->ReflectTriangle(edge_number);

				reflectedTriangle->SetBoundary(true);
				reflectedTriangle->SetParentIndex(m_triangles[triangle_counter]->Index());
				reflectedTriangle->SetIndex(m_boundingTriangles.size());

				m_boundingTriangles.push_back(reflectedTriangle);
			}
		}

	}

}


Vec4 LaxFriedrichSolver::CalculateFlux(Vec4 const &qVec, int triangleNumber, int edgeNumber) const
{
	//preparing q_minus and q_plus vectors

	//****** Calculating gaussian points ********//

	static auto const gaussian_weight = 1.0 / 2.0 + sqrt(3.0) / 6.0;

	Point2  gaussian_p_1, gaussian_p_2;
	auto const firstVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 1) % 3);
	auto const secondVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 2) % 3);

	gaussian_p_1.x = gaussian_weight * firstVertex->x() + (1 - gaussian_weight) * secondVertex->x();
	gaussian_p_1.y = gaussian_weight * firstVertex->y() + (1 - gaussian_weight) * secondVertex->y();

	gaussian_p_2.x = gaussian_weight * secondVertex->x() + (1 - gaussian_weight) * firstVertex->x();
	gaussian_p_2.y = gaussian_weight * secondVertex->y() + (1 - gaussian_weight) * firstVertex->y();

	Vec4 flux{0.0, 0.0, 0.0, 0.0};

	auto const normal = CalculateNormal(m_triangles[triangleNumber], edgeNumber);


	for (Point2 g_point : {gaussian_p_1, gaussian_p_2})
	{
		//*****************************************//

		//********** Forming q_minus vector and _minus parameters *******//
		Vec4 const q_minus = Reconstruct(qVec, m_triangles[triangleNumber], g_point, edgeNumber);

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
		auto const v1_ind =
			neighbour_triangle->getIntraTriangleIndex(m_triangles[triangleNumber]->getCorner((edgeNumber + 2) % 3));
		int neighbourEdgeNumber = (v0_ind + 1) % 3;
		if (neighbourEdgeNumber == v1_ind)
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

		//*********************************************************************//

		//**************** Calculating nu-factor ******************************//
		//Calculating A_minus norm
//		auto const A_minus_norm = std::max(std::fabs(velocityX_minus - c_minus),
//										   std::fabs(velocityX_minus + c_minus));

		//Calculating A_plus_norm

//		auto const A_plus_norm = std::max(std::fabs(velocityX_plus - c_plus),
//										  std::fabs(velocityX_plus + c_plus));

//		auto const nu_F = std::max(A_minus_norm, A_plus_norm);

		//Calculating B_minus norm
//		auto const B_minus_norm = std::max(std::fabs(velocityY_minus - c_minus),
//										   std::fabs(velocityY_minus + c_minus));

		//Calculating B_plus_norm

//		auto const B_plus_norm = std::max(std::fabs(velocityY_plus - c_plus),
//										  std::fabs(velocityY_plus + c_plus));

//		auto const nu_G = std::max(B_minus_norm, B_plus_norm);

/*		auto const norm_minus = std::fabs(sqrt(velocity_sqr_abs_minus) + c_minus);
		auto const norm_plus = std::fabs(sqrt(velocity_sqr_abs_plus) + c_plus);
		auto const nu = std::max(norm_minus, norm_plus); */

		auto const vel_n_minus = normal[0] * velocityX_minus + normal[1] * velocityY_minus;
		auto const vel_n_plus = normal[0] * velocityX_plus + normal[1] * velocityY_plus;
		auto const norm_minus = std::max(std::fabs(vel_n_minus - c_minus), std::fabs(vel_n_minus + c_minus));
		auto const norm_plus = std::max(std::fabs(vel_n_plus - c_plus), std::fabs(vel_n_plus + c_plus));


		auto const nu = std::max(norm_minus, norm_plus); 

		//************************************************************************//


		Vec4 const x_Flux = 0.5 * (F_minus + F_plus);

		Vec4 const y_Flux = 0.5 * (G_minus + G_plus); 

		flux += normal[0] * x_Flux + normal[1] * y_Flux - 0.5 * nu * (q_plus - q_minus);
	}




	return 0.5 * flux;

//	return Vec4(0,0,0,0);
}