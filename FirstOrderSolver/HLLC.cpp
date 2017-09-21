#include "HLLC.h"

using namespace euler;


void HLLC::CreateBoundingMesh()
{

	for(int triangle_counter = 0; triangle_counter < m_triangles.size(); ++triangle_counter)
	{

		for(int edge_number = 0; edge_number < 3; ++edge_number)
		{
			if(m_triangles[triangle_counter]->GetOppTriangle(edge_number) == nullptr)
			{
				auto const reflectedTriangle = m_triangles[triangle_counter]->ReflectTriangle(edge_number);
				m_boundingTriangles.push_back(reflectedTriangle);
			}
		}

	}

}


Vec4 HLLC::CalculateFlux(Vec4 const &qVec, int triangleNumber, int edgeNumber) const
{
	//preparing q_minus and q_plus vectors

	//****** Calculating gaussian points ********//

	static auto const gaussian_weight = 1.0 / 2.0 + sqrt(3.0) / 6.0;

	Point2 gaussian_p_1, gaussian_p_2;
	auto const firstVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 1) % 3);
	auto const secondVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 2) % 3);

	gaussian_p_1.x = gaussian_weight * firstVertex->x() + (1 - gaussian_weight) * secondVertex->x();
	gaussian_p_1.y = gaussian_weight * firstVertex->y() + (1 - gaussian_weight) * secondVertex->y();

	gaussian_p_2.x = gaussian_weight * secondVertex->x() + (1 - gaussian_weight) * firstVertex->x();
	gaussian_p_2.y = gaussian_weight * secondVertex->y() + (1 - gaussian_weight) * firstVertex->y();

	auto const normal = CalculateNormal(m_triangles[triangleNumber], edgeNumber);

	Vec4 flux(0.0, 0.0, 0.0, 0.0);

	for (Point2 g_point: {gaussian_p_1, gaussian_p_2})
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


		// ********* pressure estimate **************/

		double const velX_local_minus = normal[0] * velocityX_minus + normal[1] * velocityY_minus;
		double const velX_local_plus = normal[0] * velocityX_plus + normal[1] * velocityY_plus;

		double const velY_local_minus = -normal[1] * velocityX_minus + normal[0] * velocityY_minus;
		double const velY_local_plus = -normal[1] * velocityX_plus + normal[0] * velocityY_plus;



		auto const ro_bar = 0.5 * (density_minus + density_plus);
		auto const c_bar = 0.5 * (c_minus + c_plus);

		auto const P_pvrs = 0.5 * (pressure_minus + pressure_plus)
							- 0.5 * (velX_local_plus - velX_local_minus) * ro_bar * c_bar;

		auto const P_star = std::fmax(0, P_pvrs);
		//********************************************//

		// ********** wave speed estimates ************** //
		auto const q_l = P_pvrs <= pressure_minus ? 1 :
						 std::sqrt(1 + (m_gamma + 1) / (2 * m_gamma) * (P_star / pressure_minus - 1));
		auto const q_r = P_pvrs <= pressure_plus ? 1 :
						 std::sqrt(1 + (m_gamma + 1) / (2 * m_gamma) * (P_star / pressure_plus - 1));

		auto const S_l = velX_local_minus - c_minus * q_l;
		auto const S_r = velX_local_plus + c_plus * q_r;

		auto const S_star = ((pressure_plus - pressure_minus) + density_minus * velX_local_minus
																* (S_l - velX_local_minus)
							- density_plus * velX_local_plus * (S_r - velX_local_plus))
							/ (density_minus * (S_l - velX_local_minus) - density_plus * (S_r - velX_local_plus));

		// ********************************************* //

		// *************** HLLC flux ******************** //



		Vec4 const F_l(
				density_minus * velX_local_minus,
				density_minus * velX_local_minus * velX_local_minus + pressure_minus,
				density_minus * velX_local_minus * velY_local_minus,
				density_minus * velX_local_minus * H_minus
		);

		Vec4 const F_r(
				density_plus * velX_local_plus,
				density_plus * velX_local_plus * velX_local_plus + pressure_plus,
				density_plus * velX_local_plus * velY_local_plus,
				density_plus * velX_local_plus * H_plus
		);


		if(0 <= S_l)
		{
			flux += F_l;
		}
		else if((S_l <= 0) && (0 <= S_star))
		{
			Vec4 U_star_l(
					1.0,
					S_star,
					velY_local_minus,
					E_minus +
					(S_star - velX_local_minus) * (S_star + pressure_minus / (density_minus * (S_l - velX_local_minus)))
			);

			U_star_l *= density_minus * (S_l - velX_local_minus) / (S_l - S_star);
			Vec4 U_l(
					density_minus,
					density_minus * velX_local_minus,
					density_minus * velY_local_minus,
					density_minus * E_minus
			);

			flux += F_l + S_l * (U_star_l - U_l);
		}
		else if((S_star <=0) && (0 <= S_r))
		{

			Vec4 U_star_r(1.0,
						  S_star,
						  velY_local_plus,
						  E_plus + (S_star - velX_local_plus) *
								   (S_star + pressure_plus / (density_plus * (S_r - velX_local_plus)))
			);


			U_star_r *= density_plus * (S_r - velX_local_plus) / (S_r - S_star);

			Vec4 U_r(
					density_plus,
					density_plus * velX_local_plus,
					density_plus * velY_local_plus,
					density_plus * E_plus
			);


			flux += F_r + S_r * (U_star_r - U_r);
		}
		else if(S_r <= 0)
			flux += F_r;
		else
			throw 7;



	}

	return 0.5 * flux;
}