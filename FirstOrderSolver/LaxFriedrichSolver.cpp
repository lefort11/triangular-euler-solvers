#include "LaxFriedrichSolver.h"

using namespace euler;

Vec4 LaxFriedrichSolver::CalculateFlux(Vec4 const &qVec, int triangleNumber, int edgeNumber) const
{
	//preparing q_minus and q_plus vectors

	//****** Calculating gaussian point ********//

	double x_g, y_g;
	auto const firstVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 1) % 3);
	auto const secondVertex = m_triangles[triangleNumber]->getCorner((edgeNumber + 2) % 3);

	x_g = (firstVertex->x() + secondVertex->x()) / 2;
	y_g = (firstVertex->y() + secondVertex->y() / 2);

	//*****************************************//

	//********** Forming q_minus vector and _minus parameters *******//
	auto const q_minus = Reconstruct(qVec, m_triangles[triangleNumber], x_g, y_g);

//	double density_minus, velocityX_minus, velocityY_minus, pressure_minus;

//	GetGasParamsFromQ(q_minus, density_minus, velocityX_minus, velocityY_minus, pressure_minus);

//	auto const eps_minus = pressure_minus / (density_minus * (m_gamma - 1.0));

	auto const density_minus = q_minus[0];
	auto const velocityX_minus = q_minus[1] / density_minus;
	auto const velocityY_minus = q_minus[2] / density_minus;
	auto const E_minus = q_minus[3] / density_minus;
	auto const velocity_sqr_abs_minus = sqr(velocityX_minus) + sqr(velocityY_minus);
	auto const eps_minus = E_minus - 0.5 * velocity_sqr_abs_minus;
	auto const pressure_minus= (m_gamma - 1.0) * eps_minus * density_minus;
	auto const H_minus = eps_minus + pressure_minus / density_minus + 0.5 * velocity_sqr_abs_minus;
	auto const c_minus = std::sqrt(m_gamma * pressure_minus / density_minus);


	//****************************************************************//


	//********* Forming q_plus vector and _plus parameteres **********//


	Vec4 q_plus;
	Triangle* neighbour_triangle = m_triangles[triangleNumber]->GetOppTriangle(edgeNumber);
	if(neighbour_triangle == nullptr)
		neighbour_triangle = m_triangles[triangleNumber];

	auto const density_plus = neighbour_triangle->density;
	auto const velocityX_plus = neighbour_triangle->velocityX;
	auto const velocityY_plus = neighbour_triangle->velocityY;
	auto const velocity_sqr_abs_plus = sqr(velocityX_plus) + sqr(velocityY_plus);
	auto const pressure_plus = neighbour_triangle->pressure;
	auto const eps_plus = pressure_plus / (density_plus * (m_gamma - 1.0));
//	auto const E_plus = eps_plus + 0.5 * (velocity_sqr_abs_plus); // E = eps + 1/2 * |v|^2
	auto const H_plus = eps_plus + pressure_plus / density_plus + 0.5 * velocity_sqr_abs_plus;
	auto const c_plus = std::sqrt(m_gamma * pressure_plus / density_plus);

	FormQVector(q_plus, density_plus, velocityX_plus, velocityY_plus, pressure_plus);

	q_plus = Reconstruct(q_plus, neighbour_triangle, x_g, y_g);

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

	//*********************************************************************//

	//**************** Calculating nu-factor ******************************//
	//Calculating A_minus norm
//	auto const A_minus_norm = std::max(std::fabs(velocityX_minus - c_minus),
//									   std::fabs(velocityX_minus + c_minus));

	//Calculating A_plus_norm

//	auto const A_plus_norm = std::max( std::fabs(velocityX_plus - c_plus),
//									   std::fabs(velocityX_plus + c_plus));

//	auto const nu_F = std::max(A_minus_norm, A_plus_norm);

	//Calculating B_minus norm
//	auto const B_minus_norm = std::max(std::fabs(velocityY_minus - c_minus),
//									   std::fabs(velocityY_minus + c_minus));

	//Calculating B_plus_norm

//	auto const B_plus_norm = std::max(std::fabs(velocityY_plus - c_plus),
//									  std::fabs(velocityY_plus + c_plus));

//	auto const nu_G = std::max(B_minus_norm, B_plus_norm);

	auto const norm_minus = std::fabs(sqrt(velocity_sqr_abs_minus) + c_minus);
	auto const norm_plus = std::fabs(sqrt(velocity_sqr_abs_plus) + c_plus);
	auto const nu = std::max(norm_minus, norm_plus);

	//************************************************************************//


	auto const normal = CalculateNormal(triangleNumber, edgeNumber);


	auto const x_Flux = 0.5 * (F_minus + F_plus); // - nu_F * (q_plus - q_minus));

	auto const y_Flux = 0.5 * (G_minus + G_plus); //- nu_G * (q_plus - q_minus));

	auto const flux = normal[0] * x_Flux + normal[1] * y_Flux - 0.5 * nu * (q_plus - q_minus);

	return flux;

//	return Vec4(0,0,0,0);
}