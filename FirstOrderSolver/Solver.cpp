#include "Solver.h"
#include "../Maths/Algebra.h"


using namespace euler;

void Solver::Calculate(double time) const
{

	/****** Runge-Kutta method parameters ******/
	int const N_time = 1000;
	double const delta_t = time / N_time;
	/*******************************************/



	for(int time_counter = 0; time_counter < N_time; ++time_counter)
	{
		std::vector<Vec4> currentQs(m_triangles.size());

		//initializinq Q vector
		for(int i = 0; i < currentQs.size(); ++i)
		{
			auto const density = m_triangles[i]->density;
			auto const velocityX = m_triangles[i]->velocityX;
			auto const velocityY = m_triangles[i]->velocityY;
			auto const pressure = m_triangles[i]->pressure;
			auto const eps = pressure / (density * (m_gamma - 1.0));
			auto const E = eps + 0.5 * (sqr(velocityX) + sqr(velocityY)); // E = eps + 1/2 * |v|^2

			currentQs[i][0] = density;
			currentQs[i][1] = density * velocityX;
			currentQs[i][2] = density * velocityY;
			currentQs[i][3] = density * E;

		}

		std::vector<Vec4> nextQs(m_triangles.size());


		//trngl_cntr - triangle counter
		for(int trngl_cntr = 0; trngl_cntr < m_triangles.size(); ++trngl_cntr)
		{


			if(trngl_cntr == 510)
				double kek;
			auto const next_q = RungeKuttaTVDStep(currentQs[trngl_cntr], delta_t, [this, trngl_cntr](Vec4 qVec)
			{
				auto const area = CalculateTriangleArea(trngl_cntr);

				Vec4 sum = {0.0, 0.0, 0.0, 0.0};

				for(int edge_cntr = 0; edge_cntr < 3; ++edge_cntr)
				{
					auto const edge_length = CalculateTriangleEdgeLength(trngl_cntr, edge_cntr);
					sum += edge_length * CalculateFlux(qVec, trngl_cntr, edge_cntr);
				}

				sum = (-1.0)/area * sum;

				return sum;
			});

			nextQs[trngl_cntr] = next_q;
		}

		//filling next time layer data
		for(int trngl_cntr = 0; trngl_cntr < m_triangles.size(); ++trngl_cntr)
		{
			auto const density_new = nextQs[trngl_cntr][0];
			auto const velocityX_new = nextQs[trngl_cntr][1] / density_new;
			auto const velocityY_new = nextQs[trngl_cntr][2] / density_new;
			auto const E_new = nextQs[trngl_cntr][3] / density_new;

			auto const velocity_sqr_abs = sqr(velocityX_new) + sqr(velocityY_new);

			auto const eps = E_new - 0.5 * velocity_sqr_abs;

			auto const pressure_new= (m_gamma - 1.0) * eps * density_new;

			m_triangles[trngl_cntr]->density = density_new;
			m_triangles[trngl_cntr]->velocityX = velocityX_new;
			m_triangles[trngl_cntr]->velocityY = velocityY_new;
			m_triangles[trngl_cntr]->pressure = pressure_new;

		}


		currentQs = nextQs;
	}


}


Vec4 Solver::RungeKuttaTVDStep(Vec4 const&current_q, double delta_t, std::function<Vec4(Vec4)> const& f) const
{
	auto const q_0 = current_q;

	auto const q_1 = q_0 + delta_t * f(q_0);

	auto const q_2 = 3.0/4.0 * q_0 + 1.0/4.0 * q_1 + 1.0/4.0 * delta_t * f(q_1);

	auto const next_q = 1.0/3.0 * q_0 + 2.0/3.0 * q_1 + 2.0/3.0 * delta_t * f(q_2);

	return next_q;

}