#include "Solver.h"
#include <omp.h>
#include "../Maths/Algebra.h"


using namespace euler;

void Solver::Calculate(double time) const
{

	std::cout << m_triangles.size() << std::endl;

	auto currentTime = 0.0;
	auto delta_t = CalculateTimeStep();

	if(delta_t > time)
		delta_t = time;
	currentTime += delta_t;

	std::vector<Vec4> currentQs(m_triangles.size());

	//initializinq Q vector
	for(int i = 0; i < currentQs.size(); ++i)
	{
		FormQVector(currentQs[i], m_triangles[i]);
	}


	std::vector<Vec4> nextQs(m_triangles.size());


	while(delta_t != 0)
	{


		//trngl_cntr - triangle counter
		for(int trngl_cntr = 0; trngl_cntr < m_triangles.size(); ++trngl_cntr)
		{

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
#pragma omp parallel for
		for(int trngl_cntr = 0; trngl_cntr < m_triangles.size(); ++trngl_cntr)
		{

			GetGasParamsFromQ(nextQs[trngl_cntr], m_triangles[trngl_cntr]->density, m_triangles[trngl_cntr]->velocityX,
							  m_triangles[trngl_cntr]->velocityY, m_triangles[trngl_cntr]->pressure);

		}

		delta_t = CalculateTimeStep();
		if(currentTime + delta_t > time)
			delta_t = time - currentTime;
		currentTime += delta_t;


		currentQs = nextQs;

	}

}


Vec4 Solver::RungeKuttaTVDStep(Vec4 const&current_q, double delta_t, std::function<Vec4(Vec4)> const& f) const
{
	auto const q_0 = current_q;

	auto const q_1 = q_0 + delta_t * f(q_0);

	auto const next_q = q_1;

//	auto const q_2 = 3.0/4.0 * q_0 + 1.0/4.0 * q_1 + 1.0/4.0 * delta_t * f(q_1);

//	auto const next_q = 1.0/3.0 * q_0 + 2.0/3.0 * q_1 + 2.0/3.0 * delta_t * f(q_2);

	return next_q;

}