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
#pragma omp parallel for
		for(int trngl_number = 0; trngl_number < m_triangles.size(); ++trngl_number)
		{

			auto const next_q = RungeKuttaTVDStep(currentQs[trngl_number], delta_t, [this, trngl_number](Vec4 qVec)
			{
				auto const area = CalculateTriangleArea(trngl_number);

				Vec4 sum = {0.0, 0.0, 0.0, 0.0};

				for(int edge_number = 0; edge_number < 3; ++edge_number)
				{
					auto const edge_length = CalculateTriangleEdgeLength(trngl_number, edge_number);
					sum += edge_length * CalculateFlux(qVec, trngl_number, edge_number);
				}

				sum = (-1.0)/area * sum;

				return sum;
			});

			double a;
			nextQs[trngl_number] = next_q;

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


Vec4 Solver::RungeKuttaTVDStep(Vec4 const& current_q, double delta_t, std::function<Vec4(Vec4)> const& f) const
{
	auto const q_0 = current_q;

	auto const q_1 = q_0 + delta_t * f(q_0);

	auto const next_q = q_1;

//	auto const q_2 = 3.0/4.0 * q_0 + 1.0/4.0 * q_1 + 1.0/4.0 * delta_t * f(q_1);

//	auto const next_q = 1.0/3.0 * q_0 + 2.0/3.0 * q_1 + 2.0/3.0 * delta_t * f(q_2);

	return next_q;

}



double Solver::GaussianQuadrilateralIntegration(std::function<double(double, double)> const &g) const
{
	double ksi[3] = {
			- sqrt(3)/sqrt(5),
			0.0,
			sqrt(3)/sqrt(5)
	};

	double omega[3] = {
			5.0 / 9.0,
			8.0 / 9.0,
			5.0 / 9.0
	};
	double result = 0.0;
	for(int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			result += omega[i] * omega[j] * g(ksi[i], ksi[j]);
		}
	}

	return result;


}



double Solver::GaussianIntegration(std::function<double(double, double)> const &func, Triangle const *pTriangle) const
{

	auto const result = GaussianQuadrilateralIntegration([&func, pTriangle](double ksi, double eta)
											  {

												  auto const x_0 = pTriangle->getCorner(0)->x();
												  auto const x_1 = pTriangle->getCorner(1)->x();
												  auto const x_2 = pTriangle->getCorner(2)->x();

												  auto const y_0 = pTriangle->getCorner(0)->y();
												  auto const y_1 = pTriangle->getCorner(1)->y();
												  auto const y_2 = pTriangle->getCorner(2)->y();

												  auto const factor = (x_1 - x_0)*(y_2 - y_0) - (x_2 - x_0)*(y_1 - y_0);

												  auto const x = x_0 +
														  0.5 * (x_1 - x_0) * (1 + ksi)
																 + 0.25 * (x_2 - x_0) * (1 - ksi) * (1 + eta);
												  auto const y = y_0 +
																 0.5 * (y_1 - y_0) * (1 + ksi)
																 + 0.25 * (y_2 - y_0) * (1 - ksi) * (1 + eta);

												  return func(x,y) * (1 - ksi) / 8 * factor;



											  });

	return result;

}
