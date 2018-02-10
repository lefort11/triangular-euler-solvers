#include "Solver.h"
#include <omp.h>
#include "../Maths/Algebra.h"
#include <sstream>
#include <iomanip>


using namespace euler;

void Solver::Calculate(double time)
{

	std::cout << m_triangles.size() << std::endl;

	auto currentTime = 0.0;
	m_delta_t = CalculateTimeStep();

	if(m_delta_t > time)
		m_delta_t = time;
	//currentTime += delta_t;

	std::vector<Vec4> currentQs(m_triangles.size());

	//initializinq Q vector
	for(int i = 0; i < currentQs.size(); ++i)
	{
		FormQVector(currentQs[i], m_triangles[i]);
	}


	std::vector<Vec4> nextQs(m_triangles.size());

	unsigned timeLayersNumber = 0;


	while(m_delta_t > 0)
	{

		UpdateBoundingMesh(currentTime);

		std::cout << currentTime << std::endl;



		//trngl_cntr - triangle counter
#pragma omp parallel for
		for(int trngl_number = 0; trngl_number < m_triangles.size(); ++trngl_number)
		{

			auto const next_q = RungeKuttaTVDStep(currentQs[trngl_number], [this, trngl_number](Vec4 const &qVec)
			{
				auto const area = m_triangles[trngl_number]->getArea2D();

				Vec4 sum{0.0, 0.0, 0.0, 0.0};

				for(int edge_number = 0; edge_number < 3; ++edge_number)
				{
					auto const edge_length = std::sqrt(m_triangles[trngl_number]->getSquaredEdgeLength(edge_number));;
					sum += edge_length * CalculateFlux(qVec, trngl_number, edge_number);

				}

				sum = (-1.0)/area * sum;

				return sum;
			});



			nextQs[trngl_number] = next_q;

		}



		//filling next time layer data
#pragma omp parallel for
		for(int trngl_cntr = 0; trngl_cntr < m_triangles.size(); ++trngl_cntr)
		{

			GetGasParamsFromQ(nextQs[trngl_cntr], m_triangles[trngl_cntr]->density, m_triangles[trngl_cntr]->velocityX,
							  m_triangles[trngl_cntr]->velocityY, m_triangles[trngl_cntr]->pressure);

		}

		currentTime += m_delta_t;

		if(timeLayersNumber % 20 == 0)
		{
			std::stringstream stln;
			stln << std::setw(10) << std::setfill('0') << timeLayersNumber;
			std::string clcPath("results/clc/" + stln.str() + ".clc");
			ClcOutput(clcPath, m_delta_t, currentTime, timeLayersNumber);

		}

		m_delta_t = CalculateTimeStep();
		if(currentTime + m_delta_t > time)
			m_delta_t = time - currentTime;
		//currentTime += m_delta_t;

		currentQs = nextQs;


		++timeLayersNumber;


	}

	std::cout << timeLayersNumber << std::endl;

}





double Solver::GaussianQuadrilateralIntegration(std::function<double(double, double)> const &g) const
{
	double ksi[3] = {
			- sqrt(3.0)/sqrt(5.0),
			0.0,
			sqrt(3.0)/sqrt(5.0)
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
