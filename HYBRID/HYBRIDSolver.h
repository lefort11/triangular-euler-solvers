#ifndef TRIANGULAR_SOLVERS_HYBRIDSOLVER_H
#define TRIANGULAR_SOLVERS_HYBRIDSOLVER_H

#include "../WENO/WENOSolver.h"
#include "../FirstOrderSolver/RoeSolver.h"

namespace euler
{

	class HYBRIDSolver: public WENOSolver<RoeSolver>
	{

	public:

		explicit HYBRIDSolver(std::vector<Zone> const &constraints,
				std::function<void(TriangularMesh const&, TriangularMesh const&, double)> const &bcFunc,
				MeshParams const &triangleProp,
				double gamma = 5.0 / 3.0) : WENOSolver<RoeSolver>(constraints, bcFunc, triangleProp, gamma)
		{}


	protected:

		double Rotor(Triangle* const pTriangle) const
		{
			auto const area = pTriangle->getArea2D();
			auto const h = std::sqrt(area);
			double rotor = 0.0;
			for(int i = 0; i < 3; ++i)
			{
				auto const l = std::sqrt(pTriangle->getSquaredEdgeLength(i));
				auto const n = pTriangle->CalculateNormal(i);



				auto const h_n = std::sqrt(pTriangle->GetOppTriangle(i)->getArea2D());

				auto const u = (h_n * pTriangle->GetOppTriangle(i)->velocityX + h * pTriangle->velocityX) / (h + h_n);
				auto const v = (h_n * pTriangle->GetOppTriangle(i)->velocityY + h * pTriangle->velocityY) / (h + h_n);

				rotor += l * (n[0] * v - n[1] * u);

			}
			rotor /= area;

			return rotor;

		}

		double AbsGradp(Triangle* const pTriangle) const
		{
			arma::vec3 grad{0.0, 0.0, 0.0};
			auto const area = pTriangle->getArea2D();
			auto const h = std::sqrt(area);

			for(int i = 0; i < 3; ++i)
			{
				auto const l = std::sqrt(pTriangle->getSquaredEdgeLength(i));
				auto const n = pTriangle->CalculateNormal(i);
				auto const h_n = std::sqrt(pTriangle->GetOppTriangle(i)->getArea2D());

				auto const p = (h_n * pTriangle->GetOppTriangle(i)->pressure + h * pTriangle->pressure) / (h + h_n);

				grad += arma::vec3{n[0] * p, n[1] * p};
			}
			grad /= area;
			return arma::norm(grad);
		}

		double Alpha(Triangle* const pTriangle) const
		{
			double const alpha_max = 0.3;

			double const s = pTriangle->getArea2D();
			double const h = std::sqrt(s);
			double const rho = pTriangle->density;

			double const phi = 1e-2;
			double const a = -phi / alpha_max;

			double R = rho * sqr(Rotor(pTriangle)) * h / (1e-40 + AbsGradp(pTriangle));

			return -phi / ( R - a ) + alpha_max;
		}

	public:

		Vec4 CalculateFlux(Vec4 const& qVec, int triangleNumber, int edgeNumber) const override;




	};
}

#endif //TRIANGULAR_SOLVERS_HYBRIDSOLVER_H
