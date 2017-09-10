#include <iostream>
#include <cmath>

#include "SpaceMesh/Area.h"

#include "FirstOrderSolver/LaxFriedrichSolver.h"
#include "WENO/WENOLF.h"



int main()
{


	euler::ConstraintFunction circle1([](double t)
	{
		auto x = 0.5 * cos(2.0 *M_PI * t);
		auto y = 0.5 * sin(2.0 *M_PI * t);
		return GEOM_FADE2D::Point2(x, y);
	}, 50);

	euler::ConstraintFunction circle2([](double t)
	{
		auto x = 1.1 * cos(2.0 *M_PI * t);
		auto y = 1.1 * sin(2.0 *M_PI * t);
		return GEOM_FADE2D::Point2(x, y);
	}, 100);

	euler::ConstraintFunction square([](double t)
									 {
										 double x, y;
										 if((t >= 0) && (t < 0.25))
										 {
											 x = -1;
											 y = -4 * t + 1;
										 }
										 else if((t >= 0.25) && (t < 0.5))
										 {
											 y = -1;
											 x = 4*(t - 0.25) - 1;
										 }
										 else if((t >= 0.5) && (t < 0.75))
										 {
											 x = 1;
											 y = 4*(t - 0.5) - 1;
										 }
										 else
										 {
											 y = 1;
											 x = -4*(t - 0.75) + 1;
										 }

										 return GEOM_FADE2D::Point2(x, y);


									 }, 4);

	euler::Zone zone(circle1, false),
				zone1(circle2, true),
				zone3(square, true);

	std::vector<euler::Zone> vZone;
	vZone.push_back(zone);
//	vZone.push_back(zone3);


	std::array<double, 3> trProp = {27, 0.01, 1.0};



	euler::WENOLF solver(vZone, trProp);

/*	solver.Init([](GEOM_FADE2D::Point2 point)
				{
					if(point.x() < 0.0)
						return std::array<double, 4>{{2.0,0.0,0.0,5.0}};
					return std::array<double, 4>{{1.0,0.0,0.0,1.0}};

				}); */

	solver.Init([](GEOM_FADE2D::Point2 point)
				{
					return std::array<double, 4>{{1.0, 1.0, 0.0, 1.0}};
				});

	solver.Calculate(0.1);

	solver.DebugOutput("results/densityLF.txt");
	solver.Output("results/density2D.txt", "results/velocity2D.txt");



	return 0;
}