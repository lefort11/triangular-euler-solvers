#include <iostream>
#include <cmath>

#include "SpaceMesh/Area.h"

#include "FirstOrderSolver/LaxFriedrichSolver.h"



int main()
{


	euler::ConstraintFunction circle1([](double t)
	{
		auto x =cos(2*M_PI * t);
		auto y = sin(2*M_PI * t);
		return GEOM_FADE2D::Point2(x, y);
	});

	euler::ConstraintFunction circle2([](double t)
	{
		auto x = 1.1 * cos(2*M_PI * t);
		auto y = 1.1 * sin(2*M_PI * t);
		return GEOM_FADE2D::Point2(x, y);
	});

	euler::Zone zone(circle1, false, {0.0, 0.0, 0.0, 0.0}),
				zone1(circle2, true, {0.0, 0.0, 0.0, 0.0});

	std::vector<euler::Zone> vZone;
//	vZone.push_back(zone);
//	vZone.push_back(zone1);


	std::array<double, 3> trProp = {27, 0.001, 0.2};



	euler::LaxFriedrichSolver solver(vZone, 1000, trProp);

	solver.Init([](GEOM_FADE2D::Point2 point)
				{
					if(point.x() < 1.0/2.0)
						return euler::Vec4({2.0,0.0,0.0,5.0});
					return euler::Vec4({1.0,0.0,0.0,1.0});

				});

	solver.Calculate(0.177);

	solver.DebugOutput("results/density.txt");

	return 0;
}