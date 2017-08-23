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

//	euler::Area area(vZone);

	std::array<double, 3> trProp = {27, 0.01, 10};
//	area.Triangulate(1000, trProp, [](GEOM_FADE2D::Point2){
//		return std::array<double,4>{1.0,0.0,0.0,0.0};
//	});


	euler::LaxFriedrichSolver solver(vZone, 1000, trProp);

	solver.Init([](GEOM_FADE2D::Point2)
				{
					return std::array<double,4>{1.0,0.0,0.0,0.0};
				});

	solver.Calculate(0.10);


	return 0;
}