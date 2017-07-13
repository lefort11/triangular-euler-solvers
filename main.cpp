#include <iostream>
#include <cmath>

#include "SpaceMesh/Area.h"



int main()
{


	euler::ConstraintFunction kek([](double t)
	{
		auto x =cos(2*M_PI * t);
		auto y = sin(2*M_PI * t);
		return GEOM_FADE2D::Point2(x, y);
	});

	euler::ConstraintFunction keks([](double t)
	{
		auto x = 1.1 * cos(2*M_PI * t);
		auto y = 1.1 * sin(2*M_PI * t);
		return GEOM_FADE2D::Point2(x, y);
	});

	euler::Zone zone(kek, false, {0.0, 0.0, 0.0, 0.0}),
				zone1(keks, true, {0.0, 0.0, 0.0, 0.0});

	std::vector<euler::Zone> vZone;
	vZone.push_back(zone);
	//vZone.push_back(zone1);

	euler::Area area(vZone);

	std::array<double, 3> trProp = {10, 0.01, 1};
	area.Triangulate(1000, trProp, [](GEOM_FADE2D::Point2){
		return std::array<double,4>{0.0,0.0,0.0,0.0};
	});

	return 0;
}