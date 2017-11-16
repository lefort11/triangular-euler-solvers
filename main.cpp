#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "SpaceMesh/Area.h"

#include "FirstOrderSolver/LaxFriedrichSolver.h"
#include "WENO/WENOSolver.h"


int main()
{

	euler::ConstraintFunction circle1([](double t)
									  {
										  auto x = 0.2 * cos(2.0 * M_PI * t);
										  auto y = 0.2 * sin(2.0 * M_PI * t);
										  return GEOM_FADE2D::Point2(x, y);
									  }, 60);

	euler::ConstraintFunction circle2([](double t)
									  {
										  auto x = 6.0 * cos(2.0 * M_PI * t);
										  auto y = 6.0 * sin(2.0 * M_PI * t);
										  return GEOM_FADE2D::Point2(x, y);
									  }, 1000);


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
//				zone2(circle2, true),
				zone3(square, true);

	std::vector<euler::Zone> vZone;
	vZone.push_back(zone);
//	vZone.push_back(zone2);
//	vZone.push_back(zone3);


	std::array<double, 3> trProp = {30, 0.001, 0.093};



	euler::WENOSolver<euler::LaxFriedrichSolver>
			solver(vZone, [](euler::TriangularMesh const& bcmesh, euler::TriangularMesh const& mainMesh)
	{

		for (int triangle_counter = 0; triangle_counter < bcmesh.size(); ++triangle_counter)
		{
			auto const index = bcmesh[triangle_counter]->Index(); //index of the original triangle
/*			bcmesh[triangle_counter]->density = mainMesh[index]->density;
			bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
			bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
			bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure; */

			if ( bcmesh[triangle_counter]->getBarycenter().x() < -2.0 )//left boundary
			{
				bcmesh[triangle_counter]->density = 1.4;
				bcmesh[triangle_counter]->velocityX = 0.9;
				bcmesh[triangle_counter]->velocityY = 0.0;
				bcmesh[triangle_counter]->pressure = 1.0;

			}
			else if ( (bcmesh[triangle_counter]->getBarycenter().x() > 8.0)  )//right, upper and lower boundaries
			{
				bcmesh[triangle_counter]->density = mainMesh[index]->density;
				bcmesh[triangle_counter]->velocityX = 1.01 * std::fabs(mainMesh[triangle_counter]->velocityX);
//				bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
				bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
				bcmesh[triangle_counter]->pressure = 0.9 * mainMesh[index]->pressure;
			}
			else if( bcmesh[triangle_counter]->getBarycenter().y() > 4.0 )
			{
				bcmesh[triangle_counter]->density = mainMesh[index]->density;
				bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
//				bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
				bcmesh[triangle_counter]->velocityY = 1.01 * std::fabs(mainMesh[triangle_counter]->velocityY);
				bcmesh[triangle_counter]->pressure = 0.9 * mainMesh[index]->pressure;

			}
			else if ( bcmesh[triangle_counter]->getBarycenter().y() < -4.0)
			{
				bcmesh[triangle_counter]->density = mainMesh[index]->density;
				bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
//				bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
				bcmesh[triangle_counter]->velocityY = -1.01 * std::fabs(mainMesh[triangle_counter]->velocityY);
				bcmesh[triangle_counter]->pressure = 0.9 * mainMesh[index]->pressure;

			}
			else  // circle ~ wall
			{
				bcmesh[triangle_counter]->density = mainMesh[index]->density;
				bcmesh[triangle_counter]->velocityX = -mainMesh[index]->velocityX;
				bcmesh[triangle_counter]->velocityY = -mainMesh[index]->velocityY;
				bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure;
			} 

/*			if ( bcmesh[triangle_counter]->getBarycenter().x() < -3.0 )//left boundary
			{
				bcmesh[triangle_counter]->density = 1.4;
				bcmesh[triangle_counter]->velocityX = 0.9;
				bcmesh[triangle_counter]->velocityY = 0.0;
				bcmesh[triangle_counter]->pressure = 1.0;

			}
			else if(euler::sqr(bcmesh[triangle_counter]->getBarycenter().x())
					+ euler::sqr(bcmesh[triangle_counter]->getBarycenter().y()) >= 35.0)
			{
				bcmesh[triangle_counter]->density = mainMesh[index]->density;
				bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
				bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
				bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure;
			}
			else  // circle ~ wall
			{
				bcmesh[triangle_counter]->density = mainMesh[index]->density;
				bcmesh[triangle_counter]->velocityX = -mainMesh[index]->velocityX;
				bcmesh[triangle_counter]->velocityY = -mainMesh[index]->velocityY;
				bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure;
			} */


		}



	}, trProp);

	solver.Init([](GEOM_FADE2D::Point2 point)
				{

					return std::array<double, 4>{{1.4, 0.9, 0.0, 1.0}};
//					return std::array<double, 4>{{1.0, 0.9 * sqrt(5.0 / 3.0), 0.0, 1.0}};

//					if(point.x() < 0.0)
//						return std::array<double, 4>{{2.0, 0.0, 0.0, 5.0}};
//					return std::array<double, 4>{{1.0, 0.0, 0.0, 1.0}};
				});

	solver.Calculate(100.0);

	solver.Output("results/density2D.txt", "results/velocity2D.txt", "results/pressure2D.txt");
	solver.ClcOutput("results/test.clc", 0, 0.5, 1);

	double a;
	std::cin >> a;

	return 0;
}
