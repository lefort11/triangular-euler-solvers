#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "SpaceMesh/Area.h"

#include "FirstOrderSolver/LaxFriedrichSolver.h"
#include "WENO/WENOSolver.h"
#include "FirstOrderSolver/RoeSolver.h"
#include "HYBRID/HYBRIDSolver.h"


int main()
{

	euler::ConstraintFunction circle1([](double t)
									  {
										  auto x = 0.2 * cos(2.0 * M_PI * t);
										  auto y = 0.2 * sin(2.0 * M_PI * t);
										  return GEOM_FADE2D::Point2(x, y);
									  }, 50); //0,2 50

	euler::ConstraintFunction circle2([](double t)
									  {
										  auto x = 6.0 * cos(2.0 * M_PI * t);
										  auto y = 6.0 * sin(2.0 * M_PI * t);
										  return GEOM_FADE2D::Point2(x, y);
									  }, 200);


	euler::ConstraintFunction square([](double t)
									 {
										 double x, y;
										 if((t >= 0) && (t < 0.25))
										 {
											 x = X_LEFT;
											 //y = -4 * t + 1;
                                             y = Y_TOP + (Y_BOT - Y_TOP) * 4 * t;
                                            // y = -Y_TOP * (t - 0.25) / 0.25 + Y_BOT * t / 0.25;
										 }
										 else if((t >= 0.25) && (t < 0.5))
										 {
											 y = Y_BOT;
											 //x = 4*(t - 0.25) - 1;
                                             x = X_LEFT + (X_RIGHT - X_LEFT) * 4 * (t - 0.25);
										 }
										 else if((t >= 0.5) && (t < 0.75))
										 {
											 x = X_RIGHT;
											// y = 4*(t - 0.5) - 1;
                                             y = Y_BOT + (Y_TOP - Y_BOT) * 4 * (t - 0.5);
										 }
										 else
										 {
											 y = Y_TOP;
											 //x = -4*(t - 0.75) + 1;
                                             x = X_RIGHT + (X_LEFT - X_RIGHT) * 4 * (t - 0.75);
										 }

										 return GEOM_FADE2D::Point2(x, y);


									 }, 4);
    euler::ConstraintFunction smallsquare([](double t)
                                     {
                                         double x, y;
                                         if((t >= 0) && (t < 0.25))
                                         {
                                             x = -0.2;
                                             //y = -4 * t + 1;
                                             y = 0.2 + (-0.2 - 0.2) * 4 * t;
                                             // y = -Y_TOP * (t - 0.25) / 0.25 + Y_BOT * t / 0.25;
                                         }
                                         else if((t >= 0.25) && (t < 0.5))
                                         {
                                             y = -0.2;
                                             //x = 4*(t - 0.25) - 1;
                                             x = -0.2 + (0.2 + 0.2) * 4 * (t - 0.25);
                                         }
                                         else if((t >= 0.5) && (t < 0.75))
                                         {
                                             x = 0.2;
                                             // y = 4*(t - 0.5) - 1;
                                             y = -0.2 + (0.2 + 0.2) * 4 * (t - 0.5);
                                         }
                                         else
                                         {
                                             y = 0.2;
                                             //x = -4*(t - 0.75) + 1;
                                             x = 0.2 + (-0.2 - 0.2) * 4 * (t - 0.75);
                                         }

                                         return GEOM_FADE2D::Point2(x, y);


                                     }, 4);

    euler::ConstraintFunction wing([](double t)
                                   {
                                       double const x0 = -0.05;
                                       double const y0 = 0.05;
                                       double const r = 1.06;
                                       t *= 2*M_PI;
                                       double abs = euler::sqr(x0) + euler::sqr(y0) + euler::sqr(r) + 2*x0*r*cos(t) + 2*y0*r*sin(t);
                                       double x = 0.5 *(x0 + r*cos(t) + (x0 + r*cos(t))/abs) + 0.8;
                                       double y = 0.5 *(y0 + r*sin(t) - (y0 + r*sin(t))/abs) * 1.1;

                                       double const theta = -M_PI / 90;
                                       return GEOM_FADE2D::Point2(x*cos(theta) - y*sin(theta), x*sin(theta) + y * cos(theta));
                                   }, 35);

    euler::ConstraintFunction wing2([](double t)
                                    {
                                        double const x0 = -0.05;
                                        double const y0 = 0.00;
                                        double const r = 1.07;
                                        t *= 2*M_PI;
                                        double abs = euler::sqr(x0) + euler::sqr(y0) + euler::sqr(r) + 2*x0*r*cos(t) + 2*y0*r*sin(t);
                                        double x = 0.5 *(x0 + r*cos(t) + (x0 + r*cos(t))/abs) + 0.8;
                                        double y = 0.5 *(y0 + r*sin(t) - (y0 + r*sin(t))/abs) * 1.7;

                                        double const theta = -M_PI / 25;
                                        return GEOM_FADE2D::Point2(x*cos(theta) - y*sin(theta), x*sin(theta) + y * cos(theta));
                                    }, 35);



    euler::Zone smallCircle(circle1, false),
//				zone2(circle2, true),
				hugeSquare(square, true),
                smallSquare(smallsquare, false),
                wingZone(wing, false);

	std::vector<euler::Zone> vZone;
    vZone.push_back(hugeSquare);
    vZone.push_back(smallCircle);

    euler::MeshParams meshParams;
    meshParams.minAngleDegree = 30.0;
    meshParams.minEdgeLength = 0.0008;
    meshParams.maxEdgeLength = 0.15;
	meshParams.growFactor = 2.0;
    meshParams.capAspectLimit = 2.0;
	meshParams.gridVector = GEOM_FADE2D::Vector2(1.0, 0.0);
//    meshParams.gridLength = 1.5 * meshParams.maxEdgeLength;
    meshParams.gridLength = 0.0423121;

/*
	euler::MeshParams meshParams;
	meshParams.minAngleDegree = 30.0;
	meshParams.minEdgeLength = 0.001;
	meshParams.maxEdgeLength = 5.0;
//	meshParams.growFactor = 2.0;
    meshParams.capAspectLimit = 2.0;
//	meshParams.gridVector = GEOM_FADE2D::Vector2(1.0, 1.0);
	meshParams.gridLength = 0.12;
*/



	euler::HYBRIDSolver
			solver(vZone, [](euler::TriangularMesh const& bcmesh, euler::TriangularMesh const& mainMesh, double time)
	{

		for (int triangle_counter = 0; triangle_counter < bcmesh.size(); ++triangle_counter)
		{

            if(bcmesh[triangle_counter] != nullptr)
            {
                auto const index = bcmesh[triangle_counter]->ParentIndex(); //index of the original triangle

                if (bcmesh[triangle_counter]->getBarycenter().x() <= X_LEFT)//left boundary
                {
                    bcmesh[triangle_counter]->density = 1.4;
                    bcmesh[triangle_counter]->velocityX = 0.9;
                    bcmesh[triangle_counter]->velocityY = 0.0;
                    bcmesh[triangle_counter]->pressure = 1.0;

                } else if ((bcmesh[triangle_counter]->getBarycenter().x() >= X_RIGHT))//right, upper and lower boundaries
                {
                    bcmesh[triangle_counter]->density = mainMesh[index]->density;
                    bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
                    bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
                    bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure;

                } else if (bcmesh[triangle_counter]->getBarycenter().y() >= Y_TOP)
                {
                    bcmesh[triangle_counter]->density = mainMesh[index]->density;
                    bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
                    bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
                    bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure;


                } else if (bcmesh[triangle_counter]->getBarycenter().y() <= Y_BOT)
                {
                    bcmesh[triangle_counter]->density = mainMesh[index]->density;
                    bcmesh[triangle_counter]->velocityX = mainMesh[index]->velocityX;
                    bcmesh[triangle_counter]->velocityY = mainMesh[index]->velocityY;
                    bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure;
                } else  // circle ~ wall
                {
                    bcmesh[triangle_counter]->density = mainMesh[index]->density;
                    bcmesh[triangle_counter]->velocityX = -mainMesh[index]->velocityX;
                    bcmesh[triangle_counter]->velocityY = -mainMesh[index]->velocityY;
                    bcmesh[triangle_counter]->pressure = mainMesh[index]->pressure;
                }
            }




		}



	}, meshParams, 1.4);

	solver.Init([](GEOM_FADE2D::Point2 point)
				{

					return std::array<double, 4>{{1.4, 0.9, 0.0, 1.0}};
//					return std::array<double, 4>{{1.0, 0.9 * sqrt(5.0 / 3.0), 0.0, 1.0}};

//					if(point.x() < 0.0)
//						return std::array<double, 4>{{2.0, 0.0, 0.0, 5.0}};
//					return std::array<double, 4>{{1.0, 0.0, 0.0, 1.0}};
				});

	solver.Calculate(60.0);

	solver.Output("results/density2D.txt", "results/velocity2D.txt", "results/pressure2D.txt");
	solver.ClcOutput("results/test.clc", 0, 0.5, 1);


	return 0;
}
