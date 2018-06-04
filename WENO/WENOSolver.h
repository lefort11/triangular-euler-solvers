#ifndef TRIANGULAR_SOLVERS_WENOLF_H
#define TRIANGULAR_SOLVERS_WENOLF_H

#include <array>
#include "../FirstOrderSolver/LaxFriedrichSolver.h"

//#define MY_STABILITY_FIX 10.0 //100.0, 1e-6

#define CHARACTERISTIC_WISE


namespace euler
{
	template <class T>
	class WENOSolver : public T
	{
	private:

		double const m_eps = 1e-3;

		static int const gaussian_points_number = 6;

		struct FOPolynomialCoeffs
		{
			std::array<double, 3> c = {0.0, 0.0, 0.0};
		};

		struct SOPolynomialCoeffs
		{
			std::array<double, 9> gammas = {0, 0, 0, 0, 0, 0, 0, 0, 0};

			bool weights_to_be_treated = false;
			std::array<double, 9> gammas_plus = {0, 0, 0, 0, 0, 0, 0, 0, 0};
			std::array<double, 9> gammas_minus = {0, 0, 0, 0, 0, 0, 0, 0, 0};

			double sigma_plus = 0.0;
			double sigma_minus = 0.0;

			double const theta = 3.0;


		};

		struct FOReconstructionPolynomial
		{
			std::array<FOPolynomialCoeffs, gaussian_points_number> coeffsAtPoints;
			std::array<int, 3> stencil = {0, 0, 0};
		};

		struct SOReconstructionPolynomial
		{
			std::array<SOPolynomialCoeffs, gaussian_points_number> coeffsAtPoints;

		};

		struct SmoothIndicatorReconstructionData
		{
			std::array<double, 3> alpha = {0.0, 0.0, 0.0};
			std::array<double, 3> beta = {0.0, 0.0, 0.0};
		};

		struct TriangleReconstructionData
		{
			std::array<FOReconstructionPolynomial, 9> fo_polynomial;
			std::array<Point2, gaussian_points_number> gaussian_points;
			SOReconstructionPolynomial so_polynomial;
			std::array<SmoothIndicatorReconstructionData, 9> smoothIndicatorData;
		};

		std::vector<TriangleReconstructionData> m_vReconstructionData;
        std::vector<TriangleReconstructionData> m_vBoundaryReconstructionData;






	public:

		explicit WENOSolver(std::vector<Zone> const &constraints,
						std::function<void(TriangularMesh const&, TriangularMesh const&, double)> const &bcFunc,
						MeshParams const &triangleProp,
						double gamma = 5.0 / 3.0) : T(constraints, bcFunc, triangleProp, gamma)
		{}

		void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunction) override;


	protected:

		Vec4 Reconstruct(Vec4 const &qVec, Triangle const *pTriangle, int edgeNumber,  int gPointNumber) const override;

#ifdef CHARACTERISTIC_WISE
		void FormL(double density, double velX, double velY, double H, arma::mat44& L,
				   std::array<double, 2> const& normal) const;

		void FormR(double density, double velX, double velY, double H, arma::mat44& R,
				   std::array<double, 2> const& normal) const;
#endif

		void GetStencil(Triangle const* pTriangle, std::array<Triangle const*, 10> &stencil) const;

		void CreateBoundingMesh() override;

		void GetSmoothIndicatorData(TriangleReconstructionData& trRecData, Triangle const* pTriangle) const;

		void GetTriangleReconstructionData(TriangleReconstructionData &trRecData, Triangle const* triangle);

		double CalculateKsiAverage(Triangle const* pTriangle, double x_0, double y_0, double h) const;
		double CalculateEtaAverage(Triangle const* pTriangle, double x_0, double y_0, double h) const;
		double CalculateKsiSquareAverage(Triangle const* pTriangle, double x_0, double y_0, double h) const;
		double CalculateEtaSquareAverage(Triangle const* pTriangle, double x_0, double y_0, double h) const;
		double CalculateKsiEtaAverage(Triangle const* pTriangle, double x_0, double y_0, double h) const;



	};


	template<class T>
	inline double WENOSolver<T>::CalculateKsiAverage(Triangle const *pTriangle, double x_0, double y_0, double h) const
	{
		double area = pTriangle->getArea2D();
		auto const x_1 = pTriangle->getCorner(0)->x();
		auto const y_1 = pTriangle->getCorner(0)->y();
		auto const x_2 = pTriangle->getCorner(1)->x();
		auto const y_2 = pTriangle->getCorner(1)->y();
		auto const x_3 = pTriangle->getCorner(2)->x();
		auto const y_3 = pTriangle->getCorner(2)->y();


		return 1.0 / (6.0 * area * h) * ( (y_2 - y_1) * (sqr(x_2 - x_0) + (x_2 - x_0) * (x_1 - x_0) + sqr(x_1 - x_0)) +
				(y_3 - y_2) * (sqr(x_3 - x_0) + (x_3 - x_0) * (x_2 - x_0) + sqr(x_2 - x_0)) +
				(y_1 - y_3) * (sqr(x_1 - x_0) + (x_1 - x_0) * (x_3 - x_0) + sqr(x_3 - x_0)));


	}

	template<class T>
	inline double WENOSolver<T>::CalculateEtaAverage(Triangle const *pTriangle, double x_0, double y_0, double h) const
	{
		double area = pTriangle->getArea2D();
		auto const x_1 = pTriangle->getCorner(0)->x();
		auto const y_1 = pTriangle->getCorner(0)->y();
		auto const x_2 = pTriangle->getCorner(1)->x();
		auto const y_2 = pTriangle->getCorner(1)->y();
		auto const x_3 = pTriangle->getCorner(2)->x();
		auto const y_3 = pTriangle->getCorner(2)->y();


		return -1.0 / (6.0 * area * h) * ( (x_2 - x_1) * (sqr(y_2 - y_0) + (y_2 - y_0) * (y_1 - y_0) + sqr(y_1 - y_0)) +
										  (x_3 - x_2) * (sqr(y_3 - y_0) + (y_3 - y_0) * (y_2 - y_0) + sqr(y_2 - y_0)) +
										  (x_1 - x_3) * (sqr(y_1 - y_0) + (y_1 - y_0) * (y_3 - y_0) + sqr(y_3 - y_0)));

	}

	template<class T>
	inline double WENOSolver<T>::CalculateKsiSquareAverage(Triangle const *pTriangle,
														   double x_0, double y_0, double h) const
	{
		double area = pTriangle->getArea2D();
		auto const x_1 = pTriangle->getCorner(0)->x();
		auto const y_1 = pTriangle->getCorner(0)->y();
		auto const x_2 = pTriangle->getCorner(1)->x();
		auto const y_2 = pTriangle->getCorner(1)->y();
		auto const x_3 = pTriangle->getCorner(2)->x();
		auto const y_3 = pTriangle->getCorner(2)->y();

		return 1.0 / (12.0 * area * sqr(h)) * ( (y_2 - y_1) * (x_2 + x_1 - 2 * x_0) * (sqr(x_2 - x_0) + sqr(x_1 - x_0)) +
				(y_3 - y_2) * (x_3 + x_2 - 2 * x_0) * (sqr(x_3 - x_0) + sqr(x_2 - x_0)) +
				(y_1 - y_3) * (x_1 + x_3 - 2 * x_0) * (sqr(x_1 - x_0) + sqr(x_3 - x_0)));

	}


	template<class T>
	inline double WENOSolver<T>::CalculateEtaSquareAverage(Triangle const *pTriangle,
														   double x_0, double y_0, double h) const
	{
		double area = pTriangle->getArea2D();
		auto const x_1 = pTriangle->getCorner(0)->x();
		auto const y_1 = pTriangle->getCorner(0)->y();
		auto const x_2 = pTriangle->getCorner(1)->x();
		auto const y_2 = pTriangle->getCorner(1)->y();
		auto const x_3 = pTriangle->getCorner(2)->x();
		auto const y_3 = pTriangle->getCorner(2)->y();

		return -1.0 / (12.0 * area * sqr(h)) * ( (x_2 - x_1) * (y_2 + y_1 - 2 * y_0) * (sqr(y_2 - y_0) + sqr(y_1 - y_0)) +
												(x_3 - x_2) * (y_3 + y_2 - 2 * y_0) * (sqr(y_3 - y_0) + sqr(y_2 - y_0)) +
												(x_1 - x_3) * (y_1 + y_3 - 2 * y_0) * (sqr(y_1 - y_0) + sqr(y_3 - y_0)));
	}


	template<class T>
	inline double WENOSolver<T>::CalculateKsiEtaAverage(Triangle const *pTriangle, double x_0, double y_0, double h) const
	{

		double area = pTriangle->getArea2D();
		auto const x_1 = pTriangle->getCorner(0)->x();
		auto const y_1 = pTriangle->getCorner(0)->y();
		auto const x_2 = pTriangle->getCorner(1)->x();
		auto const y_2 = pTriangle->getCorner(1)->y();
		auto const x_3 = pTriangle->getCorner(2)->x();
		auto const y_3 = pTriangle->getCorner(2)->y();

		double result = 0.0;

		std::array<double, 3> delta_x = {
				x_2 - x_1,
				x_3 - x_2,
				x_1 - x_3
		};

		std::array<double, 3> delta_y = {
				y_2 - y_1,
				y_3 - y_2,
				y_1 - y_3
		};

		std::array<double, 3> delta_x_waved = {
				x_1 - x_0,
				x_2 - x_0,
				x_3 - x_0
		};

		std::array<double, 3> delta_y_waved = {
				y_1 - y_0,
				y_2 - y_0,
				y_3 - y_0
		};

		for(int i = 0; i < 3; ++i)
		{
			result += (0.25 * sqr(delta_x[i]) * delta_y[i] + 1.0 / 3.0 * sqr(delta_x[i]) * delta_y_waved[i]
					  + 2.0 / 3.0 * delta_x[i] * delta_y[i] * delta_x_waved[i]
					  + delta_x[i] * delta_x_waved[i] * delta_y_waved[i]
					  + 0.5 * delta_y[i] * sqr(delta_x_waved[i]) + delta_y_waved[i] * sqr(delta_x_waved[i])) * delta_y[i];
		}

		return result * 1.0 / (2.0 * sqr(h) * area);

	}

	template<class T>
	inline void WENOSolver<T>::CreateBoundingMesh()
	{

		for(int triangle_counter = 0; triangle_counter < T::m_triangles.size(); ++triangle_counter)
		{

			for(int edge_number = 0; edge_number < 3; ++edge_number)
			{
				if(T::m_triangles[triangle_counter]->GetOppTriangle(edge_number) == nullptr)
				{


					Triangle* pTriangle = T::m_triangles[triangle_counter];
					pTriangle->SetBoundary(true);

					std::vector<Triangle*> virtualTriangles(7);

					pTriangle->CreateWENOVirtualTriangles(edge_number, virtualTriangles);

					for(int i = 0; i < virtualTriangles.size(); ++i)
					{
						T::m_boundingTriangles.push_back(virtualTriangles[i]);

                        if((virtualTriangles[i] != nullptr) && virtualTriangles[i]->ToBeReconstructed())
                        {
                            virtualTriangles[i]->SetIndex(m_vBoundaryReconstructionData.size());
                            m_vBoundaryReconstructionData.push_back(TriangleReconstructionData());
                        }

					}

				}
			}

		}

		std::vector<GEOM_FADE2D::Triangle2*> kek(T::m_triangles.size());
		std::vector<GEOM_FADE2D::Triangle2*> kekas;

		for(int i = 0; i < kek.size(); ++i)
		{
			kek[i] = dynamic_cast<GEOM_FADE2D::Triangle2*>(T::m_triangles[i]);
		}
		for(int i = 0; i < T::m_boundingTriangles.size(); ++i)
		{
            if(T::m_boundingTriangles[i] != nullptr)
                kekas.push_back(dynamic_cast<GEOM_FADE2D::Triangle2 *>(T::m_boundingTriangles[i]));
		}

		GEOM_FADE2D::Visualizer2 vis("kekas.ps");
		vis.addObject(kek, GEOM_FADE2D::Color(GEOM_FADE2D::CBLACK));
		vis.addObject(kekas, GEOM_FADE2D::Color(GEOM_FADE2D::CRED));

		vis.writeFile();





	}

	template<class T>
	inline void WENOSolver<T>::GetStencil(Triangle const* pTriangle, std::array<Triangle const*, 10> & stencil) const
	{
		stencil[0] = pTriangle;

		stencil[1] = stencil[0]->GetOppTriangle(0);
		stencil[2] = stencil[0]->GetOppTriangle(1);
		stencil[3] = stencil[0]->GetOppTriangle(2);



		stencil[4] = stencil[1]->GetOppTriangle(stencil[1]->getIntraTriangleIndex(stencil[0]->getCorner(2)));
		stencil[5] = stencil[1]->GetOppTriangle(stencil[1]->getIntraTriangleIndex(stencil[0]->getCorner(1)));



		stencil[6] = stencil[3]->GetOppTriangle(stencil[3]->getIntraTriangleIndex(stencil[0]->getCorner(1)));
		stencil[7] = stencil[3]->GetOppTriangle(stencil[3]->getIntraTriangleIndex(stencil[0]->getCorner(0)));


		stencil[8] = stencil[2]->GetOppTriangle(stencil[2]->getIntraTriangleIndex(stencil[0]->getCorner(0)));
		stencil[9] = stencil[2]->GetOppTriangle(stencil[2]->getIntraTriangleIndex(stencil[0]->getCorner(2)));
/*
		if(stencil[5] == stencil[8])
		{
			stencil[8] = stencil[5]->GetOppTriangle(stencil[5]->getIntraTriangleIndex(stencil[0]->getCorner(2)));
		}
		if(stencil[5] == stencil[7])
		{
			stencil[7] = stencil[4]->GetOppTriangle(stencil[4]->getIntraTriangleIndex(stencil[0]->getCorner(1)));
		}
		if(stencil[9] == stencil[6])
		{
			stencil[9] = stencil[6]->GetOppTriangle(stencil[6]->getIntraTriangleIndex(stencil[0]->getCorner(0)));
		}

*/

	}

	template<class T>
	inline void WENOSolver<T>::Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const &initStateFunction)
	{
		T::Init(initStateFunction);

		m_vReconstructionData.resize(T::m_triangles.size());
 //       m_vBoundaryReconstructionData.resize(T::m_boundingTriangles.size());

#pragma omp parallel for
		for(int triangle_number = 0; triangle_number < T::m_triangles.size(); ++triangle_number)
		{
			GetTriangleReconstructionData(m_vReconstructionData[triangle_number], T::m_triangles[triangle_number]);

		}
#pragma omp parallel for
		for(int triangle_number = 0; triangle_number < T::m_boundingTriangles.size(); triangle_number++)
		{

            if((T::m_boundingTriangles[triangle_number] != nullptr) &&
					T::m_boundingTriangles[triangle_number]->ToBeReconstructed())
            {
                GetTriangleReconstructionData(m_vBoundaryReconstructionData[T::m_boundingTriangles[triangle_number]->Index()],
                                              T::m_boundingTriangles[triangle_number]);
            }

		}
	}

	template<class T>
	inline void WENOSolver<T>::GetTriangleReconstructionData(TriangleReconstructionData &trRecData, Triangle const *pTriangle)
	{

		trRecData.fo_polynomial[0].stencil = {0, 1, 2};
		trRecData.fo_polynomial[1].stencil = {0, 2, 3};
		trRecData.fo_polynomial[2].stencil = {0, 3, 1};
		trRecData.fo_polynomial[3].stencil = {0, 3, 6};
		trRecData.fo_polynomial[4].stencil = {0, 3, 7};
		trRecData.fo_polynomial[5].stencil = {0, 1, 4};
		trRecData.fo_polynomial[6].stencil = {0, 1, 5};
		trRecData.fo_polynomial[7].stencil = {0, 2, 8};
		trRecData.fo_polynomial[8].stencil = {0, 2, 9};

		static auto const gaussian_weight = 1.0 / 2.0 + sqrt(3.0) / 6.0;

		auto const first_vertex = pTriangle->getCorner(0);
		auto const second_vertex = pTriangle->getCorner(1);
		auto const third_vertex = pTriangle->getCorner(2);

		std::vector<Point2> gaussian_points(gaussian_points_number);

		gaussian_points[0].x = gaussian_weight * first_vertex->x() + (1 - gaussian_weight) * second_vertex->x();
		gaussian_points[0].y = gaussian_weight * first_vertex->y() + (1 - gaussian_weight) * second_vertex->y();

		gaussian_points[1].x = gaussian_weight * second_vertex->x() + (1 - gaussian_weight) * first_vertex->x();
		gaussian_points[1].y = gaussian_weight * second_vertex->y() + (1 - gaussian_weight) * first_vertex->y();

		gaussian_points[2].x = gaussian_weight * second_vertex->x() + (1 - gaussian_weight) * third_vertex->x();
		gaussian_points[2].y = gaussian_weight * second_vertex->y() + (1 - gaussian_weight) * third_vertex->y();

		gaussian_points[3].x = gaussian_weight * third_vertex->x() + (1 - gaussian_weight) * second_vertex->x();
		gaussian_points[3].y = gaussian_weight * third_vertex->y() + (1 - gaussian_weight) * second_vertex->y();

		gaussian_points[4].x = gaussian_weight * third_vertex->x() + (1 - gaussian_weight) * first_vertex->x();
		gaussian_points[4].y = gaussian_weight * third_vertex->y() + (1 - gaussian_weight) * first_vertex->y();

		gaussian_points[5].x = gaussian_weight * first_vertex->x() + (1 - gaussian_weight) * third_vertex->x();
		gaussian_points[5].y = gaussian_weight * first_vertex->y() + (1 - gaussian_weight) * third_vertex->y();

		std::array<Triangle const*, 10> stencil;

		GetStencil(pTriangle, stencil);


		auto const h = std::sqrt(stencil[0]->getArea2D());
		auto const x_0 = stencil[0]->getBarycenter().x();
		auto const y_0 = stencil[0]->getBarycenter().y();

		double ksi_average[10];
		double eta_average[10];

		double ksi_square_average[10];
		double eta_square_average[10];
		double ksi_eta_average[10];

		for(int i = 0; i < 10; ++i)
		{

			ksi_average[i] = CalculateKsiAverage(stencil[i], x_0, y_0, h);
			eta_average[i] = CalculateEtaAverage(stencil[i], x_0, y_0, h);
			ksi_square_average[i] = CalculateKsiSquareAverage(stencil[i], x_0, y_0, h);
			eta_square_average[i] = CalculateEtaSquareAverage(stencil[i], x_0, y_0, h);
			ksi_eta_average[i] = CalculateKsiEtaAverage(stencil[i], x_0, y_0, h);



		}

		for(int g_point_number = 0; g_point_number < gaussian_points_number; ++g_point_number)
		{
			auto& currGPoint = gaussian_points[g_point_number];
			trRecData.gaussian_points[g_point_number] = currGPoint;

			arma::vec3 b;
			b << 1.0 << (currGPoint.x - x_0) / h << (currGPoint.y - y_0) / h;

			arma::vec4 d;
			d << 1.0 << sqr( (currGPoint.x - x_0) / h ) << sqr( (currGPoint.y - y_0) / h )
			  << (currGPoint.x - x_0) * (currGPoint.y - y_0) / sqr(h);


			arma::mat M(4, 9);
			M.fill(0.0);

			//getting all first order polynomial's coeffs
			for(int polynomial_number = 0; polynomial_number < 9; ++polynomial_number)
			{
				auto& currPolynomial = trRecData.fo_polynomial[polynomial_number];

				auto const ind_0 = currPolynomial.stencil[0];
				auto const ind_1 = currPolynomial.stencil[1];
				auto const ind_2 = currPolynomial.stencil[2];


				arma::mat33 A;
				A << 1.0 << 1.0 << 1.0 << arma::endr
				  << ksi_average[ind_0] << ksi_average[ind_1] << ksi_average[ind_2] << arma::endr
				  << eta_average[ind_0] << eta_average[ind_1] << eta_average[ind_2] << arma::endr;


				arma::vec3 coeffs = arma::solve(A, b);

				currPolynomial.coeffsAtPoints[g_point_number].c[0] = coeffs[0];
				currPolynomial.coeffsAtPoints[g_point_number].c[1] = coeffs[1];
				currPolynomial.coeffsAtPoints[g_point_number].c[2] = coeffs[2];

				M(0, polynomial_number) = 1.0;
				M(1, polynomial_number) = coeffs[0] * ksi_square_average[ind_0] + coeffs[1] * ksi_square_average[ind_1]
										  + coeffs[2] * ksi_square_average[ind_2];

				M(2, polynomial_number) = coeffs[0] * eta_square_average[ind_0] + coeffs[1] * eta_square_average[ind_1]
										  + coeffs[2] * eta_square_average[ind_2];
				M(3, polynomial_number) = coeffs[0] * ksi_eta_average[ind_0] + coeffs[1] * ksi_eta_average[ind_1]
										  + coeffs[2] * ksi_eta_average[ind_2];


			}


/*
            arma::mat B(6, 10);
            B.fill(0.0);
            arma::vec::fixed<6> f;
            for (int j = 0; j < 10; ++j)
            {
                B(0, j) = 1.0;
                B(1, j) = ksi_average[j];
                B(2, j) = eta_average[j];
                B(3, j) = ksi_square_average[j];
                B(4, j) = eta_square_average[j];
                B(5, j) = ksi_eta_average[j];
            }
            f << 1.0 << (currGPoint.x - x_0) / h << (currGPoint.y - y_0) / h
              << sqr((currGPoint.x - x_0) / h) << sqr((currGPoint.y - y_0) / h)
              << (currGPoint.x - x_0) * (currGPoint.y - y_0) / sqr(h);

            arma::vec soCoeffs = arma::solve(B, f);

            //getting gammas

            arma::mat G(10, 9);
            G.fill(0.0);
            for (int j = 0; j < 9; ++j)
                G(0, j) = trRecData.fo_polynomial[j].coeffsAtPoints[g_point_number].c[0];

            G(1, 0) = trRecData.fo_polynomial[0].coeffsAtPoints[g_point_number].c[1];
            G(1, 2) = trRecData.fo_polynomial[2].coeffsAtPoints[g_point_number].c[2];
            G(1, 5) = trRecData.fo_polynomial[5].coeffsAtPoints[g_point_number].c[1];
            G(1, 6) = trRecData.fo_polynomial[6].coeffsAtPoints[g_point_number].c[1];

            G(2, 0) = trRecData.fo_polynomial[0].coeffsAtPoints[g_point_number].c[2];
            G(2, 1) = trRecData.fo_polynomial[1].coeffsAtPoints[g_point_number].c[1];
            G(2, 7) = trRecData.fo_polynomial[7].coeffsAtPoints[g_point_number].c[1];
            G(2, 8) = trRecData.fo_polynomial[8].coeffsAtPoints[g_point_number].c[1];

            G(3, 1) = trRecData.fo_polynomial[1].coeffsAtPoints[g_point_number].c[2];
            G(3, 2) = trRecData.fo_polynomial[2].coeffsAtPoints[g_point_number].c[1];
            G(3, 3) = trRecData.fo_polynomial[3].coeffsAtPoints[g_point_number].c[1];
            G(3, 4) = trRecData.fo_polynomial[4].coeffsAtPoints[g_point_number].c[1];


            G(4, 5) = trRecData.fo_polynomial[5].coeffsAtPoints[g_point_number].c[2];
            G(5, 6) = trRecData.fo_polynomial[6].coeffsAtPoints[g_point_number].c[2];
            G(6, 3) = trRecData.fo_polynomial[3].coeffsAtPoints[g_point_number].c[2];
            G(7, 4) = trRecData.fo_polynomial[4].coeffsAtPoints[g_point_number].c[2];
            G(8, 7) = trRecData.fo_polynomial[7].coeffsAtPoints[g_point_number].c[2];
            G(9, 8) = trRecData.fo_polynomial[8].coeffsAtPoints[g_point_number].c[2];

*/
			//getting gammas
            M *= 10000;
            d *= 10000;
			arma::vec9 gammas = arma::solve(M, d);


			for(int i = 0; i < 9; ++i)
			{
				trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas[i] = gammas[i];
				if(gammas[i] < 0)
					trRecData.so_polynomial.coeffsAtPoints[g_point_number].weights_to_be_treated = true;

			}

			if(trRecData.so_polynomial.coeffsAtPoints[g_point_number].weights_to_be_treated)
			{

				for(int i = 0; i < 9; ++i)
				{
					trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas_plus[i] =
							0.5 * (gammas[i] + trRecData.so_polynomial.coeffsAtPoints[g_point_number].theta *
                                                       std::fabs(gammas[i]));
					trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas_minus[i] =
							trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas_plus[i] - gammas[i];

					trRecData.so_polynomial.coeffsAtPoints[g_point_number].sigma_plus
							+= trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas_plus[i];

					trRecData.so_polynomial.coeffsAtPoints[g_point_number].sigma_minus
							+= trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas_minus[i];


				}

				for(int i = 0; i < 9; ++i)
				{

					trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas_plus[i]
							/= trRecData.so_polynomial.coeffsAtPoints[g_point_number].sigma_plus;

					trRecData.so_polynomial.coeffsAtPoints[g_point_number].gammas_minus[i]
							/= trRecData.so_polynomial.coeffsAtPoints[g_point_number].sigma_minus;


				}

			}


		}

		GetSmoothIndicatorData(trRecData, pTriangle);



	}


	template<class T>
	inline void WENOSolver<T>::GetSmoothIndicatorData(TriangleReconstructionData &trRecData, Triangle const *pTriangle) const
	{
		std::array<Triangle const*, 10> stencil;
		GetStencil(pTriangle, stencil);

		auto const h = pTriangle->getArea2D();
		auto const x_0 = pTriangle->getBarycenter().x();
		auto const y_0 = pTriangle->getBarycenter().y();

		double ksi_average[10];
		double eta_average[10];


		for(int i = 0; i < 10; ++i)
		{
			ksi_average[i] = CalculateKsiAverage(stencil[i], x_0, y_0, h);
			eta_average[i] = CalculateEtaAverage(stencil[i], x_0, y_0, h);

		}


		for(int polynomial_number = 0; polynomial_number < 9; ++polynomial_number)
		{
			arma::mat99 A;
			arma::vec9 c;

			auto const ind_0 = trRecData.fo_polynomial[polynomial_number].stencil[0];
			auto const ind_1 = trRecData.fo_polynomial[polynomial_number].stencil[1];
			auto const ind_2 = trRecData.fo_polynomial[polynomial_number].stencil[2];

			int i = 0;
			for (int g_point_number = 0; g_point_number < gaussian_points_number; g_point_number+=2)
			{
				//@todo check +=2!
				auto const c_0 = trRecData.fo_polynomial[polynomial_number].coeffsAtPoints[g_point_number].c[0];
				auto const c_1 = trRecData.fo_polynomial[polynomial_number].coeffsAtPoints[g_point_number].c[1];
				auto const c_2 = trRecData.fo_polynomial[polynomial_number].coeffsAtPoints[g_point_number].c[2];

				auto const ksi_G = (trRecData.gaussian_points[g_point_number].x - x_0) / h;
				auto const eta_G = (trRecData.gaussian_points[g_point_number].y - y_0) / h;


				//first row, u = 1;
				A(i, 0) = ksi_G;
				A(i, 1) = ksi_G;
				A(i, 2) = ksi_G;

				A(i, 3) = eta_G;
				A(i, 4) = eta_G;
				A(i, 5) = eta_G;

				A(i, 6) = 1.0;
				A(i, 7) = 1.0;
				A(i, 8) = 1.0;

				c(i) = 1.0;

				//second row, u = ksi
				A(i + 1, 0) = ksi_average[ind_0] * ksi_G;
				A(i + 1, 1) = ksi_average[ind_1] * ksi_G;
				A(i + 1, 2) = ksi_average[ind_2] * ksi_G;

				A(i + 1, 3) = ksi_average[ind_0] * eta_G;
				A(i + 1, 4) = ksi_average[ind_1] * eta_G;
				A(i + 1, 5) = ksi_average[ind_2] * eta_G;

				A(i + 1, 6) = ksi_average[ind_0];
				A(i + 1, 7) = ksi_average[ind_1];
				A(i + 1, 8) = ksi_average[ind_2];
				c(i + 1) = ksi_G;

				//third row, u = eta;
				A(i + 2, 0) = eta_average[ind_0] * ksi_G;
				A(i + 2, 1) = eta_average[ind_1] * ksi_G;
				A(i + 2, 2) = eta_average[ind_2] * ksi_G;

				A(i + 2, 3) = eta_average[ind_0] * eta_G;
				A(i + 2, 4) = eta_average[ind_1] * eta_G;
				A(i + 2, 5) = eta_average[ind_2] * eta_G;

				A(i + 2, 6) = eta_average[ind_0];
				A(i + 2, 7) = eta_average[ind_1];
				A(i + 2, 8) = eta_average[ind_2];
				c(i + 2) = eta_G;



				i+=3;
			}


 //           assert(arma::rank(A) == 9);
			arma::vec9 solution = arma::solve(A, c);

			trRecData.smoothIndicatorData[polynomial_number].alpha[0] = solution[0];
			trRecData.smoothIndicatorData[polynomial_number].alpha[1] = solution[1];
			trRecData.smoothIndicatorData[polynomial_number].alpha[2] = solution[2];

			trRecData.smoothIndicatorData[polynomial_number].beta[0] = solution[3];
			trRecData.smoothIndicatorData[polynomial_number].beta[1] = solution[4];
			trRecData.smoothIndicatorData[polynomial_number].beta[2] = solution[5];



		}
	}

	template<class T>
	inline Vec4 WENOSolver<T>::Reconstruct(Vec4 const& qVec, Triangle const *pTriangle, int edgeNumber,
										   int gPointNumber) const
	{
		std::array<Triangle const*, 10> stencil;
		GetStencil(pTriangle, stencil);


		std::array<Vec4, 10> q;

		q[0] = qVec;

#ifdef MY_STABILITY_FIX
        auto max_norm = arma::norm(q[0], 2);
#endif


		for(int i = 1; i < 10; ++i)
		{
			T::FormQVector(q[i], stencil[i]);
#ifdef MY_STABILITY_FIX
            auto const norm = arma::norm(q[i], 2);
			if(norm > max_norm)
				max_norm = norm;
#endif


		}

#ifdef MY_STABILITY_FIX
		for(int i = 0; i < 10; ++i)
		{
			q[i] /= MY_STABILITY_FIX * max_norm;
		}
#endif



		//Reconstruction!

		Vec4 q_reconstructed{0.0, 0.0, 0.0, 0.0};

#ifdef CHARACTERISTIC_WISE


		auto const density_minus = qVec[0];
		auto const velocityX_minus = qVec[1] / density_minus;
		auto const velocityY_minus = qVec[2] / density_minus;
		auto const E_minus = qVec[3] / density_minus;
		auto const velocity_sqr_abs_minus = sqr(velocityX_minus) + sqr(velocityY_minus);
		auto const eps_minus = E_minus - 0.5 * velocity_sqr_abs_minus;
		auto const pressure_minus = (T::m_gamma - 1.0) * eps_minus * density_minus;
		auto const H_minus = eps_minus + pressure_minus / density_minus + 0.5 * velocity_sqr_abs_minus;

		auto const neighbour_triangle = pTriangle->GetOppTriangle(edgeNumber);
		auto const density_plus = neighbour_triangle->density;
		auto const velocityX_plus = neighbour_triangle->velocityX;
		auto const velocityY_plus = neighbour_triangle->velocityY;
		auto const pressure_plus = neighbour_triangle->pressure;
		auto const velocity_sqr_abs_plus = sqr(velocityX_plus) + sqr(velocityY_plus);
		auto const eps_plus = pressure_plus / (density_plus * (T::m_gamma - 1.0));
		auto const H_plus = eps_plus + pressure_plus / density_plus + 0.5 * velocity_sqr_abs_plus;

		auto const density_star = std::sqrt(density_minus * density_plus);

		auto const velocityX_star = (std::sqrt(density_minus) * velocityX_minus
									 + std::sqrt(density_plus) * velocityX_plus)
									/ (std::sqrt(density_minus) + std::sqrt(density_plus));

		auto const velocityY_star = (std::sqrt(density_minus) * velocityY_minus
									 + std::sqrt(density_plus) * velocityY_plus)
									/ (std::sqrt(density_minus) + std::sqrt(density_plus));

		auto const H_star = (std::sqrt(density_minus) * H_minus
							 + std::sqrt(density_plus) * H_plus)
							/ (std::sqrt(density_minus) + std::sqrt(density_plus));

        arma::mat44 R;
		arma::mat44 L;

		auto const normal = pTriangle->CalculateNormal(edgeNumber);
		FormL(density_star, velocityX_star, velocityY_star, H_star, L, normal);
		FormR(density_star, velocityX_star, velocityY_star, H_star, R, normal);

		for(int i = 0; i < 10; ++i)
		{
			q[i] = L * q[i];
		}


#endif

		TriangleReconstructionData const& triangleRecData =
				pTriangle->IsVirtual()? m_vBoundaryReconstructionData[pTriangle->Index()] :
									 m_vReconstructionData[pTriangle->Index()];




		int current_g_n = (2 * (edgeNumber + 1) + gPointNumber) % 6;


		std::array<Vec4, 9> omega_waved;
		Vec4 omega_waved_sum{0.0, 0.0, 0.0, 0.0};

		std::array<Vec4, 9> omega_waved_plus;
		Vec4 omega_waved_plus_sum{0.0, 0.0, 0.0, 0.0};
		std::array<Vec4, 9> omega_waved_minus;
		Vec4 omega_waved_minus_sum{0.0, 0.0, 0.0, 0.0};

//        double area = pTriangle->getArea2D();
        //double eps0 = 4 * sqr(std::min(T::m_meshProperties.gridLength, T::m_meshProperties.maxEdgeLength)) / 144;
        //double eps0 = std::sqrt(2 * pTriangle->getArea2D()) / 144;


        //double eps0 = std::max(std::min(T::m_max_area / T::m_total_area, 1e-2), 1e-6);

        double eps0 = 2.0 * pTriangle->getArea2D() / T::m_total_area;
/*        auto const barycenter = pTriangle->getBarycenter();
		if(!((barycenter.x() > -2.0) && (barycenter.x() < 7.3) && (std::fabs(barycenter.y()) < 5.0)))
			eps0 = 1e-3; */



        Vec4 const eps{eps0, eps0, eps0, eps0};

		auto const weights_to_be_treated =
                triangleRecData.so_polynomial.coeffsAtPoints[current_g_n].weights_to_be_treated;



		for(int polynom_num = 0; polynom_num < 9; ++polynom_num)
		{

			auto const ind_0 = triangleRecData.fo_polynomial[polynom_num].stencil[0];
			auto const ind_1 = triangleRecData.fo_polynomial[polynom_num].stencil[1];
			auto const ind_2 = triangleRecData.fo_polynomial[polynom_num].stencil[2];

			SmoothIndicatorReconstructionData const& smIndData = triangleRecData.smoothIndicatorData[polynom_num];

			Vec4 smoothIndicator = 1.0 / pTriangle->getArea2D() * (arma::square(smIndData.alpha[0] * q[ind_0]
										   + smIndData.alpha[1] * q[ind_1]
										   + smIndData.alpha[2] * q[ind_2]) +
							  arma::square(smIndData.beta[0] * q[ind_0]
										   + smIndData.beta[1] * q[ind_1]
										   + smIndData.beta[2] * q[ind_2]));


			if(!weights_to_be_treated)
			{
                auto const gamma = triangleRecData.so_polynomial.coeffsAtPoints[current_g_n].gammas[polynom_num];
				omega_waved[polynom_num] = Vec4{gamma, gamma, gamma, gamma} / arma::square(eps + smoothIndicator); //changed 1.5


				omega_waved_sum += omega_waved[polynom_num];
			}
			else
			{
                auto const gamma_plus = triangleRecData.so_polynomial.coeffsAtPoints[current_g_n].gammas_plus[polynom_num];
                auto const gamma_minus = triangleRecData.so_polynomial.coeffsAtPoints[current_g_n].gammas_minus[polynom_num];

				omega_waved_plus[polynom_num] = Vec4{gamma_plus, gamma_plus, gamma_plus, gamma_plus}
												/ arma::square(eps + smoothIndicator);


				omega_waved_plus_sum += omega_waved_plus[polynom_num];

				omega_waved_minus[polynom_num] = Vec4{gamma_minus, gamma_minus, gamma_minus, gamma_minus}
                                                    / arma::square(eps + smoothIndicator);

				omega_waved_minus_sum += omega_waved_minus[polynom_num];

			}

		}


		for(int polynom_num = 0; polynom_num < 9; ++polynom_num)
		{
			FOReconstructionPolynomial const& current_polynomial = triangleRecData.fo_polynomial[polynom_num];

			auto const ind_0 = current_polynomial.stencil[0];
			auto const ind_1 = current_polynomial.stencil[1];
			auto const ind_2 = current_polynomial.stencil[2];

			auto const c_0 = current_polynomial.coeffsAtPoints[current_g_n].c[0];
			auto const c_1 = current_polynomial.coeffsAtPoints[current_g_n].c[1];
			auto const c_2 = current_polynomial.coeffsAtPoints[current_g_n].c[2];


			if(!weights_to_be_treated)
			{
				Vec4 const omega = omega_waved[polynom_num] / omega_waved_sum;

				q_reconstructed += omega % (c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2]);



            }
			else
			{
				Vec4 const omega_plus = omega_waved_plus[polynom_num] / omega_waved_plus_sum;
				Vec4 const omega_minus = omega_waved_minus[polynom_num] / omega_waved_minus_sum;


				Vec4 const p = c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2];


				q_reconstructed += triangleRecData.so_polynomial.coeffsAtPoints[current_g_n].sigma_plus * omega_plus % p -
								   triangleRecData.so_polynomial.coeffsAtPoints[current_g_n].sigma_minus * omega_minus % p;
			}


		}




#ifdef MY_STABILITY_FIX
        q_reconstructed *= MY_STABILITY_FIX * max_norm;
#endif


#ifdef CHARACTERISTIC_WISE

		q_reconstructed = R * q_reconstructed;
#endif

		if(!((q_reconstructed[0] > 0) && (q_reconstructed[3] > 0)))
		{
			std::cout << "kek" << std::endl;
            T::ClcOutput("results/fail.clc", 0, 0.5, 1);
			throw 1;
		}


		return q_reconstructed;



	}




#ifdef CHARACTERISTIC_WISE

	template<class T>
	inline void WENOSolver<T>::FormL(double density, double velocityX, double velocityY,
									   double H, arma::mat44 &L,  std::array<double, 2> const& normal) const
	{

		auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);
		auto const c = std::sqrt((T::m_gamma - 1) * (H - 0.5 * velocity_sqr_abs));

		auto const vel_n_x = normal[0] * velocityX + normal[1] * velocityY;

		L.row(0) = arma::rowvec::fixed<4>{ ( (T::m_gamma - 1) * 0.5 * velocity_sqr_abs + c * vel_n_x ) / (2 * sqr(c)),
						 ( (1 - T::m_gamma) * velocityX - c * normal[0] ) / (2 * sqr(c)),
						 ( (1 - T::m_gamma) * velocityY - c * normal[1] ) / (2 * sqr(c)),
						 (T::m_gamma - 1) / (2 * sqr(c))
		};
		L.row(1) = arma::rowvec::fixed<4>{ ( sqr(c) - (T::m_gamma - 1) * 0.5 * velocity_sqr_abs ) / sqr(c),
						 (T::m_gamma - 1) * velocityX / sqr(c),
						 (T::m_gamma - 1) * velocityY / sqr(c),
						 (1 - T::m_gamma) / sqr(c)
		};
		L.row(2) = arma::rowvec::fixed<4>{ (velocityY * normal[0] - velocityX * normal[1]),
						 normal[1],
						 -normal[0],
						 0.0
		};
		L.row(3) = arma::rowvec::fixed<4>{ ( (T::m_gamma - 1) * 0.5 * velocity_sqr_abs - c * vel_n_x ) / (2 * sqr(c)),
						 ( (1 - T::m_gamma) * velocityX + c * normal[0] ) / (2 * sqr(c)),
						 ( (1 - T::m_gamma) * velocityY + c * normal[1] ) / (2 * sqr(c)),
						 (T::m_gamma - 1) / (2 * sqr(c))
		};




	}



	template<class T>
	inline void WENOSolver<T>::FormR(double density, double velocityX, double velocityY,
									   double H, arma::mat44 &R, std::array<double, 2> const& normal) const
	{

		auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);
		auto const c = std::sqrt((T::m_gamma - 1) * (H - 0.5 * velocity_sqr_abs));

		auto const vel_n_x = normal[0] * velocityX + normal[1] * velocityY;

		R.col(0) = Vec4{1.0,
						velocityX - c * normal[0],
						velocityY - c * normal[1],
						H - c * vel_n_x
		};
		R.col(1) = Vec4{1.0,
						velocityX,
						velocityY,
						0.5 * velocity_sqr_abs
		};
		R.col(2) = Vec4{0.0,
						normal[1],
						-normal[0],
						velocityX * normal[1] - velocityY * normal[0]
		};
		R.col(3) = Vec4{1.0,
						velocityX + c * normal[0],
						velocityY + c * normal[1],
						H + c * vel_n_x
		};

	}



#endif

}








#endif //TRIANGULAR_SOLVERS_WENOLF_H
