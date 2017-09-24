#ifndef TRIANGULAR_SOLVERS_WENOLF_H
#define TRIANGULAR_SOLVERS_WENOLF_H

#include "../FirstOrderSolver/LaxFriedrichSolver.h"

//#define CHARACTERISTIC_WISE

#define MY_STABILITY_FIX 100.0 //100.0, 1e-6

namespace euler
{
	template <class T>
	class WENOSolver : public T
	{
	private:

		double const m_eps = 1e-6;

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






	public:

		explicit WENOSolver(std::vector<Zone> const &constraints,
						std::function<void(TriangularMesh const&, TriangularMesh const&)> const &bcFunc,
						std::array<double, 3> const &triangleProp = {0.0, 0.0, 0.0},
						double gamma = 5.0 / 3.0) : T(constraints, bcFunc, triangleProp, gamma)
		{}

		void Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunction) override;


	protected:

		Vec4 Reconstruct(Vec4 const &qVec, Triangle const *pTriangle,
						 Point2 const& gaussianPoint, int edgeNumber) const override;


		void FormL_A(double density, double velX, double velY, double H, arma::mat44& L_A) const;
		void FormL_B(double density, double velX, double velY, double H, arma::mat44& L_B) const;

		void FormR_A(double density, double velX, double velY, double H, arma::mat44& R_A) const;
		void FormR_B(double density, double velX, double velY, double H, arma::mat44& R_B) const;

		void GetStencil(Triangle const* pTriangle, std::array<Triangle const*, 10> &stencil) const;

		void CreateBoundingMesh() override;

		void GetSmoothIndicatorData(TriangleReconstructionData& trRecData, Triangle const* pTriangle) const;

		void GetTriangleReconstructionData(TriangleReconstructionData &trRecData, Triangle const* triangle);

	};

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

					auto const reflectedTriangle0 = pTriangle->ReflectTriangle(edge_number);
					T::m_boundingTriangles.push_back(reflectedTriangle0);

					auto const ind_1 =
							reflectedTriangle0->getIntraTriangleIndex(pTriangle->getCorner((edge_number + 1) % 3));
					auto const ind_2 =
							reflectedTriangle0->getIntraTriangleIndex(pTriangle->getCorner((edge_number + 2) % 3));

					auto const reflectedTriangle1 = reflectedTriangle0->ReflectTriangle(ind_1);
					auto const reflectedTriangle2 = reflectedTriangle0->ReflectTriangle(ind_2);
					T::m_boundingTriangles.push_back(reflectedTriangle1);
					T::m_boundingTriangles.push_back(reflectedTriangle2);


/*					auto const p1_1 = reflectedTriangle0->getCorner((ind_1 + 1) % 3);
					auto const p2_1 = reflectedTriangle0->getCorner((ind_1 + 2) % 3);
					auto const ind_1_1 = reflectedTriangle1->getIntraTriangleIndex(p1_1);
					auto const ind_2_1 = reflectedTriangle1->getIntraTriangleIndex(p2_1);
					auto const reflectedTriangle3 = reflectedTriangle1->ReflectTriangle(ind_1_1);
					auto const reflectedTriangle4 = reflectedTriangle1->ReflectTriangle(ind_2_1);
					T::m_boundingTriangles.push_back(reflectedTriangle3);
					T::m_boundingTriangles.push_back(reflectedTriangle4);

					auto const p1_2 = reflectedTriangle0->getCorner((ind_2 + 1) % 3);
					auto const p2_2 = reflectedTriangle0->getCorner((ind_2 + 2) % 3);
					auto const ind_1_2 = reflectedTriangle2->getIntraTriangleIndex(p1_2);
					auto const ind_2_2 = reflectedTriangle2->getIntraTriangleIndex(p2_2);
					auto const reflectedTriangle5 = reflectedTriangle2->ReflectTriangle(ind_1_2);
					auto const reflectedTriangle6 = reflectedTriangle2->ReflectTriangle(ind_2_2);
					T::m_boundingTriangles.push_back(reflectedTriangle5);
					T::m_boundingTriangles.push_back(reflectedTriangle6);

*/

				}
			}


		}

	}

	template<class T>
	inline void WENOSolver<T>::GetStencil(Triangle const* pTriangle, std::array<Triangle const*, 10> & stencil) const
	{
		stencil[0] = pTriangle;

		stencil[1] = stencil[0]->GetOppTriangle(0);
		stencil[2] = stencil[0]->GetOppTriangle(1);
		stencil[3] = stencil[0]->GetOppTriangle(2);
//		if((stencil[0] == nullptr) || (stencil[2] == nullptr) || stencil[3] == nullptr)
//			return;


		stencil[4] = stencil[1]->GetOppTriangle(stencil[1]->getIntraTriangleIndex(stencil[0]->getCorner(2)));

		stencil[5] = stencil[1]->GetOppTriangle(stencil[1]->getIntraTriangleIndex(stencil[0]->getCorner(1)));


		stencil[6] = stencil[3]->GetOppTriangle(stencil[3]->getIntraTriangleIndex(stencil[0]->getCorner(1)));

		stencil[7] = stencil[3]->GetOppTriangle(stencil[3]->getIntraTriangleIndex(stencil[0]->getCorner(0)));


		stencil[8] = stencil[2]->GetOppTriangle(stencil[2]->getIntraTriangleIndex(stencil[0]->getCorner(0)));

		stencil[9] = stencil[2]->GetOppTriangle(stencil[2]->getIntraTriangleIndex(stencil[0]->getCorner(2)));


	}

	template<class T>
	inline void WENOSolver<T>::Init(std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const &initStateFunction)
	{
		T::Init(initStateFunction);

		m_vReconstructionData.resize(T::m_triangles.size());

#pragma omp parallel for
		for(int triangle_number = 0; triangle_number < T::m_triangles.size(); ++triangle_number)
		{
			GetTriangleReconstructionData(m_vReconstructionData[triangle_number], T::m_triangles[triangle_number]);

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
		for(int i = 0; i < 10; ++i)
			if(stencil[i] == nullptr)
				return;

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
			auto const area = stencil[i]->getArea2D();

			auto const x = stencil[i]->getBarycenter().x();
			auto const y = stencil[i]->getBarycenter().y();

			ksi_average[i] = (x - x_0) / h;
			eta_average[i] = (y - y_0) / h;

			ksi_square_average[i] = sqr(x - x_0) / sqr(h);
			eta_square_average[i] = sqr(y - y_0) / sqr(h);
			ksi_eta_average[i] = (x - x_0) * (y - y_0) / sqr(h);

			/*ksi_average[i] = 1.0 / area * T::GaussianIntegration([h, x_0](double x, double y)
															{
																return (x - x_0) / h;
															}, stencil[i]);

			eta_average[i] = 1.0 / area * T::GaussianIntegration([h, y_0](double x, double y)
															{
																return (y - y_0) / h;
															}, stencil[i]);

			ksi_square_average[i] = 1.0 / area * T::GaussianIntegration([h, x_0](double x, double y)
																   {
																	   return sqr(x-x_0) / sqr(h);
																   }, stencil[i]);
			eta_square_average[i] = 1.0 / area * T::GaussianIntegration([h, y_0](double x, double y)
																   {
																	   return sqr(y-y_0) / sqr(h);
																   }, stencil[i]);
			ksi_eta_average[i] = 1.0 / area * T::GaussianIntegration([h, x_0, y_0](double x, double y)
																{
																	return (x - x_0) * (y - y_0)  / sqr(h);
																}, stencil[i]); */




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
			M.fill(1.0);



			for(int polynomial_number = 0; polynomial_number < 9; ++polynomial_number)
			{
				auto& currPolynomial = trRecData.fo_polynomial[polynomial_number];

				auto const ind_0 = currPolynomial.stencil[0];
				auto const ind_1 = currPolynomial.stencil[1];
				auto const ind_2 = currPolynomial.stencil[2];


				arma::mat33 A;
				A << 1.0 << 1.0 << 1.0 << arma::endr
				  << ksi_average[ind_0] << ksi_average[ind_1]
				  << ksi_average[ind_2] << arma::endr
				  << eta_average[ind_0] << eta_average[ind_1]
				  << eta_average[ind_2] << arma::endr;

				arma::vec3 coeffs = arma::solve(A, b);

				currPolynomial.coeffsAtPoints[g_point_number].c[0] = coeffs[0];
				currPolynomial.coeffsAtPoints[g_point_number].c[1] = coeffs[1];
				currPolynomial.coeffsAtPoints[g_point_number].c[2] = coeffs[2];

				M(1, polynomial_number) = coeffs[0] * ksi_square_average[ind_0] + coeffs[1] * ksi_square_average[ind_1]
										  + coeffs[2] * ksi_square_average[ind_2];

				M(2, polynomial_number) = coeffs[0] * eta_square_average[ind_0] + coeffs[1] * eta_square_average[ind_1]
										  + coeffs[2] * eta_square_average[ind_2];
				M(3, polynomial_number) = coeffs[0] * ksi_eta_average[ind_0] + coeffs[1] * ksi_eta_average[ind_1]
										  + coeffs[2] * ksi_eta_average[ind_2];


			}


			M *= 1000;
			d *= 1000;
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
							0.5 * (gammas[i] + trRecData.so_polynomial.coeffsAtPoints[g_point_number].theta * gammas[i]);
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
			auto const area = stencil[i]->getArea2D();
			auto const x = stencil[i]->getBarycenter().x();
			auto const y = stencil[i]->getBarycenter().y();
/*			ksi_average[i] = 1.0 / area * T::GaussianIntegration([h, x_0](double x, double y)
															{
																return (x - x_0) / h;
															}, stencil[i]);

			eta_average[i] = 1.0 / area * T::GaussianIntegration([h, y_0](double x, double y)
															{
																return (y - y_0) / h;
															}, stencil[i]); */
			ksi_average[i] = (x - x_0) / h;
			eta_average[i] = (y - y_0) / h;

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

				auto const ksi = (trRecData.gaussian_points[g_point_number].x - x_0) / h;
				auto const eta = (trRecData.gaussian_points[g_point_number].y - y_0) / h;


				//first row, u = 1;
				A(i, 0) = ksi;
				A(i, 1) = ksi;
				A(i, 2) = ksi;

				A(i, 3) = eta;
				A(i, 4) = eta;
				A(i, 5) = eta;

				A(i, 6) = 1.0;
				A(i, 7) = 1.0;
				A(i, 8) = 1.0;

				c(i) = c_0 + c_1 + c_2;

				//second row, u = ksi
				A(i + 1, 0) = ksi_average[ind_0] * ksi;
				A(i + 1, 1) = ksi_average[ind_1] * ksi;
				A(i + 1, 2) = ksi_average[ind_2] * ksi;

				A(i + 1, 3) = ksi_average[ind_0] * eta;
				A(i + 1, 4) = ksi_average[ind_1] * eta;
				A(i + 1, 5) = ksi_average[ind_2] * eta;

				A(i + 1, 6) = ksi_average[ind_0];
				A(i + 1, 7) = ksi_average[ind_1];
				A(i + 1, 8) = ksi_average[ind_2];
				c(i + 1) = c_0 * ksi_average[ind_0] + c_1 * ksi_average[ind_1] + c_2 * ksi_average[ind_2];

				//third row, u = eta;
				A(i + 2, 0) = eta_average[ind_0] * ksi;
				A(i + 2, 1) = eta_average[ind_1] * ksi;
				A(i + 2, 2) = eta_average[ind_2] * ksi;

				A(i + 2, 3) = eta_average[ind_0] * eta;
				A(i + 2, 4) = eta_average[ind_1] * eta;
				A(i + 2, 5) = eta_average[ind_2] * eta;

				A(i + 2, 6) = eta_average[ind_0];
				A(i + 2, 7) = eta_average[ind_1];
				A(i + 2, 8) = eta_average[ind_2];
				c(i + 2) = c_0 * eta_average[ind_0] + c_1 * eta_average[ind_1] + c_2 * eta_average[ind_2];



				i+=3;
			}


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
	inline Vec4 WENOSolver<T>::Reconstruct(Vec4 const& qVec, Triangle const *pTriangle,
										   Point2 const &gaussianPoint, int edgeNumber) const
	{
		std::array<Triangle const*, 10> stencil;
		GetStencil(pTriangle, stencil);

		for(int i = 0; i < 10; ++i)
		{
			if(stencil[i] == nullptr)
				return qVec;
		}

		std::array<Vec4, 10> q;
		q[0] = qVec;
		double max_norm = std::sqrt(sqr(q[0][0]) + sqr(q[0][1]) + sqr(q[0][2]) + sqr(q[0][3]));
		for(int i = 1; i < 10; ++i)
		{
			T::FormQVector(q[i], stencil[i]);
#ifdef MY_STABILITY_FIX
			if(max_norm < std::sqrt(sqr(q[i][0]) + sqr(q[i][1]) + sqr(q[i][2]) + sqr(q[i][3])))
				max_norm = std::sqrt(sqr(q[i][0]) + sqr(q[i][1]) + sqr(q[i][2]) + sqr(q[i][3]));
#endif

		}
#ifdef MY_STABILITY_FIX
		for(int i = 0; i < 10; ++i)
		{
			q[i] *= 1.0 / (MY_STABILITY_FIX * max_norm);
		}
#endif


#ifdef CHARACTERISTIC_WISE
		//Forming Riemann invariants
		arma::mat44 L_A, L_B, R_A, R_B;
		//Vec4 q_avrg = 0.5 * (q[0] + q[1 + edgeNumber]);

		auto const density_minus = q[0][0];
		auto const velocityX_minus = q[0][1] / density_minus;
		auto const velocityY_minus = q[0][2] / density_minus;
		auto const E_minus = q[0][3] / density_minus;
		auto const velocity_sqr_abs_minus = sqr(velocityX_minus) + sqr(velocityY_minus);
		auto const eps_minus = E_minus - 0.5 * velocity_sqr_abs_minus;
		auto const pressure_minus = (T::m_gamma - 1.0) * eps_minus * density_minus;
		auto const H_minus = eps_minus + pressure_minus / density_minus + 0.5 * velocity_sqr_abs_minus;
		auto const c_minus = std::sqrt(T::m_gamma * pressure_minus / density_minus);

		auto const density_plus = q[1 + edgeNumber][0];
		auto const velocityX_plus = q[1 + edgeNumber][1] / density_plus;
		auto const velocityY_plus = q[1 + edgeNumber][2] / density_plus;
		auto const E_plus = q[1 + edgeNumber][3] / density_plus;
		auto const velocity_sqr_abs_plus = sqr(velocityX_plus) + sqr(velocityY_plus);
		auto const eps_plus = E_plus - 0.5 * velocity_sqr_abs_plus;
		auto const pressure_plus = (T::m_gamma - 1.0) * eps_plus * density_plus;
		auto const H_plus = eps_plus + pressure_plus / density_plus + 0.5 * velocity_sqr_abs_plus;
		auto const c_plus = std::sqrt(T::m_gamma * pressure_plus / density_plus);

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

		auto const velocity_sqr_abs_star = sqr(velocityX_star) + sqr(velocityY_star);



//		FormL_A(density_star, velocityX_star, velocityY_star, H_star, L_A);
//		FormL_B(density_star, velocityX_star, velocityY_star, H_star, L_B);
		FormR_A(density_star, velocityX_star, velocityY_star, H_star, R_A);
		FormR_B(density_star, velocityX_star, velocityY_star, H_star, R_B);

		bool char_wise = true;

		auto const normal = T::CalculateNormal(pTriangle, edgeNumber);
		arma::mat44 R = normal[0] * R_A + normal[1] * R_B;
		arma::mat44 L;
		if(std::fabs(det(R)) > 0.0001)
			 L = R.i();
		else
			char_wise = false; // then we reconstruct in component-wise way

		std::array<Vec4, 10> w;
		Vec4 w_reconstructed(0.0, 0.0, 0.0, 0.0);

		double max_norm_w = 0.0;
		if(char_wise)
		{
			for (int i = 0; i < 10; ++i)
			{
				w[i] = L * q[i];
#ifdef MY_STABILITY_FIX
				if(max_norm_w < std::sqrt(sqr(w[i][0]) + sqr(w[i][1]) + sqr(w[i][2]) + sqr(w[i][3])))
					max_norm_w = std::sqrt(sqr(w[i][0]) + sqr(w[i][1]) + sqr(w[i][2]) + sqr(w[i][3]));
#endif
			}
#ifdef MY_STABILITY_FIX
			for(int i = 0; i < 10; ++i)
				w[i] *= 1.0 / (MY_STABILITY_FIX * max_norm_w);
#endif
		}

#endif
		//Reconstruction!
		Vec4 q_reconstructed(0.0, 0.0, 0.0, 0.0);

		auto const triangleReconstructionData = m_vReconstructionData[pTriangle->Index()];

		//Searching for coressponding gaussian point
		int curr_g_point_n = 0;
		for(curr_g_point_n = 0; curr_g_point_n < gaussian_points_number; ++curr_g_point_n)
		{
			if(triangleReconstructionData.gaussian_points[curr_g_point_n] == gaussianPoint)
				break;
		}
		if(curr_g_point_n == gaussian_points_number)
			throw 1;


		std::array<Vec4, 9> omega, omega_waved;
		Vec4 o_wave_sum(0.0, 0.0, 0.0, 0.0);


		std::array<Vec4, 9> omega_plus, omega_waved_plus;
		Vec4 o_wave_sum_plus(0.0, 0.0, 0.0, 0.0);

		std::array<Vec4, 9> omega_minus, omega_waved_minus;
		Vec4 o_wave_sum_minus(0.0, 0.0, 0.0, 0.0);


		bool weights_to_be_treated =
				triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].weights_to_be_treated;




		Vec4 smoothIndicator(0.0, 0.0, 0.0, 0.0);


		static Vec4 const eps(1e-3, m_eps, m_eps, 1e-3);
		for(int i = 0; i < 9; ++i)
		{


			auto const ind_0 = triangleReconstructionData.fo_polynomial[i].stencil[0];
			auto const ind_1 = triangleReconstructionData.fo_polynomial[i].stencil[1];
			auto const ind_2= triangleReconstructionData.fo_polynomial[i].stencil[2];

#ifndef CHARACTERISTIC_WISE

			smoothIndicator = (sqr(triangleReconstructionData.smoothIndicatorData[i].alpha[0] * q[ind_0]
								   + triangleReconstructionData.smoothIndicatorData[i].alpha[1] * q[ind_1]
								   + triangleReconstructionData.smoothIndicatorData[i].alpha[2] * q[ind_2])
							   + sqr(triangleReconstructionData.smoothIndicatorData[i].beta[0] * q[ind_0]
									 + triangleReconstructionData.smoothIndicatorData[i].beta[1] * q[ind_1]
									 + triangleReconstructionData.smoothIndicatorData[i].beta[2] * q[ind_2]));

#else
			if(char_wise)
			{
				smoothIndicator = (sqr(triangleReconstructionData.smoothIndicatorData[i].alpha[0] * w[ind_0]
									   + triangleReconstructionData.smoothIndicatorData[i].alpha[1] * w[ind_1]
									   + triangleReconstructionData.smoothIndicatorData[i].alpha[2] * w[ind_2])
								   + sqr(triangleReconstructionData.smoothIndicatorData[i].beta[0] * w[ind_0]
										 + triangleReconstructionData.smoothIndicatorData[i].beta[1] * w[ind_1]
										 + triangleReconstructionData.smoothIndicatorData[i].beta[2] * w[ind_2]));

			}
			else
			{
				smoothIndicator =  (sqr(triangleReconstructionData.smoothIndicatorData[i].alpha[0] * q[ind_0]
									   + triangleReconstructionData.smoothIndicatorData[i].alpha[1] * q[ind_1]
									   + triangleReconstructionData.smoothIndicatorData[i].alpha[2] * q[ind_2])
								   + sqr(triangleReconstructionData.smoothIndicatorData[i].beta[0] * q[ind_0]
										 + triangleReconstructionData.smoothIndicatorData[i].beta[1] * q[ind_1]
										 + triangleReconstructionData.smoothIndicatorData[i].beta[2] * q[ind_2]));
			}
#endif


			if(!weights_to_be_treated)
			{
				omega_waved[i] = Vec4(triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i],
									  triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i],
									  triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i],
									  triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas[i])
								 / sqr(eps + smoothIndicator);

				o_wave_sum += omega_waved[i];

			}
			else // weights to be treated
			{
				omega_waved_plus[i] =
						Vec4(triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_plus[i],
							 triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_plus[i],
							 triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_plus[i],
							 triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_plus[i])
						/ sqr(eps + smoothIndicator);

				omega_waved_minus[i] =
						Vec4(triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_minus[i],
							 triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_minus[i],
							 triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_minus[i],
							 triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].gammas_minus[i])
						/ sqr(eps + smoothIndicator);

				o_wave_sum_plus += omega_waved_plus[i];
				o_wave_sum_minus += omega_waved_minus[i];

			}

		}



		for(int i = 0; i < 9; ++i)
		{
			auto const c_0 = triangleReconstructionData.fo_polynomial[i].coeffsAtPoints[curr_g_point_n].c[0];
			auto const c_1= triangleReconstructionData.fo_polynomial[i].coeffsAtPoints[curr_g_point_n].c[1];
			auto const c_2 = triangleReconstructionData.fo_polynomial[i].coeffsAtPoints[curr_g_point_n].c[2];

			auto const ind_0 = triangleReconstructionData.fo_polynomial[i].stencil[0];
			auto const ind_1 = triangleReconstructionData.fo_polynomial[i].stencil[1];
			auto const ind_2 = triangleReconstructionData.fo_polynomial[i].stencil[2];


			if(!weights_to_be_treated)
			{
				omega[i] = omega_waved[i] / o_wave_sum;

#ifndef CHARACTERISTIC_WISE
				q_reconstructed += omega[i] *
								   (c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2]);
#else

				if(char_wise)
				{
					w_reconstructed += omega[i] * (c_0 * w[ind_0] + c_1 * w[ind_1] + c_2 * w[ind_2]);

				}
				else
					q_reconstructed += omega[i] *
									   (c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2]);
#endif
			}
			else //weights to be treated
			{


				omega_plus[i] = omega_waved_plus[i] / o_wave_sum_plus;
				omega_minus[i] = omega_waved_minus[i] / o_wave_sum_minus;

#ifndef CHARACTERISTIC_WISE

				q_reconstructed +=
						(triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].sigma_plus
						 * omega_plus[i]
						 - triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].sigma_minus
						   * omega_minus[i]) * (c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2]);
#else
				if(char_wise)
				{

					w_reconstructed +=
							(triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].sigma_plus
							 * omega_plus[i]
							 - triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].sigma_minus
							   * omega_minus[i]) * (c_0 * w[ind_0] + c_1 * w[ind_1] + c_2 * w[ind_2]);
				}
				else
				{
					q_reconstructed +=
							(triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].sigma_plus
							 * omega_plus[i]
							 - triangleReconstructionData.so_polynomial.coeffsAtPoints[curr_g_point_n].sigma_minus
							   * omega_minus[i]) * (c_0 * q[ind_0] + c_1 * q[ind_1] + c_2 * q[ind_2]);
				}
#endif
			}


		}

#ifdef CHARACTERISTIC_WISE

#ifdef MY_STABILITY_FIX
		w_reconstructed *= (MY_STABILITY_FIX * max_norm_w);
#endif

		if(char_wise)
		{
			q_reconstructed = R * w_reconstructed;
		}
#endif

#ifdef MY_STABILITY_FIX
		q_reconstructed *= (MY_STABILITY_FIX * max_norm);
#endif
//		if(!((q_reconstructed(3) > 0) && (q_reconstructed(0) > 0)))
//			throw 1;


		return q_reconstructed;


	}

	template<class T>
	inline void WENOSolver<T>::FormL_A(double density, double velocityX, double velocityY,
									   double H, arma::mat44 &L_A) const
	{

		auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);
		auto const sound_speed = std::sqrt((T::m_gamma - 1) * (H - 0.5 * velocity_sqr_abs));


		L_A(0, 0) = velocity_sqr_abs / 2 + velocityX * sound_speed / (T::m_gamma - 1);
		L_A(0, 1) = -velocityX - sound_speed / (T::m_gamma - 1);
		L_A(0, 2) = -velocityY;
		L_A(0, 3) = 1.0;

		L_A(1, 0) = sound_speed * sound_speed / (T::m_gamma - 1) - velocity_sqr_abs / 2 - velocityY * sound_speed /(T::m_gamma - 1);
		L_A(1, 1) = velocityX;
		L_A(1, 2) = velocityY + sound_speed / (T::m_gamma - 1);
		L_A(1, 3) = -1.0;

		L_A(2, 0) = sound_speed * sound_speed / (T::m_gamma - 1) - velocity_sqr_abs / 2 + velocityY * sound_speed /(T::m_gamma - 1);
		L_A(2, 1) = velocityX;
		L_A(2, 2) = velocityY - sound_speed / (T::m_gamma - 1);
		L_A(2, 3) = -1.0;

		L_A(3, 0) = velocity_sqr_abs / 2 - velocityX * sound_speed / (T::m_gamma - 1);
		L_A(3, 1) = -velocityX + sound_speed / (T::m_gamma - 1);
		L_A(3, 2) = -velocityY;
		L_A(3, 3) = 1.0;

		L_A *= (T::m_gamma - 1) / (2 * sound_speed * sound_speed);

	}


	template<class T>
	inline void WENOSolver<T>::FormL_B(double density, double velocityX, double velocityY,
									   double H, arma::mat44 &L_B) const
	{

		auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);
		auto const sound_speed = std::sqrt((T::m_gamma - 1) * (H - 0.5 * velocity_sqr_abs));


		L_B(0, 0) = velocity_sqr_abs / 2 + velocityY * sound_speed / (T::m_gamma - 1);
		L_B(0, 1) = -velocityX;
		L_B(0, 2) = -velocityY - sound_speed / (T::m_gamma - 1);
		L_B(0, 3) = 1.0;

		L_B(1, 0) = sound_speed * sound_speed / (T::m_gamma - 1) - velocity_sqr_abs / 2 - velocityX * sound_speed / (T::m_gamma - 1);
		L_B(1, 1) = velocityX + sound_speed / (T::m_gamma - 1);
		L_B(1, 2) = velocityY;
		L_B(1, 3) = -1.0;

		L_B(2, 0) = sound_speed * sound_speed / (T::m_gamma - 1) - velocity_sqr_abs / 2 + velocityX * sound_speed / (T::m_gamma - 1);
		L_B(2, 1) = velocityX - sound_speed / (T::m_gamma - 1);
		L_B(2, 2) = velocityY;
		L_B(2, 3) = -1.0;

		L_B(3, 0) = velocity_sqr_abs / 2 - velocityY * sound_speed / (T::m_gamma - 1);
		L_B(3, 1) = -velocityX;
		L_B(3, 2) = -velocityY + sound_speed / (T::m_gamma - 1);
		L_B(3, 3) = 1.0;

		L_B *= (T::m_gamma - 1) / (2 * sound_speed * sound_speed);

	}


	template<class T>
	inline void WENOSolver<T>::FormR_A(double density, double velocityX, double velocityY,
									   double H, arma::mat44 &R_A) const
	{

		auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);
		auto const sound_speed = std::sqrt((T::m_gamma - 1) * (H - 0.5 * velocity_sqr_abs));



		R_A(0, 0) = 1;
		R_A(0, 1) = 1;
		R_A(0, 2) = 1;
		R_A(0, 3) = 1;

		R_A(1, 0) = velocityX - sound_speed;
		R_A(1, 1) = velocityX;
		R_A(1, 2) = velocityX;
		R_A(1, 3) = velocityX + sound_speed;

		R_A(2, 0) = velocityY;
		R_A(2, 1) = velocityY + sound_speed;
		R_A(2, 2) = velocityY - sound_speed;
		R_A(2, 3) = velocityY;

		R_A(3, 0) = H - velocityX * sound_speed;
		R_A(3, 1) = velocity_sqr_abs / 2 + velocityY * sound_speed;
		R_A(3, 2) = velocity_sqr_abs / 2 - velocityY * sound_speed;
		R_A(3, 3) = H + velocityX * sound_speed;


	}

	template<class T>
	inline void WENOSolver<T>::FormR_B(double density, double velocityX, double velocityY,
									   double H, arma::mat44 &R_B) const
	{

		auto const velocity_sqr_abs = sqr(velocityX) + sqr(velocityY);
		auto const sound_speed = std::sqrt((T::m_gamma - 1) * (H - 0.5 * velocity_sqr_abs));


		R_B(0, 0) = 1;
		R_B(0, 1) = 1;
		R_B(0, 2) = 1;
		R_B(0, 3) = 1;

		R_B(1, 0) = velocityX;
		R_B(1, 1) = velocityX + sound_speed;
		R_B(1, 2) = velocityX - sound_speed;
		R_B(1, 3) = velocityX;

		R_B(2, 0) = velocityY - sound_speed;
		R_B(2, 1) = velocityY;
		R_B(2, 2) = velocityY;
		R_B(2, 3) = velocityY + sound_speed;

		R_B(3, 0) = H - velocityY * sound_speed;
		R_B(3, 1) = velocity_sqr_abs / 2 + velocityX * sound_speed;
		R_B(3, 2) = velocity_sqr_abs / 2 - velocityX * sound_speed;
		R_B(3, 3) = H + velocityY * sound_speed;

	}



}







#endif //TRIANGULAR_SOLVERS_WENOLF_H
