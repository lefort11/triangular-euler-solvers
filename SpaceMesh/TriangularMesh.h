#ifndef TRIANGULAR_SOLVERS_TRIANGULARMESH_H
#define TRIANGULAR_SOLVERS_TRIANGULARMESH_H

#include <vector>
#include <Fade_2D.h>
#include <array>
#include "../Maths/Algebra.h"

namespace euler
{

	class Triangle: public GEOM_FADE2D::Triangle2
	{
	private:
		std::array<Triangle*, 3> m_pOppTriangles;

		int m_index;

		struct
		{
			bool isVirtual = false;
			bool isBoundary = false;
			int parentIndex = -1;

			bool toBeReconstructed = true;

		} m_boundaryProperties;


	public:
		double density = 0, velocityX = 0, velocityY = 0, pressure = 0;

	public:

		Triangle():m_index(-1), m_boundaryProperties()
		{
			m_pOppTriangles[0] = nullptr;
			m_pOppTriangles[1] = nullptr;
			m_pOppTriangles[2] = nullptr;
		}

		Triangle(GEOM_FADE2D::Triangle2 const& tr2): GEOM_FADE2D::Triangle2(tr2)
		{
			m_pOppTriangles[0] = nullptr;
			m_pOppTriangles[1] = nullptr;
			m_pOppTriangles[2] = nullptr;
		}

		void SetState(std::array<double, 4> const& state)
		{
			density = state[0];
			velocityX = state[1];
			velocityY = state[2];
			pressure = state[3];
		}

		void SetIndex(int index)
		{
			m_index = index;
		}

		int Index() const
		{
			return m_index;
		}

		bool IsVirtual() const
		{
			return m_boundaryProperties.isVirtual;
		}

		bool IsBoundary() const
		{
			return m_boundaryProperties.isBoundary;
		}

		bool ToBeReconstructed() const
		{
			return m_boundaryProperties.toBeReconstructed;
		}


		int ParentIndex() const
		{
			return m_boundaryProperties.parentIndex;
		}

		void SetParentIndex(int const ind)
		{
			m_boundaryProperties.parentIndex = ind;
		}



		void SetVirtual(bool isvirtual)
		{
			m_boundaryProperties.isVirtual = isvirtual;
		}

		void SetBoundary(bool isbndry)
		{
			m_boundaryProperties.isBoundary = isbndry;
		}

		void SetReconstructionNecessity(bool tobereconstructed)
		{
			m_boundaryProperties.toBeReconstructed = tobereconstructed;
		}

		void SetOppTriangle(const int ith, Triangle* pTriangle)
		{
			m_pOppTriangles[ith] = pTriangle;
		}


		Triangle* GetOppTriangle(const int ith) const
		{
			return m_pOppTriangles[ith];
		}

		std::array<double, 2> CalculateNormal(int edgeNumber) const
		{
			auto const vertex1 = getCorner((edgeNumber + 1) % 3);
			auto const vertex2 = getCorner((edgeNumber + 2) % 3);

			auto const p_x = vertex2->x() - vertex1->x();
			auto const p_y = vertex2->y() - vertex1->y();

			auto const p_abs = std::sqrt(p_x * p_x + p_y * p_y);

			return {p_y/p_abs, -p_x/p_abs};
		};

        GEOM_FADE2D::Point2 ReflectPoint(GEOM_FADE2D::Point2* const pointToBeReflected, GEOM_FADE2D::Point2* const p1,
                                         GEOM_FADE2D::Point2* const p2) const
        {
            auto const n_x = p2->y() - p1->y();
            auto const n_y = -(p2->x() - p1->x());

            arma::mat22 A;
            A << n_x << n_y << arma::endr
              << -n_y << n_x << arma::endr;

            auto const c = -(n_x * p2->x() + n_y * p2->y());

            auto const d = n_y * pointToBeReflected->x() - n_x * pointToBeReflected->y();

            arma::vec2 rightSide;
            rightSide[0] = -c;
            rightSide[1] = -d;

            arma::vec2 middlePoint = arma::solve(A, rightSide);

            return GEOM_FADE2D::Point2(2 * middlePoint(0) - pointToBeReflected->x(),
                                       2 * middlePoint(1) - pointToBeReflected->y());


        }

		Triangle* ReflectTriangle(const int edgeNumber)
		{




			auto const p1 = getCorner((edgeNumber + 1) % 3);
			auto const p2 = getCorner((edgeNumber + 2) % 3);

			auto const p0 = getCorner(edgeNumber); //point to be reflected



			auto const reflectedTriangle = new Triangle();
			reflectedTriangle->SetIndex(-1);

            auto const reflectedPoint = new GEOM_FADE2D::Point2(ReflectPoint(p0, p1, p2));
			//reflectedTriangle->setProperties(p1, reflectedPoint, p2);
			reflectedTriangle->setVertexPointer(0, p1);
			reflectedTriangle->setVertexPointer(1, reflectedPoint);
			reflectedTriangle->setVertexPointer(2, p2);

			reflectedTriangle->SetOppTriangle(0, nullptr);
			reflectedTriangle->SetOppTriangle(1, nullptr);
			reflectedTriangle->SetOppTriangle(2, nullptr);


			reflectedTriangle->SetOppTriangle(reflectedTriangle->getIntraTriangleIndex(reflectedPoint), this);

			SetOppTriangle(edgeNumber, reflectedTriangle);

			return reflectedTriangle;

		}



		void CreateWENOVirtualTriangles(int const edgeNumber, std::array<Triangle*, 7> & vtriangles)
		{
			vtriangles[0] = ReflectTriangle(edgeNumber);
			vtriangles[0]->SetVirtual(true);
			vtriangles[0]->SetParentIndex(m_index);

			auto const p0 = getCorner(edgeNumber);
			auto const p1 = getCorner((edgeNumber + 1) % 3);
			auto const p2 = getCorner((edgeNumber + 2) % 3);
			auto const refPoint = vtriangles[0]->getCorner(1);


			auto const ind1 = vtriangles[0]->getIntraTriangleIndex(getCorner((edgeNumber + 1) % 3));
			auto const ind2 = vtriangles[0]->getIntraTriangleIndex(getCorner((edgeNumber + 2) % 3));


            if((GetOppTriangle((edgeNumber + 1) % 3) != nullptr) &&
                    !(GetOppTriangle((edgeNumber + 1) % 3)->IsVirtual()))
			{

				auto const tr1 = GetOppTriangle((edgeNumber + 1) % 3);
				auto const ind1_0 = tr1->getIntraTriangleIndex(p0);

                vtriangles[1] = new Triangle();
                vtriangles[1]->setVertexPointer(1, new GEOM_FADE2D::Point2(
                        ReflectPoint(tr1->getCorner((ind1_0 + 2) % 3), p1, p2)));
				vtriangles[1]->setVertexPointer(0, refPoint);
				vtriangles[1]->setVertexPointer(2, p2);

				vtriangles[1]->SetVirtual(true);
				assert(!tr1->IsVirtual());
                vtriangles[1]->SetParentIndex(tr1->Index());

				vtriangles[0]->SetOppTriangle(0, vtriangles[1]);


				auto const rp1 = vtriangles[1]->getCorner(1);

				if((tr1->GetOppTriangle((ind1_0 + 1) % 3) != nullptr) &&
				   !(tr1->GetOppTriangle((ind1_0 + 1) % 3)->IsVirtual()))
				{
					auto const tr3 = tr1->GetOppTriangle((ind1_0 + 1) % 3);
					auto const ind3_0 = tr3->getIntraTriangleIndex(p0);

					vtriangles[2] = new Triangle();
					vtriangles[2]->setVertexPointer(1, new GEOM_FADE2D::Point2(
							ReflectPoint(tr3->getCorner((ind3_0 + 2) % 3), p1, p2)));
					vtriangles[2]->setVertexPointer(0, refPoint);
					vtriangles[2]->setVertexPointer(2, rp1);
					vtriangles[2]->SetVirtual(true);
					vtriangles[2]->SetParentIndex(tr3->Index());

					vtriangles[1]->SetOppTriangle(2, vtriangles[2]);
					vtriangles[2]->SetOppTriangle(1, vtriangles[1]);

				}
				else
				{

					vtriangles[2] = vtriangles[1]->ReflectTriangle(2);
                    vtriangles[2]->SetVirtual(true);
                    vtriangles[2]->SetParentIndex(vtriangles[1]->ParentIndex());

				}

				if((tr1->GetOppTriangle(ind1_0)!= nullptr) &&
				   !(tr1->GetOppTriangle(ind1_0)->IsVirtual()))
				{

					auto const tr5 = tr1->GetOppTriangle(ind1_0);
					auto const ind5_2 = tr5->getIntraTriangleIndex(p2);

					vtriangles[3] = new Triangle();
					vtriangles[3]->setVertexPointer(1, new GEOM_FADE2D::Point2(
							ReflectPoint(tr5->getCorner((ind5_2 + 1) % 3), p1, p2)));
					vtriangles[3]->setVertexPointer(0, rp1);
					vtriangles[3]->setVertexPointer(2, p2);
					vtriangles[3]->SetVirtual(true);
					vtriangles[3]->SetParentIndex(tr5->Index());

					vtriangles[1]->SetOppTriangle(0, vtriangles[3]);
					vtriangles[3]->SetOppTriangle(1, vtriangles[1]);

				}
				else
				{

					vtriangles[3] = vtriangles[1]->ReflectTriangle(0);
                    if((std::fabs(vtriangles[3]->getCorner(1)->x() - tr1->getCorner(ind1_0)->x()) < 1e-4) &&
                            (std::fabs((vtriangles[3]->getCorner(1)->y() - tr1->getCorner(ind1_0)->y())) < 1e-4))
                    {
                        delete vtriangles[3];
                        vtriangles[3] = nullptr;
                        vtriangles[1]->SetOppTriangle(0, tr1);
                    }
                    else
                    {
                        vtriangles[3]->SetVirtual(true);
                        vtriangles[3]->SetParentIndex(vtriangles[1]->ParentIndex());
                    }

				}


			}
			else
			{
				auto const vtr1 = vtriangles[0]->SummonThreeTriangles(ind1);
				for(int i = 0; i < 3; ++i)
				{
					vtr1[i]->SetVirtual(true);
					vtr1[i]->SetParentIndex(m_index);
					vtriangles[i + 1] = vtr1[i];
				}

			}



			if((GetOppTriangle((edgeNumber + 2) % 3) != nullptr) &&
			   !(GetOppTriangle((edgeNumber + 2) % 3)->IsVirtual()))
			{

				auto const tr2 = GetOppTriangle((edgeNumber + 2) % 3);
				auto const ind2_0 = tr2->getIntraTriangleIndex(p0);
				vtriangles[4] = new Triangle();
				vtriangles[4]->setVertexPointer(1, new GEOM_FADE2D::Point2(
						ReflectPoint(tr2->getCorner((ind2_0 + 1) % 3), p1, p2)));
				vtriangles[4]->setVertexPointer(0, p1);
				vtriangles[4]->setVertexPointer(2, refPoint);
				assert(!tr2->IsVirtual());


				vtriangles[0]->SetOppTriangle(2, vtriangles[4]);

				vtriangles[4]->SetVirtual(true);

				vtriangles[4]->SetParentIndex(tr2->Index());

				auto const rp2 = vtriangles[4]->getCorner(1);

				if((tr2->GetOppTriangle((ind2_0 + 2) % 3) != nullptr) &&
				   !(tr2->GetOppTriangle((ind2_0 + 2) % 3)->IsVirtual()))
				{
					auto const tr4 = tr2->GetOppTriangle((ind2_0 + 2) % 3);
					auto const ind4_0 = tr4->getIntraTriangleIndex(p0);

					vtriangles[5] = new Triangle();
					vtriangles[5]->setVertexPointer(1, new GEOM_FADE2D::Point2(
							ReflectPoint(tr4->getCorner((ind4_0 + 1) % 3), p1, p2)));
					vtriangles[5]->setVertexPointer(0, rp2);
					vtriangles[5]->setVertexPointer(2, refPoint);
					vtriangles[5]->SetVirtual(true);
					vtriangles[5]->SetParentIndex(tr4->Index());

					vtriangles[4]->SetOppTriangle(0, vtriangles[5]);
					vtriangles[5]->SetOppTriangle(1, vtriangles[4]);

				}
				else
				{

					vtriangles[5] = vtriangles[4]->ReflectTriangle(0);
                    vtriangles[5]->SetVirtual(true);
                    vtriangles[5]->SetParentIndex(vtriangles[4]->ParentIndex());

				}

				if((tr2->GetOppTriangle(ind2_0)!= nullptr) &&
				   !(tr2->GetOppTriangle(ind2_0)->IsVirtual()))
				{

					auto const tr6 = tr2->GetOppTriangle(ind2_0);
					auto const ind6_1 = tr6->getIntraTriangleIndex(p1);

					vtriangles[6] = new Triangle();
					vtriangles[6]->setVertexPointer(1, new GEOM_FADE2D::Point2(
							ReflectPoint(tr6->getCorner((ind6_1 + 2) % 3), p1, p2)));
					vtriangles[6]->setVertexPointer(0, p1);
					vtriangles[6]->setVertexPointer(2, rp2);
					vtriangles[6]->SetVirtual(true);
					vtriangles[6]->SetParentIndex(tr6->Index());

					vtriangles[4]->SetOppTriangle(2, vtriangles[6]);
					vtriangles[6]->SetOppTriangle(1, vtriangles[4]);

				}
				else
				{

					vtriangles[6] = vtriangles[4]->ReflectTriangle(2);
                    if((std::fabs(vtriangles[6]->getCorner(1)->x() - tr2->getCorner(ind2_0)->x()) < 1e-4) &&
                            (std::fabs(vtriangles[6]->getCorner(1)->y() -  tr2->getCorner(ind2_0)->y()) < 1e-4))
                    {
                        delete vtriangles[6];
                        vtriangles[6] = nullptr;
                        vtriangles[4]->SetOppTriangle(2, tr2);
                    } else
                    {
                        vtriangles[6]->SetVirtual(true);
                        vtriangles[6]->SetParentIndex(vtriangles[4]->ParentIndex());
                    }

				}


			}
			else
			{
				auto const vtr2 = vtriangles[0]->SummonThreeTriangles(ind2);
				for(int i = 0; i < 3; ++i)
				{
					vtr2[i]->SetVirtual(true);
					vtr2[i]->SetParentIndex(m_index);
					vtriangles[i + 4] = vtr2[i];
				}

			}


			for(int i = 0; i < 7; ++i)
            {
                if(vtriangles[i] != nullptr)
                    vtriangles[i]->SetReconstructionNecessity(false);

            }
            vtriangles[0]->SetReconstructionNecessity(true);
		}




		std::array<Triangle*, 3> SummonThreeTriangles(const int edgeNumber)
		{
			auto triangle0 = ReflectTriangle(edgeNumber);
			triangle0->SetIndex(-1);
			auto const barycenter = new GEOM_FADE2D::Point2(triangle0->getBarycenter());
			auto const p0 = triangle0->getCorner(0);
			auto const p1 = triangle0->getCorner(1);
			auto const p2 = triangle0->getCorner(2);

			auto triangle1 = new Triangle();
			auto triangle2 = new Triangle();
			triangle1->SetIndex(-1);
			triangle2->SetIndex(-1);

			triangle1->setVertexPointer(0, p0);
			triangle1->setVertexPointer(1, p1);
			triangle1->setVertexPointer(2, barycenter);

			triangle2->setVertexPointer(0, p1);
			triangle2->setVertexPointer(1, p2);
			triangle2->setVertexPointer(2, barycenter);

			triangle0->setVertexPointer(1, barycenter);


			triangle0->SetOppTriangle(0, triangle2);
			triangle0->SetOppTriangle(1, this);
			triangle0->SetOppTriangle(2, triangle1);


			triangle1->SetOppTriangle(0, triangle2);
			triangle1->SetOppTriangle(1, triangle0);
			triangle1->SetOppTriangle(2, nullptr);

			triangle2->SetOppTriangle(0, triangle0);
			triangle2->SetOppTriangle(1, triangle1);
			triangle2->SetOppTriangle(2, nullptr);

			return std::array<Triangle*, 3>{triangle0, triangle1, triangle2};
		}

	};




	class TriangularMesh:public std::vector<Triangle*>
	{

	public:

		TriangularMesh(): std::vector<Triangle*>()
		{}

		explicit TriangularMesh(std::vector<Triangle*> const& mesh):
				std::vector<Triangle*>(mesh)
		{}

		TriangularMesh(std::vector<Triangle*> const& mesh, std::vector<std::array<double, 4>> const& initialState):
				std::vector<Triangle*>(mesh)
		{
			assert(this->size() == initialState.size());
			for(int i = 0; i < this->size(); ++i)
				(*this)[i]->SetState(initialState[i]);
		}

	};
}


#endif //TRIANGULAR_SOLVERS_TRIANGULARMESH_H
