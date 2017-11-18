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
			bool isBoundary = false;
			int parentIndex = -1;

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

		bool IsBoundary() const
		{
			return m_boundaryProperties.isBoundary;
		}

		int ParentIndex() const
		{
			return m_boundaryProperties.parentIndex;
		}

		void SetParentIndex(int const ind)
		{
			m_boundaryProperties.parentIndex = ind;
		}

		void SetBoundary(bool isbndry)
		{
			m_boundaryProperties.isBoundary = isbndry;
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

		Triangle* ReflectTriangle(const int edgeNumber)
		{

			auto const normal = CalculateNormal(edgeNumber);

			arma::mat22 A;
			A << normal[0] << normal[1] << arma::endr
			  << -normal[1] << normal[0] << arma::endr;

			auto const p1 = getCorner((edgeNumber + 1) % 3);
			auto const p2 = getCorner((edgeNumber + 2) % 3);

			auto const p0 = getCorner(edgeNumber); //point to be reflected

			auto const c = -(normal[0] * p2->x() + normal[1] * p2->y());


			auto const d = -(-normal[1] * p0->x() + normal[0] * p0->y());

			arma::vec rightSide;
			rightSide << -c << -d;

			arma::vec middlePoint = arma::solve(A, rightSide);

			auto const reflectedPoint = new GEOM_FADE2D::Point2(2 * middlePoint(0) - p0->x(), 2 * middlePoint(1) - p0->y());

			auto const reflectedTriangle = new Triangle();
			reflectedTriangle->SetIndex(-1);

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
