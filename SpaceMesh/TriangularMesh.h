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


	public:
		double density = 0, velocityX = 0, velocityY = 0, pressure = 0;

	public:

		Triangle():m_index(-1)
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

			auto reflectedPoint = new GEOM_FADE2D::Point2(2 * middlePoint(0) - p0->x(), 2 * middlePoint(1) - p0->y());

			auto const reflectedTriangle = new Triangle();
			reflectedTriangle->SetIndex(Index());

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
			auto reflectedTriangle = ReflectTriangle(edgeNumber);
			auto const barycenter = new GEOM_FADE2D::Point2(reflectedTriangle->getBarycenter());
			auto const p0 = reflectedTriangle->getCorner(0);
			auto const p1 = reflectedTriangle->getCorner(1);
			auto const p2 = reflectedTriangle->getCorner(2);

			reflectedTriangle->setVertexPointer(1, barycenter);

			auto triangle0 = new Triangle();
			auto triangle1 = new Triangle();
			triangle0->SetIndex(m_index);
			triangle1->SetIndex(m_index);

			triangle0->setVertexPointer(0, p2);
			triangle0->setVertexPointer(1, barycenter);
			triangle0->setVertexPointer(2, p1);

			triangle1->setVertexPointer(0, p0);
			triangle1->setVertexPointer(1, p1);
			triangle1->setVertexPointer(2, barycenter);

			triangle0->SetOppTriangle(0, triangle1);
			triangle0->SetOppTriangle(1, nullptr);
			triangle0->SetOppTriangle(2, reflectedTriangle);

			triangle1->SetOppTriangle(0, triangle0);
			triangle1->SetOppTriangle(1, reflectedTriangle);
			triangle1->SetOppTriangle(2, nullptr);

			reflectedTriangle->SetOppTriangle(0, triangle0);
			reflectedTriangle->SetOppTriangle(1, this);
			reflectedTriangle->SetOppTriangle(2, triangle1);


			std::array<Triangle*, 3> triangles = {triangle0, reflectedTriangle, triangle1};

			return triangles;
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
