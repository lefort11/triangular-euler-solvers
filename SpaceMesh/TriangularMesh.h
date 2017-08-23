#ifndef TRIANGULAR_SOLVERS_TRIANGULARMESH_H
#define TRIANGULAR_SOLVERS_TRIANGULARMESH_H

#include <vector>
#include <Fade_2D.h>
#include <array>

namespace euler
{

	class Triangle: public GEOM_FADE2D::Triangle2
	{
	private:
		Triangle* m_pOppTriangles[3];


	public:
		double density = 0, velocityX = 0, velocityY = 0, pressure = 0;

	public:

		Triangle()
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

		void SetOppTriangle(const int ith, Triangle* const pTriangle)
		{
			m_pOppTriangles[ith] = pTriangle;
		}

		Triangle* GetOppTriangle(const int ith)
		{
			return m_pOppTriangles[ith];
		}

	};




	class TriangularMesh:public std::vector<Triangle*>
	{


	public:

		TriangularMesh(): std::vector<Triangle*>()
		{}

		explicit TriangularMesh(std::vector<Triangle*> const& mesh): std::vector<Triangle*>(mesh)
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
