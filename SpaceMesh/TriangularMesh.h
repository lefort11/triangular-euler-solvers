#ifndef TRIANGULAR_SOLVERS_TRIANGULARMESH_H
#define TRIANGULAR_SOLVERS_TRIANGULARMESH_H

#include <vector>
#include <Fade_2D.h>
#include <array>

namespace euler
{

	class Triangle //~triangle in fade_2d triangulation;
	{


		GEOM_FADE2D::Triangle2* triangle;
		std::array<double, 4> m_currentState; // {ro, v, w, P} at current time in the centroid of the triangle

	public:
		Triangle(std::array<double, 4> const& currentState): m_currentState(currentState),
															 triangle(nullptr)
		{}


		void SetState(std::array<double, 4> const& currentState)
		{
			m_currentState = currentState;
		}

		void SetTriangle2(GEOM_FADE2D::Triangle2* tr2)
		{
			triangle = tr2;
		}

		GEOM_FADE2D::Point2 GetBarycenter() const
		{
			return triangle->getBarycenter();
		}

		GEOM_FADE2D::Point2* GetCorner(int const i) const
		{
			return triangle->getCorner(i);
		}

/*		GEOM_FADE2D::Triangle2* GetOppositeTriangle(int const i) const
		{
			return triangle->getOppositeTriangle(i);
		} */

	};

	class TriangularMesh
	{
		std::vector<Triangle> m_Mesh;

	public:
		TriangularMesh(std::vector<Triangle> const& vec): m_Mesh(vec)
		{}

		//!@param i-th elem of initialStates corresponds to i-th triangle
		TriangularMesh(std::vector<GEOM_FADE2D::Triangle2*> const& vTr2,
					   std::vector<std::array<double, 4>> const& initialStates):
				m_Mesh()
		{
			assert(vTr2.size() == initialStates.size());
			for(int i = 0; i < vTr2.size(); ++i)
			{
				auto triangle = Triangle(initialStates[i]);
				triangle.SetTriangle2(vTr2[i]);

				m_Mesh.push_back(triangle);
			}
		}
	};
}


#endif //TRIANGULAR_SOLVERS_TRIANGULARMESH_H
