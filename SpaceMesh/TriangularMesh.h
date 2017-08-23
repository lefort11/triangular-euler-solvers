#ifndef TRIANGULAR_SOLVERS_TRIANGULARMESH_H
#define TRIANGULAR_SOLVERS_TRIANGULARMESH_H

#include <vector>
#include <Fade_2D.h>
#include <array>

namespace euler
{
/*
	class Triangle: public GEOM_FADE2D::Triangle2
	{


	public:


		void SetState(std::array<double, 4> const& currentState)
		{
			m_CurrentState = currentState;
		}





		double& Density()
		{
			return m_CurrentState[0];
		}

		double& VelocityX() // global x coord of velocity
		{
			return m_CurrentState[1];
		}

		double& VelocityY() // global y coord of velocyty
		{
			return m_CurrentState[2];
		}

		double& Pressure()
		{
			return m_CurrentState[3];
		}



	}; */

	//fade2d library doesn't allow to modernize Triangle2 class saving the graph. Library was fixed in "Triangle2.h"
	typedef GEOM_FADE2D::Triangle2 Triangle;



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

/*		void SortByX()
		{
			std::sort(m_mesh.begin(), m_mesh.end(), [](Triangle const& first, Triangle const& scnd)
			{
				return (first.GetBarycenter().x() <  scnd.GetBarycenter().x());
			});
		}

		void SortByY()
		{
			std::sort(m_mesh.begin(), m_mesh.end(), [](Triangle const& first, Triangle const& scnd)
			{
				return (first.GetBarycenter().y() <  scnd.GetBarycenter().y());
			});
		} */
	};
}


#endif //TRIANGULAR_SOLVERS_TRIANGULARMESH_H
