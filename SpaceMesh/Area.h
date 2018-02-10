#ifndef TRIANGULAR_SOLVERS_AREA_H
#define TRIANGULAR_SOLVERS_AREA_H

#include <vector>
#include <functional>

#include "TriangularMesh.h"

#define X_LEFT (-4.2)
#define X_RIGHT 7.8
#define Y_BOT (-6.0)
#define Y_TOP 6.0


namespace euler
{

	class ConstraintFunction: public std::function<GEOM_FADE2D::Point2(double)> // function of argument t from 0 to 1
	{
		unsigned m_discrPointNumber;

	public:

		explicit ConstraintFunction(std::function<GEOM_FADE2D::Point2(double)> const& func, unsigned discrPointNumber):
				std::function<GEOM_FADE2D::Point2(double)>(func), m_discrPointNumber(discrPointNumber)
		{}


		std::vector<GEOM_FADE2D::Point2> Discretize() const;


	};

    class MeshParams
    {
    public:
        double minAngleDegree = 20.0;
        double minEdgeLength = 0.001;
        double maxEdgeLength = DBL_MAX;
        GEOM_FADE2D::Vector2 gridVector = {1.0, 0.0};
        double gridLength = DBL_MAX;
        double growFactor = 5.0;
        double growFactorMinArea = 0.001;

    };



    //!@brief Zone class describes a zone bounded by some functions.
	class Zone
	{
		const ConstraintFunction m_function;
		const bool m_isInside;


//		const std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> m_boundaryConditions; //returns {ro, u, v, P}

	public:
		Zone(ConstraintFunction const& func, bool isInside): m_function(func),
															 m_isInside(isInside)
		{}


		std::vector<GEOM_FADE2D::Point2> Discretize() const
		{
			return m_function.Discretize();
		}

		bool IsInside() const
		{
			return m_isInside;
		}


	};

	//!@brief Area class describes an area inside the bounding rectangle
	class Area
	{
		std::vector<Zone> m_zones;

		GEOM_FADE2D::Fade_2D m_globalArea;

	public:

		Area(std::vector<Zone> const& constraints): m_zones(constraints)
		{}


		/**@brief trianglulate area consisted of zones
		 * @param triangleProperties = {minAngleDegree, minEdgeLength, maxEdgeLength}
		 * @param discrPointNumber - point number of constraint function discretezation
		 * @param initStateFunc - initial function of the problem
		**/
		TriangularMesh Triangulate(MeshParams const& triangleProperties,
								   std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> const& initStateFunc);

	private:

		/**@brief Finds a position of the Triangle in the triangle vector
		 *
		 * @param vpTriangle2 Triangle vector
		 * @param pTr triangle pointer
		 * @return an index
		 */
		int FindTriangleIndex(std::vector<GEOM_FADE2D::Triangle2*> const& vpTriangle2,
							  GEOM_FADE2D::Triangle2* const pTr) const
		{
			int i = 0;
			bool found = false;
			while((!found) && (i < vpTriangle2.size()))
			{
				if(vpTriangle2[i] == pTr)
					found = true;
				else
					++i;
			}
			if(found)
				return i;
			return -1;
		}

		/**@brief Makes a Triangle graph from a GEOM_FADE2D::Triangle2 graph. Needed to generate a mesh
		 *
		 * @param destination
		 * @param source
		**/
		void CopyGraph(std::vector<Triangle*> & destination,
					   std::vector<GEOM_FADE2D::Triangle2*> const& source) const
		{

			destination.resize(source.size());
			for(int trngl_cntr = 0; trngl_cntr < source.size(); ++trngl_cntr)
			{
				destination[trngl_cntr] = new Triangle((*source[trngl_cntr]));
				destination[trngl_cntr]->SetIndex(trngl_cntr);
			}
			for(int trngl_cntr = 0; trngl_cntr < source.size(); ++trngl_cntr)
			{
				for(int ith = 0; ith < 3; ++ith)
				{
					auto const oppTriangle = source[trngl_cntr]->getOppositeTriangle(ith);
					auto const oppTrianglePosition = FindTriangleIndex(source, oppTriangle);
					if(oppTrianglePosition != -1)
						destination[trngl_cntr]->SetOppTriangle(ith, destination[oppTrianglePosition]);
					else
						destination[trngl_cntr]->SetOppTriangle(ith, nullptr);
				}
			}

		}
	};


}


#endif //TRIANGULAR_SOLVERS_AREA_H
