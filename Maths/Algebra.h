#ifndef TRIANGULAR_SOLVERS_ALGEBRA_H
#define TRIANGULAR_SOLVERS_ALGEBRA_H

#include <array>
#include <cmath>


namespace euler
{

	class Point2
	{
		double x, y;
	};



	class Vec4: public std::array<double, 4>
	{

	public:

		Vec4(double x = 0.0, double y = 0.0, double z = 0.0, double w = 0.0): std::array<double, 4>{{x, y, z, w}}
		{}


		friend const Vec4 operator+(Vec4 const& left, Vec4 const& right);

		friend const Vec4 operator-(Vec4 const& left, Vec4 const& right);

		friend const Vec4 operator*(double alpha, Vec4 const& vec);

		friend const Vec4 operator*(Vec4 const& vec, double alpha);

		friend euler::Vec4& operator+=(euler::Vec4& first, euler::Vec4 const& right);

	};


	class Matrix4x4
	{
		class Row
		{
			std::array<double, 4> elements;
		public:
			double& operator[] (const int j)
			{
				if(j >= 4)
					throw 1;
				return elements[j];
			}
			double const& operator[] (const int j) const
			{
				if(j >= 4)
					throw 1;
				return elements[j];
			}
		};
		std::array<Row, 5> rows;

	public:
		Row& operator[] (const int i)
		{
			if(i >= 4)
				throw 1;
			return rows[i];
		}
		Row const& operator[] (const int i) const
		{
			if(i >= 4)
				throw 1;
			return rows[i];
		}


		friend const Vec4 operator*(Matrix4x4 const& mat, euler::Vec4 const& vec3);

		friend Matrix4x4& operator*=(Matrix4x4& mat, double alpha);

	};

	Matrix4x4& operator*=(Matrix4x4& mat, double alpha);

	const euler::Vec4 operator*(Matrix4x4 const& mat, euler::Vec4 const& vec3);

	euler::Vec4& operator+=(euler::Vec4& first, euler::Vec4 const& right);

	const Vec4 operator+(Vec4 const& left, Vec4 const& right);

	const Vec4 operator-(Vec4 const& left, Vec4 const& right);

	const Vec4 operator*(double alpha, Vec4 const& vec);

	const Vec4 operator*(Vec4 const& vec, double alpha);

	double dot(Vec4 const &first, Vec4 const &second);

	double sgn(double x);

	double minmod(double x, double y);

	Vec4 minmod(euler::Vec4 const& x, euler::Vec4 const& y);

	double sqr(double x);

	Vec4 sqr(euler::Vec4 const& vec);



}


#endif //TRIANGULAR_SOLVERS_ALGEBRA_H
