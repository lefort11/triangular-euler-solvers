#ifndef TRIANGULAR_SOLVERS_ALGEBRA_H
#define TRIANGULAR_SOLVERS_ALGEBRA_H

#include <array>
#include <cmath>

#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


namespace euler
{
 /*	class Vec4: public arma::vec::fixed<4>
	{
	public:
		Vec4(std::array<double, 4> const& vec)
		{
			(*this)(0) = vec[0];
			(*this)(1) = vec[1];
			(*this)(2) = vec[2];
			(*this)(3) = vec[3];
		}


	}; */


	class Vec4: public arma::vec::fixed<4>
	{

	public:

/*		Vec4(double x = 0.0, double y = 0.0, double z = 0.0, double w = 0.0): std::array<double, 4>({x,y,z,w})
		{
		} */

		Vec4(std::array<double, 4> const &vec)
		{
			(*this)(0) = vec[0];
			(*this)(1) = vec[1];
			(*this)(2) = vec[2];
			(*this)(3) = vec[3];
		}

		Vec4& operator=(const std::array<double, 4>& vec) {
			(*this)(0) = vec[0];
			(*this)(1) = vec[1];
			(*this)(2) = vec[2];
			(*this)(3) = vec[3];
			return (*this);
		}

		Vec4(double x = 0.0, double y = 0.0, double z = 0.0, double w = 0.0)
		{
			(*this)(0) = x;
			(*this)(1) = y;
			(*this)(2) = z;
			(*this)(3) = w;
		}

/*		double &operator[](const int i)
		{
			return (*this)(i);
		}

		double const& operator[](const int i) const
		{
			return (*this)(i);
		} */


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
