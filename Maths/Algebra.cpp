#include "Algebra.h"
#include "assert.h"


using namespace euler;

Point2 operator+(Point2 const& left, Point2 const& right)
{
	return Point2(left.x + right.x, left.y + right.y);
}

Point2 operator-(Point2 const& left, Point2 const& right)
{
	return Point2(left.x - right.x, left.y - right.y);
}


Point2 operator*(double alpha, Point2 const& point)
{
	return Point2(alpha * point.x,  alpha * point.y);
}

bool operator==(Point2 const& left, Point2 const& right)
{
	return ((left.x == right.x) && (left.y == right.y));
}


double euler::sgn(double x)
{
	if( x > 0 )
		return 1;
	if( x < 0 )
		return -1;
	return 0;
}

double euler::minmod(double x, double y)
{
	return 0.5 * (sgn(x) + sgn(y)) * std::min(fabs(x), fabs(y));
}


euler::Vec4 euler::minmod(euler::Vec4 const &x, euler::Vec4 const &y)
{
	return euler::Vec4(euler::minmod(x[0], y[0]), euler::minmod(x[1], y[1]), euler::minmod(x[2], y[2]),
					   euler::minmod(x[3], y[3]));
}


double euler::dot(Vec4 const &first, Vec4 const &second)
{
	return first[0] * second[0] + first[1] * second[1] + first[2] * second[2] + first[3] * second[3];
}


euler::Vec4 euler::operator+(Vec4 const &left, Vec4 const &right)
{
	return Vec4(left[0] + right[0], left[1] + right[1], left[2] + right[2], left[3] + right[3]);
}

euler::Vec4 euler::operator-(Vec4 const &left, Vec4 const &right)
{
	return Vec4(left[0] - right[0], left[1] - right[1], left[2] - right[2], left[3] - right[3]);
}

euler::Vec4 euler::operator*(Vec4 const &vec, double alpha)
{
	return Vec4(alpha * vec[0], alpha * vec[1], alpha * vec[2], alpha * vec[3]);
}


euler::Vec4 euler::operator*(double alpha, Vec4 const &vec)
{
	return Vec4(alpha * vec[0], alpha * vec[1], alpha * vec[2], alpha * vec[3]);
}

euler::Vec4& euler::operator+=(euler::Vec4 &left, euler::Vec4 const &right)
{
	left[0] += right[0];
	left[1] += right[1];
	left[2] += right[2];
	left[3] += right[3];
	return left;
}
euler::Vec4& euler::operator-=(euler::Vec4 &left, euler::Vec4 const &right)
{
	left[0] -= right[0];
	left[1] -= right[1];
	left[2] -= right[2];
	left[3] -= right[3];
	return left;
}

euler::Matrix4x4& euler::operator*=(Matrix4x4 &mat, double alpha)
{
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
			mat[i][j] *= alpha;
	return mat;
}

double euler::sqr(double x)
{
	return x*x;
}


euler::Vec4 euler::sqr(euler::Vec4 const& vec)
{
	return Vec4(sqr(vec[0]), sqr(vec[1]), sqr(vec[2]), sqr(vec[3]));
}


euler::Vec4 euler::operator*(euler::Matrix4x4 const &mat, euler::Vec4 const &vec4)
{
	euler::Vec4 result;
	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j < 4; ++j)
			result[i] += mat[i][j] * vec4[j];
	}
	return result;
}

euler::Vec4 euler::absv(euler::Vec4 const &vec)
{
	return Vec4(std::fabs(vec[0]), std::fabs(vec[1]), std::fabs(vec[2]), std::fabs(vec[3]));
}
