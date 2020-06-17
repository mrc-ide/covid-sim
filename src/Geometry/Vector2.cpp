#include "Vector2.h"

using namespace Geometry;

Vector2<double> Geometry::operator*(const Vector2<double> &left, const Vector2<float> &right){
	return left * Vector2<double>(right);
}
Vector2<double> Geometry::operator*(const Vector2<float> &left, const Vector2<double> &right){
	return Vector2<double>(left) * right;
}

Vector2<double> Geometry::operator-(const Vector2<double> &left, const Vector2<int> &right){
	return left - Vector2<double>(right);
}
Vector2<double> Geometry::operator-(const Vector2<int> &left, const Vector2<double> &right){
	return Vector2<double>(left) - right;
}
