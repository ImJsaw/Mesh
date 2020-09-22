#pragma once
#include "GUA_OM.h"

using namespace OMT;

inline double Dot(const Point& a, const Point& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Point Cross(const Point& a, const Point& b) {
	return Point(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

inline Point Cross(const Point& a, const Point& b, const Point& P) {
	return Cross(P - a, b - a);
}

inline double _length(const Point p) {
	return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

inline bool InTriangle(const Point& p1, const Point& p2, const Point& p3, const Point& target) {
	double s = _length(Cross(p1 - p2, p1 - p3));
	double s1 = _length(Cross(p2 - target, p3 - target));
	double s2 = _length(Cross(target - p3, target - p1));
	double s3 = _length(Cross(target - p1, target - p2));

	double a = s1 / s;
	double b = s2 / s;
	double c = s3 / s;

	return 0 <= b && 0 <= a && 0 <= c && c <= 1 && b <= 1 && a <= 1 && (abs(b + a + c - 1) <= 0.0001f);
}