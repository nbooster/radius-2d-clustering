#pragma once

#include <compare>
#include <iostream>

class Point
{
	double x_{ 0.0 };
	double y_{ 0.0 };

public:

	Point() = default;

	Point(const double x, const double y) :
		x_{ x }, y_{ y } {}

	double x() const noexcept
	{
		return x_;
	}

	double y() const noexcept
	{
		return y_;
	}

	auto operator<=>(const Point&) const = default;

	friend std::ostream& operator << (std::ostream& os, const Point& point)
	{
		os << "( " << point.x_ << ", " << point.y_ << " )";
		return os;
	}
};
