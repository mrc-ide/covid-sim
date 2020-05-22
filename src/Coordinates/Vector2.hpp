#pragma once

template<class T>
struct Vector2 {
	T x;
	T y;

	Vector2(T x, T y);

	Vector2<T> abs() const;
	T length() const;
	T length_squared() const;

	T distance_to(Vector2<T> other) const;
	T distance_squared_to(Vector2<T> other) const;

	Vector2<T> operator+(Vector2<T> other) const;
	Vector2<T> operator-(Vector2<T> other) const;
	Vector2<T> operator*(Vector2<T> other) const;
	Vector2<T> operator/(Vector2<T> other) const;

	template<class U>
	explicit operator Vector2<U>() const;
};

