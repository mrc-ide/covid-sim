#pragma once

template<class T>
struct Vector2 {
	T x;
	T y;

	Vector2();
	Vector2(T s);
	Vector2(T x, T y);

	T length() const;
	T length_squared() const;

	T distance_to(const Vector2<T> &other) const;
	T distance_squared_to(const Vector2<T> &other) const;

	Vector2<T> abs() const;
	Vector2<T> floor() const;
	Vector2<T> ceil() const;

	Vector2<T> operator+(const Vector2<T> &other) const;
	Vector2<T>& operator+=(const Vector2<T> &other);
	Vector2<T> operator-(const Vector2<T> &other) const;
	Vector2<T>& operator-=(const Vector2<T> &other);
	Vector2<T> operator*(const Vector2<T> &other) const;
	Vector2<T> operator/(const Vector2<T> &other) const;

	Vector2<T> operator+(const T &other) const;
	Vector2<T> operator-(const T &other) const;
	Vector2<T> operator*(const T &other) const;
	Vector2<T> operator/(const T &other) const;

	template<class U>
	explicit operator Vector2<U>() const;
};

