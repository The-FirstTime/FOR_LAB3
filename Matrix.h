#pragma once
#include <initializer_list>
#include <iostream>
namespace linalg {
	template<typename T>
	class Matrix {
	public:
		template<typename T2> friend class Matrix;
		Matrix() : m_ptr(nullptr), m_rows(0), m_columns(0) {}//дефолтный конструктор
		Matrix(int rows, int cols);//конструктор заполняет нулями
		Matrix(int rows);

		template<typename T2>
		Matrix(std::initializer_list<std::initializer_list<T2>> lst);
		template <typename T2>
		Matrix(std::initializer_list<T2> lst);
		Matrix(const Matrix& object);//конструктор копирования

		template<typename T2>
		Matrix(const Matrix<T2>& obj);
		Matrix(Matrix&& object)noexcept;

		int rows() const noexcept { return m_rows; }
		int columns() const noexcept { return m_columns; }
		bool empty() const { return m_ptr == nullptr; }
		void reshape(int rows, int cols);

		size_t capacity() const { return m_capacity; }
		size_t size() const noexcept { return m_rows*m_columns; }
		void reserve(size_t n);
		void shrink_to_fit();
		void clear();

		T norm() const;
		T trace() const;
		void gauss_forward();
		void gauss_backward();
		T determinant() const;
		int rank() const;
		double getMinor(int firstRow, int firstCol, int n) const;

		Matrix& operator=(const Matrix& object);
		Matrix& operator=(Matrix&& object)noexcept;

		template<typename T2>
		Matrix& operator=(const Matrix<T2>& obj);

		T& operator()(size_t row, size_t col);
		T operator()(size_t row, size_t col) const;

		template<typename T2>
		Matrix& operator +=(const Matrix<T2>& object);

		template<typename T2>
		Matrix& operator -=(const Matrix<T2>& object);

		template<typename T2>
		Matrix& operator*=(const Matrix<T2>& object);

		Matrix& operator*=(double digit);
		~Matrix();
		Matrix operator -() const;
	private:
		template<typename T2>
		void copy_constr(const Matrix<T2>& obj);
	private:
		void swapRows(size_t row_index, size_t max_row, size_t col_index);
		/*void swapRows(size_t row1, size_t row2);*/
		T* m_ptr=nullptr;
		int m_rows=0;
		int m_columns=0;
		int m_capacity = 0;

	};

	template<typename T1, typename  T2>
	auto operator*(const Matrix<T1> object, T2 value);

	template<typename T1, typename  T2>
	auto operator*(T1 value, const Matrix<T2> object);

	template<typename T1, typename  T2>
	auto operator*(Matrix<T1>& object1, Matrix<T2>& object2);

	template<typename T1, typename  T2>
	auto operator+(const Matrix<T1>& object1, const Matrix<T2>& object2);

	template<typename T1, typename  T2>
	auto operator-(const Matrix<T1>& object1, const Matrix<T2>& object2);

	template<typename T1, typename  T2>
	bool operator==(const Matrix<T1>& object1, const Matrix<T2>& object2) noexcept;

	template<typename T1, typename T2>
	bool operator!=(const Matrix<T1>& object1, const Matrix<T2>& object2) noexcept;

	template<typename T1>
	std::ostream& operator<<(std::ostream& out, const Matrix<T1>& matrix);

	template<typename T1, typename T2>
	auto concatenate(const Matrix<T1>& left, const Matrix<T2>& right);//??

	template<typename T1>
	Matrix<T1> transpose(Matrix<T1> matr);

	template<typename T1>
	Matrix<T1> invert(Matrix<T1> matr);

	template<typename T1>
	Matrix<T1> power(Matrix<T1>& matr, int pow);//??

	template<typename T1>
	Matrix<T1> solve(const Matrix<T1>& matr_a, const Matrix<T1>& vek_f);//??
}
#include "Matrix.hpp"
