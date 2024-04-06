#include "Matrix.h"
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <sstream>

template<typename T>
template<typename T2>
void linalg::Matrix<T>::copy_constr(const Matrix<T2>& obj) {
    m_ptr = reinterpret_cast<T*>(operator new(obj.m_rows * obj.m_cols * sizeof(T)));
    T* ptr_el = m_ptr;
    try {
        for (; ptr_el!=m_ptr+obj.m_columns*obj.m_rows; ++ptr_el) {
                new (ptr_el) T(obj.m_ptr[ptr_el - m_ptr]);
            }
    }
    catch (...) {
        while (--ptr_el >= m_ptr)
            ptr_el->~T();
        operator delete (reinterpret_cast<void*>(m_ptr));
        throw;
    }
    m_rows = obj.m_rows;
    m_cols = obj.m_cols;
    m_capacity = m_rows*m_columns;

}
template<typename T>
template<typename T2>
linalg::Matrix<T>::Matrix(const Matrix<T2>& obj) {
    copy_constr(obj);
}
template<typename T>
linalg::Matrix<T>::Matrix(const Matrix& obj) {
    copy_constr(obj);
}

template<typename T>
linalg::Matrix<T>::Matrix(int rows) :m_rows(rows) {//???? НАДО ИЛИ НЕТ??? ЕСЛИ ДА-ИСПРАВИТЬ КАПАСИТИ
    if (m_rows <= 0) {
        throw std::runtime_error("Размеры матрицы неккоректны");
    }
    m_columns = 1;
    m_ptr = reinterpret_cast<T*>(operator new (sizeof(T) * rows));
    T* ptr = m_ptr;
    try {
        for(; ptr!=m_ptr+m_rows; ++ptr) {
            new(ptr) T();
        }
    }
    catch (...) {
        while (--ptr >= m_ptr)
            ptr->~T();
        operator delete (reinterpret_cast<void*>(m_ptr));
        throw;
    }
    m_rows = rows;
    m_columns = 1;
    m_ptr = ptr;
    m_capacity = m_rows;
}

template<typename T>
linalg::Matrix<T>::Matrix(int rows, int cols) : m_rows(rows), m_columns(cols) {//???? НАДО ИЛИ НЕТ??? ЕСЛИ ДА-ИСПРАВИТЬ КАПАСИТИ
    if (m_rows <= 0 || m_columns <= 0) {
        throw std::runtime_error("Размеры матрицы некорректны");
    }
    m_ptr = reinterpret_cast<T*>(operator new (sizeof(T) * rows * columns));
    T* ptr = m_ptr;
    try {
        for (; ptr != m_ptr + m_rows*m_columns; ++ptr) {
            new(ptr) T();
        }
    }
    catch (...) {
        while (--ptr >= m_ptr)
            ptr->~T();
        operator delete (reinterpret_cast<void*>(m_ptr));
        throw;
    }
    m_rows = rows;
    m_columns = columns;
    m_ptr = ptr;
    m_capacity = m_rows*m_columns;
}
template<typename T>
linalg::Matrix<T>::Matrix(Matrix&& object) noexcept {// конструктор перемещения
    std::swap(m_ptr, object.m_ptr);
    std::swap(m_rows, object.m_rows);
    std::swap(m_columns, object.m_columns);
    std::swap(m_capacity, object.m_capacity);
}
template<typename T>
template<typename T2>
linalg::Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T2>> list) {
    size_t columns = list.begin()->size();
    for (std::initializer_list<T2> row : list) {
        if (row.size() != columns) {
            throw std::runtime_error("Incorrect sizes of rows");
        }
    }
    if (list.begin()->size() == 0) {
        return;
    }
    size_t i = 0;
    m_ptr = reinterpret_cast<T*>(operator new (sizeof(T) * list.size() * list.begin()->size()));
    T* ptr_el = m_ptr;
    try {
        for (std::initializer_list<T2> row : list) {
            for (const T2& element : row) {
                new(ptr_el) T(element);
                ++ptr_el;
            }
        }
    }
    catch (...) {
        while (--ptr_el >= m_ptr)
            ptr_el->~T();
        operator delete (reinterpret_cast<void*>(m_ptr));
        throw;
    }
    m_rows = list.size();
    m_columns = columns;
    m_capacity = m_rows*m_columns;
}

template <typename T>
template <typename T2>
linalg::Matrix<T>::Matrix(std::initializer_list<T2> lst) {
    if (lst.size() == 0) {
        return;
    }
    m_ptr = reinterpret_cast<T*>(operator new(lst.size() * sizeof(T)));
    T* ptr_el = m_ptr;
    try {
        for (const T2& el : lst) {
            new (ptr_el) T(el);
            ++ptr_el;
        }
    }
    catch (...) {
        while (--ptr_el >= m_ptr)
            ptr_el->~T();
        operator delete (reinterpret_cast<void*>(m_ptr));
        throw;
    }
    m_rows = lst.size();
    m_columns = 1;
    m_capacity = m_rows * m_columns;
}


template<typename T>
linalg::Matrix<T>::~Matrix() {
    for (size_t i = 0; i < m_rows*m_columns; ++i) {
        (m_ptr + i)->~T();
    }
    operator delete (reinterpret_cast<void*>(m_ptr));
}

template<typename T>
void linalg::Matrix<T>::reshape(int rows, int cols) {
    if (rows * cols != m_rows * m_columns) {
        throw std::runtime_error("Нельзя изменить размерность сохраняя количество элементов");
    }
    m_rows = rows;
    m_columns = cols;
}

template<typename T>
template<typename T2>
linalg::Matrix<T>& linalg::Matrix<T>::operator=(const Matrix<T2>& object) {
    if (m_capacity < object.size())
        return *this = Vector<T2>{ object };
    for (size_t i = 0; i < std::min(m_rows*m_columns, object.size()); ++i)
        m_ptr[i] = object.m_ptr[i];
    T* ptr_el = m_ptr + m_rows*m_columns;
    if (m_rows*m_columns < object.size()) {
        try {
            for (T2* ptr_obj = obj.m_ptr + m_rows*m_columns;
                ptr_obj < object.m_ptr + object.size();
                ++ptr_obj, ++ptr_el) new(ptr_el) T(*ptr_obj);
        }
        catch (...) {
            while (--ptr_el >= m_ptr + m_rows*m_columns)
                ptr_el->~T();
            throw;
        }
    }
    else {
        while (--ptr_el >= m_ptr + object.size())
            ptr_el->~T();;
    }
    m_rows = object.m_rows;
    m_columns = object.m_columns;
    return *this;

}

template<typename T>
linalg::Matrix<T>& linalg::Matrix<T>::operator=(const Matrix<T>& obj) {
    if (this == &obj)
        return *this;
    return operator=<T>(obj);
}


template <typename T>
void linalg::Matrix<T>::reserve(size_t n) {
    if (n <= m_capacity) {
        return;
    }
    T* new_ptr = reinterpret_cast<T*>(operator new(n * sizeof(T)));
    for (int i = 0; i < m_rows * m_columns; ++i) {
        new (&new_ptr[i]) T(std::move(m_ptr[i]));
    }
    for (size_t i = 0; i < m_rows * m_columns; ++i)
        (m_ptr + i)->~T();
    operator delete (reinterpret_cast<void*>(m_ptr));
    m_ptr = new_ptr;
    m_capacity = n;
}


template<typename T>
void linalg::Matrix<T>::shrink_to_fit() {
    if (m_rows * m_columns == m_capacity) {
        return;
    }
    T* ptr_m = m_ptr;
    T* new_ptr = reinterpret_cast<T*>(operator new(sizeof(T) * m_rows * m_columns));
    for (T* ptr = new_ptr; ptr_m < m_ptr + m_rows * m_columns; ++ptr, ++ptr_m) {
        new(ptr) T(*ptr_m);
    }
    for (T* ptr = m_ptr; ptr != m_ptr + m_capacity; ++ptr) {
        ptr->~T();
    }
    operator delete (reinterpret_cast<void*>(m_ptr));
    m_capacity = m_rows * m_columns;
    m_ptr = new_ptr;

}

template <typename T>
void linalg::Matrix<T>::clear() {
    for (T* ptr_new = m_ptr; ptr_new != m_ptr + m_rows * m_columns; ++ptr_new) {
        ptr_new->~T();
    }

    m_rows = 0;
    m_columns = 0;
}


template<typename T>
 linalg::Matrix<T>& linalg::Matrix<T>::operator=(Matrix&& object) noexcept {//оператор присваивания с перемещением
     std::swap(m_ptr, object.m_ptr);
     std::swap(m_rows, object.m_rows);
     std::swap(m_columns, object.m_columns);
     std::swap(m_capacity, object.m_capacity)
     return (*this);
 }
 template<typename T>
 T& linalg::Matrix<T>::operator()(size_t row, size_t col) {// оператор вызова функции для изменения значения
     if (row >= m_rows || col >= m_columns) {
         throw std::runtime_error("Неверный номер элемента");
     }

     return m_ptr[row * m_columns + col];
 }
 template<typename T>
T linalg::Matrix<T>::operator()(size_t row, size_t col) const {
    if (row >= m_rows || col >= m_columns) {
        throw std::runtime_error("Неверный номер элементаЫ");
    }

    return m_ptr[row * m_columns + col];
 }
template<typename T1>
static size_t number_of_digits(T1 element, std::ios_base::fmtflags flags) {//набор форматных флагов (flags), которые будут применены при преобразовании числа в строку.
    std::ostringstream stream;//создается объект класса std::ostringstream, который представляет собой поток, используемый для форматированного вывода данных в строку.
    stream.flags(flags); //Здесь устанавливаются форматные флаги для объекта stream.
    //Форматные флаги влияют на способ преобразования числа в строку, включая количество знаков после запятой, ширину поля и другие параметры форматирования.
    stream << element;//Значение element выводится в поток stream. Применяются форматные флаги, установленные ранее, к преобразованию числа в строку.
    return stream.str().size();//Результат преобразования числа в строку извлекается с помощью вызова stream.str()
    
}
template<typename T1>
std::ostream& linalg:: operator << (std::ostream& out, const Matrix<T1>& obj) {
    if (obj.empty()) // Отдельно обработаем случай пустых массивов
        return out << "|empty|\n";
    size_t* max_digits = new size_t[obj.columns()]();
    size_t num_digits = 0;
    size_t digits_max = number_of_digits(obj(0, 0), out.flags());
    for (size_t j = 0; j < obj.columns(); ++j) {
        digits_max = 0;
        for (size_t i = 0; i < obj.rows(); ++i) {
            digits_max = std::max(digits_max, number_of_digits(obj(i, j), out.flags()));;
        }
        max_digits[j] = digits_max;
    }
  
    // Вывод матрицы с выравниванием по столбцам
    for (size_t i = 0; i < obj.rows(); ++i) {
        out << '|';
        for (size_t j = 0; j < obj.columns(); ++j) {
            out << std::setw(max_digits[j]) << obj(i, j);

            if (j != obj.columns() - 1) {
                out << ' ';
            }
        }
        out << "|\n";
    }
    delete[] max_digits;
    return out;
  
}

template<typename T>
template<typename T2>
linalg::Matrix<T>& linalg::Matrix<T>:: operator+=(const Matrix<T2>& object) {
     if (m_rows != object.m_rows || m_columns != object.m_columns) {
         throw std::runtime_error("Матрицы нельзя сложить");
     }

     for (size_t i = 0; i < m_rows; ++i) {
         for (size_t j = 0; j < m_columns; ++j) {
             (*this)(i, j) += object(i, j);
         }
     }

     return *this;
 }

template<typename T1, typename T2>
auto linalg::operator+(const Matrix<T1>& object1, const Matrix<T2>& object2) {
    Matrix<decltype(T1{} + T2{}) > result(object1);

    if (object1.rows() != object2.rows() || object1.columns() != object2.columns() || (object1.empty() && object2.empty())) {
        throw std::runtime_error("Матрицы нельзя сложить");
    }

    if (object1.empty() && !object2.empty()) {
        return object2;
    }
    if (object2.empty() && !object1.empty()) {
        return object1;
    }

    result += object2;

    return result;
}

template<typename T>
template<typename T2>
linalg::Matrix<T>&linalg::Matrix<T>:: operator-=(const Matrix<T2>& object) {
    if (m_rows != object.m_rows || m_columns != object.m_columns || object.empty() || (*this).empty()) {
        throw std::runtime_error("Матрицы нельзя вычитать");
    }

    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            (*this)(i, j) -= object(i, j);
        }
    }

    return *this;
}

template<typename T1, typename T2>
auto linalg:: operator-(const Matrix<T1>& object1, const Matrix<T2>& object2){
    Matrix<decltype(T1{} + T2{}) > result(object1);
    if ( (object1.rows() != object2.rows()) || (object1.columns() != object2.columns()) || (object1.empty() && object2.empty()) ) {
        throw std::runtime_error("Матрицы нельзя вычесть");
    }
    if (object2.empty() && !object1.empty()) {
        return object1;
    }

    result -= object2;

    return result;
}

template<typename T>
template<typename T2>
linalg::Matrix<T>& linalg::Matrix<T>:: operator*=(const Matrix<T2>& object) {
    if (m_columns != object.m_rows || object.empty() || (*this).empty()) {
        throw std::runtime_error("Матрицы нельзя перемножить");
    }

    Matrix result(m_rows, object.m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < object.m_columns; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < m_columns; ++k) {
                sum += (*this)(i, k) * object(k, j);
            }
            result(i, j) = sum;
        }
    }
    *this = result;
    return *this;
}

template<typename T1, typename  T2>
auto linalg:: operator*(Matrix<T1>& object1, Matrix<T2>& object2) {
    Matrix<decltype(T1{} + T2{}) > result(object1);
    if ((object1.columns() != object2.rows()) || object1.empty()||object2.empty()) {
        throw std::runtime_error("Невозможно умножить матрицы");
    }

    return result*=object2;
}

template<typename T1, typename  T2>
auto linalg::operator*(T1 value, const Matrix<T2> object) {
    Matrix<decltype(T1{} + T2{}) > result(object);
    return result *= value;
}


template<typename T1, typename  T2>
auto linalg:: operator*(const Matrix<T1> object, T2 value) {
    Matrix<decltype(T1{} + T2{}) > result(object);
    return result *= value;
}

template<typename T>
linalg::Matrix<T>& linalg::Matrix<T>:: operator*=(double value) {
    if ((*this).empty() == true) {
        throw std::runtime_error("Матрица пустая,нет смысла умножать на число");
    }
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            (*this)(i, j) *= value;
        }
    }
    return *this;
}

template<typename T1,typename T2>
bool linalg::operator==(const Matrix<T1>& object1, const Matrix<T2>& object2) noexcept{

   int n1 = object1.rows();
   int m1 = object1.columns();
   int n2 = object2.rows();
   int m2 = object2.columns();
    if (n1 != n2 || m1 != m2) {
        return false; // т.к. размерности не совпадают
    }

    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < m1; ++j) {
            if ((object1)(i, j) != object2(i, j)) {
                return false; //не совпадает хотя бы 1 элемент
            }
        }
    }

    return true; 
}
template<typename T1, typename T2>
bool linalg:: operator!=(const Matrix<T1>& object1, const Matrix<T2>& object2)noexcept{
    return !(object1==object2);// вернули не (оператор совпадения)=>вернули несовпадение, чтобы не писать лишний код
}
template<typename T>
T linalg::Matrix<T>:: norm() const {
    if (m_rows == 0 || m_columns == 0 || (*this).empty()) {
        throw std::invalid_argument("Размеры матрицы неккоректны");
    }

    double sum = 0.0;
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            sum += m_ptr[i * m_columns + j] * m_ptr[i * m_columns + j];
        }
    }
    return sqrt(sum);
}
template<typename T>
T linalg::Matrix<T>:: trace() const {
    if (m_rows != m_columns){
        throw std::runtime_error("Матрица не квадратная, посчитать след нельзя");
    }
    if (m_rows == 0 || m_columns == 0 || (*this).empty()) {
        throw std::invalid_argument("Размеры матрицы неккоректны");
    }
    double tracesum = 0.0;

    for (size_t i = 0; i < m_rows; ++i) {
        tracesum += m_ptr[i * m_columns + i];
    }
    return tracesum;
}
template<typename T>
T linalg::Matrix<T>::determinant() const {
    if (m_rows != m_columns) {
        throw std::runtime_error("Матрица не квадратная");
    }

    if (m_rows == 1) {
        return m_ptr[0]; // Для матрицы 1x1 определитель равен её единственному элементу
    }

    if (m_rows == 2) {
        return m_ptr[0] * m_ptr[3] - m_ptr[1] * m_ptr[2];
    }

    Matrix<T> copy(*this);
    copy.gauss_forward();
    double det = 1.0;
    for (size_t i = 0; i < copy.m_rows; ++i)
    {
        for (size_t j = 0; j < copy.m_columns; ++j) {
            if (i == j) {
                det *= copy(i, j);
            }
        }
    }
    return det;
}

template<typename T>
void linalg::Matrix<T>::swapRows(size_t row_index, size_t max_row, size_t col_index) {
    // Меняем местами текущую строку и строку с максимальным элементом
    for (size_t i = col_index; i < m_columns; ++i) {
        double temp = (*this)(row_index, i);
        (*this)(row_index, i) = (*this)(max_row, i);
        (*this)(max_row, i) = -temp; // Меняем знак строки
    }
}
/*void linalg::Matrix::swapRows(size_t row1, size_t row2) {
    // Меняем местами текущую строку и строку с максимальным элементом
    for (size_t i = row1; i < m_columns; ++i) {
        double temp = (*this)(row1, i);
        (*this)(row1, i) = (*this)(row2, i);
        (*this)(row2, i) = temp; // Меняем знак строки
    }
}*/

/*void linalg::Matrix::gauss_forward() {
    if ((*this).empty()) {
        throw std::runtime_error("Матрица пуста, не получится выполнить прямой ход Гаусса.");
    }
    for (size_t nowrow = 0; nowrow < m_rows; ++nowrow) {
        // Найти первый ненулевой элемент в текущем столбце
        size_t nozerorow = nowrow;//nozerow отвечает за строку, итерируется по столбцам

        while (nozerorow < m_rows && (*this)(nozerorow, nowrow) == 0.0) {
            ++nozerorow;
        }

        if (nozerorow == m_rows) {// не найдено ненулевых элементов в столбце=> переходим к следующему столбцу
            continue;
        }

        if (nozerorow != nowrow) {
            swapRows(nozerorow, nowrow);// поменять местами строки, чтобы перенести строку с ненулевым элементом  на диагональ
        }
        for (size_t j = nowrow + 1; j < m_rows; ++j) {
            double divide = (*this)(j, nowrow) / (*this)(nowrow, nowrow);
            for (size_t i = nowrow; i < m_columns; ++i) {//итерация по столбцам(идём по элементам одной строки)
                (*this)(j, i) -= divide * (*this)(nowrow, i);
            }
        }
    }
}*/
 template<typename T>
 void linalg::Matrix<T>::gauss_forward() {
        const double mistake = 1e-6;
        if ((*this).empty()) {
            throw std::runtime_error("Матрица пустая. Невозможно выполнить прямой ход метода Гаусса");
        }
        size_t row_index = 0, col_index = 0;

        while (row_index < m_rows && col_index < m_columns) {
            // Находим максимальный элемент в текущем столбце
            size_t max_row = row_index;
            double max_element = fabs((*this)(row_index, col_index));

            for (size_t i = row_index; i < m_rows; ++i) {
                double current_element = fabs((*this)(i, col_index));
                if (current_element > max_element) {
                    max_element = current_element;
                    max_row = i;
                }
            }

            if (max_element <=mistake) {
                for (size_t i = row_index; i < m_rows; ++i) {
                    (*this)(i, col_index) = 0.0;
                }
                ++col_index; // Переходим к следующему столбцу
            }
            else {
                if (max_row != row_index) {
                    // Меняем местами текущую строку и строку с максимальным элементом
                    swapRows(row_index, max_row, col_index);
                }

                for (size_t i = row_index + 1; i < m_rows; ++i) {
                    double factor = -((*this)(i, col_index) / (*this)(row_index, col_index));
                    (*this)(i, col_index) = 0.0;

                    for (size_t j = col_index + 1; j < m_columns; ++j) {
                        (*this)(i, j) += factor * (*this)(row_index, j);
                    }
                }

                ++row_index;
                ++col_index;
            }
        }

}
 template<typename T>
void linalg::Matrix<T>::gauss_backward() {
    if (empty()) {
        throw std::runtime_error("Матрица пустая. Невозможно выполнить обратный ход метода Гаусса");
    }

    for (int row = m_rows - 1; row >= 0; --row) {
        int nowColumn = -1; // Инициализируем столбец как -1

        // Находим первый ненулевой элемент в текущей строке
        for (size_t col = 0; col < m_columns; ++col) {
            if ((*this)(row, col)!=0) {
                nowColumn = col;
                break;
            }
        }

        if (nowColumn == -1) {
            continue;
        }

        // Нормализуем строку, чтобы текущий элемент стал равен 1
        double pivotValue = (*this)(row, nowColumn);
        for (size_t col = 0; col < m_columns; ++col) {
                (*this)(row, col) /= pivotValue;
        }

        // Обнуляем элементы выше текущего элемента в столбце
        for (int upper = row - 1; upper >= 0; --upper) {
            double value = - (*this)(upper, nowColumn);
            for (size_t col = 0; col < m_columns; ++col) {
                if (col == nowColumn) {
                    (*this)(upper, col) = 0;
                }
                else {
                    (*this)(upper, col) += value * (*this)(row, col);
                }
            }
        }
    }
}
template<typename T>
int linalg::Matrix<T>::rank() const {
    if (m_rows == 0 || m_columns == 0 || (*this).empty()) {
        throw std::invalid_argument("Размеры матрицы неккоректны");
    }
    int rank = 0;
    Matrix temp(*this);//копия,потому что метод Гаусса меняет матрицу
    temp.gauss_forward();
    for (size_t i = 0; i < temp.m_rows; ++i) {
        bool nonullrow = false;
        for (size_t j = 0; j < temp.m_columns; ++j) {
            if (temp(i, j) != 0.0) {
                nonullrow = true;
                break;
            }
        }
        if (nonullrow) {
            ++rank;
        }
    }
    return rank;
}

template<typename T1, typename T2>
auto linalg::concatenate(const Matrix<T1> &left, const Matrix<T2> &right) {

    if (left.empty()||right.empty()) {
        throw std::runtime_error("Матрица пустая. Невозможно выполнить обратный ход метода Гаусса");
    }
    if (left.rows() != right.rows()) {
        throw std::runtime_error("Невозможно объединить матрицы с разным количеством строк.");
    }

    int commonCols = left.columns() + right.columns();
    Matrix<decltype(T1{} + T2{}) > result(left.rows(), commonCols);
    for (int i = 0; i < left.rows(); ++i) {
        for (int j = 0; j < left.columns(); ++j) {
            result(i, j) = left(i, j);
        }
        for (int j = 0; j < right.columns(); ++j) {
            result(i, left.columns() + j) = right(i, j);
        }
    }

    return result;
}

template<typename T1>
linalg::Matrix<T1> linalg::transpose(Matrix<T1> matr) {
    if (matr.empty()) {
        throw std::runtime_error("Невозможно транспонировать");
    }
    int n = matr.rows();
    int m= matr.columns();
    Matrix result(m, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result(j, i) = matr(i, j);
        }
    }

    return result;
}

template<typename T1>
linalg::Matrix<T1> linalg:: invert(Matrix<T1> matr) {
    const double mistake = 1e-6;
    if (matr.empty()) {
        throw std::runtime_error("Матрица пустая");
    }
    if (matr.columns() != matr.rows()){
        throw std::runtime_error("Матрица не квадратная, невозможно найти обратную");
    }
    if (matr.determinant() == 0) {
        throw std::runtime_error("Невозможно найти обратную матрицу");
    }

    size_t n = matr.rows();
    size_t m = matr.columns();
    Matrix backMatr(n, 2 * n); 


    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < 2 * n; ++j) {
            if (j < n) {
                backMatr(i, j) = matr(i, j);
            }
            else if (j == (i + n)) {
                backMatr(i, j) = 1.0;
            }
            else {
                backMatr(i, j) = 0.0;
            }
        }
    }

    backMatr.gauss_forward();
    backMatr.gauss_backward();




    Matrix result(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result(i, j) = backMatr(i, j + n);
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (fabs(result(i, j)) <= mistake) {
                result(i, j) = 0;
            }
        }
    }

    return result;
}

template<typename T1>
linalg::Matrix<T1> linalg::power(Matrix<T1>& matr, int pow) {
    if (matr.empty()) {
        throw std::runtime_error("Матрица пустая");
    }
    if (matr.columns() != matr.rows()) {
        throw std::runtime_error("Матрица не квадратная, невозможно возвести в степень");
    }
    Matrix MatrPow(matr.columns(), matr.rows());
    int n = matr.columns();
    if (pow == 0) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    MatrPow(i, j) = 1.0;
                }
                else {
                    MatrPow(i, j) = 0.0;
                }
    
            }
        }
        return MatrPow;
    }

    if (pow == 1) {
        return matr;
    }

    if (pow > 0) {
        Matrix MatrPow = matr;
        for (int i = 0; i < pow - 1; ++i) {
            MatrPow = MatrPow * matr;
        }
        return MatrPow;
    }

    if (pow < 0) {
        MatrPow = invert(matr);
        for (int i = 0; i < -pow - 1; ++i) {
            MatrPow = MatrPow * matr;
        }
        return MatrPow;

    }
}

template<typename T1>
linalg::Matrix<T1> linalg::solve(const Matrix<T1>& matr_a, const Matrix<T1>& vec_f) {
    
    Matrix solution = concatenate(matr_a, vec_f);
    std::cout << solution.rank();
    std::cout << matr_a.rank();
    if (matr_a.empty() || vec_f.empty()){
        throw std::runtime_error("Ошибка ввода данных");
    }
    if (solution.rank() != matr_a.rank()){
        throw std::runtime_error("Система не имеет решений");
    }
    if ((solution.rank() == matr_a.rank()) && (solution.rank()!=matr_a.columns())){
        throw std::runtime_error("Система имеет бесконечно много решений");
    }
    solution.gauss_forward();
    solution.gauss_backward();
    linalg::Matrix solution_vector(vec_f.rows(), 1);
    for (size_t i = 0; i < solution_vector.rows(); i++) {
        solution_vector(i, 0) = solution(i, solution.columns() - 1);
    }

    return solution_vector;
}

/*Дополнительно*/
template<typename T>
linalg::Matrix<T> linalg::Matrix<T>::operator-() const { //Оператор унарного минуса
    linalg::Matrix result(m_rows, m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            result(i, j) = -(*this)(i, j);
        }
    }
    return result;
}
template<typename T>
double linalg::Matrix<T>::getMinor(int firstRow, int firstCol, int n) const {
    if (firstRow < 0 || firstCol < 0 || firstRow + n > m_rows || firstCol + n > m_columns) {
        throw std::runtime_error("Параметры некорректны");
    }

    linalg::Matrix minorMatrix(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            minorMatrix(i, j) = (*this)(firstRow + i, firstCol + j);
        }
    }
    return minorMatrix.determinant();
}
