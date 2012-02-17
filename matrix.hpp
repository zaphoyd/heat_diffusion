#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T>
std::ostream& operator<< (std::ostream& out, std::vector<T> a) {
    out << "[";
    
    std::string col_separator = "";
    
    out.precision(6);
    out.setf(std::ios::fixed,std::ios::floatfield);

    for (int i = 0; i < a.size(); i++) {
        out << col_separator << a[i];
        col_separator = ", ";
    }
    
    out << "]";

    return out;
}

template <typename T>
class matrix {
public:
    matrix(size_t nx, size_t ny) : m_nx(nx), m_ny(ny), m_data(nx) {
        for (size_t i = 0; i < nx; i++) {
            m_data[i].resize(ny);
        }
    }

    std::vector<T>& operator[] (size_t x) {
        return m_data[x];
    }

    const std::vector<T>& operator[] (size_t x) const {
        return m_data[x];
    }
    
    size_t nx() const {return m_nx;}
    size_t ny() const {return m_ny;}


    friend std::ostream& operator<< (std::ostream &out, matrix<T>& m) {
        out << "[";
        
        bool separator = false;
        
        for (int i = 0; i < m.nx(); i++) {
            if (separator) {
                out << std::endl;
            } else {
                separator = true;
            }

            out << m[i];
        }
        
        out << "]";

        return out;
    }
private:
    size_t m_nx;
    size_t m_ny;
    std::vector<std::vector<T>> m_data;
};


template <typename T>
void upper_triangulate(matrix<T>& A,std::vector<T>& b) {
    T scale_factor;
    size_t X = A.nx();

    if (A.nx() != A.ny() || A.nx() != b.size()) {
        throw std::invalid_argument("upper triangulate: matrix must be square and the same size as vector b");
    }

    // loop over columns
    for (int j = 0; j < X; ++j) {
        // for each row under the diagonal
        
        // find largest value in column under diagonal
        size_t max_col = 0;
        for (int k = j; k < X; ++k) {
            if (max_col == 0 || std::abs(A[k][j]) > std::abs(A[max_col][j])) {
                max_col = k;
            }
        }
        
        if (max_col == 0) {
            // there are no rows under the pivot
            continue;
        }
        
        if (max_col > 0) {
            // comment these out to not use largest abs for pivot
            std::swap(A[j],A[max_col]);
            std::swap(b[j],b[max_col]);
        }
            
        for (int i = j+1; i < X; ++i) {
            // eliminate
            scale_factor = A[i][j]/A[j][j];
            for (int k = 0; k < X; k++) {
                A[i][k] -= scale_factor*A[j][k];
            }
            b[i] -= scale_factor*b[j];

        }
    }
}

template <typename T>
void gauss_jordan(matrix<T>& A) {
    T scale_factor;
    
    size_t X = A.nx();
    size_t Y = A.ny();
    size_t size = std::min(X,Y);
    
    // loop over columns
    for (int j = 0; j < size; ++j) {
        // for each row under the diagonal
        for (int i = j+1; i < size; ++i) {
            // find largest value in column under diagonal
            size_t max_col = 0;
            for (int k = i; k < size; ++k) {
                if (max_col == 0 || std::abs(A[k][j]) > std::abs(A[max_col][j])) {
                    max_col = k;
                }
            }

            if (max_col == 0) {
                // there are no rows under the pivot
                continue;
            }
            

            if (max_col > i) {
                std::swap(A[i],A[max_col]);
            } else if (max_col == 0) {
                std::cout << "row finished early" << std::endl;
                continue;
            }
            
            // eliminate
            scale_factor = A[i][j]/A[j][j];
            for (int k = 0; k < Y; k++) {
                A[i][k] -= scale_factor*A[j][k];
            }
        }
    }
    
    // eliminate above
    for (int j = size-1; j >= 0; --j) {
        // for each row under the diagonal
        for (int i = j-1; i >= 0; --i) {
            // TODO find largest value in column above diagonal
            
            // eliminate
            scale_factor = A[i][j]/A[j][j];
            for (int k = Y-1; k >= 0; --k) {
                A[i][k] -= scale_factor*A[j][k];
            }
        }
    }
    
    // scale all to 1.0
    for (int i = 0; i < X; ++i) {
        scale_factor = 1.0/A[i][i];
        for (int j = 0; j < Y; ++j) {
            A[i][j] *= scale_factor;
        }
    }
}

template <typename T>
void back_sub(const matrix<T>& A,const std::vector<T>& b, std::vector<T>& x) {
    size_t X = A.nx();

    if (A.ny() != X || b.size() != X || x.size() != X) {
        throw std::invalid_argument("back_sub: matrix must be square and the same size as vectors b and x.");
    }

    x[X-1] = b[X-1]/A[X-1][X-1];
    
    for (int i = X-2; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i+1; j < X; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

