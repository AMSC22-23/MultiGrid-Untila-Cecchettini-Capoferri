#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <functional>

// Prodotto scalare di due vettori
double dot(const std::vector<double>& v1, const std::vector<double>& v2) {
    assert(v1.size() == v2.size());
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i)
        sum += v1[i] * v2[i];
    return sum;
}

// Norm (norma euclidea) di un vettore
double norm(const std::vector<double>& v) {
    return std::sqrt(dot(v, v));
}

// Moltiplicazione matrice-vettore
std::vector<double> matVecMul(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    size_t n = A.size();
    std::vector<double> result(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            result[i] += A[i][j] * x[j];
    return result;
}

// Metodo del gradiente coniugato
std::vector<double> conjugateGradient(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int max_iter = 1000, double tol = 1e-6) {
    int n = b.size();
    std::vector<double> x(n, 0.0); // Stima iniziale x = 0

    std::vector<double> r = b; // Residuo r = b - A*x (ma x = 0 inizialmente)
    std::vector<double> p = r; // Direzione di discesa iniziale p = r
    double rs_old = dot(r, r);

    for (int i = 0; i < max_iter; ++i) {
        std::vector<double> Ap = matVecMul(A, p);
        double alpha = rs_old / dot(p, Ap);

        // Aggiorna la soluzione x
        for (int j = 0; j < n; ++j)
            x[j] += alpha * p[j];

        // Aggiorna il residuo r
        for (int j = 0; j < n; ++j)
            r[j] -= alpha * Ap[j];

        double rs_new = dot(r, r);

        // Verifica la convergenza
        if (std::sqrt(rs_new) < tol)
            break;

        // Aggiorna la direzione di discesa p
        for (int j = 0; j < n; ++j)
            p[j] = r[j] + (rs_new / rs_old) * p[j];

        rs_old = rs_new;
    }

    return x;
}



int main(int argc, char* argv[]) {
    // Esempio di utilizzo

    int x1 = 10;
    int y1 = 20;
    double (*funcPtr)(double, double);

    funcPtr = createMultiplier();

    int result = funcPtr(x1, y1);
    std::cout << "Result: " << result << std::endl;

    std::vector<std::vector<double>> A = {
        {4, 1, 0},
        {1, 3, 1},
        {0, 1, 2}
    };
    std::vector<double> b = {4, 5, 3};

    std::vector<double> x = conjugateGradient(A, b);

    std::cout << "Soluzione x: ";
    for (const auto& val : x) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}
