#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 10000
#define TOLERANCE 1e-6

double **allocate_matrix(int rows, int cols) {
    double **matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));
    }
    return matrix;
}

void free_matrix(double **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void initialize_temperature(double **T, int rows, int cols, double T_left, double T_right, double T_top, double T_bottom) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == 0) {
                T[i][j] = T_top;
            } else if (i == rows - 1) {
                T[i][j] = T_bottom;
            } else if (j == 0) {
                T[i][j] = T_left;
            } else if (j == cols - 1) {
                T[i][j] = T_right;
            } else {
                T[i][j] = 0.0;
            }
        }
    }
}

double heat_generation(double x, double y) {
    // Heat generation function (w/m³)
    // You can modify this function as per your requirements
    return 1000.0;
}

void solve_heat_equation(double **T, int rows, int cols, double dx, double dy, double alpha, double q) {
    double **T_new = allocate_matrix(rows, cols);
    double convergence_factor = 1.0;
    int iter = 0;

    while (iter < MAX_ITER) {
        double max_diff = 0.0;

        for (int i = 1; i < rows - 1; i++) {
            for (int j = 1; j < cols - 1; j++) {
                double T_curr = T[i][j];
                double T_left = T[i][j - 1];
                double T_right = T[i][j + 1];
                double T_top = T[i - 1][j];
                double T_bottom = T[i + 1][j];
                double heat_gen = heat_generation((j - 1) * dx, (rows - i - 1) * dy);

                T_new[i][j] = convergence_factor * (
                    (T_left + T_right) / (dx * dx) +
                    (T_top + T_bottom) / (dy * dy) +
                    heat_gen / (alpha * alpha)
                ) / (2.0 / (dx * dx) + 2.0 / (dy * dy));

                double diff = fabs(T_new[i][j] - T_curr);
                if (diff > max_diff) {
                    max_diff = diff;
                }
            }
        }

        double **temp = T;
        T = T_new;
        T_new = temp;

        iter++;

        if (max_diff < TOLERANCE) {
            break;
        }
    }

    free_matrix(T_new, rows);
}

int main() {
    int rows = 11, cols = 11;
    double dx = 0.1, dy = 0.1;
    double alpha = 1.0;
    double T_left = 100.0, T_right = 200.0, T_top = 300.0, T_bottom = 400.0;

    double **T = allocate_matrix(rows, cols);
    initialize_temperature(T, rows, cols, T_left, T_right, T_top, T_bottom);

    solve_heat_equation(T, rows, cols, dx, dy, alpha, 1000.0);

    // Print the temperature distribution
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.2f ", T[i][j]);
        }
        printf("\n");
    }

    free_matrix(T, rows);

    return 0;
}