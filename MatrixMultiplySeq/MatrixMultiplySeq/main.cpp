#include "/opt/openmpi/include/mpi.h"
#include <iostream>
#include <fstream>

void Serial(int matrix_size, double *matrix_a, double *matrix_b, double *matrix_res)
{
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            matrix_res[matrix_size * i + j] = 0;
            for (int k = 0; k < matrix_size; k++) {
                matrix_res[matrix_size * i + j] += matrix_a[i * matrix_size + k] * matrix_b[k * matrix_size + j];
            }
        }
    }
}

int main(int argc, char **argv)
{
    int procRank, procNum;
    MPI_Init(&argc, &argv);
    
    double seq_start, seq_stop;
    
    int matrix_size = 3000;
    double *matrix_a = new double[matrix_size*matrix_size];
    double *matrix_b = new double[matrix_size*matrix_size];
    double *matrix_res = new double[matrix_size*matrix_size];
    
    std::ifstream in_a("/Users/kitaev/Documents/Desktop/Тех/7 семестр/Практика/matrix3000_a.txt");
    std::ifstream in_b("/Users/kitaev/Documents/Desktop/Тех/7 семестр/Практика/matrix3000_b.txt");
    
    for(int i = 0; i < matrix_size; i++) {
        in_a >> matrix_a[i];
        in_b >> matrix_b[i];
    }
    
    in_a.close();
    in_b.close();
    
    seq_start = MPI_Wtime();
    Serial(matrix_size, matrix_a, matrix_b, matrix_res);
    seq_stop = MPI_Wtime();
    
    std::cout << "Seq time: " << seq_stop - seq_start << std::endl;
    
    
    delete[] matrix_a;
    delete[] matrix_b;
    delete[] matrix_res;
    
    
    MPI_Finalize();
    return 0;
}
