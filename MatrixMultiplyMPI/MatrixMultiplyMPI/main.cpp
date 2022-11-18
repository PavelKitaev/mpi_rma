#include "/opt/openmpi/include/mpi.h"
#include <iostream>
#include <fstream>

void PrintMatrixAB(int matrix_size, double *matrix_a, double *matrix_b)
{
    for (int i = 0; i < matrix_size; i++){
        for (int j = 0; j < matrix_size; j++) {
            std::cout << matrix_a[matrix_size * i + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
    for (int i = 0; i < matrix_size; i++){
        for (int j = 0; j < matrix_size; j++) {
            std::cout << matrix_b[matrix_size * i + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void PrintMatrixRes(int matrix_size, double *matrix)
{
    for (int i = 0; i < matrix_size; i++){
        for (int j = 0; j < matrix_size; j++) {
            std::cout << matrix[matrix_size * i + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void GetData(int *sendcounts, int *displs, int matrix_size, int procNum)
{
    if (matrix_size % procNum != 0) {
        int temp = matrix_size - 1;
        for (int i = 1; i < matrix_size; i++) {
            if (temp % procNum == 0) {
                break;
            }
            else {
                temp = temp - 1;
            }
        }
        
        for (int i = 0; i < procNum; i++) {
            sendcounts[i] = (temp / procNum) * matrix_size;
        }
        
        int residue = matrix_size - temp;
        int i = 0;
        
        while (residue != 0) {
            sendcounts[i] += matrix_size;
            residue--;
            i++;
            
            if (i == procNum) {
                i = 0;
            }
        }
    }
    else {
        for (int i = 0 ; i < procNum; i++) {
            sendcounts[i] = (matrix_size / procNum) * matrix_size;
        }
    }
    
    displs[0] = 0;
    for (int i = 1; i < procNum; i++) {
        displs[i] = sendcounts[i-1] + displs[i-1];
    }
}

bool Comparison(double *matrix_a, double *matrix_b, int size)
{
    bool flag = true;
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (matrix_a[size * i + j] != matrix_b[size * i + j]) {
                return false;
            }
        }
    }
    
    return flag;
}

int main(int argc, char **argv)
{
    int procRank, procNum;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    int matrix_size = 3000;
    
    double *matrix_a = nullptr;
    double *matrix_b = new double[matrix_size*matrix_size];
    double *matrix_res_parallel = nullptr;
    double *local_matrix = nullptr;
    
    int *sendcounts = nullptr;
    int *displs = nullptr;
    
    int sendcount = -1;
    int displ = -1;
    int local_size = -1;
    
    double par_start = -1, par_stop = -1;
    
    if (procRank == 0) {
        matrix_a = new double[matrix_size*matrix_size];
        matrix_res_parallel = new double[matrix_size*matrix_size];
        
        sendcounts = new int[procNum];
        displs = new int[procNum];
        
        std::ifstream in_a("/Users/kitaev/Documents/Desktop/Тех/7 семестр/Практика/matrix3000_a.txt");
        std::ifstream in_b("/Users/kitaev/Documents/Desktop/Тех/7 семестр/Практика/matrix3000_b.txt");
        
        for(int i = 0; i < matrix_size; i++) {
            in_a >> matrix_a[i];
            in_b >> matrix_b[i];
        }
        
        in_a.close();
        in_b.close();
        
        par_start = MPI_Wtime();
        
        if (matrix_size < procNum) {
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
        GetData(sendcounts, displs, matrix_size, procNum);
        
        sendcount = sendcounts[0];
        displ = displs[0];
        
        char buff[100];
        for (int i = 1; i < procNum; i++) {
            int pos = 0;
            MPI_Pack(&sendcounts[i], 1, MPI_INT, buff, 100, &pos, MPI_COMM_WORLD);
            MPI_Pack(&displs[i], 1, MPI_INT, buff, 100, &pos, MPI_COMM_WORLD);
            
            MPI_Send(buff, pos, MPI_PACKED, i, 0, MPI_COMM_WORLD);
        }
    }
    
    if (procRank != 0) {
        int temp[2];
        MPI_Recv(&temp, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
        
        sendcount = temp[0];
        displ = temp[1];
    }
    
    local_matrix = new double[sendcount];
    MPI_Bcast(&matrix_b[0], matrix_size*matrix_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&matrix_a[0], sendcounts, displs, MPI_DOUBLE, &local_matrix[0], sendcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    double *local_res = new double[sendcount];
    int count = 0;
    
    for (int i = 0; i < sendcount/matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            local_res[count] = 0;
            for (int k = 0; k < matrix_size; k++) {
                local_res[count] += local_matrix[i * matrix_size + k] * matrix_b[k * matrix_size + j];
            }
            count++;
        }
    }
    
    MPI_Gatherv(&local_res[0], sendcount, MPI_DOUBLE, &matrix_res_parallel[0], sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (procRank == 0) {
        par_stop = MPI_Wtime();

        std::cout << "Par time: " << par_stop - par_start << std::endl;
        
        //PrintMatrixAB(matrix_size, matrix_a, matrix_b);
        //PrintMatrixRes(matrix_size, matrix_res_parallel);
        
        delete[] matrix_a;
        delete[] matrix_res_parallel;
        delete[] sendcounts;
        delete[] displs;
    }

    delete[] matrix_b;
    delete[] local_res;
    delete[] local_matrix;
    
    MPI_Finalize();
    return 0;
}
