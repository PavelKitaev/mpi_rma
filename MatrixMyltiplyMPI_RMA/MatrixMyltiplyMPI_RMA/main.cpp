#include "/opt/openmpi/include/mpi.h"
#include <iostream>
#include <fstream>

void PrintMatrixAB(int matrix_size, double *matrix_a, double *matrix_b)
{
  std::cout << std::endl << std::endl;
  
  for (int i = 0; i < matrix_size; i++){
    for (int j = 0; j < matrix_size; j++) {
      std::cout << matrix_a[matrix_size * i + j] << " ";
    }
    std::cout << std::endl;
  }
  
  std::cout << std::endl << std::endl;
  
  for (int i = 0; i < matrix_size; i++){
    for (int j = 0; j < matrix_size; j++) {
      std::cout << matrix_b[matrix_size * i + j] << " ";
    }
    std::cout << std::endl;
  }
  
  std::cout << std::endl << std::endl;
}

void PrintMatrixRes(int matrix_size, double *matrix)
{
  std::cout << std::endl << std::endl;
  
  for (int i = 0; i < matrix_size; i++){
    for (int j = 0; j < matrix_size; j++) {
      std::cout << matrix[matrix_size * i + j] << " ";
    }
    std::cout << std::endl;
  }
  
  std::cout << std::endl << std::endl;
}

bool Comparison(double *matrix_a, double *matrix_b, int size)
{
    bool flag = true;
    
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++) {
            if (matrix_a[size * i + j] != matrix_b[size * i + j]) {
                return false;
            }
        }
    }
    
    return flag;
}

void GetSendcounts(int procNum, int size, int *sendcounts, int *displs)
{
  if (size % procNum != 0) {
    int temp = size - 1;
    for (int i = 1; i < size; i++) {
      if (temp % procNum == 0) {
        break;
      }
      else {
        temp = temp - 1;
      }
    }
    
    for (int i = 0; i < procNum; i++) {
      sendcounts[i] = (temp / procNum) * size;
    }
    
    int residue = size - temp;
    int i = 0;
    
    while (residue != 0) {
      sendcounts[i] += size;
      residue--;
      i++;
      
      if (i == procNum) {
        i = 0;
      }
    }
  }
  else {
    for (int i = 0 ; i < procNum; i++) {
      sendcounts[i] = (size / procNum) * size;
    }
  }
  
  displs[0] = 0;
  for (int i = 1; i < procNum; i++) {
    displs[i] = sendcounts[i-1] + displs[i-1];
  }
}

int main(int argc, char *argv[])
{
    int procRank, procNum;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);            //Количество потоков
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);           //Номер потока
    
    int matrix_size = 3000;
    
    double *matrix_a = nullptr;
    double *matrix_b = new double[matrix_size*matrix_size];
    double *matrix_res_parallel;
    double *local_matrix = nullptr;
    double *local_res = nullptr;
    int *sendcounts = nullptr;
    int *displs = nullptr;
    
    int local_size = -1;
    int local_displ = -1;
    
    double par_start = -1, par_stop = -1;
    
    if (procRank == 0) {
        
        
        matrix_a = new double[matrix_size*matrix_size];
        
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
        
        sendcounts = new int[procNum];
        displs = new int[procNum];
        
        GetSendcounts(procNum, matrix_size, sendcounts, displs);
        
        local_size = sendcounts[0];
        local_displ = displs[0];
    }
    
    MPI_Win win_res;          //Создаем окно доступа к матрице с результатом
    MPI_Win_allocate((MPI_Aint)(matrix_size * matrix_size * sizeof(double)), sizeof(double),
                     MPI_INFO_NULL, MPI_COMM_WORLD, &matrix_res_parallel, &win_res);
    MPI_Win_fence(0, win_res); //Открываем доступ к матрице с результатом
    
    MPI_Win win_a;            //Создаем окно доступа к матрице А
    MPI_Win_create(&matrix_a[0], (matrix_size*matrix_size) * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_a);
    MPI_Win_fence(0, win_a);    //Открываем доступ к матрице А
    
    MPI_Win win_b;            //Создаем окно доступа к матрице А
    MPI_Win_create(&matrix_b[0], (matrix_size*matrix_size) * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_b);
    MPI_Win_fence(0, win_b);    //Открываем доступ к матрице В
    
    MPI_Win win_sc;           //Создаем окно доступа к sendcounts
    MPI_Win_create(&sendcounts[0], procNum * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_sc);
    MPI_Win_fence(0, win_sc);   //Открываем доступ к sendcounts
    
    MPI_Win win_dl;           //Создаем окно доступа к displs
    MPI_Win_create(&displs[0], procNum * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_dl);
    MPI_Win_fence(0, win_dl);   //Открываем доступ к displs
    
    if (procRank != 0) {
        MPI_Get(&local_size, 1, MPI_INT, 0, procRank, 1, MPI_INT, win_sc); //Берем данные из sendcounts нулевого потока
        MPI_Get(&local_displ, 1, MPI_INT, 0, procRank, 1, MPI_INT, win_dl); //Берем данные из sendcounts нулевого потока
        MPI_Get(&matrix_b[0], matrix_size * matrix_size, MPI_DOUBLE, 0, 0,
                matrix_size * matrix_size, MPI_DOUBLE, win_b); //Берем матрицу В из нулевого потока
    }
    
    MPI_Win_fence(0, win_sc); //Закрываем окно к sendcounts
    MPI_Win_fence(0, win_dl); //Закрываем окно к displs
    MPI_Win_fence(0, win_b); //Закрываем окно к матрице А
    
    local_matrix = new double[local_size];
    MPI_Get(&local_matrix[0], local_size, MPI_DOUBLE, 0, local_displ, local_size, MPI_DOUBLE, win_a);
    
    MPI_Win_fence(0, win_a); //Закрываем окно к матрице B
    
    local_res = new double[local_size];
    int count = 0;
    
    for (int i = 0; i < local_size/matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            local_res[count] = 0;
            for (int k = 0; k < matrix_size; k++) {
                local_res[count] += local_matrix[i * matrix_size + k] * matrix_b[k * matrix_size + j];
            }
            count++;
        }
    }
    
    //Кладем локальные результаты в основную матрицу с результатом
    MPI_Put(&local_res[0], local_size, MPI_DOUBLE, 0, local_displ, local_size, MPI_DOUBLE, win_res);
    
    MPI_Win_fence(0, win_res); //Закрывам доступ к матрице с результатом
    
    if (procRank == 0) {
        par_stop = MPI_Wtime();
        
        std::cout << "Par time: " << par_stop - par_start << std::endl;
        
        delete[] sendcounts;
        delete[] displs;
        delete[] matrix_a;
    }
    
    delete[] matrix_b;
    delete[] local_matrix;
    delete[] local_res;
    
    //Удаляем окна
    MPI_Win_free(&win_a);
    MPI_Win_free(&win_b);
    MPI_Win_free(&win_sc);
    MPI_Win_free(&win_dl);
    MPI_Win_free(&win_res); //--> delete[] matrix_res_parallel;
    
    MPI_Finalize();
    return 0;
}
