// Copyright 2019 Dubovskoy Andrey

#include <omp.h>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <random>

unsigned int* get_randomized_array(int size) {
  unsigned int* array = new unsigned int[size];
  std::default_random_engine generator;
  std::uniform_int_distribution<unsigned int> distribution(0, 10000);
  for (int i = 0; i < size; i++) {
    array[i] = distribution(generator);
  }
  return array;
}

void print_array(unsigned int* array, int size) {
  if (size > 25) {
    return;
  }
  printf(" Array: ");
  for (int i = 0; i < size; i++) {
    printf(" %u", array[i]);
  }
  printf("\n");
}
void copy_array(unsigned int* array1, unsigned int* array2, int size) {
  for (int i = 0; i < size; i++) {
    array2[i] = array1[i];
  }
}

void countSort(unsigned int* arr, int n, int exp) {
  unsigned int *output = new unsigned int[n];
  int i;
  unsigned int count[10] = { 0 };

  for (i = 0; i < n; i++)
    count[(arr[i] / exp) % 10]++;
  for (i = 1; i < 10; i++)
    count[i] += count[i - 1];
  for (i = n - 1; i >= 0; i--) {
    output[count[(arr[i] / exp) % 10] - 1] = arr[i];
    count[(arr[i] / exp) % 10]--;
  }

  for (i = 0; i < n; i++)
    arr[i] = output[i];
}

void radix_sort(unsigned int *arr, int n) {
  unsigned int max = arr[0];
  for (int i = 1; i < n; i++)
    if (arr[i] > max)
      max = arr[i];
  for (int exp = 1; max / exp > 0; exp *= 10)
    countSort(arr, n, exp);
}

void check_result(unsigned int * array, int size) {
  for (int i = 1; i < size; i++)
    if (array[i] < array[i - 1]) {
      printf("Array not sorted \n");
      return;
    }
  printf("Array sorted\n");
}

unsigned int * Splitter(unsigned int * arr1, unsigned int* arr2, int size1, int size2) {
  unsigned int * tmp = new unsigned int[size1 + size2];
  int a = 0, b = 0, i = 0;
  while (a < size1 && b < size2) {
    if (arr1[a] < arr2[b]) {
      tmp[i++] = arr1[a++];
    } else {
      tmp[i++] = arr2[b++];
    }
  }
  while (a < size1) {
    tmp[i++] = arr1[a++];
  }
  while (b < size2) {
    tmp[i++] = arr2[b++];
  }
  return tmp;
}
int BinSearch(unsigned int *mas, int l, int r, unsigned int x) {
  if (l == r) return l;
  if (l + 1 == r) {
    if (x < mas[l]) return l;
    else
      return r;
  }

  int m = (l + r) / 2;
  if (x < mas[m]) {
    r = m;
  } else {
    if (x > mas[m]) l = m;
    else
      return m;
  }

  return BinSearch(mas, l, r, x);
}

void dac_sort(unsigned int * array, int size, int threads) {
  printf("point 1\n");
  int * piece_mas = new int[threads];
  int piece = size / threads;
  int remainder = size % threads;
  for (int i = 0; i < threads; i++) {
    piece_mas[i] = size / threads;
  }
  if (size / threads != 0) {
    piece_mas[threads - 1] = piece_mas[threads - 1] + remainder;
  }
  printf("point 2\n");
  omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < threads; i++)
    radix_sort(array + i * piece, piece_mas[i]);
  printf("point 3\n");
  print_array(array, size);


  int counter = static_cast<int>(std::log(threads) / std::log(2));
  printf(" Counter = %d \n", counter);

  int size_j = 0;
  for (int c = 0; c < counter; c++) {
    if (c != 0) {
      piece = piece * 2;
    }
    size_j = static_cast<int>(threads / pow(2, c) / 2);
    printf(" Piece = %d \n", piece);
    printf(" size_j = %d \n", size_j);
    int* r1 = new int[size_j];
    int* l1 = new int[size_j];
    int* r2 = new int[size_j];
    int* l2 = new int[size_j];
    for (int j = 0; j < size_j; j++) {
      l1[j] = j * piece * 2;
      r1[j] = j * piece * 2 + piece - 1;
      l2[j] = j * piece * 2 + piece;
      r2[j] = BinSearch(array, l2[j], j * piece * 2 + piece * 2 - 1, array[r1[j]]);
    }
    r2[size_j - 1] = BinSearch(array, l2[size_j - 1], size - 1, array[r1[size_j - 1]]);
    for (int j = 0; j < size_j; j++) {
      printf(" l1 = %d r1 = %d \n", l1[j], r1[j]);
      printf(" l2 = %d r2 = %d \n", l2[j], r2[j]);
    }
    int tid = 0;
#pragma omp parallel for schedule(dynamic, 1) private(tid)
    for (int i = 0; i < size_j; i++) {
      tid = omp_get_thread_num();
      printf("Thread: %d \n", tid);
      unsigned int * tmp = new unsigned int[r1[i] - l1[i] + 1 + r2[i] - l2[i] + 1];
      tmp = Splitter(array + l1[i], array + l2[i], r1[i] - l1[i] + 1, r2[i] - l2[i] + 1);
      int j = l1[i], g = 0;
      while (j <= r2[i]) {
        array[j] = tmp[g++];
        j++;
      }
      delete[] tmp;
    }
    print_array(array, size);
    printf("\n");
  }
}
int main() {
  int size = 20, threads = 4;
  printf(" Size = %d Threads = %d\n", size, threads);
  double time1, time2;
  unsigned int *a = new unsigned int[size];
  a = get_randomized_array(size);
  unsigned int *b = new unsigned int[size];
  copy_array(a, b, size);
  print_array(a, size);
  printf("\n\n");

  time1 = omp_get_wtime();
  dac_sort(a, size, threads);
  time2 = omp_get_wtime();
  printf("\n\n");
  print_array(a, size);
  printf(" Parralel time %f\n", time2 - time1);
  check_result(a, size);
  printf("\n\n");

  time1 = omp_get_wtime();
  radix_sort(b, size);
  time2 = omp_get_wtime();
  print_array(b, size);
  printf(" Linear time %f\n", time2 - time1);
  check_result(b, size);
  delete[] a;
  delete[] b;
  return 0;
}
