// Copyright 2019 Batova Darya
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <random>
#include <ctime>
#include <algorithm>

void FillingArray(int* array, int n, int min, int max) {
    static std::random_device random_device;
    static std::mt19937 generator(random_device());
    static std::uniform_int_distribution<> distribution(min, max);
    int i = 0;
    while (i < n) {
        array[i++] = distribution(generator);
    }
}

void quickSort(int arr[], int left, int right) {
    int i = left, j = right;
    int tmp;
    int pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }

    /* recursion */
    if (left < j)
        quickSort(arr, left, j);
    if (i < right)
        quickSort(arr, i, right);
}

void ShowArray(int* arr, int n) {
    if (n < 19) {
        int i = 0;
        while (i < n) {
            std::cout << arr[i++] << " ";
        }
        std::cout << std::endl;
    }
}

void Merge_Two_Arrays_Into_Tmp(int* arr1, int* arr2,
        int* tmp, int size1, int size2) {
    int index1 = 0;
    int index2 = 0;
    int i = 0;

    for (; (index1 != size1) && (index2 != size2); i++) {
        if (arr1[index1] <= arr2[index2]) {
            tmp[i] = arr1[index1];
            index1++;
        } else {
            tmp[i] = arr2[index2];
            index2++;
        }
    }

    if (index1 == size1) {
        int j = index2;
        while (j < size2) {
            tmp[i++] = arr2[j++];
        }
    } else {
        int j = index1;
        while (j < size1) {
            tmp[i++] = arr1[j++];
        }
    }
}

void Merge_Two_Arrays_Into_Res(int* arr1, int* arr2,
                     int size1, int size2, int* res) {
    int x = arr1[size1/2];
    auto it = std::lower_bound(arr2, arr2 + size2, x);
    int x1 = std::distance(arr2, it);
    int nSize1 = size1 / 2 + x1;
    Merge_Two_Arrays_Into_Tmp(arr1, arr2, res, size1 / 2, x1);
    Merge_Two_Arrays_Into_Tmp(arr1 + size1 / 2, arr2 + x1,
                 res + nSize1, size1 - size1 / 2, size2 - x1);
}

void QuickSortMerge(int* array, int* res, int size) {
    int size1 = size / 2;
    int size2 = size - size1;
    quickSort(array, 0, size1-1);
    quickSort(array + size1, 0, size2 - 1);
    Merge_Two_Arrays_Into_Res(array, array + size1, size1, size2, res);
}

bool AreArraysEqual(int* arr1, int* arr2, int size) {
    for (int i = 0; i < size; ++i) {
        if (arr1[i] != arr2[i]) {
            return false;
        }
    }
    return true;
}

void CopyArray(int* src, int* dst, int size) {
    for (int i = 0; i < size; ++i) {
        dst[i] = src[i];
    }
}

int main() {
    int* array;
    int* copy;
    int size = 15;
    array = new int[size];
    copy = new int[size];
    FillingArray(array, size, 0, 30);
    std::cout << "Initially the array looks like this: " << std::endl;
    ShowArray(array, size);
    CopyArray(array, copy, size);
    quickSort(copy, 0, size-1);
    std::cout << "Sorted array using quick sort: " << std::endl;
    ShowArray(copy, size);
    int* res = new int[size];
    QuickSortMerge(array, res, size);
    std::cout << "Sorted array using divide and conquer: " << std::endl;
    ShowArray(res, size);
    if (AreArraysEqual(copy, res, size)) {
        std::cout << "Excellent! Arrays are matching" << std::endl;
    } else {
        std::cout << "Error! Arrays are not matching" << std::endl;
    }
    delete[] array;
    delete[] copy;
    delete[] res;
    return 0;
}
