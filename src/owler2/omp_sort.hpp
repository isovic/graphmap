#ifndef OMP_SORT_HPP
#define OMP_SORT_HPP

#include <thread>
#include <algorithm>

static size_t K = 19;
//static size_t num_threads = std::max((size_t) std::thread::hardware_concurrency(), (size_t) 1U);

template<typename T>
void swap(T* a, T* b) {
    T temp = *a;
    *a = *b;
    *b = temp;
}

template<typename T>
void insertionSort(T* a, size_t length) {

    size_t i = 1;
    for (; i < length; ++i) {
        T x = a[i];
        size_t j = i;

        while (j > 0 && x < a[j - 1]) {
            a[j] = a[j - 1];
            --j;
        }
        a[j] = x;
    }
}

template<typename T>
void quickSort(T* a, size_t length) {

    if (length < 2) {
        return;
    }

    if (length < K) {
        insertionSort(a, length);
        return;
    }

    T pivot = a[length / 2];
    size_t left = 0;
    size_t right = length - 1;

    for (;; ++left, --right) {
        while (a[left] < pivot) {
            ++left;
        }
        while (pivot < a[right]) {
            --right;
        }
        if (left >= right) {
            break;
        }
        swap(&a[left], &a[right]);
    }

    #pragma omp task
    quickSort(a, left);

    #pragma omp task
    quickSort(a + left, length - left);
}

template<typename T>
void pquickSort(T* a, size_t length, int32_t num_threads) {
    #pragma omp parallel num_threads(num_threads)
    #pragma omp single nowait
    quickSort(a, length);
}

#endif
