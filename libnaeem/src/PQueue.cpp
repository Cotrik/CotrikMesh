#include "PQueue.h"

template <class T>
PQueue<T>::PQueue() {}

template <class T>
PQueue<T>::~PQueue() {}

template <class T>
bool PQueue<T>::empty() {
    return q.empty();
}

template <class T>
int PQueue<T>::size() {
    return q.size();
}

template <class T>
int PQueue<T>::parent(int index) {
    return floor(0.5 * (index - 1));
}

template <class T>
int PQueue<T>::left(int index) {
    return (2 * index) + 1;
}

template <class T>
int PQueue<T>::right(int index) {
    return (2 * index) + 2;
}

template <class T>
bool PQueue<T>::comp(double a, double b) {
    if (min_queue) return a < b;
    return a > b;
}

template <class T>
bool PQueue<T>::compSp(double a, double b, int id1, int id2) {
    if (min_queue) return a < b || (a == b && id1 < id2);
    return a > b || (a == b && id1 < id2);
}

template <class T>
void PQueue<T>::swap(int index_a, int index_b) {
    key_ref.at(q.at(index_a).key) = index_b;
    key_ref.at(q.at(index_b).key) = index_a;
    qElement<T> temp = q.at(index_a);
    q.at(index_a) = q.at(index_b);
    q.at(index_b) = temp;
}

template <class T>
void PQueue<T>::moveDown(int index) {
    int left_index = left(index);
    int right_index = right(index);
    int final_index;
    if (specialCmp) {
        final_index = left_index <= size() - 1 && compSp(q.at(left_index).priority, q.at(index).priority, q.at(left_index).data->GetCenterId(), q.at(index).data->GetCenterId()) ? left_index : index;
        final_index = right_index <= size() - 1 && compSp(q.at(right_index).priority, q.at(left_index).priority, q.at(right_index).data->GetCenterId(), q.at(left_index).data->GetCenterId()) && compSp(q.at(right_index).priority, q.at(index).priority, q.at(right_index).data->GetCenterId(), q.at(index).data->GetCenterId()) ? right_index: final_index;
    } else {
        final_index = left_index <= size() - 1 && comp(q.at(left_index).priority, q.at(index).priority) ? left_index : index;
        final_index = right_index <= size() - 1 && comp(q.at(right_index).priority, q.at(left_index).priority) && comp(q.at(right_index).priority, q.at(index).priority) ? right_index: final_index;
    }
    if (final_index == index) return;
    swap(final_index, index);
    moveDown(final_index);
}

template <class T>
void PQueue<T>::moveUp(int index) {
    while (true) {
        int parent_index = parent(index);
        if (specialCmp) {
            if (parent_index < 0 || !compSp(q.at(index).priority, q.at(parent_index).priority, q.at(index).data->GetCenterId(), q.at(parent_index).data->GetCenterId())) break;
        } else {
            if (parent_index < 0 || !comp(q.at(index).priority, q.at(parent_index).priority)) break;
        }
        swap(index, parent_index);
        index = parent_index;
    }
}

template <class T>
void PQueue<T>::setMinQueueOn() {
    min_queue = true;
}

template <class T>
void PQueue<T>::setMaxQueueOn() {
    min_queue = false;
}

template <class T>
void PQueue<T>::setSpecialComparisonOn() {
    specialCmp = true;
}

template <class T>
void PQueue<T>::setSpecialComparisonOff() {
    specialCmp = false;
}

template <class T>
void PQueue<T>::clear() {
    q.clear();
    key_ref.clear();
}

template <class T>
void PQueue<T>::insert(double priority, int key, T data) {
    qElement<T> el;
    el.priority = priority;
    el.key = key;
    el.data = data;
    int index = size();
    key_ref.insert({key, index});
    q.push_back(el);
    moveUp(index);
}

template <class T>
T PQueue<T>::pop() {
    qElement<T> res = q.front();
    q.at(0) = q.back();
    key_ref.at(q.front().key) = 0;
    q.pop_back();
    // key_ref.erase(res.key);
    if (size() > 1) moveDown(0);
    return res.data;
}

template <class T>
T PQueue<T>::getByKey(int key) {
    if (key_ref.find(key) == key_ref.end() || empty()) return NULL;
    return q.at(key_ref.at(key)).data;
}


template <class T>
void PQueue<T>::update(double priority, int key) {
    int index = key_ref.at(key);
    double old_priority = q.at(index).priority;
    q.at(index).priority = priority;
    if (comp(priority, old_priority)) {
        moveUp(index);
    } else {
        moveDown(index);
    }
}

template <class T>
void PQueue<T>::print() {
    for (auto el: q) {
        // std::cout << el.key << " ";
        std::cout << el.priority << " ";
    }
    std::cout << "*********************" << std::endl;
}