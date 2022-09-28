/*
* PQueue.h
*
*  Created on: November 18, 2021
*      Author: https://github.com/naeem014
*/

#ifndef P_QUEUE_H_
#define P_QUEUE_H_

#include <iostream>
#include <math.h>
#include <vector>
#include <unordered_map>

template <class T>
struct qElement {
    int key;
    double priority;
    T data;
};

template <class T>
class PQueue {
    public:
        PQueue();
        ~PQueue();
        
        bool empty();
        int size();
        void setMinQueueOn();
        void setMaxQueueOn();
        void setSpecialComparisonOn();
        void setSpecialComparisonOff();
        
        void clear();
        void insert(double priority, int key, T data);
        T pop();
        T getByKey(int key);
        void update(double priority, int key);

        void print();

    private:
        int parent(int index);
        int left(int index);
        int right(int index);
        bool comp(double a, double b);
        bool compSp(double a, double b, int id1, int id2);
        void swap(int index_a, int index_b);
        void moveDown(int index);
        void moveUp(int index);
        
        std::vector<qElement<T>> q;
        std::unordered_map<int, int> key_ref;
        bool min_queue = true;
        bool specialCmp = false;
};

#endif