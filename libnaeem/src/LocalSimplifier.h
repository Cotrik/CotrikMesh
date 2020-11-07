/*
* Simplifier.h
*
*  Created on: October 7, 2020
*      Author: https://github.com/naeem014
*/

#ifndef LOCAL_SIMPLIFIER_H
#define LOCAL_SIMPLIFIER_H

// #include <bits/stdc++.h>

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "Simplifier.h"

struct LocalOperation {
    std::string type;
    double profitability = 0;
    std::set<size_t> canceledFids;
    std::vector<Face> newFaces;
};

class LocalSimplifier : public Simplifier {
    public:
        LocalSimplifier(Mesh& mesh);
        virtual ~LocalSimplifier();
    private:
        LocalSimplifier();
        LocalSimplifier(const LocalSimplifier&);
        LocalSimplifier& operator = (const LocalSimplifier&);
    public:
        void Simplify();
        void getOperationsPriorities(std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& OptimizationOps, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& CollapseOps);
        void sortPriorities(std::multiset<LocalOperation>& OptimizationOps, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& CollapseOps);
        void VertexRotate(Vertex& v, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops);
        void EdgeRotate(Edge& e, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops);
        void EdgeCollapse(Edge& e, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops);
        void DiagonalCollapse(Face& f, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops);
        void RemoveDoublets(std::set<size_t>& canceledFids);
        double getLength(size_t vid1, size_t vid2);
};

#endif // !LOCAL_SIMPLIFIER_H