
#ifndef MEMENTO_H_
#define MEMENTO_H_

#include "Mesh.h"

class MeshMemento {
    public:
        MeshMemento(Mesh* mesh_) {
            std::copy(mesh_->V.begin(), mesh_->V.end(), std::back_inserter(vertices_));
            std::copy(mesh_->E.begin(), mesh_->E.end(), std::back_inserter(edges_));
            std::copy(mesh_->F.begin(), mesh_->F.end(), std::back_inserter(faces_));
            std::copy(mesh_->C.begin(), mesh_->C.end(), std::back_inserter(cells_));
        }
        
        int getVerticesSize() const { return vertices_.size(); }
        int getEdgesSize() const { return edges_.size(); }
        int getFacesSize() const { return faces_.size(); }
        int getCellsSize() const { return cells_.size(); }
        
        std::vector<Vertex> vertices_;
        std::vector<Edge> edges_;
        std::vector<Face> faces_;
        std::vector<Cell> cells_;
};

class MeshCaretaker {
public:
    void saveState(Mesh* mesh) {
        mementos_.push_back(std::make_unique<MeshMemento>(mesh));
    }
    void restoreState(Mesh* mesh) {
        if (!mementos_.empty()) {
            auto memento = std::move(mementos_.back());
            mementos_.pop_back();
            mesh->V.clear(); std::copy(memento->vertices_.begin(), memento->vertices_.end(), std::back_inserter(mesh->V));
            mesh->E.clear(); std::copy(memento->edges_.begin(), memento->edges_.end(), std::back_inserter(mesh->E));
            mesh->F.clear(); std::copy(memento->faces_.begin(), memento->faces_.end(), std::back_inserter(mesh->F));
            mesh->C.clear(); std::copy(memento->cells_.begin(), memento->cells_.end(), std::back_inserter(mesh->C));
        }
    }

    int vSize() {
        return mementos_.back()->getVerticesSize();
    }
private:
    std::vector<std::unique_ptr<MeshMemento>> mementos_;
};

#endif