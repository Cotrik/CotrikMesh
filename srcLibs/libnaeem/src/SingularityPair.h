#ifndef SINGULARITY_FIVE_PAIR_H
#define SINGULARITY_FIVE_PAIR_H

#include <glm/glm.hpp>

#include "SimplificationOperation.h"
#include "DiagonalCollapse.h"
#include "EdgeRotation.h"
#include "QuadSplit.h"
#include "VertexSplit.h"

#include "Mesh.h"
#include "MeshUtil.h"
#include "Smooth.h"

class SingularityPair {
    public:
        SingularityPair() {}
        SingularityPair(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t threeId_, size_t fiveId_) {
            mesh = &mesh_;
            mu = &mu_;
            smoother = &smoother_;
            threeId = threeId_;
            fiveId = fiveId_;
        }
        ~SingularityPair() {}

        // void MoveUpperLeft();
        // void MoveUpperRight();
        // void MoveLeft();
        // void MoveRight();
        // void MoveLowerLeft();
        // void MoveLowerRight();
        // void SplitSixSingularity();
        virtual int Move(size_t dest, double& delta, bool skipCheck = true) = 0;
        // void SetResolvedSingularity(size_t dest, int ref);
        // int GetResolvedSingularity();
        // std::vector<size_t> GetPairIds();
        
    protected:
        void CheckValidity() {
            if (mesh == NULL) {
                std::cout << "No mesh to use for Three Five Pair." << std::endl;
                exit(0);
            }
            if (mu == NULL) {
                std::cout << "MeshUtil is not initialized for Three Five Pair." << std::endl;
                exit(0);
            }
            if (smoother == NULL) {
                std::cout << "Smoother is not intitalized for Three Five Pair." << std::endl;
                exit(0);
            }
            if (mesh->V.size() == 0 || mesh->F.size() == 0 || mesh->C.size() == 0) {
                std::cout << "No mesh to use for Three Five Pair." << std::endl;
                exit(0);
            }
        }
        bool IsOperationValid() {
            if (mesh->V.at(threeId).N_Vids.size() != 3 || mesh->V.at(fiveId).N_Vids.size() != 5) return false;
            return true;
        }
        // void SetEdgeId();
        void AddContents(std::vector<size_t>& a, std::vector<size_t> b) {mu->AddContents(a, b);}
        void UpdateContents(std::vector<size_t>& a, std::vector<size_t> b) {mu->UpdateContents(a, b);}
        std::vector<size_t> GetDifference(std::vector<size_t> a, std::vector<size_t> b) {return mu->GetDifference(a, b);}
        std::vector<size_t> GetUnion(std::vector<size_t> a, std::vector<size_t> b) {return mu->GetUnion(a, b);}
        std::vector<size_t> GetIntersection(std::vector<size_t> a, std::vector<size_t> b) {return mu->GetIntersection(a, b);}
    
        Mesh* mesh;
        MeshUtil* mu;
        Smoother* smoother;
        size_t threeId, fiveId;
};

#endif