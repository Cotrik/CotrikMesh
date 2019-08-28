/*
 * OpenQuadMeshOpt.cpp
 *
 *  Created on: June 20, 2019
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "SurfaceMeshOpt.h"
#include "ArgumentManager.h"
#include "Util.h"
#include "MeshQuality.h"
#include "GetSet.h"
#include <math.h>
#include <random>
#include <time.h>

class MeshQualityImprover {
public:
    MeshQualityImprover(Mesh& mesh, const Mesh* refMesh = nullptr): mesh(mesh), refMesh(refMesh){

    }
    virtual ~MeshQualityImprover() {

    }
private:
    MeshQualityImprover();
    MeshQualityImprover(const MeshQualityImprover&);
    MeshQualityImprover& operator = (const MeshQualityImprover&);
public:
    virtual int Run(int iters = 20) = 0;
    virtual size_t UpdateQuality(const double invalidThreshold = 0.0) = 0; // return #invalidElements
    virtual double GetQuality() const {
        return 1.0;
    }
protected:
    Mesh& mesh;
    const Mesh* refMesh = NULL;
};

class MeshMSJImprover : public MeshQualityImprover {
public:
    MeshMSJImprover(Mesh& mesh, const Mesh* refMesh = nullptr): MeshQualityImprover(mesh, refMesh) {

    }
    virtual ~MeshMSJImprover() {

    }
private:
    MeshMSJImprover();
    MeshMSJImprover(const MeshMSJImprover&);
    MeshMSJImprover& operator = (const MeshMSJImprover&);
public:
    virtual void AddOrthogonalTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& B, size_t& row) = 0;   // * alpha
    virtual void AddStraightTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& B, size_t& row) = 0;     // * alpha
    virtual void AddSingularTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& B, size_t& row) = 0;     // * alpha
    virtual size_t AddBoundaryTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& B, size_t& row) = 0;     // * alpha

    GETSET(double, quality)
    GETSET(double, alpha)
    GETSET(double, beta)
    GETSET(double, gamma)
    GETSET(double, anisotropy)
    GETSET(double, targetMSJ)

    GETSET(double, stepSize)
    GETSET(double, avgMeshEdgeLength)

    GETSET(bool, useAverageTargetLength)
    GETSET(bool, recoverable)
    GETSET(bool, allowBigStep)
    GETSET(bool, changeBoundary)
    GETSET(bool, projectToTargetSurface)
    GETSET(bool, useProjection)

private:
    double quality = INT_MIN;
    double alpha = 1000000.0;
    double beta = 1000000.0;
    double gamma = 1.0;
    double anisotropy = 0.05;
    double targetMSJ = 0.0;

    double stepSize = 0.9;
    double avgMeshEdgeLength = 1.0;
    bool useAverageTargetLength = false;
    bool recoverable = true;
    bool allowBigStep = false;
    bool changeBoundary = false;
    bool projectToTargetSurface = false;
    bool useProjection = true;
};

class QuadMeshMSJImprover : public MeshMSJImprover {
public:
    QuadMeshMSJImprover(Mesh& mesh, const Mesh* refMesh = nullptr): MeshMSJImprover(mesh, refMesh) {

    }
    virtual ~QuadMeshMSJImprover() {

    }
private:
    QuadMeshMSJImprover();
    QuadMeshMSJImprover(const QuadMeshMSJImprover&);
    QuadMeshMSJImprover& operator = (const QuadMeshMSJImprover&);
public:
    double GetQuality() const {
        return Getquality();
    }
    size_t UpdateQuality(const double invalidThreshold = 0.0) {
        double minMSJ = 1.0;
        auto numOfInvalidElements = GetMinScaledJacobianQuad(mesh, minMSJ, invalidThreshold);
        Setquality(minMSJ);
        return numOfInvalidElements;
    }
    void UpdateAvgEdgeLength() {
        Eigen::initParallel();
        if (GetstepSize() < 0.1) {
            std::cout << "ERROR! stepSize < 0.1\n";
            return;
        }
        double sumEdgeLength = 0.0;
        // size_t numOfBoundaryEdges = 0;
        for (auto& e : mesh.E) {
            // if (e.isBoundary)
            {
                const Vertex& v1 = mesh.V.at(e.Vids.at(0));
                const Vertex& v2 = mesh.V.at(e.Vids.at(1));
                e.length = glm::length(v1 - v2);
                sumEdgeLength += e.length;
                // numOfBoundaryEdges++;
            }
        }
        SetavgMeshEdgeLength(sumEdgeLength / mesh.E.size());
    }
    virtual int Run(int iters = 20) { // return number of iterations that is needed to achieve target mesh
        UpdateAvgEdgeLength();
        std::cout << "Average Edge Length = " << GetavgMeshEdgeLength() << std::endl;
        auto numOfInvertdElements = UpdateQuality(this->GettargetMSJ());
        std::cout << "iter = " << 0 << " #inverted = " << numOfInvertdElements << " MSJ = " << this->GetQuality() << std::endl;
        //----------------------------------------------
        int iter = 0;
        double prevMinimumScaledJacobian = -1.0;
        bool converged = false;
        double initStepSize = GetstepSize();
        bool initUseAverageTargetLength = GetuseAverageTargetLength();
        bool untangled = false;
        ComputeMeshTargetLength();
        static int global_count = 0;
        global_count++;
        auto bestMesh = mesh;
        while (!converged && iter++ < iters) {
            if (!initUseAverageTargetLength  && iter == 1) SetuseAverageTargetLength(true);
            mesh.ClearLabelOfSurface();
            mesh.LabelSurface();
            std::string filename = std::string("Faces.") + std::to_string(iter) + ".vtk";
            converged = UpdateVertices();
            if (!initUseAverageTargetLength && iter == 1) SetuseAverageTargetLength(false);
            SetstepSize(GetstepSize() * initStepSize);

            auto numOfInvertdElements = UpdateQuality(this->GettargetMSJ());
            std::cout << "iter = " << iter << " #inverted = " << numOfInvertdElements << " MSJ = " << this->GetQuality() << std::endl;
            if (prevMinimumScaledJacobian > this->GettargetMSJ() && this->GetQuality() <= prevMinimumScaledJacobian) {
                filename = std::string("MeshOpt.") + std::to_string(iter - 1) + ".vtk";
                std::cout << "*************************" << std::endl;
                std::cout << "\033[1;32mBest Mesh is " << filename << "\033[0m" << std::endl;
                std::cout << "*************************" << std::endl;
                MeshFileWriter optwriter(bestMesh, "opt.vtk");
                optwriter.WriteFile();
                untangled = true;
                break;
            } else if (converged) {
                std::cout << "*************************" << std::endl;
                std::cout << "Converged at " << filename << std::endl;
                std::cout << "*************************" << std::endl;
            }
            prevMinimumScaledJacobian = this->GetQuality();
            bestMesh = mesh;
        }
        MeshFileWriter optwriter(bestMesh, "output.vtk");
        optwriter.WriteFile();
        return iter - 1;
    }
    void ComputeMeshTargetLength() {
        std::cout << "=============================\n";
        std::cout << "   Computing TargetLength    \n";

        std::vector<Eigen::Triplet<double>> coefficients;
        std::vector<double> b;
        coefficients.reserve(mesh.E.size() * 5);
        b.reserve(mesh.E.size() * 5);
        size_t row  = 0;
        for (size_t i = 0; i < mesh.E.size(); i++) {
            auto& e = mesh.E.at(i);
            coefficients.push_back(Eigen::Triplet<double>(row, i, 1));
            if (e.isBoundary) b.push_back(e.length);
            else b.push_back(GetavgMeshEdgeLength());
            ++row;

            for (auto parallelEid : e.parallelEids) {
                coefficients.push_back(Eigen::Triplet<double>(row, i, 1));
                coefficients.push_back(Eigen::Triplet<double>(row++, parallelEid, -1));
                b.push_back(0);
            }
            for (auto consecutiveEid : e.consecutiveEids) {
                coefficients.push_back(Eigen::Triplet<double>(row, i, 1));
                coefficients.push_back(Eigen::Triplet<double>(row++, consecutiveEid, -1));
                b.push_back(0);
            }
//            for (auto orthogonalEid : e.orthogonalEids) {
//                coefficients.push_back(Trip(row, i, 1));
//                coefficients.push_back(Trip(row, orthogonalEid, -1));
//                b.push_back(0);
//            }
//            row += e.orthogonalEids.size();
        }

//        for (const Edge& e : mesh.E) {
//            if (!e.isSingularity) continue;
//            if (e.N_Cids.size() != 5 && e.N_Cids.size() != 3 && e.N_Cids.size() != 6) continue;
//            auto combs = Util::combine(e.orthogonalEids.size(), 2);
//            for (auto& comb : combs) {
//                auto eid1 = e.orthogonalEids[comb[0]];
//                auto eid2 = e.orthogonalEids[comb[1]];
//                coefficients.push_back(Eigen::Triplet<double>(row, eid1, 1));
//                coefficients.push_back(Eigen::Triplet<double>(row++, eid2, -1));
//                b.push_back(0);
//            }
//        }

        size_t col = mesh.E.size();
        Eigen::SparseMatrix<double> A(row, col);
        A.setFromTriplets(coefficients.begin(), coefficients.end());
        Eigen::SparseMatrix<double> A_T = A.transpose();
        Eigen::SparseMatrix<double> ATA = A_T * A;
        VectorXd X(col), B(row);
        for (size_t i = 0; i < row; i++)
            B(i) = b[i];

        VectorXd ATB = A_T * B;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(ATA);
        X = chol.solve(ATB);

        for (auto& e : mesh.E) {
            e.length = X[e.id];
            //std::cout << X[i] << std::endl;
        }
        std::cout << "=============================\n";
    }
    VectorXd Solve() {
        std::vector<Eigen::Triplet<double>> A_Entries;
        std::vector<double> b;
        size_t row = 0;
        size_t currentRow = row;
        size_t currentNumOfA_Entries = A_Entries.size();
        //AddSingularTerm(A_Entries, b, row);
        //std::cout << "AddSingularTerm--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
        currentRow = row;
        currentNumOfA_Entries = A_Entries.size();
        AddOrthogonalTerm(A_Entries, b, row);
        //std::cout << "AddOrthogonalTerm\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
        currentRow = row;
        currentNumOfA_Entries = A_Entries.size();
        AddStraightTerm(A_Entries, b, row);
        //std::cout << "AddStraightTerm-\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
        currentRow = row;
        currentNumOfA_Entries = A_Entries.size();
        size_t a_id = AddBoundaryTerm(A_Entries, b, row);
        //std::cout << "AddBoundaryTerm--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

        const size_t col = 3 * mesh.V.size() + a_id;

    //    std::cout << "----------- A Info ------------\n";
    //    std::cout << "Entries = " << A_Entries.size() << " row = " << row << " col = " << col << std::endl;

        Eigen::SparseMatrix<double> A(row, col);
        A.setFromTriplets(A_Entries.begin(), A_Entries.end());
        Eigen::SparseMatrix<double> AT = A.transpose();
        Eigen::SparseMatrix<double> ATA = AT * A;

        VectorXd B(row);
        for (size_t i = 0; i < row; i++)
            B(i) = b[i];

        VectorXd ATB = AT * B;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(ATA);
        VectorXd X(col);
        X = chol.solve(ATB);

        return X;
    }

protected:
    bool UpdateVertices() {
        auto X = Solve();
        size_t lastNumOfInvertedElements = UpdateQuality(GettargetMSJ());
        bool converged = true;
        if (Getrecoverable()) {
            Mesh targetMesh(mesh);
            std::vector<glm::dvec3> oldV(mesh.V.size());
            for (auto& v : mesh.V)
                oldV[v.id] = v.xyz();

//            for (size_t i = 0; i < mesh.V.size(); i++) {
//                if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 1e-6
//                    || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 1e-6
//                    //|| fabs(double(mesh.V[i].z - X[3 * i + 2])) > 1e-6
//                    )
//                {
//                    converged = false;
//                }
//                if (!GetallowBigStep())
//                    if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * GetavgMeshEdgeLength()
//                        || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * GetavgMeshEdgeLength()
//                        //|| fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * GetavgMeshEdgeLength()
//                )
//                    {
//                        continue;
//                    }
//                mesh.V[i].x = GetstepSize() * X[3 * i + 0] + (1.0 - GetstepSize()) * mesh.V[i].x;
//                mesh.V[i].y = GetstepSize() * X[3 * i + 1] + (1.0 - GetstepSize()) * mesh.V[i].y;
//                //mesh.V[i].z = GetstepSize() * X[3 * i + 2] + (1.0 - GetstepSize()) * mesh.V[i].z;
//            }
            for (auto& v : mesh.V) {
                if (fabs(v.x - X[3 * v.id + 0]) > 1e-6 || fabs(v.y - X[3 * v.id + 1]) > 1e-6 || fabs(v.z - X[3 * v.id + 2]) > 1e-6)
                    converged = false;
                if (fabs(v.x - X[3 * v.id + 0]) < 2.0 * GetavgMeshEdgeLength())
                    v.x = GetstepSize() * X[3 * v.id + 0] + (1.0 - GetstepSize()) * v.x;
                if (fabs(v.y - X[3 * v.id + 1]) < 2.0 * GetavgMeshEdgeLength())
                    v.y = GetstepSize() * X[3 * v.id + 1] + (1.0 - GetstepSize()) * v.y;
                if (fabs(v.z - X[3 * v.id + 2]) < 2.0 * GetavgMeshEdgeLength())
                    v.z = GetstepSize() * X[3 * v.id + 2] + (1.0 - GetstepSize()) * v.z;
            }
            if (GetuseProjection()) {
                if (!GetprojectToTargetSurface()) mesh.ProjectTo(targetMesh);
                else mesh.FastProjectTo(targetMesh);
            }

            size_t numOfInvertedElements = UpdateQuality(GettargetMSJ());
            if (numOfInvertedElements > lastNumOfInvertedElements) {
                std::cout << "Recover previous mesh\n";
                for (size_t i = 0; i < mesh.V.size(); i++)
                    mesh.V[i] = oldV[i];
            }
        }
        else {
//            for (size_t i = 0; i < mesh.V.size(); i++) {
//                if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 1e-6
//                    || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 1e-6
//                    || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 1e-6
//                    )
//                {
//                    converged = false;
//                }
//                if (!GetallowBigStep())
//                if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * GetavgMeshEdgeLength()
//                    || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * GetavgMeshEdgeLength()
//                    || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * GetavgMeshEdgeLength()
//                )
//                {
//                    continue;
//                }
//                mesh.V[i].x = GetstepSize() * X[3 * i + 0] + (1.0 - GetstepSize()) * mesh.V[i].x;
//                mesh.V[i].y = GetstepSize() * X[3 * i + 1] + (1.0 - GetstepSize()) * mesh.V[i].y;
//                mesh.V[i].z = GetstepSize() * X[3 * i + 2] + (1.0 - GetstepSize()) * mesh.V[i].z;
//            }
            for (auto& v : mesh.V) {
                if (fabs(v.x - X[3 * v.id + 0]) > 1e-6 || fabs(v.y - X[3 * v.id + 1]) > 1e-6 || fabs(v.z - X[3 * v.id + 2]) > 1e-6)
                    converged = false;
                if (fabs(v.x - X[3 * v.id + 0]) < 2.0 * GetavgMeshEdgeLength())
                    v.x = GetstepSize() * X[3 * v.id + 0] + (1.0 - GetstepSize()) * v.x;
                if (fabs(v.y - X[3 * v.id + 1]) < 2.0 * GetavgMeshEdgeLength())
                    v.y = GetstepSize() * X[3 * v.id + 1] + (1.0 - GetstepSize()) * v.y;
                if (fabs(v.z - X[3 * v.id + 2]) < 2.0 * GetavgMeshEdgeLength())
                    v.z = GetstepSize() * X[3 * v.id + 2] + (1.0 - GetstepSize()) * v.z;
            }
        }
        return converged;
    }
    virtual void AddOrthogonalTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& b, size_t& row) {
        for (auto& e: mesh.E)
            e.energyOrthogonality = 0;
        double energy = 0.0;
        for (auto& edge1 : mesh.E) {
            for (auto edgeId2 : edge1.orthogonalEids) {
                auto& edge2 = mesh.E.at(edgeId2);
                const size_t vId1_1 = edge1.Vids.at(0);
                const size_t vId1_2 = edge1.Vids.at(1);
                const size_t vId2_1 = edge2.Vids.at(0);
                const size_t vId2_2 = edge2.Vids.at(1);
                size_t shareVId = vId1_1;
                if (vId1_2 == vId2_1 || vId1_2 == vId2_2)
                    shareVId = vId1_2;
                const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
                const Vertex& vetexC = mesh.V.at(shareVId);
                const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
                const glm::dvec3 v2 = glm::dvec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
                const glm::dvec3 v2n = glm::normalize(v2);
                double length_v1 = edge1.length;
                if (!GetuseAverageTargetLength()) {
                    length_v1 = glm::length(glm::dvec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                    // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                    length_v1 = length_v1 > GetavgMeshEdgeLength() * Getanisotropy() ? length_v1 : GetavgMeshEdgeLength() * Getanisotropy();
                    length_v1 = std::isnan(length_v1) ? GetavgMeshEdgeLength() * Getanisotropy() : length_v1;
                }
                ///////////////////////////////////////////////
                const double length_v1_ = 1.0/length_v1;
                const glm::dvec3 v1n = glm::dvec3((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_, (vetexC.z - vetex1.z)*length_v1_);
                const double energy_orthognality = (glm::dot(v1n, v2n) - 0.0) * (glm::dot(v1n, v2n) - 0.0) * Getgamma();
                energy += energy_orthognality;
                edge1.energyOrthogonality += energy_orthognality;

                double weight = Getgamma();
                for (size_t n = 0; n < 3; n++) {
                    double a = weight * v2n[n] / length_v1;
                    a = std::isnan(a) ? 1.0 : a;
                    A_Entries.push_back(Eigen::Triplet<double>(row, 3 * vetexC.id + n, a));
                    A_Entries.push_back(Eigen::Triplet<double>(row, 3 * vetex1.id + n, -a));
                }
                b.push_back(0.0 * weight);
                row++;
            }
        }
        EOrthogonality.push_back(energy);
    }
    virtual void AddStraightTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& b, size_t& row) {
        for (auto& e: mesh.E)
            e.energyStraightness = 0;
        double initgamma = Getgamma();
        Setgamma(1.0);
        double energy = 0.0;
        for (auto& edge1 : mesh.E) {
            for (auto edgeId2 : edge1.consecutiveEids) {
                /*const */Edge& edge2 = mesh.E.at(edgeId2);
                const size_t vId1_1 = edge1.Vids.at(0);
                const size_t vId1_2 = edge1.Vids.at(1);
                const size_t vId2_1 = edge2.Vids.at(0);
                const size_t vId2_2 = edge2.Vids.at(1);
                size_t shareVId = vId1_1;
                if (vId1_2 == vId2_1 || vId1_2 == vId2_2) shareVId = vId1_2;
                const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
                const Vertex& vetexC = mesh.V.at(shareVId);
                const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
                const glm::dvec3 v2 = glm::dvec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
                const glm::dvec3 v2n = glm::normalize(v2);
                double length_v1 = edge1.length;
                if (!GetuseAverageTargetLength()) {
                    length_v1 = glm::length(glm::dvec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                    // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                    length_v1 = length_v1 > GetavgMeshEdgeLength() * Getanisotropy() ? length_v1 : GetavgMeshEdgeLength() * Getanisotropy();
                    length_v1 = std::isnan(length_v1) ? GetavgMeshEdgeLength() * Getanisotropy() : length_v1;
                }
                ///////////////////////////////////////////////
                const double length_v1_ = 1.0 / length_v1;
                const glm::dvec3 v1n = glm::dvec3((vetexC.x - vetex1.x) * length_v1_, (vetexC.y - vetex1.y) * length_v1_,
                        (vetexC.z - vetex1.z) * length_v1_);
                const double energy_straightness = (glm::dot(v1n, v2n) + 1.0) * (glm::dot(v1n, v2n) + 1.0) * Getgamma();
                energy += energy_straightness;
                edge1.energyStraightness += energy_straightness;
                ///////////////////////////////////////////////
                // v1 is variable, and use v2 as constant
                // <v1/length_v1, v2n> = -1
                // double weight = energy_straightness * 10 > gamma ? energy_straightness  * 10: gamma;
                // if (stepSize > 0.9)
                double weight = Getgamma();
                for (size_t n = 0; n < 3; n++) {
                    double a = weight * v2n[n] / length_v1;
                    a = std::isnan(a) ? 1.0 : a;
                    A_Entries.push_back(Eigen::Triplet<double>(row, 3 * vetexC.id + n, a));
                    A_Entries.push_back(Eigen::Triplet<double>(row, 3 * vetex1.id + n, -a));
                }
                b.push_back(-1.0 * weight);
                row++;
            }
        }
        EStraightness.push_back(energy);
        Setgamma(initgamma);
    }
    virtual void AddSingularTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& B, size_t& row) {
        for (auto& e: mesh.E)
            e.energySingularity = 0;
        double energy = 0.0;
    }
    virtual size_t AddBoundaryTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& b, size_t& row) {
        auto alpha = Getalpha();
        auto beta = Getbeta();
        size_t a_id = 0;
        for (const auto& v : mesh.V) {
            if (!v.isBoundary) continue;
            for (size_t n = 0; n < 3; n++) {
//                int a_index = a_id++ + 3 * mesh.V.size();
//                A_Entries.push_back(Eigen::Triplet<double>(row, a_index, alpha));
//                b.push_back(0);
//                row++;
                A_Entries.push_back(Eigen::Triplet<double>(row, 3 * v.id + n, alpha));
                b.push_back(alpha * v[n]);
                row++;
            }
        }
//        return a_id;

        for (size_t i = 0; i < mesh.V.size(); i++) {
            const Vertex& v = mesh.V.at(i);
            if (!v.isBoundary)
                continue;
            if (v.type == REGULAR) {
                //beta(nv + d)
                const glm::dvec3& implicit_n = v.normal;
                float implicit_d = glm::dot(implicit_n, v.xyz());
                for (size_t k = 0; k < 3; k++) {
                    A_Entries.push_back(Eigen::Triplet<double>(row, 3 * i + k, beta * implicit_n[k]));
                }
                b.push_back(beta * implicit_d);
                row++;
            } else if (v.type == FEATURE) {
                int a_index = a_id + 3 * mesh.V.size();
                for (size_t k = 0; k < 3; k++) {
                    A_Entries.push_back(Eigen::Triplet<double>(row, 3 * i + k, alpha));
                    A_Entries.push_back(Eigen::Triplet<double>(row, a_index, -alpha * v.tangent[k]));
                    b.push_back(alpha * v[k]);
                    row++;
                }
            } else if (v.type == CORNER) {
                for (size_t k = 0; k < 3; k++) {
                    A_Entries.push_back(Eigen::Triplet<double>(row, 3 * i + k, alpha));
                    b.push_back(alpha * v[k]);
                    row++;
                }
            }
        }
        a_id = 0;
        for (size_t i = 0; i < mesh.V.size(); i++) {
            const Vertex& v = mesh.V.at(i);
            if (v.isBoundary) {
                int a_index = a_id++ + 3 * mesh.V.size();
                A_Entries.push_back(Eigen::Triplet<double>(row, a_index, alpha));
                b.push_back(0);
                row++;
            }
        }

        return a_id;
    }
private:
    std::vector<double> ESingularity;
    std::vector<double> EOrthogonality;
    std::vector<double> EStraightness;
};

class OpenQuadMeshMSJImprover : public MeshMSJImprover {
public:
	OpenQuadMeshMSJImprover(Mesh& mesh, const Mesh* refMesh = nullptr) : MeshMSJImprover(mesh, refMesh) {

	}
	virtual ~OpenQuadMeshMSJImprover() {

	}
private:
	OpenQuadMeshMSJImprover();
	OpenQuadMeshMSJImprover(const OpenQuadMeshMSJImprover&);
	OpenQuadMeshMSJImprover& operator = (const OpenQuadMeshMSJImprover&);
public:
	double GetQuality() const {
		return Getquality();
	}
	size_t UpdateQuality(const double invalidThreshold = 0.0) {
		double minMSJ = 1.0;
		auto numOfInvalidElements = GetMinScaledJacobianQuad(mesh, minMSJ, invalidThreshold);
		Setquality(minMSJ);
		return numOfInvalidElements;
	}
	void UpdateAvgEdgeLength() {
		Eigen::initParallel();
		if (GetstepSize() < 0.1) {
			std::cout << "ERROR! stepSize < 0.1\n";
			return;
		}
		SetavgMeshEdgeLength(mesh.GetAvgEdgeLength());
	}
	virtual int Run(int iters = 20) { // return number of iterations that is needed to achieve target mesh
		UpdateAvgEdgeLength();
		std::cout << "Average Edge Length = " << GetavgMeshEdgeLength() << std::endl;
		auto numOfInvertdElements = UpdateQuality(this->GettargetMSJ());
		std::cout << "iter = " << 0 << " #inverted = " << numOfInvertdElements << " MSJ = " << this->GetQuality() << std::endl;
		//----------------------------------------------
		int iter = 0;
		double prevMinimumScaledJacobian = -1.0;
		bool converged = false;
		double initStepSize = GetstepSize();
		bool initUseAverageTargetLength = GetuseAverageTargetLength();
		bool untangled = false;
		ComputeMeshTargetLength();
		static int global_count = 0;
		global_count++;
		auto bestMesh = mesh;
		while (!converged && iter++ < iters) {
			if (!initUseAverageTargetLength  && iter == 1) SetuseAverageTargetLength(true);
			mesh.ClearLabelOfSurface();
			mesh.LabelSurface();
			std::string filename = std::string("Faces.") + std::to_string(iter) + ".vtk";
			converged = UpdateVertices();
			if (!initUseAverageTargetLength && iter == 1) SetuseAverageTargetLength(false);
			SetstepSize(GetstepSize() * initStepSize);

			auto numOfInvertdElements = UpdateQuality(this->GettargetMSJ());
			std::cout << "iter = " << iter << " #inverted = " << numOfInvertdElements << " MSJ = " << this->GetQuality() << std::endl;
			if (prevMinimumScaledJacobian > this->GettargetMSJ() && this->GetQuality() <= prevMinimumScaledJacobian) {
				filename = std::string("MeshOpt.") + std::to_string(iter - 1) + ".vtk";
				std::cout << "*************************" << std::endl;
				std::cout << "\033[1;32mBest Mesh is " << filename << "\033[0m" << std::endl;
				std::cout << "*************************" << std::endl;
				MeshFileWriter optwriter(bestMesh, "opt.vtk");
				optwriter.WriteFile();
				untangled = true;
				break;
			} else if (converged) {
				std::cout << "*************************" << std::endl;
				std::cout << "Converged at " << filename << std::endl;
				std::cout << "*************************" << std::endl;
			}
			prevMinimumScaledJacobian = this->GetQuality();
			bestMesh = mesh;
		}
		MeshFileWriter optwriter(bestMesh, "output.vtk");
		optwriter.WriteFile();
		return iter - 1;
	}
	void ComputeMeshTargetLength() {
		std::cout << "=============================\n";
		std::cout << "   Computing TargetLength    \n";

		std::vector<Eigen::Triplet<double>> coefficients;
		std::vector<double> b;
		coefficients.reserve(mesh.E.size() * 5);
		b.reserve(mesh.E.size() * 5);
		size_t row = 0;
		for (size_t i = 0; i < mesh.E.size(); i++) {
			auto& e = mesh.E.at(i);
			coefficients.push_back(Eigen::Triplet<double>(row, i, 1));
			if (e.isBoundary) b.push_back(e.length);
			else b.push_back(GetavgMeshEdgeLength());
			++row;

			for (auto parallelEid : e.parallelEids) {
				coefficients.push_back(Eigen::Triplet<double>(row, i, 1));
				coefficients.push_back(Eigen::Triplet<double>(row++, parallelEid, -1));
				b.push_back(0);
			}
			for (auto consecutiveEid : e.consecutiveEids) {
				coefficients.push_back(Eigen::Triplet<double>(row, i, 1));
				coefficients.push_back(Eigen::Triplet<double>(row++, consecutiveEid, -1));
				b.push_back(0);
			}
			//            for (auto orthogonalEid : e.orthogonalEids) {
			//                coefficients.push_back(Trip(row, i, 1));
			//                coefficients.push_back(Trip(row, orthogonalEid, -1));
			//                b.push_back(0);
			//            }
			//            row += e.orthogonalEids.size();
		}

		//        for (const Edge& e : mesh.E) {
		//            if (!e.isSingularity) continue;
		//            if (e.N_Cids.size() != 5 && e.N_Cids.size() != 3 && e.N_Cids.size() != 6) continue;
		//            auto combs = Util::combine(e.orthogonalEids.size(), 2);
		//            for (auto& comb : combs) {
		//                auto eid1 = e.orthogonalEids[comb[0]];
		//                auto eid2 = e.orthogonalEids[comb[1]];
		//                coefficients.push_back(Eigen::Triplet<double>(row, eid1, 1));
		//                coefficients.push_back(Eigen::Triplet<double>(row++, eid2, -1));
		//                b.push_back(0);
		//            }
		//        }

		size_t col = mesh.E.size();
		Eigen::SparseMatrix<double> A(row, col);
		A.setFromTriplets(coefficients.begin(), coefficients.end());
		Eigen::SparseMatrix<double> A_T = A.transpose();
		Eigen::SparseMatrix<double> ATA = A_T * A;
		VectorXd X(col), B(row);
		for (size_t i = 0; i < row; i++)
			B(i) = b[i];

		VectorXd ATB = A_T * B;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(ATA);
		X = chol.solve(ATB);

		for (auto& e : mesh.E) {
			e.length = X[e.id];
			//std::cout << X[i] << std::endl;
		}
		std::cout << "=============================\n";
	}
	VectorXd Solve() {
		std::vector<Eigen::Triplet<double>> A_Entries;
		std::vector<double> b;
		size_t row = 0;
		size_t currentRow = row;
		size_t currentNumOfA_Entries = A_Entries.size();
		//AddSingularTerm(A_Entries, b, row);
		//std::cout << "AddSingularTerm--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
		currentRow = row;
		currentNumOfA_Entries = A_Entries.size();
		AddOrthogonalTerm(A_Entries, b, row);
		//std::cout << "AddOrthogonalTerm\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
		currentRow = row;
		currentNumOfA_Entries = A_Entries.size();
		AddStraightTerm(A_Entries, b, row);
		//std::cout << "AddStraightTerm-\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
		currentRow = row;
		currentNumOfA_Entries = A_Entries.size();
		size_t a_id = AddBoundaryTerm(A_Entries, b, row);
		//std::cout << "AddBoundaryTerm--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

		const size_t col = 3 * mesh.V.size() + a_id;

		//    std::cout << "----------- A Info ------------\n";
		//    std::cout << "Entries = " << A_Entries.size() << " row = " << row << " col = " << col << std::endl;

		Eigen::SparseMatrix<double> A(row, col);
		A.setFromTriplets(A_Entries.begin(), A_Entries.end());
		Eigen::SparseMatrix<double> AT = A.transpose();
		Eigen::SparseMatrix<double> ATA = AT * A;

		VectorXd B(row);
		for (size_t i = 0; i < row; i++)
			B(i) = b[i];

		VectorXd ATB = AT * B;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(ATA);
		VectorXd X(col);
		X = chol.solve(ATB);

		return X;
	}

protected:
	bool UpdateVertices() {
		auto X = Solve();
		size_t lastNumOfInvertedElements = UpdateQuality(GettargetMSJ());
		bool converged = true;
		if (Getrecoverable()) {
			Mesh targetMesh(mesh);
			std::vector<glm::dvec3> oldV(mesh.V.size());
			for (auto& v : mesh.V)
				oldV[v.id] = v.xyz();

			for (auto& v : mesh.V) {
				if (fabs(v.x - X[2 * v.id + 0]) > 1e-6 || fabs(v.y - X[2 * v.id + 1]) > 1e-6)
					converged = false;
				if (fabs(v.x - X[2 * v.id + 0]) < 2.0 * GetavgMeshEdgeLength())
					v.x = GetstepSize() * X[2 * v.id + 0] + (1.0 - GetstepSize()) * v.x;
				if (fabs(v.y - X[2 * v.id + 1]) < 2.0 * GetavgMeshEdgeLength())
					v.y = GetstepSize() * X[2 * v.id + 1] + (1.0 - GetstepSize()) * v.y;
			}
			if (GetuseProjection()) {
				if (!GetprojectToTargetSurface()) mesh.ProjectTo(targetMesh);
				else mesh.FastProjectTo(targetMesh);
			}

			size_t numOfInvertedElements = UpdateQuality(GettargetMSJ());
			if (numOfInvertedElements > lastNumOfInvertedElements) {
				std::cout << "Recover previous mesh\n";
				for (size_t i = 0; i < mesh.V.size(); i++)
					mesh.V[i] = oldV[i];
			}
		} else {
			for (auto& v : mesh.V) {
				if (fabs(v.x - X[2 * v.id + 0]) > 1e-6 || fabs(v.y - X[2 * v.id + 1]) > 1e-6)
					converged = false;
				if (fabs(v.x - X[2 * v.id + 0]) < 2.0 * GetavgMeshEdgeLength())
					v.x = GetstepSize() * X[2 * v.id + 0] + (1.0 - GetstepSize()) * v.x;
				if (fabs(v.y - X[2 * v.id + 1]) < 2.0 * GetavgMeshEdgeLength())
					v.y = GetstepSize() * X[2 * v.id + 1] + (1.0 - GetstepSize()) * v.y;
			}
		}
		return converged;
	}
	virtual void AddOrthogonalTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& b, size_t& row) {
		for (auto& e : mesh.E)
			e.energyOrthogonality = 0;
		double energy = 0.0;
		for (auto& edge1 : mesh.E) {
			for (auto edgeId2 : edge1.orthogonalEids) {
				auto& edge2 = mesh.E.at(edgeId2);
				const size_t vId1_1 = edge1.Vids.at(0);
				const size_t vId1_2 = edge1.Vids.at(1);
				const size_t vId2_1 = edge2.Vids.at(0);
				const size_t vId2_2 = edge2.Vids.at(1);
				size_t shareVId = vId1_1;
				if (vId1_2 == vId2_1 || vId1_2 == vId2_2)
					shareVId = vId1_2;
				const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
				const Vertex& vetexC = mesh.V.at(shareVId);
				const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
				const glm::dvec2 v2 = glm::dvec2(vetexC.x - vetex2.x, vetexC.y - vetex2.y);
				const glm::dvec2 v2n = glm::normalize(v2);
				double length_v1 = edge1.length;
				if (!GetuseAverageTargetLength()) {
					length_v1 = glm::length(glm::dvec2(vetexC.x - vetex1.x, vetexC.y - vetex1.y));
					// length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
					length_v1 = length_v1 > GetavgMeshEdgeLength() * Getanisotropy() ? length_v1 : GetavgMeshEdgeLength() * Getanisotropy();
					length_v1 = std::isnan(length_v1) ? GetavgMeshEdgeLength() * Getanisotropy() : length_v1;
				}
				///////////////////////////////////////////////
				const double length_v1_ = 1.0 / length_v1;
				const glm::dvec2 v1n = glm::dvec2((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_);
				const double energy_orthognality = (glm::dot(v1n, v2n) - 0.0) * (glm::dot(v1n, v2n) - 0.0) * Getgamma();
				energy += energy_orthognality;
				edge1.energyOrthogonality += energy_orthognality;

				double weight = Getgamma();
				for (size_t n = 0; n < 2; n++) {
					double a = weight * v2n[n] / length_v1;
					a = std::isnan(a) ? 1.0 : a;
					A_Entries.push_back(Eigen::Triplet<double>(row, 2 * vetexC.id + n, a));
					A_Entries.push_back(Eigen::Triplet<double>(row, 2 * vetex1.id + n, -a));
				}
				b.push_back(0.0 * weight);
				row++;
			}
		}
		EOrthogonality.push_back(energy);
	}
	virtual void AddStraightTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& b, size_t& row) {
		for (auto& e : mesh.E)
			e.energyStraightness = 0;
		double initgamma = Getgamma();
		Setgamma(1.0);
		double energy = 0.0;
		for (auto& edge1 : mesh.E) {
			for (auto edgeId2 : edge1.consecutiveEids) {
				/*const */Edge& edge2 = mesh.E.at(edgeId2);
				const size_t vId1_1 = edge1.Vids.at(0);
				const size_t vId1_2 = edge1.Vids.at(1);
				const size_t vId2_1 = edge2.Vids.at(0);
				const size_t vId2_2 = edge2.Vids.at(1);
				size_t shareVId = vId1_1;
				if (vId1_2 == vId2_1 || vId1_2 == vId2_2) shareVId = vId1_2;
				const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
				const Vertex& vetexC = mesh.V.at(shareVId);
				const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
				const glm::dvec2 v2 = glm::dvec2(vetexC.x - vetex2.x, vetexC.y - vetex2.y);
				const glm::dvec2 v2n = glm::normalize(v2);
				double length_v1 = edge1.length;
				if (!GetuseAverageTargetLength()) {
					length_v1 = glm::length(glm::dvec2(vetexC.x - vetex1.x, vetexC.y - vetex1.y));
					// length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
					length_v1 = length_v1 > GetavgMeshEdgeLength() * Getanisotropy() ? length_v1 : GetavgMeshEdgeLength() * Getanisotropy();
					length_v1 = std::isnan(length_v1) ? GetavgMeshEdgeLength() * Getanisotropy() : length_v1;
				}
				///////////////////////////////////////////////
				const double length_v1_ = 1.0 / length_v1;
				const glm::dvec2 v1n = glm::dvec2((vetexC.x - vetex1.x) * length_v1_, (vetexC.y - vetex1.y) * length_v1_);
				const double energy_straightness = (glm::dot(v1n, v2n) + 1.0) * (glm::dot(v1n, v2n) + 1.0) * Getgamma();
				energy += energy_straightness;
				edge1.energyStraightness += energy_straightness;
				///////////////////////////////////////////////
				// v1 is variable, and use v2 as constant
				// <v1/length_v1, v2n> = -1
				// double weight = energy_straightness * 10 > gamma ? energy_straightness  * 10: gamma;
				// if (stepSize > 0.9)
				double weight = Getgamma();
				for (size_t n = 0; n < 2; n++) {
					double a = weight * v2n[n] / length_v1;
					a = std::isnan(a) ? 1.0 : a;
					A_Entries.push_back(Eigen::Triplet<double>(row, 2 * vetexC.id + n, a));
					A_Entries.push_back(Eigen::Triplet<double>(row, 2 * vetex1.id + n, -a));
				}
				b.push_back(-1.0 * weight);
				row++;
			}
		}
		EStraightness.push_back(energy);
		Setgamma(initgamma);
	}
	virtual void AddSingularTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& B, size_t& row) {
		for (auto& e : mesh.E)
			e.energySingularity = 0;
		double energy = 0.0;
	}
	virtual size_t AddBoundaryTerm(std::vector<Eigen::Triplet<double>>& A_Entries, std::vector<double>& b, size_t& row) {
		auto alpha = Getalpha();
		auto beta = Getbeta();
		size_t a_id = 0;
		for (const auto& v : mesh.V) {
			if (!v.isBoundary) continue;
			for (size_t n = 0; n < 3; n++) {
				//                int a_index = a_id++ + 3 * mesh.V.size();
				//                A_Entries.push_back(Eigen::Triplet<double>(row, a_index, alpha));
				//                b.push_back(0);
				//                row++;
				A_Entries.push_back(Eigen::Triplet<double>(row, 2 * v.id + n, alpha));
				b.push_back(alpha * v[n]);
				row++;
			}
		}
		//        return a_id;

		for (size_t i = 0; i < mesh.V.size(); i++) {
			const Vertex& v = mesh.V.at(i);
			if (!v.isBoundary)
				continue;
			if (v.type == REGULAR) {
				////beta(nv + d)
				//const glm::dvec2& implicit_n = v.normal;
				//float implicit_d = glm::dot(implicit_n, v.xyz());
				//for (size_t k = 0; k < 2; k++) {
				//	A_Entries.push_back(Eigen::Triplet<double>(row, 2 * i + k, beta * implicit_n[k]));
				//}
				//b.push_back(beta * implicit_d);
				//row++;
			} else if (v.type == FEATURE) {
				int a_index = a_id + 2 * mesh.V.size();
				for (size_t k = 0; k < 2; k++) {
					A_Entries.push_back(Eigen::Triplet<double>(row, 2 * i + k, alpha));
					A_Entries.push_back(Eigen::Triplet<double>(row, a_index, -alpha * v.tangent[k]));
					b.push_back(alpha * v[k]);
					row++;
				}
			} else if (v.type == CORNER) {
				for (size_t k = 0; k < 2; k++) {
					A_Entries.push_back(Eigen::Triplet<double>(row, 2 * i + k, alpha));
					b.push_back(alpha * v[k]);
					row++;
				}
			}
		}
		a_id = 0;
		for (size_t i = 0; i < mesh.V.size(); i++) {
			const Vertex& v = mesh.V.at(i);
			if (v.isBoundary) {
				int a_index = a_id++ + 3 * mesh.V.size();
				A_Entries.push_back(Eigen::Triplet<double>(row, a_index, alpha));
				b.push_back(0);
				row++;
			}
		}

		return a_id;
	}
private:
	std::vector<double> ESingularity;
	std::vector<double> EOrthogonality;
	std::vector<double> EStraightness;
};

//class OpenQuadMeshMSJImprover : public QuadMeshMSJImprover {
//public:
//    OpenQuadMeshMSJImprover(Mesh& mesh, const Mesh* refMesh = nullptr): QuadMeshMSJImprover(mesh, refMesh) {
//
//    }
//    virtual ~OpenQuadMeshMSJImprover() {
//
//    }
//private:
//    OpenQuadMeshMSJImprover();
//    OpenQuadMeshMSJImprover(const OpenQuadMeshMSJImprover&);
//    OpenQuadMeshMSJImprover& operator = (const OpenQuadMeshMSJImprover&);
//public:
//    int Run(int iters = 20) {
//        int iter = 0;
//        while (iter++ < iters) {
//            UpdateQuality(GettargetMSJ());
//            std::cout << "MSJ = " << GetQuality() << std::endl;
//            auto invertedVids = GetInvertedVids();
//            for (auto& vid : invertedVids) {
//                auto& v = mesh.V.at(vid);
//                if (v.isBoundary) continue;
//                auto sampleVertices = GetSampleVertices(v.id, 20);
//                auto bestSampleVertexId = GetBestSampleVertexId(v.id, sampleVertices);
//                v = sampleVertices.at(bestSampleVertexId);
//            }
//        }
//        MeshFileWriter optwriter(mesh, "best.vtk");
//        optwriter.WriteFile();
//        return iter;
//    }
//private:
//    std::vector<size_t> GetInvertedVids() const {
//        std::map<double, size_t> sj_f;
////        int minMSJ = 0;
////        for (auto& f : mesh.F)
////            if (GetScaledJacobianQuad(mesh, f.id) < 0)
////                res.insert(f.Vids.begin(), f.Vids.end());
//        return {};
//    }
//    std::vector<glm::dvec3> GetSampleVertices(size_t vid, int sampleResolution = 10) {
//        auto res = GetSampleVerticesInFaces(vid, sampleResolution);
//        auto vids = GetSampleVerticesInEdges(vid, sampleResolution);
//        std::copy(vids.begin(), vids.end(), std::back_inserter(res));
//        vids = GetSampleVerticesInVertices(vid, sampleResolution);
//        std::copy(vids.begin(), vids.end(), std::back_inserter(res));
//        return res;
//    }
//    size_t GetBestSampleVertexId(size_t vid, const std::vector<glm::dvec3>& samples) {
//        size_t res = 0;
//        double maxMSJ = -1.0;
//        double avgMSJ = -1.0;
//        auto& v = mesh.V.at(vid);
//        const auto orig = v.xyz();
//        int id = 0;
//        for (auto& sample : samples) {
//            v = sample;
//            double localMSJ = 1.0;
//            for (auto fid : v.N_Fids)
//                localMSJ = std::min(localMSJ, GetScaledJacobianQuad(mesh, fid));
//            double totalSJ = 0.0;
//            for (auto fid : v.N_Fids)
//                totalSJ += GetScaledJacobianQuad(mesh, fid);
//            if (/*fabs(localMSJ) > 0.05 && */localMSJ > maxMSJ && totalSJ / v.N_Fids.size() > avgMSJ) {
//                maxMSJ = localMSJ;
//                avgMSJ = totalSJ / v.N_Fids.size();
//                res = id;
//            }
//            ++id;
//        }
//        v = orig;
//        return res;
//    }
//
//    std::vector<glm::dvec3> GetSampleVerticesInFaces(size_t vid, int resolution = 10) {
//        std::vector<glm::dvec3> res;
//        auto& v = mesh.V.at(vid);
//        glm::dvec3 vertex;
//        for (auto fid : v.N_Fids) {
//            auto& f = mesh.F.at(fid);
//            for (double u = 0; u < resolution; ++u)
//                for (double v = 0; v < resolution; ++v) {
//                    auto base = 1.0 / (resolution + 1);
//                    auto v01 = (1.0 + u) * base * mesh.V.at(f.Vids[0]).xyz() + (resolution - u) * base * mesh.V.at(f.Vids[1]).xyz();
//                    auto v32 = (1.0 + u) * base * mesh.V.at(f.Vids[3]).xyz() + (resolution - u) * base * mesh.V.at(f.Vids[2]).xyz();
//                    vertex = (1.0 + v) * base * v01 + (resolution - v) * base * v32;
//                    res.push_back(vertex);
//                }
//        }
//        return res;
//    }
//
//    std::vector<glm::dvec3> GetSampleVerticesInEdges(size_t vid, int resolution = 10) {
//        std::vector<glm::dvec3> res;
//        auto& v = mesh.V.at(vid);
//        glm::dvec3 vertex;
//        for (auto eid : v.N_Eids) {
//            auto& e = mesh.E.at(eid);
//            for (double u = 0; u < resolution; ++u) {
//                auto base = 1.0 / (resolution + 1);
//                vertex = (1.0 + u) * base * mesh.V.at(e.Vids[0]).xyz() + (resolution - u) * base * mesh.V.at(e.Vids[1]).xyz();
//                res.push_back(vertex);
//            }
//        }
//        return res;
//    }
//
//    std::vector<glm::dvec3> GetSampleVerticesInVertices(size_t vid, int resolution = 10) {
//        std::vector<glm::dvec3> res;
//        auto& v = mesh.V.at(vid);
////        glm::dvec3 vertex;
////        for (auto vid : v.N_Vids)
////            res.push_back(mesh.V.at(vid).xyz());
//        res.push_back(v.xyz());
//        return res;
//    }
//};

const double PI = 3.1415926536;
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: MeshOpt input.vtk output.vtk iters=50 alpha=1000000 beta=1000000 gamma=1.0 cosangle=0.984807753 "
                  << "stepSize=0.9 anisotropy=0.5 useAvgLength=false converge=true allowBigStep=false\n targetMSJ=0.0";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = "input.vtk";
    int iters = 50;
    double alpha = 1000000.0;
    double beta = 1000000.0;
    double gamma = 1.0;
    double stepSize = 0.9;
    double anisotropy = 0.05;
    double angle = 170;
    double targetMSJ = 0.0;
    bool useAvgLength = false;
    bool allowBigStep = false;
    bool changeBoundary = false;
    bool converge = true;
    bool useProjection = false;
    {
        argumentManager.get("iters", iters);
        argumentManager.get("alpha", alpha);
        argumentManager.get("beta", beta);
        argumentManager.get("gamma", gamma);
        argumentManager.get("stepSize", stepSize);
        argumentManager.get("anisotropy", anisotropy);
        argumentManager.get("targetMSJ", targetMSJ);
        argumentManager.get("angle", angle);
        argumentManager.get("useAvgLength", useAvgLength);
        argumentManager.get("converge", converge);
        argumentManager.get("allowBigStep", allowBigStep);
        argumentManager.get("changeBoundary", changeBoundary);
        argumentManager.get("useProjection", useProjection);

        std::cout << "-----------------------------------\n";
        std::cout << "input = " << filename << std::endl;
        std::cout << "iters = " << iters << std::endl;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "gamma = " << gamma << std::endl;
        std::cout << "stepSize = " << stepSize << std::endl;
        std::cout << "anisotropy = " << anisotropy << std::endl;
        std::cout << "targetMSJ = " << targetMSJ << std::endl;
        std::cout << "angle = " << angle << std::endl;
        std::cout << "useAvgLength = " << useAvgLength << std::endl;
        std::cout << "converge = " << converge << std::endl;
        std::cout << "allowBigStep = " << allowBigStep << std::endl;
        std::cout << "changeBoundary = " << changeBoundary << std::endl;
        std::cout << "useProjection = " << useProjection << std::endl;
        std::cout << "-----------------------------------\n";
    }
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    mesh.SetCosAngleThreshold(cos((180.0 - angle) * PI / 180.0));
    mesh.LabelSurface();
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();
    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;


    size_t nF = 0;
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary) ++nF;
    std::vector<Face> F(nF);
    nF = 0;
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary) {
            F[nF].id = nF;
            F[nF++].Vids = mesh.F.at(i).Vids;
        }
    std::vector<Cell> C(nF);
    nF = 0;
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary) {
            C[nF].id = nF;
            C[nF++].Vids = mesh.F.at(i).Vids;
        }
    //Mesh surfaceMesh(mesh.V, C, QUAD);

//    srand (time(NULL));
//    for (auto& v : mesh.V)
//        v.z += (rand() % 10) * 1e-4 ;
    MeshFileWriter surfaceWriter(mesh.V, C, "surface.quad.off", QUAD);
    surfaceWriter.WriteFile();
    MeshFileReader surfaceReader("surface.quad.off");
    Mesh& surfaceMesh = (Mesh&)surfaceReader.GetMesh();

    surfaceMesh.RemoveUselessVertices();
    surfaceMesh.BuildAllConnectivities();
    surfaceMesh.ExtractBoundary();
    surfaceMesh.ExtractSingularities();
    surfaceMesh.SetCosAngleThreshold(cos((180.0 - angle) * PI / 180.0));
	surfaceMesh.SetFeatureAngleThreshold(angle);
    //surfaceMesh.LabelSurface();
	surfaceMesh.Label2DSurfaceVertices();
    surfaceMesh.BuildParallelE();
    surfaceMesh.BuildConsecutiveE();
    surfaceMesh.BuildOrthogonalE();
    surfaceMesh.GetNormalOfSurfaceFaces();
    surfaceMesh.GetNormalOfSurfaceVertices();

    OpenQuadMeshMSJImprover meshOpt(surfaceMesh, &surfaceMesh);
    meshOpt.Setalpha(alpha);
    meshOpt.Setbeta(beta);
    meshOpt.Setgamma(gamma);
    meshOpt.SetstepSize(stepSize);
    meshOpt.SetuseAverageTargetLength(useAvgLength);
    meshOpt.Setanisotropy(anisotropy);
    meshOpt.SettargetMSJ(targetMSJ);
    meshOpt.Setrecoverable(converge);
    meshOpt.SetallowBigStep(allowBigStep);
    meshOpt.SetuseProjection(useProjection);
//    meshOpt.SetChangeBoundary(changeBoundary);

    size_t iter = meshOpt.Run(iters);
//    MeshFileWriter optwriter(mesh, "out.vtk");
//    optwriter.WriteFile();
    return 0;
}

