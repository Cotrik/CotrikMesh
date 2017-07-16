/*
 * CheckMesh..cpp
 *
 *  Created on: Dec 9, 2016
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshQuality.h"
#include <iostream>
bool IsAllInvertedCellOnBoundary(const Mesh& mesh, const std::vector<size_t>& badCellIds, size_t& innerCount);
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: CheckMesh vtkFile\n";
        return -1;
    }
    std::string filename(argv[1]);
    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();

    double minimumScaledJacobian = 0.0;
    double averageScaledJacobian = 0.0;
    double maximumScaledJacobian = 0.0;
    std::vector<size_t> badCellIds;
    size_t numOfInvertdElements = GetQuality(filename.c_str(), minimumScaledJacobian, averageScaledJacobian, maximumScaledJacobian, badCellIds);

    size_t innerCount = 0;
    if (!IsAllInvertedCellOnBoundary(mesh, badCellIds, innerCount))
        std::cout << "\033[1;31m" <<  innerCount << " inverted cell NOT on boundary!!!\033[0m\n";
    else
        std::cout << "\033[1;32mAll inverted cell on boundary\033[0m\n";

    return 0;
}

bool IsAllInvertedCellOnBoundary(const Mesh& mesh, const std::vector<size_t>& badCellIds, size_t& innerCount)
{
    for (size_t i = 0; i < badCellIds.size(); i++)
        if (!mesh.C.at(badCellIds.at(i)).isBoundary)
            innerCount++;

    if (innerCount == 0)
        return true;
    return false;
}
