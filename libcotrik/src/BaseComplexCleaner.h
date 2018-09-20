/*
 * BaseComplexCleaner.h
 *
 *  Created on: Aug 28, 2018
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXCLEANER_H_
#define LIBCOTRIK_SRC_BASECOMPLEXCLEANER_H_

#include "BaseComplex.h"

class BaseComplexCleaner {
public:
    BaseComplexCleaner(BaseComplex& baseComplex);
    virtual ~BaseComplexCleaner();
private:
    BaseComplexCleaner(const BaseComplexCleaner&);
    BaseComplexCleaner();
    BaseComplexCleaner& operator = (const BaseComplexCleaner&);
public:
    void Run();

private:
    bool NeighborOnBoundary(const SingularV& sv) const;
    bool AllNeighborsValence4(const SingularV& sv) const;
    bool CubeStructureLabeled(const SingularV& sv, const std::vector<bool>& labeled) const;
    bool IsCubeStructureCorrect(const SingularV& sv) const;
    void LabelCubeStructure(const SingularV& sv, std::vector<bool>& labeled, int& cubelabel, std::vector<int>& celldata) const;
//private:
    BaseComplex& baseComplex;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEXCLEANER_H_ */
