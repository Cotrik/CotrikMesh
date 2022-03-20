#include "EdgeCollapse.h"

EdgeCollapse::EdgeCollapse() {}
EdgeCollapse::EdgeCollapse(Mesh& mesh_, MeshUtil& mu_, std::shared_ptr<SimplificationOperation> vr_, std::shared_ptr<SimplificationOperation> dc_) : SimplificationOperation(mesh_, mu_) {
    vr = std::move(vr_);
    dc = std::move(dc_);
}
EdgeCollapse::~EdgeCollapse() {}

void EdgeCollapse::SetRanking(glm::dvec3 d) {
    CheckValidity();
}

bool EdgeCollapse::IsOperationValid() {
    CheckValidity();

    return true;
}

void EdgeCollapse::PerformOperation() {
    CheckValidity();
    if (!IsOperationValid()) return;

    vr->PerformOperation();
    dc->PerformOperation();
}
