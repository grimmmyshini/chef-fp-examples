#include "clad/Differentiator/EstimationModel.h"

/// This is our dummy estimation model class.
// We will be using this to override the virtual function in the
// FPErrorEstimationModel class.
class PrintModel : public clad::FPErrorEstimationModel {
public:
  PrintModel(clad::DerivativeBuilder& builder)
      : FPErrorEstimationModel(builder) {}
  // Return an expression of the following kind:
  //  dfdx * delta_x
  clang::Expr* AssignError(clad::StmtDiff refExpr, std::string name) override;
};