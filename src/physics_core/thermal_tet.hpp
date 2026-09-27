#pragma once

#include "chrono/core/ChMatrix33.h"
#include "chrono/fea/ChElementGeneric.h"
#include "chrono/fea/ChNodeFEAxyz.h"

#include <memory>

namespace rcs {

/// Four-node corotational tet. Chrono's element is elastic only, so thermal
/// strain and J2 return live here. The tangent stays the temperature-scaled
/// elastic stiffness; the internal force is what sits on the yield surface.
class ThermalTet4 : public chrono::fea::ChElementGeneric {
public:
    ThermalTet4();

    void SetNodes(std::shared_ptr<chrono::fea::ChNodeFEAxyz> a,
                  std::shared_ptr<chrono::fea::ChNodeFEAxyz> b,
                  std::shared_ptr<chrono::fea::ChNodeFEAxyz> c,
                  std::shared_ptr<chrono::fea::ChNodeFEAxyz> d);

    void configure(double youngs, double nu, double yield, double density, double beta);
    void setThermalState(double thermalStrain, double stiffnessScale, double strengthScale);
    /// Rankine limits, Pa. Zero leaves that side uncapped. Brittle solids use
    /// this the way J2 uses yield: the tangent stays elastic, the force cannot
    /// exceed the material strength.
    void setRankine(double tensionLimit, double compressionLimit);
    /// Equivalent plastic strain at which the tet tears. Zero disables the
    /// damage ramp. Past half of this the stress and tangent fade toward a
    /// small residual so the element sheds load before it is removed.
    void setFailureStrain(double strain);

    double vonMises() const { return _vonMises; }
    double plasticStrain() const { return _plastic; }
    double tension() const { return _tension; }
    /// True when the last stress evaluation had to project onto Rankine.
    bool overstrength() const { return _overstrength; }
    const chrono::ChVector3d& tensionAxis() const { return _tensionAxis; }
    double volume() const { return _volume; }

    unsigned int GetNumNodes() override { return 4; }
    unsigned int GetNumCoordsPosLevel() override { return 12; }
    unsigned int GetNodeNumCoordsPosLevel(unsigned int) override { return 3; }
    std::shared_ptr<chrono::fea::ChNodeFEAbase> GetNode(unsigned int n) override { return _nodes[n]; }

    void GetStateBlock(chrono::ChVectorDynamic<>& state) override;
    void ComputeKRMmatricesGlobal(chrono::ChMatrixRef H, double Kfactor, double Rfactor, double Mfactor) override;
    void ComputeInternalForces(chrono::ChVectorDynamic<>& Fi) override;
    void SetupInitial(chrono::ChSystem* system) override;
    void Update() override;

    /// Keep the plastic strain from the last converged step. Called once after
    /// the implicit solve, so Newton iterations do not ratchet it.
    void commitPlastic();
    /// Throw away a rejected Newton iterate so a failed impact cannot store
    /// a 1e8 plastic strain or a Rankine flag that the positions do not have.
    void discardTrial();
    void copyPlasticFrom(const ThermalTet4& other);

private:
    void noteTension(const chrono::ChVectorDynamic<>& stress);
    void refreshRotation();
    void rebuildReference();
    void elasticLaw(chrono::ChMatrixDynamic<>& C) const;
    void trialStress(chrono::ChVectorDynamic<>& stress, chrono::ChVectorDynamic<>& strain);
    /// Fraction of the elastic law still carried. Damage fades it near tear.
    /// The yield return caps the force; the tangent stays elastic so a soft
    /// punch can still converge. Fast impacts use contact substeps instead.
    double tangentScale() const;

    std::shared_ptr<chrono::fea::ChNodeFEAxyz> _nodes[4];
    chrono::ChMatrix33<> _rotation;
    chrono::ChMatrixNM<double, 4, 4> _shape;
    chrono::ChMatrixDynamic<> _B;
    chrono::ChMatrixDynamic<> _K;
    double _volume = 0.0;
    double _youngs = 1.0;
    double _nu = 0.33;
    double _yield = 0.0;
    double _failPlastic = 0.0;
    double _rankineTension = 0.0;
    double _rankineCompression = 0.0;
    double _density = 1.0;
    double _beta = 0.0;
    double _thermal = 0.0;
    double _stiffness = 1.0;
    double _strength = 1.0;
    double _vonMises = 0.0;
    double _tension = 0.0;
    chrono::ChVector3d _tensionAxis{1, 0, 0};
    double _plastic = 0.0;
    bool _overstrength = false;
    double _plasticStrain[6] = {};
    double _plasticPending[6] = {};
    bool _ready = false;
};

}  // namespace rcs
