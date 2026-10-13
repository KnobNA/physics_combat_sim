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
    /// Handbook UTS. With fail plastic this sets isotropic hardening so the
    /// 0.2% offset sits at `yieldStress`.
    void setHardening(double ultimate);

    double vonMises() const { return _vonMises; }
    double plasticStrain() const { return _plastic; }
    double tension() const { return _tension; }
    /// Principal tension before Rankine projection. Peel uses this with
    /// Griffith energy so sitting on the cap is not the same as a crack.
    double rawTension() const { return _rawTension; }
    /// Magnitude of the most compressive principal, before Rankine. Used to
    /// tell a hoop-tension spike under crush from a tension-dominated crack.
    double rawCompression() const { return _rawCompression; }
    /// Hydrostatic mean of the trial stress, before Rankine.
    double rawMean() const { return _rawMean; }
    double trialEnergy() const { return _trialEnergy; }
    double tensileStrain() const { return _tensileStrain; }
    /// 1 when an accepted step met the Griffith / tensile-strain gate.
    double damage() const { return _damage; }
    /// True when the last stress evaluation had to project onto Rankine.
    bool overstrength() const { return _overstrength; }
    const chrono::ChVector3d& tensionAxis() const { return _tensionAxis; }
    double volume() const { return _volume; }
    /// True when the last deformation gradient was singular, so polar
    /// decomposition was skipped and this element contributed no force.
    bool deformationSingular() const { return _collapsed; }
    /// True when the last F had a healthy negative determinant.
    bool inverted() const { return _inverted; }

    unsigned int GetNumNodes() override { return 4; }
    unsigned int GetNumCoordsPosLevel() override { return 12; }
    unsigned int GetNodeNumCoordsPosLevel(unsigned int) override { return 3; }
    std::shared_ptr<chrono::fea::ChNodeFEAbase> GetNode(unsigned int n) override { return _nodes[n]; }

    void GetStateBlock(chrono::ChVectorDynamic<>& state) override;
    void ComputeKRMmatricesGlobal(chrono::ChMatrixRef H, double Kfactor, double Rfactor, double Mfactor) override;
    void ComputeInternalForces(chrono::ChVectorDynamic<>& Fi) override;
    void SetupInitial(chrono::ChSystem* system) override;
    void Update() override;

    /// Rebuild B and rest volume from the current X0. Placement moves X0
    /// after SetupInitial; without this the corotated strain is in the
    /// authoring frame while the nodes live in world coordinates.
    void refreshRest() { rebuildReference(); }
    /// Length of the current HHT slice. Trial strain can include the
    /// deformational velocity over this window so a stiff Newton iterate
    /// that has not yet turned closing into displacement still yields.
    void setImpactStep(double h) { _impactStep = std::max(0.0, h); }
    /// Keep the plastic strain from the last converged step. Called once after
    /// the implicit solve, so Newton iterations do not ratchet it.
    void commitPlastic();
    /// Brittle crack variable. Set on accepted steps only; a rejected iterate
    /// must not peel from a leftover Rankine clip.
    void commitDamage(double value);
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
    void deformationalSpeed(chrono::ChVectorDynamic<>& speed) const;
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
    double _ultimate = 0.0;
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
    double _rawTension = 0.0;
    double _rawCompression = 0.0;
    double _rawMean = 0.0;
    double _trialEnergy = 0.0;
    double _tensileStrain = 0.0;
    double _damage = 0.0;
    chrono::ChVector3d _tensionAxis{1, 0, 0};
    double _plastic = 0.0;
    bool _overstrength = false;
    double _plasticStrain[6] = {};
    double _plasticPending[6] = {};
    double _impactStep = 0.0;
    chrono::ChMatrix33<> _stretch;
    bool _ready = false;
    bool _collapsed = false;
    bool _inverted = false;
};

}  // namespace rcs
