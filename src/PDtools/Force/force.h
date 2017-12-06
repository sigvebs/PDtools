#ifndef FORCE_H
#define FORCE_H

#include <unordered_map>
#include <armadillo>

#include "Elements/pd_element.h"

using namespace std;
namespace PDtools
{
// Forward declerations
class PD_Particles;
//class PD_element;
//class PD_quadElement;

//------------------------------------------------------------------------------
class Force
{
protected:
    PD_Particles &m_particles;
    bool m_numericalInitialization = false;
    int m_dim;
    double m_E;
    double m_nu;
    double m_h;
    double m_delta;
    double m_lc;
    bool m_hasSurfaceCorrection = false;
    int m_counter = 0;

    arma::mat & m_r;
    arma::mat & m_v;
    arma::mat & m_r0;
    arma::mat & m_F;
    arma::mat & m_data;
    unordered_map<int, int> & m_idToCol;
    arma::ivec & m_colToId;
    int m_calulateStress = false;

    // Element based Pd (EPD)
    unordered_map<int, int> &m_idToElement;
    vector<PD_quadElement> &m_quadElements;
    vector<PD_triElement> &m_triElements;

    enum enum_coordinates{X, Y, Z};
    static constexpr double THRESHOLD = 2.2204e-016;

    // For force scaling
    const std::vector<std::string> forceScalingStringIds = {"g_x", "g_y", "g_z"};
    arma::ivec3 m_g;

    // Ghost parameters
    vector<string> m_initialGhostParameters;
    vector<string> m_ghostParameters;

    vector<pair<string, int>> neededProperties; // name and updatefrquency

    // Modifiers
    bool m_continueState = false;
    bool m_hasStaticModifier = false;
    bool m_hasStepOneModifier = false;
    bool m_hasStepTwoModifier = false;
    bool m_hasUpdateState = false;
public:
    const string name;
    Force(PD_Particles &particles, string _type="none");

    virtual ~Force();

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    virtual void
    calculateForces(const int id, const int i) = 0;

    virtual double
    calculatePotentialEnergyDensity(const int id_i, const int i);

    virtual void
    calculatePotentialEnergy(const int id_i, const int i,
                                          int indexPotential);
    virtual double
    calculateBondEnergy(const std::pair<int, int> & idCol,
                                          std::pair<int, std::vector<double>> & con);

    virtual void
    calculateStress(const int id_i, const int i,
                                 const int (&indexStress)[6]);
    virtual void
    updateState();

    virtual void
    updateState(int id, int i);

    virtual double
    calculateStableMass(const int id_i, const int i,
                                     double dt);
    void
    numericalInitialization(bool ni);

    void
    setDim(int dim);

    virtual int
    initializeSurfaceCorrection();

    virtual void
    applySurfaceCorrectionStep1(double strain = 0.001);

    virtual void
    applySurfaceCorrectionStep2();

    virtual void
    applyStrainCorrection(double strain);

    virtual void
    applyShearCorrection(double shear=0.001);

    virtual std::vector<string>
    initalGhostDependencies();

    virtual std::vector<string>
    ghostDependencies();

    virtual std::vector<string>
    getSurfaceCorrectionGhostParameters();

    vector<pair<string, int> >
    getNeededProperties() const;
    int getCalulateStress() const;

    // Modifiers
    virtual void
    evaluateStepOne();

    virtual void
    evaluateStepOne(int id, int i);

    virtual void
    evaluateStepTwo(int id, int i);

    virtual void
    evaluateStatic(int id, int i);

    bool getHasStepOneModifier() const;
    bool getHasStepTwoModifier() const;
    bool getContinueState() const;
    bool getHasStaticModifier() const;
    bool getHasUpdateState() const;
};
//------------------------------------------------------------------------------
}
#endif // FORCE_H
