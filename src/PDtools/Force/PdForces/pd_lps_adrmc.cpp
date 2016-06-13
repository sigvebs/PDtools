#include "pd_lps_adrmc.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_LPS_adrmc::PD_LPS_adrmc(PD_Particles &particles, double phi, double C, double T, bool planeStress):
    LPS_mc(particles, 1., phi, C, T, planeStress)
{
    particles.setNeedGhostVelocity(false);
    m_hasStaticModifier = true;
}
//------------------------------------------------------------------------------
void PD_LPS_adrmc::calculateForces(const int id, const int i)
{
    const double theta_i = m_data(i, m_iTheta);
    const double m_i = m_data(i, m_iMass);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);

    const int nConnections = PDconnections.size();
    double dr0_ij[m_dim];
    double dr_ij[m_dim];
    _F.zeros();

    double thetaNew = 0;
    int nConnected = 0;

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];

        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double m_j = m_data(j, m_iMass);
        const double theta_j = m_data(j, m_iTheta);
        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling = con.second[m_iVolumeScaling];
        const double vol = vol_j*volumeScaling;
        const double w = 1./dr0;

        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        double bond = m_c*(theta_i/m_i + theta_j/m_j)*dr0;
        bond += m_alpha*(1./m_i + 1./m_j)*ds;
        bond *= w*vol/dr;
        thetaNew += w*dr0*ds*vol;

        for(int d=0; d<m_dim; d++)
        {
            m_F(i, d) += dr_ij[d]*bond;

            for(int d2=0; d2<m_dim; d2++)
            {
                _F(d, d2) += m_delta*w*dr_ij[d]*dr0_ij[d2]*vol;
            }
        }

        nConnected++;
    }

    if(nConnections <= 3)
    {
        m_data(i, m_iThetaNew) = 0;
    }
    else
        m_data(i, m_iThetaNew) = m_dim/m_i*thetaNew;

    computeStress(id, i, nConnected);
    m_continueState = false;

//    //--------------------------------------------------------------------------
//    // Stress
//    //--------------------------------------------------------------------------
//    if(m_data(i, m_indexBrokenNow))
//    {
//        computeK(id, i);
//        m_data(i, m_indexBrokenNow) = 0;
//    }
//    m_K(0, 0) = m_data(i, m_indexK[0]);
//    m_K(1, 1) = m_data(i, m_indexK[1]);
//    m_K(0, 1) = m_data(i, m_indexK[2]);
//    m_K(1, 0) = m_K(0, 1);

//    if(nConnected <= 3)
//    {
//        m_data(i, m_indexStrain[0]) = 0;
//        m_data(i, m_indexStrain[1]) = 0;
//        m_data(i, m_indexStrain[2]) = 0;
//        m_data(i, m_indexStress[0]) = 0;
//        m_data(i, m_indexStress[1]) = 0;
//        m_data(i, m_indexStress[2]) = 0;
//        return;
//    }
//    else
//    {
//        _F = _F*m_K; // K = inv(K);
//        if(m_smallStrain)
//        {
//            m_strain = 0.5*(_F.t() + _F);
//            m_strain(0,0) -= 1.;
//            m_strain(1,1) -= 1.;

//        }
//        else
//        {
//            m_strain = 0.5*_F.t()*_F;
//            m_strain(0,0) -= 0.5;
//            m_strain(1,1) -= 0.5;
//        }
//    }

//    m_data(i, m_indexStrain[0]) = m_strain(0,0);
//    m_data(i, m_indexStrain[1]) = m_strain(1,1);
//    m_data(i, m_indexStrain[2]) = m_strain(0,1);

//    // Assuming linear elasticity
//    if(m_dim == 2)
//    {
//        // Constituent model, linear elastic
//        // Computing the second PK stress
//        if(m_planeStress)
//        {
//            const double a = m_E/(1. - m_nu*m_nu);
//            m_P(0,0) = a*(m_strain(0,0) + m_nu*m_strain(1,1));
//            m_P(1,1) = a*(m_strain(1,1) + m_nu*m_strain(0,0));
//            m_P(0,1) = a*(1 - m_nu)*m_strain(0,1);
//            m_P(1,0) = m_P(0,1);
//        }
//        else // Plane strain
//        {
//            const double a = m_E/((1. + m_nu)*(1. - 2*m_nu));
//            m_P(0,0) = a*((1. - m_nu)*m_strain(0,0) + m_nu*m_strain(1,1));
//            m_P(1,1) = a*((1. - m_nu)*m_strain(1,1) + m_nu*m_strain(0,0));
//            m_P(0,1) = a*0.5*(1. - 2.*m_nu)*m_strain(0,1);
//            m_P(1,0) = m_P(0,1);
//        }
//        // Converting PK stress to Cauchy stress
//        if(m_greenStrain)
//        {
//            const double detF = 1./det(_F);
//            m_P = detF*_F*m_P*_F.t();
//        }
//        m_data(i, m_indexStress[0]) = m_P(0,0);
//        m_data(i, m_indexStress[1]) = m_P(1,1);
//        m_data(i, m_indexStress[2]) = m_P(0,1);
//    }

}
//------------------------------------------------------------------------------
void PD_LPS_adrmc::evaluateStatic(int id, int i)
{
    if(m_data(i, m_indexUnbreakable) >= 1)
        return;
#if CALCULATE_NUMMERICAL_PRINCIPAL_STRESS
    arma::vec eigval(m_dim);
#endif
    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);

    if(m_dim == 2)
    {
        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = m_idToCol[id_j];

            if(m_data(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const double sx = 0.5*(m_data(i, m_indexStress[0]) + m_data(j, m_indexStress[0]));
            const double sy = 0.5*(m_data(i, m_indexStress[1]) + m_data(j, m_indexStress[1]));
            const double sxy = 0.5*(m_data(i, m_indexStress[2]) + m_data(j, m_indexStress[2]));

            const double first = 0.5*(sx + sy);
            const double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

            double p_2 = first + second; // max
            double p_1 = first - second; // min

            const double shear = fabs(0.5*(p_1 - p_2)*m_sin_theta);
            const double normal = 0.5*(p_1 + p_2) + 0.5*(p_1 - p_2)*m_cos_theta;

            const double criticalShear = fabs(shear) - fabs(m_C - m_d*normal);
            const double criticalTensile = p_2 - m_T;

//            if(criticalShear >= 0  && normal < 0)
            if(criticalShear >= 0)
            {
                m_data(i, m_indexBrokenNow) = 1;
                con.second[m_indexConnected] = 0;
                m_continueState = true;
            }
            else if(criticalTensile >= 0)
            {
                m_data(i, m_indexBrokenNow) = 1;
                con.second[m_indexConnected] = 0;
                m_continueState = true;
            }
        }
    }
}

//------------------------------------------------------------------------------
void PD_LPS_adrmc::evaluateStepTwo(int id, int i)
{
    (void) id;
    (void) i;
}
//------------------------------------------------------------------------------
}
