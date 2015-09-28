#include "pdfunctions.h"

#include <armadillo>
#include <unordered_map>
#include "Particles/pd_particles.h"
#include "Grid/grid.h"
#include "PDtools/Force/force.h"

#define USE_EXTENDED_RANGE_RADIUS 0
#define USE_EXTENDED_RANGE_LC 0

namespace PDtools
{
//------------------------------------------------------------------------------
void setPdConnections(PD_Particles & particles,
                      Grid & grid,
                      double delta,
                      double lc)
{
    const unordered_map<int, GridPoint*> &gridpoints = grid.gridpoints();
    const vector<int> &mygridPoints = grid.myGridPoints();
    const mat & R = particles.r();
    const int indexRadius = particles.getParamId("radius");
    const mat &data  = particles.data();

    // The order is important!
    particles.registerPdParameter("dr0");
    particles.registerPdParameter("connected");

    int nFour = 0;
    int nElse = 0;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<mygridPoints.size(); i++)
    {
        double dx, dy, dz;
        int gridId = mygridPoints.at(i);
        const GridPoint & gridPoint = *gridpoints.at(gridId);

        for(const pair<int, int> & idCol_i:gridPoint.particles())
        {
            int id_i = idCol_i.first;
            int col_i = idCol_i.second;
            const vec & r_i = R.col(col_i);
            unordered_map<int, vector<double>> connections;
            vector<pair<int, vector<double>>> connectionsVector;
            const double radius_i = data(col_i, indexRadius);

            for(const pair<int, int> & idCol_j:gridPoint.particles())
            {
                const int id_j = idCol_j.first;
                const int col_j = idCol_j.second;
                if(id_i == id_j)
                    continue;

                const double radius_j = data(col_j, indexRadius);
#if USE_EXTENDED_RANGE_RADIUS
                const double l_delta = delta + 0.5*(radius_i + radius_j);
#elif USE_EXTENDED_RANGE_LC
                const double l_delta = delta + 0.5*lc;
#else
                const double l_delta = delta;
#endif
                const double deltaSquared = pow(l_delta, 2);

                const vec & r_j = R.col(col_j);
                dx = r_i(0) - r_j(0);
                dy = r_i(1) - r_j(1);
                dz = r_i(2) - r_j(2);

                double drSquared = dx*dx + dy*dy + dz*dz;

                if(drSquared <= deltaSquared)
                {
                    vector<double> connectionData;
                    const double dr = sqrt(drSquared);

                    connectionData.push_back(dr);
                    connectionData.push_back(1.0); // Connected
                    connections[id_j] = connectionData;
                    connectionsVector.push_back(pair<int, vector<double>>(id_j, connectionData) );
                }
            }

            // Neighbouring cells
            const vector<GridPoint*> & neighbours = gridPoint.neighbours();

            for(const GridPoint *neighbour:neighbours)
            {
                for(const pair<int, int> & idCol_j:neighbour->particles())
                {
                    const int col_j = idCol_j.second;

                    const double radius_j = data(col_j, indexRadius);
#if USE_EXTENDED_RANGE_RADIUS
                const double l_delta = delta + 0.5*(radius_i + radius_j);
#elif USE_EXTENDED_RANGE_LC
                const double l_delta = delta + 0.5*lc;
#else
                const double l_delta = delta;
#endif
                    const double deltaSquared = pow(l_delta, 2);

                    const vec & r_j = R.col(col_j);
                    dx = r_i(0) - r_j(0);
                    dy = r_i(1) - r_j(1);
                    dz = r_i(2) - r_j(2);

                    double drSquared = dx*dx + dy*dy + dz*dz;

                    if(drSquared <= deltaSquared)
                    {
                        int id_j = idCol_j.first;
                        vector<double> connectionData;
                        const double dr = sqrt(drSquared);
                        connectionData.push_back(dr);
                        connectionData.push_back(1.0); // Connected
                        connections[id_j] = connectionData;
                        connectionsVector.push_back(pair<int, vector<double>>(id_j, connectionData) );
                    }
                }
            }

#ifdef USE_OPENMP
#pragma omp critical
            {
                particles.setPdConnections(id_i, connectionsVector);

                if(connectionsVector.size() == 4)
                    nFour++;
                else
                    nElse++;


            }
#else
            particles.setPdConnections(id_i, connectionsVector);
#endif
        }
    }

    double scale = nFour / double(nFour + nElse);

    cout << scale << endl;

}
//------------------------------------------------------------------------------
void applyVolumeCorrection(PD_Particles &particles, double delta, double lc)
{
    arma::mat & data = particles.data();
    std::unordered_map<int, int> &pIds = particles.pIds();
    const int indexDr0 = particles.getPdParamId("dr0");
    const int indexRadius = particles.getParamId("radius");
    const int indexVolumeScaling = particles.getPdParamId("volumeScaling");
    const int indexVolume= particles.getParamId("volume");
    double avg = 0;

#ifdef USE_OPENMP
# pragma omp parallel for reduction(+:avg)
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        int pId = i;
        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
        double vol_delta = 0;
        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = pIds[id_j];
            const double dr = con.second[indexDr0];
#if USE_EXTENDED_RANGE_RADIUS
            const double radius_j = data(col_j, indexRadius);
#elif USE_EXTENDED_RANGE_LC
            const double radius_j = 0.5*lc;
#else
            const double radius_j = data(col_j, indexRadius);
#endif
            const double rc = delta - radius_j;

            double volumeCorrection = 1.0;
            if(dr > rc)
            {
                volumeCorrection = 0.5*(delta + radius_j - dr)/radius_j;
            }
            con.second[indexVolumeScaling] = volumeCorrection;
            //con.second[indexVolumeScaling] = 1;
            const double vol_i = data(i, indexVolume);
            vol_delta += vol_i*volumeCorrection;
        }
        const double v_a = 0.0005*M_PI*delta*delta;
        //cout <<  vol_delta / v_a << endl;
        avg += vol_delta / v_a;
    }

    avg /= particles.nParticles();
    cout << avg << endl;
    /*
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        int pId = i;
        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = pIds[id_j];
            const double radius_i = data(col_j, indexRadius);
            const double dr = con.second[indexDr0];
//            const double l_delta = delta + radius_i;
            const double rc = delta - radius_i;
            const double l_delta = delta;

            double volumeCorrection = 1.0;
            if(dr > rc)
            {
//                volumeCorrection = 0.5*(l_delta - dr)/radius_i;
//                volumeCorrection = (l_delta - dr)/rc + 0.5;
                volumeCorrection = 0.5*(l_delta + radius_i - dr)/radius_i;
            }
            con.second[indexVolumeScaling] = volumeCorrection;
        }
    }
    */
}
//------------------------------------------------------------------------------
void reCalculatePdMicromodulus(PD_Particles &particles, int dim)
{
    arma::mat & data = particles.data();
    std::unordered_map<int, int> &pIds = particles.pIds();
    int indexVolume = particles.getParamId("volume");
    int indexMicromodulus = particles.getParamId("micromodulus");
    int indexVolumeScaling = particles.getPdParamId("volumeScaling");
    int indexDr0 = particles.getPdParamId("dr0");
    double dimScaling = 2.*pow(dim, 2.);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        int pId = i;
        int col_i = i;
        double dRvolume = 0;

        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = pIds[id_j];
            double volumeScaling = con.second[indexVolumeScaling];
            double dr0Len = con.second[indexDr0];

            dRvolume += dr0Len*data(col_j, indexVolume);
            //            dRvolume += dr0Len*data(col_j, indexVolume)*volumeScaling;
        }
        data(col_i, indexMicromodulus) *= dimScaling/dRvolume;
    }
}
//------------------------------------------------------------------------------
void reCalculatePdFractureCriterion(PD_Particles &particles, double G0,
                                    double delta, double h)
{
    /*
    arma::mat & data = particles.data();
    std::unordered_map<int, int> & m_pIds = particles.pIds(); // !
    int indexMicromodulus = particles.getParamId("micromodulus");
    //    std::unordered_map<int, int> &pIds = particles.pIds();
    int indexS0 = particles.getParamId("s0");

    std::unordered_map<int, int> &pIds = particles.pIds();
    int indexVolume = particles.getParamId("volume");
    int indexVolumeScaling = particles.getPdParamId("volumeScaling");
    int indexDr0 = particles.getPdParamId("dr0");

    double delta5 = pow(delta, 5);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        int pId = i;
        int col_i = i;

        double c_i = data(col_i, indexMicromodulus);
        double s0 = sqrt(10.*G0/(c_i*M_PI*delta5));
        data(col_i, indexS0) = s0;

        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = pIds[id_j];
            const double c_j = data(col_j, indexMicromodulus);
            const double c_ij = 0.5*(c_i + c_j);
            const double dr0 = con.second[indexDr0];
            const double volumeScaling = con.second[indexVolumeScaling];
            const double V_j = data(col_j, indexVolume)*volumeScaling;
            const double s_ij = pow(6*G0/(V_j*c_ij*dr0*dr0), 1./3.);
            cout << s_ij << endl;
        }
//        data(col_i, indexS0) = s_i;
    }
    */
    arma::mat & data = particles.data();
    std::unordered_map<int, int> & m_pIds = particles.pIds(); // !
    int indexMicromodulus = particles.getParamId("micromodulus");
    //    std::unordered_map<int, int> &pIds = particles.pIds();
    int indexS0 = particles.getParamId("s0");
    double delta3 = delta*delta*delta;
    double delta4 = delta*delta*delta*delta;
    double delta5 = delta4*delta;


    std::unordered_map<int, int> &pIds = particles.pIds();
    int indexVolume = particles.getParamId("volume");
    int indexVolumeScaling = particles.getPdParamId("volumeScaling");
    int indexDr0 = particles.getPdParamId("dr0");

    if(h != -1)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            int pId = i;
            int col_i = i;
            double c = data(col_i, indexMicromodulus);
            double s0 = sqrt(4*G0/(c*h*delta4));
            data(col_i, indexS0) = s0;
        }
    }
    else
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            int pId = i;
            int col_i = i;

            double c = data(col_i, indexMicromodulus);
            double s0 = sqrt(10.*G0/(c*M_PI*delta5));;
            data(col_i, indexS0) = s0;

            double v = 0;

            vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
            for(auto &con:PDconnections)
            {
                int id_j = con.first;
                int col_j = pIds[id_j];
                double volumeScaling = con.second[indexVolumeScaling];
                v += data(col_j, indexVolume)*volumeScaling;
            }
            double l_delta = pow((3.*v/(4.*M_PI)), 1./3.);
            double d5 = pow(l_delta, 5);
            double s_i = sqrt(10.*G0/(c*M_PI*d5));
//            double s_i = sqrt(60.*G0/(c*M_PI*delta5));
            data(col_i, indexS0) = s_i;
        }
    }
}
//------------------------------------------------------------------------------
void calculateRadius(PD_Particles &particles, int dim, double h)
{
    arma::mat & data = particles.data();
    int indexVolume = particles.getParamId("volume");
    int indexRadius = particles.getParamId("radius");
    double scale = 0.9;

    if(dim == 3)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            //            int pId = i;
            int col_i = i;

            double V = data(col_i, indexVolume);
            double r = pow(3.*V/(4.*M_PI), 1./3.);
            data(col_i, indexRadius) = scale*r;
        }
    }
    else if (dim == 2)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            //            int pId = i;
            int col_i = i;

            double V = data(col_i, indexVolume);
            double r = sqrt(V/(M_PI * h));
            data(col_i, indexRadius) = scale*r;
        }
    }
    else if (dim == 1)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            int col_i = i;

            double V = data(col_i, indexVolume);
            double r = 0.5*V/(h * h);
            data(col_i, indexRadius) = scale*r;
        }
    }
    else
    {
        cerr << "Dimension not supported in radius calulations: " << dim << endl;
        throw;
    }
}
//------------------------------------------------------------------------------
void surfaceCorrection(PD_Particles &particles, vector<Force *> &forces,
                       double E, double nu, int dim)
{
#if 0
    double strain = 0.001;
    vec3 scaleFactor;
    arma::mat & r = particles.r();
    std::unordered_map<int, int> & m_pIds = particles.pIds();
    arma::mat g = zeros(particles.nParticles(), dim);

    int indexDr0 = particles.getPdParamId("dr0");
    int indexForceScaling = particles.getPdParamId("forceScalingBond");

    // Stretching all particle in the x-direction
    scaleFactor(0) = strain;
    scaleFactor(1) = 0;
    scaleFactor(2) = 0;
    double W_infty = 0;

    for(int a=0; a<dim; a++)
    {
        if(a == 1)
            scaleFactor.swap_rows(0,1);
        else if(a == 2)
            scaleFactor.swap_rows(1,2);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Loading the geometry
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<dim; d++)
            {
                r(d, col_i) = (1 + scaleFactor(d))*r(d, col_i);
            }
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Calculating the elastic energy density
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;
            double W = 0;

            for(Force *force:forces)
            {
                W += force->calculatePotentialEnergyDensity(idCol);
            }
            g(col_i, a) = W;
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Resetting the positions
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<dim; d++)
            {
                r(d, col_i) = r(d, col_i)/(1 + scaleFactor(d));
            }
        }
    }
    W_infty = 0.6*E*pow(strain, 2);

    // Scaling the energy with the median energy, which we assume
    // to be the bulk energy
    for(int a=0; a<dim; a++)
    {
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;
            double g_i = g(col_i, a);
            double W =  W_infty/g_i;
            g(col_i, a) = W;
//            cout << W << endl;
        }
    }

    // Calculating the scaling
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        int pId = idCol.first;
        int col_i = idCol.second;

        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = m_pIds[id_j];

            double dr0Len = con.second[indexDr0];
            vec3 n = (r.col(col_i) - r.col(col_j))/dr0Len;

            vec3 g_mean;
            double G = 0;
            for(int d=0; d<dim; d++)
            {
                g_mean(d) = 0.5*(g(col_i, d) + g(col_j, d));
                G += pow(n(d)/g_mean(d), 2);
            }

            G = pow(G, -0.5);
            con.second[indexForceScaling] *= G;
        }
    }
#else
    double strain = 0.001;
    vec3 scaleFactor;
    arma::mat & r = particles.r();
    std::unordered_map<int, int> & m_pIds = particles.pIds();
    arma::mat g = zeros(particles.nParticles(), dim);

    int indexDr0 = particles.getPdParamId("dr0");
    int indexForceScaling = particles.getPdParamId("forceScalingBond");

    // Stretching all particle in the x-direction
    scaleFactor(0) = strain;
    scaleFactor(1) = -nu*strain;
    scaleFactor(2) = -nu*strain;
    double W_infty = 0;

    for(int a=0; a<dim; a++)
    {
        if(a == 1)
            scaleFactor.swap_rows(0,1);
        else if(a == 2)
            scaleFactor.swap_rows(1,2);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Loading the geometry
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<dim; d++)
            {
                r(d, col_i) = (1 + scaleFactor(d))*r(d, col_i);
            }
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Calculating the elastic energy density
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;
            double W = 0;

            for(Force *force:forces)
            {
                W += force->calculatePotentialEnergyDensity(idCol);
            }
            g(col_i, a) = W;
        }
        double median_g = arma::median(g.col(a));
        W_infty += median_g;

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Resetting the positions
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<dim; d++)
            {
                r(d, col_i) = r(d, col_i)/(1 + scaleFactor(d));
            }
        }
    }
    W_infty /= (dim);
    W_infty = 0.5*E*pow(strain, 2);

    // Scaling the energy with the median energy, which we assume
    // to be the bulk energy
    for(int a=0; a<dim; a++)
    {
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;
            double g_i = g(col_i, a);
            double W =  W_infty/g_i;
            g(col_i, a) = W;
        }
    }

    // Calculating the scaling
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        int pId = idCol.first;
        int col_i = idCol.second;

        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = m_pIds[id_j];

            double dr0Len = con.second[indexDr0];
            vec3 n = (r.col(col_i) - r.col(col_j))/dr0Len;

            vec3 g_mean;
            double G = 0;
            for(int d=0; d<dim; d++)
            {
                g_mean(d) = 0.5*(g(col_i, d) + g(col_j, d));
                G += pow(n(d)/g_mean(d), 2);
            }

            G = pow(G, -0.5);
            con.second[indexForceScaling] *= G;
        }
    }
#endif
}
//------------------------------------------------------------------------------
void applyInitialStrainStrain(PD_Particles &particles, double strain, int axis, pair<double, double> area)
{
    arma::mat & r = particles.r();
    arma::mat & r0 = particles.r0();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        int col_i = idCol.second;
        if(area.first  <= r0(axis, col_i) && r0(axis, col_i) <= area.second)
        {
            double rOld = r0(axis, col_i);
            double rNew = (1 + strain)*r0(axis, col_i);
            r(axis, col_i) = rNew;
        }
    }
}

//------------------------------------------------------------------------------
void setPD_N3L(PD_Particles &particles)
{
    const int indexComputeId = particles.registerPdParameter("compute", -1);
    std::unordered_map<int, int> &pIds = particles.pIds();
    int n = 0;
    int nFound = 0;
    arma::mat & r = particles.r();
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        int id_i = i;
        vector<pair<int, vector<double>>> & PDconnections_i = particles.pdConnections(id_i);

        for(auto &con_i:PDconnections_i)
        {
            const int id_j = con_i.first;
            const int j = pIds[id_j];
            const double compute = con_i.second[indexComputeId];
#ifdef USE_N3L
            if(compute == -1)
            {
                con_i.second[indexComputeId] = 1;
                vector<pair<int, vector<double>>> & PDconnections_j = particles.pdConnections(id_j);
                bool notFound = true;
                for(auto &con_j:PDconnections_j)
                {
                    if(con_j.first == id_i)
                    {
                        con_j.second[indexComputeId] = 0;
                        notFound = false;
                        nFound++;
                    }
                }

                if(notFound)
                {
                    double dx, dy, dz;
                    cerr << "N3L connections settings not found for"
                         << id_i << " and " << id_j << endl;
                    const vec & r_i = r.col(i);
                    const vec & r_j = r.col(j);
                    dx = r_i(0) - r_j(0);
                    dy = r_i(1) - r_j(1);
                    dz = r_i(2) - r_j(2);
                    double dr = sqrt(dx*dx + dy*dy + dz*dz);
                    cerr << r.col(i) << r.col(j) << endl;
                    cerr << dr << endl;
                    n++;
                }
            }
#else
            con_i.second[indexComputeId] = 1;
#endif
        }
    }
}

//------------------------------------------------------------------------------

}
