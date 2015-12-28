#include "pdfunctions.h"

#include <armadillo>
#include <unordered_map>
#include "Particles/pd_particles.h"
#include "Grid/grid.h"
#include "PDtools/Force/force.h"

#define USE_EXTENDED_RANGE_RADIUS 0
#define USE_EXTENDED_RANGE_LC 1

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
//    const int indexRadius = particles.getParamId("radius");
//    const mat &data  = particles.data();

    // The order is important!
    particles.registerPdParameter("dr0");
    particles.registerPdParameter("connected");

    int nFour = 0;
    int nElse = 0;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(unsigned int i=0; i<mygridPoints.size(); i++)
    {
        double dx, dy, dz;
        int gridId = mygridPoints.at(i);
        const GridPoint & gridPoint = *gridpoints.at(gridId);

        for(const pair<int, int> & idCol_i:gridPoint.particles())
        {
            int id_i = idCol_i.first;
            int col_i = idCol_i.second;
            const vec & r_i = R.row(col_i).t();
            unordered_map<int, vector<double>> connections;
            vector<pair<int, vector<double>>> connectionsVector;

            for(const pair<int, int> & idCol_j:gridPoint.particles())
            {
                const int id_j = idCol_j.first;
                const int col_j = idCol_j.second;
                if(id_i == id_j)
                    continue;

#if USE_EXTENDED_RANGE_RADIUS
                const double radius_i = data(col_i, indexRadius);
                const double radius_j = data(col_j, indexRadius);
                const double l_delta = delta + 0.5*(radius_i + radius_j);
#elif USE_EXTENDED_RANGE_LC
                const double l_delta = delta + 0.5*lc;
#else
                const double l_delta = delta;
#endif
                const double deltaSquared = pow(l_delta, 2);

                const vec & r_j = R.row(col_j).t();
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

#if USE_EXTENDED_RANGE_RADIUS
                const double radius_j = data(col_j, indexRadius);
                const double l_delta = delta + 0.5*(radius_i + radius_j);
#elif USE_EXTENDED_RANGE_LC
                const double l_delta = delta + 0.5*lc;
#else
                const double l_delta = delta;
#endif
                    const double deltaSquared = pow(l_delta, 2);

                    const vec & r_j = R.row(col_j).t();
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
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        int pId = i;
        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
        double vol_delta = 0;
        for(auto &con:PDconnections)
        {
            const double dr = con.second[indexDr0];
            const int id_j = con.first;
            const int col_j = pIds[id_j];
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
            const double vol_i = data(i, indexVolume);
            vol_delta += vol_i*volumeCorrection;
        }
    }
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
//    int indexVolumeScaling = particles.getPdParamId("volumeScaling");
    int indexDr0 = particles.getPdParamId("dr0");
    double dimScaling = 2.*pow(dim, 2.);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        int pId = i;
        int col_i = i;
        double dRvolume = 0;

        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = pIds[id_j];
//            double volumeScaling = con.second[indexVolumeScaling];
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
//    std::unordered_map<int, int> & m_pIds = particles.pIds(); // !
    const int indexMicromodulus = particles.getParamId("micromodulus");
    //    std::unordered_map<int, int> &pIds = particles.pIds();
    int indexS0 = particles.getParamId("s0");
    const double delta4 = delta*delta*delta*delta;
    const double delta5 = delta4*delta;


    std::unordered_map<int, int> &pIds = particles.pIds();
    int indexVolume = particles.getParamId("volume");
    int indexVolumeScaling = particles.getPdParamId("volumeScaling");

    if(h != -1)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
//            const int pId = i;
            const int col_i = i;
            const double c = data(col_i, indexMicromodulus);
            const double s0 = sqrt(4*G0/(c*h*delta4));
            data(col_i, indexS0) = s0;
        }
    }
    else
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
            const int pId = i;
            const int col_i = i;

            const double c = data(col_i, indexMicromodulus);
            const double s0 = sqrt(10.*G0/(c*M_PI*delta5));;
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
    const int indexVolume = particles.getParamId("volume");
    const int indexRadius = particles.getParamId("radius");
    const double scale = 0.9;

    if(dim == 3)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
            //            int pId = i;
            const int col_i = i;

            const double V = data(col_i, indexVolume);
            const double r = pow(3.*V/(4.*M_PI), 1./3.);
            data(col_i, indexRadius) = scale*r;
        }
    }
    else if (dim == 2)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
            //            int pId = i;
            const int col_i = i;

            const double V = data(col_i, indexVolume);
            const double r = sqrt(V/(M_PI * h));
            data(col_i, indexRadius) = scale*r;
        }
    }
    else if (dim == 1)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
            const int col_i = i;
            const double V = data(col_i, indexVolume);
            const double r = 0.5*V/(h * h);
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
                r(col_i, d) = (1 + scaleFactor(d))*r(col_i, d);
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
                r(col_i, d) = r(col_i, d)/(1 + scaleFactor(d));
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
            vec3 n = (r.row(col_i).t() - r.row(col_j).t())/dr0Len;

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
    const double strain = 0.001;
    vec3 scaleFactor;
    arma::mat & r = particles.r();
    std::unordered_map<int, int> & m_pIds = particles.pIds();
    arma::mat g = zeros(particles.nParticles(), dim);

    const int indexDr0 = particles.getPdParamId("dr0");
    const int indexForceScaling = particles.getPdParamId("forceScalingBond");

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
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;

            for(int d=0; d<dim; d++)
            {
                r(col_i, d) = (1 + scaleFactor(d))*r(col_i, d);
            }
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Calculating the elastic energy density
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;
            double W = 0;

            for(Force *force:forces)
            {
                W += force->calculatePotentialEnergyDensity(idCol);
            }
            g(col_i, a) = W;
        }
        double median_g = arma::median(g.row(a));
        W_infty += median_g;

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Resetting the positions
        for(unsigned int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<dim; d++)
            {
                r(col_i, d) = r(col_i, d)/(1 + scaleFactor(d));
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
        for(unsigned int i=0; i<particles.nParticles(); i++)
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
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        const int pId = idCol.first;
        const int col_i = idCol.second;

        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = m_pIds[id_j];

            const double dr0Len = con.second[indexDr0];
            vec3 n = (r.row(col_i).t() - r.row(col_j).t())/dr0Len;

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
    const arma::mat & r0 = particles.r0();
    int counter = 0;
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        const int col_i = idCol.second;
        if(area.first  <= r0(col_i, axis) && r0(col_i, axis) <= area.second)
        {
            const double rNew = (1 + strain)*r0(col_i, axis);
            r(col_i, axis) = rNew;
            counter++;
        }
    }
}
//------------------------------------------------------------------------------
void setPD_N3L(PD_Particles &particles)
{
    const int indexComputeId = particles.registerPdParameter("compute", -1);
#ifdef USE_N3L
    std::unordered_map<int, int> &pIds = particles.pIds();
    int n = 0;
    int nFound = 0;
    arma::mat & r = particles.r();
#endif
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        int id_i = i;
        vector<pair<int, vector<double>>> & PDconnections_i = particles.pdConnections(id_i);

        for(auto &con_i:PDconnections_i)
        {
#ifdef USE_N3L
            const int id_j = con_i.first;
            const double compute = con_i.second[indexComputeId];
            const int j = pIds[id_j];
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
                    const vec & r_i = r.row(i).t();
                    const vec & r_j = r.row(j).t();
                    dx = r_i(0) - r_j(0);
                    dy = r_i(1) - r_j(1);
                    dz = r_i(2) - r_j(2);
                    double dr = sqrt(dx*dx + dy*dy + dz*dz);
                    cerr << r.row.t()(i) << r.row.t()(j) << endl;
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
void removeVoidConnections(PD_Particles &particles, Grid &grid,
                           const double delta, const double lc)
{
    (void) grid;
    (void) delta;
    (void) lc;

    arma::mat & data = particles.data();
    std::unordered_map<int, int> &pIds = particles.pIds();
    const mat & R0 = particles.r0();
    const int indexDr0 = particles.getPdParamId("dr0");
    const int indexRadius = particles.getParamId("radius");
    const int dim = 2;
    const int m = 15;
    const int indexConnected = particles.getPdParamId("connected");
    const double sf = 1.5;

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        const int pId = i;
        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
        const vec & r_a = R0.row(i).t();
        ivec filled(m);
        arma::mat filledCenters = arma::zeros(DIM, m);
        const double radius_i = sf*data(i, indexRadius);

        for(auto &con:PDconnections)
        {
            const int id_g = con.first;
            const int g = pIds[id_g];
            const double dr0 = con.second[indexDr0];
            const vec & r_g = R0.row(g).t();
            const vec & n = (r_g - r_a)/dr0;
            const double radius_g = sf*data(g, indexRadius);
            const double radius_ig = radius_i + radius_g;

            if(radius_ig >= dr0)
                continue;

            const double spacing = (dr0 - radius_ig)/(m + 1);

            for(int k=0; k<m; k++)
            {
                filled(k) = 0;
                filledCenters.col(k) = r_a + ((k+1)*spacing + radius_i)*n;
            }

            int nFilled = 0;

            for(auto &con_b:PDconnections)
            {
                const int id_b = con_b.first;
                const int b = pIds[id_b];
                if(g == b)
                    continue;

                const vec & r_b = R0.row(b).t();
                const double radius = sf*data(b, indexRadius);

                // Distance to each point
                for(int k=0; k<m; k++)
                {
                    const arma::vec& diff = filledCenters.col(k) - r_b;
                    double len_sq = 0;
                    for(int d=0;d<dim; d++)
                    {
                        len_sq += diff(d)*diff(d);
                    }

                    if(len_sq <= radius*radius)
                    {
                        filled(k) = 1;
                        nFilled++;
                    }
                }

                if(nFilled == m)
                    continue;
            }
            bool remove = false;

            for(int l=0; l<m; l++)
            {
                if(filled(l) == 0)
                {
                    remove = true;
                }
            }

            if(remove)
            {
                con.second[indexConnected] = 0;
            }
        }
    }

    // Enforcing symmetry
    int nFound = 0;

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        int id_i = i;
        vector<pair<int, vector<double>>> & PDconnections_i = particles.pdConnections(id_i);

        for(auto &con_j:PDconnections_i)
        {
            const int id_j = con_j.first;
            const bool connected = con_j.second[indexConnected];

            if(!connected)
            {
                vector<pair<int, vector<double>>> & PDconnections_j = particles.pdConnections(id_j);
                for(auto &con_k:PDconnections_j)
                {
                    if(con_k.first == id_i)
                    {
                        con_k.second[indexConnected] = 0;
                        nFound++;
                    }
                }
            }
        }
    }

    cout << "Void bonds removed due to symmetry = " << nFound << endl;
}
//------------------------------------------------------------------------------
void cleanUpPdConnections(PD_Particles &particles)
{
    const int indexConnected = particles.getPdParamId("connected");

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        const int id_i = i;
        vector<pair<int, vector<double>>> & PDconnections_i = particles.pdConnections(id_i);
        vector<pair<int, vector<double>>> removeConnections;

        for(auto &con_j:PDconnections_i)
        {
            const bool connected = con_j.second[indexConnected];
            if(!connected)
            {
                removeConnections.push_back(con_j);
            }
        }

        for(auto &con_j:removeConnections)
        {
            PDconnections_i.erase(std::remove(PDconnections_i.begin(), PDconnections_i.end(), con_j), PDconnections_i.end());
        }
    }
}
//------------------------------------------------------------------------------
void addFractures(PD_Particles &particles, const vector<pair<double, double> > &domain)
{
    std::unordered_map<int, int> &pIds = particles.pIds();
    const mat & R0 = particles.r0();
    const int indexConnected = particles.getPdParamId("connected");

    const double xFrac = 0.5*(domain[0].second + domain[0].first);
    const double yFrac = 0.5*(domain[1].second + domain[1].first);
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        const int pId = i;
        vector<pair<int, vector<double>>> & PDconnections = particles.pdConnections(pId);
        const vec & r_a = R0.row(i).t();

        for(auto &con:PDconnections)
        {
            const int id_b = con.first;
            const int b = pIds[id_b];
            const vec & r_b = R0.row(b).t();
            bool remove = false;

            if(r_a(1) > yFrac && r_b(1) < yFrac && r_a(0) < xFrac && r_b(0) < xFrac)
                remove = true;

            if(r_a(1) < yFrac && r_b(1) > yFrac && r_a(0) < xFrac && r_b(0) < xFrac)
                remove = true;

            if(remove)
            {
                con.second[indexConnected] = 0;
            }
        }
    }
}
//------------------------------------------------------------------------------
}
