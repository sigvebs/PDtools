#include "pdfunctions.h"

#include <armadillo>
#include <unordered_map>
#include "Particles/pd_particles.h"
#include "Grid/grid.h"
#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
void setPdConnections(PD_Particles & particles,
                      Grid & grid,
                      double radius,
                      double alpha)
{
    double radiusSquared = radius*radius;
    double volumeCorrectionRadius = radius*alpha;

    const unordered_map<int, GridPoint*> &gridpoints = grid.gridpoints();
    const vector<int> &mygridPoints = grid.myGridPoints();
    const mat & R = particles.r();

    // The order is important!
    particles.registerPdParameter("dr0");
    particles.registerPdParameter("volumeScaling");

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
            int pos_i = idCol_i.second;
            const vec & r_i = R.col(pos_i);
            unordered_map<int, vector<double>> connections;
            vector<pair<int, vector<double>>> connectionsVector;

            for(const pair<int, int> & idCol_j:gridPoint.particles())
            {
                int id_j = idCol_j.first;
                if(id_i == id_j)
                    continue;

                const vec & r_j = R.col(idCol_j.second);
                dx = r_i(0) - r_j(0);
                dy = r_i(1) - r_j(1);
                dz = r_i(2) - r_j(2);

                double drSquared = dx*dx + dy*dy + dz*dz;

                if(drSquared < radiusSquared)
                {
                    vector<double> connectionData;

                    double dr = sqrt(drSquared);
                    double volumeCorrection = 1.0;

                    if(dr > radius*(1 - alpha))
                        volumeCorrection = 0.5*(radius*(1 + alpha) - dr)
                                / volumeCorrectionRadius;

                    connectionData.push_back(dr);
                    connectionData.push_back(volumeCorrection);
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
                    const vec & r_j = R.col(idCol_j.second);
                    dx = r_i(0) - r_j(0);
                    dy = r_i(1) - r_j(1);
                    dz = r_i(2) - r_j(2);

                    double drSquared = dx*dx + dy*dy + dz*dz;

                    if(drSquared < radiusSquared)
                    {
                        int id_j = idCol_j.first;
                        vector<double> connectionData;

                        double dr = sqrt(drSquared);
                        double volumeCorrection = 1.0;

                        if(dr > radius*(1 - alpha))
                            volumeCorrection = 0.5*(radius*(1 + alpha) - dr)
                                    / volumeCorrectionRadius;

                        connectionData.push_back(dr);
                        connectionData.push_back(volumeCorrection);
                        connections[id_j] = connectionData;
                        connectionsVector.push_back(pair<int, vector<double>>(id_j, connectionData) );
                    }
                }
            }

#ifdef USE_OPENMP
#pragma omp critical
            {
                particles.setPdConnections(id_i, connectionsVector);
            }
#else
            particles.setPdConnectionsVector(id_i, connectionsVector);
#endif
        }
    }
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
    arma::mat & data = particles.data();
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
//            int pId = i;
            int col_i = i;
            double c = data(col_i, indexMicromodulus);
            double s0 = sqrt(4*G0/(c*h*delta4));
//            s0 = 6*G0/(c*h*delta3);
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
            //            int pId = i;
            int col_i = i;

            double c = data(col_i, indexMicromodulus);
            double s0 = sqrt(10.*G0/(c*M_PI*delta5));;
            data(col_i, indexS0) = s0;
        }
    }
}
//------------------------------------------------------------------------------
void calculateRadius(PD_Particles &particles, int dim, double h)
{
    arma::mat & data = particles.data();
    int indexVolume = particles.getParamId("volume");
    int indexRadius = particles.getParamId("radius");
    double scale = 0.8;

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
    else
    {
        cerr << "Dimension not supported in radius calulations: " << dim << endl;
        throw;
    }
}
//------------------------------------------------------------------------------
void surfaceCorrection(PD_Particles &particles, vector<Force *> &forces,
                       double k, double nu, int dim)
{
    double strain = 0.001;
    vec3 scaleFactor;
    arma::mat & r = particles.r();
    std::unordered_map<int, int> & m_pIds = particles.pIds();
    arma::mat g = zeros(particles.nParticles(), dim);

    int indexDr0 = particles.getPdParamId("dr0");
    int indexVolumeScaling = particles.getPdParamId("volumeScaling");

    // Stretching all particle sin the x-direction
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
        W_infty += arma::median(g.col(a));

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
    W_infty /= dim;
    int a  = 2;
    W_infty = 0.5*k*pow(strain, 2)*(1 + (dim-1)*pow(nu, 2));
    int b  = 2;
    a = 3;
    b = 4;
    // Scaling the energy with the median energy, which we assume
    // to be the bulk energy
    for(int a=0; a<dim; a++)
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            double W =  W_infty/g(col_i, a);
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
            con.second[indexVolumeScaling] *= G;
        }
    }
}
//------------------------------------------------------------------------------

}
