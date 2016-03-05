#include "mohrcoulombmaxfracture.h"

#include "PDtools/Particles/pd_particles.h"

#include <stdlib.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

namespace PDtools
{
//------------------------------------------------------------------------------
MohrCoulombMaxFracture::MohrCoulombMaxFracture(double mu, double C, double T, int dim):
    m_C(C), m_T(T), m_dim(dim)
{
    m_phi = mu*M_PI/180.;
    m_d = tan(m_phi);
    m_neededProperties = {pair<string, int>("stress", 1), pair<string, int>("damage", 1)};
    m_C = 0.5*C*(1./(sqrt(m_d*m_d + 1.) + m_d));
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::registerParticleParameters()
{
    m_data = &m_particles->data();
    m_indexUnbreakable =  m_particles->registerParameter("unbreakable");
    m_indexConnected = m_particles->registerPdParameter("connected");
    m_indexRadius = m_particles->registerParameter("radius");
    m_indexCompute = m_particles->registerPdParameter("compute");
    m_indexBroken = m_particles->registerParameter("broken", 0);
    m_indexDamage = m_particles->registerParameter("damage");
    m_idToCol = &m_particles->idToCol();

    switch(m_dim)
    {
    case 1:
        m_ghostParameters = {"nx"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexNormal[0] = m_particles->registerParameter("nx");
        break;
    case 2:
        m_ghostParameters = {"nx", "ny"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        m_indexNormal[0] = m_particles->registerParameter("nx");
        m_indexNormal[1] = m_particles->registerParameter("ny");
        break;
    case 3:
        m_ghostParameters = {"nx", "ny", "nz"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        m_indexStress[3] = m_particles->registerParameter("s_zz");
        m_indexStress[4] = m_particles->registerParameter("s_xz");
        m_indexStress[5] = m_particles->registerParameter("s_yz");
        m_indexNormal[0] = m_particles->registerParameter("nx");
        m_indexNormal[1] = m_particles->registerParameter("ny");
        m_indexNormal[2] = m_particles->registerParameter("nz");
        break;
    }
    m_ghostParameters.push_back("broken");
    m_ghostParameters.push_back("damage");
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::initialize()
{
    srand(time(NULL));
    m_broken = false;
    m_cosTheta = cos(M_PI/2. + m_phi);
    m_sinTheta = sin(M_PI/2. + m_phi);
    m_brokenParticles.clear();
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::evaluateStepOne(const int id_i, const int i)
{
    mat & data = *m_data;
    const mat & r = m_particles->r();

    if(data(i, m_indexUnbreakable) >= 1)
        return;

    const int broken_i = data(i, m_indexBroken);
    double n_i[m_dim];
    double n_j[m_dim];
    double dr_ij[m_dim];

    if(broken_i)
    {
        for(int d=0; d<m_dim; d++)
        {
            n_i[d] = data(i, m_indexNormal[d]);
        }
    }

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    for(auto &con:PDconnections)
    {
        const int id_j = con.first;
        const int j = (*m_idToCol)[id_j];

        if(data(j, m_indexUnbreakable) >= 1)
            continue;

        if(con.second[m_indexConnected] <= 0.5)
            continue;

        if(broken_i)
        {
            double dotproduct = 0;
            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = r(j, d) - r(i, d);
                dotproduct += n_i[d]*dr_ij[d];
            }

            if(dotproduct > 0)
                con.second[m_indexConnected] = 0;

            continue;
        }

        const int broken_j = data(j, m_indexBroken);
        if(broken_j > 0)
        {
            double dotproduct = 0;
            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = r(i, d) - r(j, d);
                n_j[d] = data(j, m_indexNormal[d]);
                dotproduct += n_j[d]*dr_ij[d];
            }

            if(dotproduct > 0)
                con.second[m_indexConnected] = 0;

            //--------------------------------------------------------------
            // Cheking all the othder bonds crossing the fracture
            //--------------------------------------------------------------
//            const double radius_j = 1.01*data(j, m_indexRadius);
//            const double radius_j = 1.1*data(j, m_indexRadius);
            const double radius_j = 1.1*data(j, m_indexRadius);
            const double r2 = radius_j*radius_j;
            // Normal in fracture direction
            n_j[0] = -data(j, m_indexNormal[1]);
            n_j[1] = data(j, m_indexNormal[0]);
            double x1, x2, x3, x4;
            double y1, y2, y3, y4;
            x1 = r(j, 0);
            y1 = r(j, 1);
            x2 = r(j, 0) + radius_j*n_j[0];
            y2 = r(j, 1) + radius_j*n_j[1];
            x3 = r(i, 0);
            y3 = r(i, 1);

            vector<pair<int, vector<double>>> & PDconnections2 = m_particles->pdConnections(id_i);
            for(auto &con2:PDconnections2)
            {
                if(con2.first == id_j)
                {
                    continue;
                }
                const int id_k = con2.first;
                const int k = (*m_idToCol).at(id_k);

                if(data(k, m_indexUnbreakable) >= 1)
                    continue;

                if(con2.second[m_indexConnected] <= 0.5)
                    continue;

                x4 = r(k, 0);
                y4 = r(k, 1);

                // If the denominator is close to zero the two lines are parallell of coincident
                const double threshold = 1e-7;
                const double det = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);
                if(fabs(det) > threshold)
                    continue;

                // Finding the intersection between the two lines
                const double px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4))/det;
                const double py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4))/det;
                const double dist = (px - x1)*(px - x1) + (py - y1)*(py - y1);

                if(dist < r2)
                {
                    double dotproduct = (px - r(i, 0))*(px - r(k, 0));
                    dotproduct += (py - r(i, 1))*(py - r(k, 1));
                    if(dotproduct < 0)
                    {
//                        cout << "#break " << id_i << "(" << id_j << ") - " << id_k << "\t " << px << "\t " << py << "\t (" << x1 << ", " << y1 << ")" << endl;
//                        cout << "px = " << px << "\npy = " << py << endl;
//                        cout << "x1 =" << x1 << "\ny1 = " << y1 << "\nx2 = " << x2 << "\ny2 = " << y2 << "\nx3 = " << x3 << "\ny3 = " << y3 << "\nx4 = " << x4 << "\ny4 = " << y4 << endl;
                        con2.second[m_indexConnected] = 0;
                        m_brokenParticles[id_k].push_back(id_i);
                    }
                }
            }
            /*
            continue;
            if(broken_j != 2)
                continue;

            //--------------------------------------------------------------
            // If shear failure we check both sides
            //--------------------------------------------------------------
            n_j[0] = data(j, m_indexNormal[0]);
            n_j[1] = data(j, m_indexNormal[1]);

            x1 = r(j, 0);
            y1 = r(j, 1);
            x2 = r(j, 0) + radius_j*n_j[0];
            y2 = r(j, 1) + radius_j*n_j[1];
            x3 = r(i, 0);
            y3 = r(i, 1);

            center_x = r(j, 0) + 0.5*radius_j*n_j[0];
            center_y = r(j, 1) + 0.5*radius_j*n_j[1];

//            vector<pair<int, vector<double>>> & PDconnections2 = m_particles->pdConnections(id_i);
            for(auto &con2:PDconnections2)
            {
                if(con2.first == id_j)
                {
                    continue;
                }
                const int id_k = con2.first;
                const int k = (*m_idToCol).at(id_k);
                x4 = r(k, 0);
                y4 = r(k, 1);

                // If the denominator is close to zero the two lines are parallell of coincident
                const double threshold = 1e-7;
                const double det = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);
                if(fabs(det) > threshold)
                    continue;

                // Finding the intersection between the two lines
                const double px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4))/det;
                const double py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4))/det;
                const double dist = (px-center_x)*(px-center_x) + (py-center_y)*(py-center_y);

                if(dist < r2)
                {
                    double dotproduct = (px - r(i, 0))*(px - r(k, 0));
                    dotproduct += (py - r(i, 1))*(py - r(k, 1));
                    if(dotproduct < 0)
                    {
//                        cout << "#break " << id_i << "(" << id_j << ") - " << id_k << "\t " << px << "\t " << py << "\t (" << x1 << ", " << y1 << ")" << endl;
//                        cout << "px = " << px << "\npy = " << py << endl;
//                        cout << "x1 =" << x1 << "\ny1 = " << y1 << "\nx2 = " << x2 << "\ny2 = " << y2 << "\nx3 = " << x3 << "\ny3 = " << y3 << "\nx4 = " << x4 << "\ny4 = " << y4 << endl;
                        con2.second[m_indexConnected] = 0;
                        m_brokenParticles[id_k].push_back(id_i);
                    }
                }
            }
            */
        }
    }
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::evaluateStepOnePost()
{
#if USE_MPI
    exchangeBrokenParticlesMPI();
#endif
    const ivec &colToId = m_particles->colToId();
    std::unordered_map<int, int> & idToCol = m_particles->idToCol();
    const int nParticles = m_particles->nParticles();
#if USE_MPI
    const int myRank = MPI::COMM_WORLD.Get_rank( );
#endif
    for(auto broken_idList:m_brokenParticles)
    {
        const int id_i = broken_idList.first;
        const vector<int> & broken = broken_idList.second;

#if USE_MPI
        if(!idToCol.count(id_i))
        {
            continue;
        }
        // Check if the particle is a ghost particle
        if(idToCol.at(id_i) >= nParticles)
        {
            continue;
        }
#endif
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

        for(auto &con:PDconnections)
        {
            for(int id:broken)
            {
                if(id == id_i)
                {
                    con.second[m_indexConnected] = 0;
                }
            }
        }
    }
    m_brokenParticles.clear();
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::evaluateStepTwo(const int id_i, const int i)
{
    if((*m_data)(i, m_indexUnbreakable) >= 1)
        return;

    mat & data = *m_data;
    mat & r = m_particles->r();

    if(m_dim == 2)
    {
        double sx, sy, sxy;
        sx = data(i, m_indexStress[0]);
        sy = data(i, m_indexStress[1]);
        sxy = data(i, m_indexStress[2]);

        double first = 0.5*(sx + sy);
        double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

        double p_2 = first + second; // max
        double p_1 = first - second; // min

        const double shear = 0.5*(p_1 - p_2)*m_sinTheta;
        const double normal = 0.5*(p_1 + p_2) + 0.5*(p_1 - p_2)*m_cosTheta;

        double s_max = max(sx, sy);
        double s_min = min(sx, sy);
        double theta = 0.5*atan(2.*sxy/(s_min - s_max));
        int broken = 0;
        const double criticalShear = fabs(shear) - fabs(m_C - m_d*normal);
        const double criticalTensile = p_2 - m_T;

//        if(fabs(shear) >= fabs(m_C - m_d*normal) && normal < 0)
//        if(fabs(shear) >= fabs(m_C - m_d*normal))
        if(criticalShear >= 0  && normal < 0)
        {
            data(i, m_indexBroken) = 2;
            m_broken = true;
            broken = 2;
            theta = -theta; // Changing the rotation to the front
        }
        //        else if(p_2 >= m_T)
        else if(criticalTensile >= 0)
        {
            data(i, m_indexBroken) = 1;
            m_broken = true;
            broken = 1;
        }
        else
        {
            data(i, m_indexBroken) = 0;
        }

        if(broken > 0)
        {

//            if(broken == 2)
//            {
//                const double shearAngle = 0.5*m_phi + 0.25*M_PI;
//                theta = shearAngle - theta;
////                if(theta > 0)
////                    theta = shearAngle - theta;
////                else
////                    theta = -shearAngle - theta;
//            }

////            if(sx < sy)
////            {
////                theta += 0.5*M_PI;
////                cout << "HAHAHAHAHAHAH" << endl;
////            }


            vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
            //------------------------------------------------------------------
            // Choosing the least damaged side for tensile fracture
            if(broken == 1)
            {
                if(sy > sx)
                {
                    theta += 0.5*M_PI;
//                    theta -= 0.5*M_PI;
                }

                double damageLeft = 0;
                double damageRight = 0;
                double n[m_dim];
                double dr_ij[m_dim];
//                n[0] = -sin(theta);
//                n[1] = cos(theta);
                n[0] = cos(theta);
                n[1] = sin(theta);

                for(auto &con:PDconnections)
                {
                    const int id_j = con.first;
                    const int j = (*m_idToCol)[id_j];

                    double dotproduct = 0;
                    for(int d=0; d<m_dim; d++)
                    {
                        dr_ij[d] = r(j, d) - r(i, d);
                        dotproduct += n[d]*dr_ij[d];
                    }
                    if(dotproduct > 0)
                        damageRight += data(j, m_indexDamage);
                    else
                        damageLeft += data(j, m_indexDamage);
                }
                if(damageLeft > damageRight)
                    theta += M_PI;
            }
            int i1, j1;
            double theta0 = -500;
            double t1, t2;
            double damageDiag_1 = 0;
            double damageDiag_2 = 0;
            //------------------------------------------------------------------
            // Choosing the shear fracture in the direction of most failure
            //------------------------------------------------------------------
            if(broken == 2)
            {
                if(sy < sx)
                {
                    theta += 0.5*M_PI;
                }

                double np1[m_dim];
                double np2[m_dim];
                theta0 = theta;

                // Principal axis
                np1[0] = cos(theta0);
                np1[1] = sin(theta0);
                np2[0] = -np1[1];
                np2[1] = np1[0];

                double positive_np1 = 0;
                double negative_np1 = 0;
                double positive_np2 = 0;
                double negative_np2 = 0;
                double n[m_dim];
                double dr_ij[m_dim];

                for(auto &con:PDconnections)
                {
                    const int id_j = con.first;
                    const int j = (*m_idToCol)[id_j];

                    double dotproduct_np1 = 0;
                    double dotproduct_np2 = 0;
                    for(int d=0; d<m_dim; d++)
                    {
                        dr_ij[d] = r(j, d) - r(i, d);
                        dotproduct_np1 += np1[d]*dr_ij[d];
                        dotproduct_np2 += np2[d]*dr_ij[d];
                    }
                    if(dotproduct_np1 > 0)
                        positive_np1 += data(j, m_indexDamage);
                    else
                        negative_np1 += data(j, m_indexDamage);

                    if(dotproduct_np2 > 0)
                        positive_np2 += data(j, m_indexDamage);
                    else
                        negative_np2 += data(j, m_indexDamage);

                    if((dotproduct_np1 > 0 && dotproduct_np2 > 0)
                            || (dotproduct_np1 < 0 && dotproduct_np2 < 0))
                        damageDiag_1 += data(j, m_indexDamage);
                    else
                        damageDiag_2 += data(j, m_indexDamage);
                }

//                m_cosTheta = cos(M_PI/2. + m_phi);
                const double shearAngle = 0.5*m_phi + 0.25*M_PI;


//                if((positive_np1 > negative_np1 && positive_np2 > negative_np2)
//                        || (positive_np1 < negative_np1 && positive_np2 < negative_np2))
//                {
//                    theta = - shearAngle + theta;
//                }else
//                {
//                    theta = shearAngle - theta;
//                }


                //                if(positive_np1 > negative_np1)
                //                    i1 = 1;
                //                else if(positive_np1 == negative_np1)
                //                    i1 = rand() % 2;
                //                else
                //                    i1 = 0;

                //                if(positive_np2 > negative_np2)
                //                    j1 = 1;
                //                else if(positive_np2 == negative_np2)
                //                    j1 = rand() % 2;
                //                else
                //                    j1 = 0;

                //                t1 = -shearAngle - theta0;
                //                t2 =  shearAngle - theta0;

                //                if(i1 == j1)
                //                    theta = -shearAngle - theta0;
                //                else
                //                    theta = shearAngle - theta0;

                if(damageDiag_1 == damageDiag_2)
                {
                    if(rand() % 2)
                        theta = -shearAngle - theta0;
                    else
                        theta = shearAngle - theta0;
                }
                else if(damageDiag_1 > damageDiag_2)
                    theta = -shearAngle - theta0;
                else
                    theta = shearAngle - theta0;

                t1 = -shearAngle - theta0;
                t2 =  shearAngle - theta0;
            }
            //------------------------------------------------------------------
//            data(i, m_indexNormal[0]) = -sin(theta);
//            data(i, m_indexNormal[1]) = cos(theta);
            data(i, m_indexNormal[0]) = cos(theta);
            data(i, m_indexNormal[1]) = sin(theta);

            double radius = data(i, m_indexRadius);

            double fx = -sin(t1);
            double fy = cos(t1);
            double nx = -fy;
            double ny = fx;
            double fx_ = -sin(t2);
            double fy_ = cos(t2);
            double nx_ = -fy_;
            double ny_ = fx_;
            cout << "# b: " << id_i
                 << "\t t1:" << t1*180./M_PI << " t2:" << t2*180./M_PI  << " tp:" << theta*180./M_PI
                 << "\ti:" << damageDiag_1 << " j:" << damageDiag_2
                 << "\t b:" << broken << endl;
            cout << "fx[0] = "     << r(i, 0) - radius*fx
                 << "\nfy[0] = "   << r(i, 1) - radius*fy
                 << "\nfx[1] = "   << r(i, 0) + radius*fx
                 << "\nfy[1] = "   << r(i, 1) + radius*fy
                 << "\nnx[0] = "  << r(i, 0) - radius*nx
                 << "\nny[0] = "  << r(i, 1) - radius*ny
                 << "\nnx[1] = "  << r(i, 0) + radius*nx
                 << "\nny[1] = "  << r(i, 1) + radius*ny
                 << "\nfx_[0] = "  << r(i, 0) - radius*fx_
                 << "\nfy_[0] = "  << r(i, 1) - radius*fy_
                 << "\nfx_[1] = "  << r(i, 0) + radius*fx_
                 << "\nfy_[1] = "  << r(i, 1) + radius*fy_
                 << "\nnx_[0] = " << r(i, 0) - radius*nx_
                 << "\nny_[0] = " << r(i, 1) - radius*ny_
                 << "\nnx_[1] = " << r(i, 0) + radius*nx_
                 << "\nny_[1] = " << r(i, 1) + radius*ny_
                 << endl;
//            cout << "# b: " << id_i << "\t theta:" << theta*180./M_PI << " theta:" << theta0*180./M_PI <<  "\t x:" << nx << "\t y:" << ny << "\ti:" << i1 << " j:" << j1 << "\t b:" << broken << endl;
//            cout << "# sx:" << sx << "\tsy:" << sy << "\tsxy:" << sxy << endl;
        }
    }
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::evaluateStepTwo()
{
    if(m_broken)
    {
        m_state = true;
    }
    else
    {
        m_state = false;
    }

    m_broken = false;

}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::exchangeBrokenParticlesMPI()
{
#if USE_MPI
    vector<int> sendData;
    const vector<int> &neighbouringCores = m_grid->neighbouringCores();

    // Collecting all the send data
    for(auto broken_idList:m_brokenParticles)
    {
        const int id_i = broken_idList.first;
        const vector<int> & broken = broken_idList.second;

        sendData.push_back(id_i);
        sendData.push_back(broken.size());

        for(int b:broken)
            sendData.push_back(b);
    }
    const int myRank = MPI::COMM_WORLD.Get_rank( );

    for(int toCore:neighbouringCores)
    {
        int nSendElements = sendData.size();
        int nRecieveElements;
        MPI_Status status;

        MPI_Sendrecv(&nSendElements, 1, MPI_INT,
                     toCore, myRank,
                     &nRecieveElements, 1,MPI_INT,
                     toCore, toCore,
                     MPI_COMM_WORLD, &status);

        int recieveData[nRecieveElements];

        MPI_Sendrecv(&sendData[0], nSendElements, MPI_INT,
                toCore, 1,
                &recieveData, nRecieveElements, MPI_INT,
                toCore, 1,
                MPI_COMM_WORLD, &status);

        int i = 0;
        while(i<nRecieveElements)
        {
            const int id = recieveData[i++];
            int nBrokenBonds = recieveData[i++];
            for(int j=0;j<nBrokenBonds; j++)
            {
                m_brokenParticles[id].push_back(recieveData[i++]);
            }
        }
    }
#endif
}
//------------------------------------------------------------------------------
}
