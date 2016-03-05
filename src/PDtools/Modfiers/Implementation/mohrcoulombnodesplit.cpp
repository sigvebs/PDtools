#include "mohrcoulombnodesplit.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
MohrCoulombNodeSplit::MohrCoulombNodeSplit(double mu, double C, double T, int dim):
    m_C(C), m_T(T), m_dim(dim)
{
    m_phi = mu*M_PI/180.;
    m_d = tan(m_phi);
    m_neededProperties= {pair<string, int>("stress",1)};
    m_C = 0.5*C*(1./(sqrt(m_d*m_d + 1.) + m_d));

    m_weight1 = 0.95;
    m_weight2 = 1. - m_weight1;
}
//------------------------------------------------------------------------------
void MohrCoulombNodeSplit::registerParticleParameters()
{
    m_data = &m_particles->data();
    m_indexUnbreakable =  m_particles->registerParameter("unbreakable");
    m_indexRadius =  m_particles->registerParameter("radius");
    m_indexVolume =  m_particles->registerParameter("volume");
    m_indexDr0 = m_particles->getPdParamId("dr0");
    m_indexConnected = m_particles->registerPdParameter("connected");
    m_indexCompute = m_particles->registerPdParameter("compute");
    m_indexNewConnectionId_1 = m_particles->registerParameter("newConnectionId1");
    m_indexNewConnectionId_2 = m_particles->registerParameter("newConnectionId2");
    m_indexBroken = m_particles->registerParameter("broken", 0);
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
    m_ghostParameters.push_back("newConnectionId1");
    m_ghostParameters.push_back("newConnectionId2");
}
//------------------------------------------------------------------------------
void MohrCoulombNodeSplit::initialize()
{

}
//------------------------------------------------------------------------------
void MohrCoulombNodeSplit::evaluateStepOne()
{
    const ivec &colToId = m_particles->colToId();
    int nParticles =  m_particles->nParticles();
    mat & data = m_particles->data();
    mat & r = m_particles->r();
    mat & r0 = m_particles->r0();
    mat & v = m_particles->v();

    int newCol= nParticles;
    m_toBeDeleted.clear();

    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        data(i, m_indexBroken) = 0;

        if(data(i, m_indexUnbreakable) >= 1)
            continue;
        if(data(i, m_indexBroken) >= 1)
            continue;

        if(m_dim == 2)
        {
            double sx, sy, sxy;
            sx = data(i, m_indexStress[0]);
            sy = data(i, m_indexStress[1]);
            sxy = data(i, m_indexStress[2]);

            const double first = 0.5*(sx + sy);
            const double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

            const double p_2 = first + second;
            if(p_2 >= m_T)
            {
                data(i, m_indexBroken) = 1;
            }

            if(data(i, m_indexBroken) > 0 && m_particles->pdConnections(id_i).size() <= 0)
            {
                data(i, m_indexBroken) = 0;
            }

            if(data(i, m_indexBroken) > 0)
            {
                m_toBeDeleted.push_back(id_i);

                double s_max = max(sx, sy);
                double s_min = min(sx, sy);
                double theta = 0.5*atan(2.*sxy/(s_max - s_min));
//                double theta = 0.5*atan2(2.*(sxy), (s_max - s_min));
                double normal[m_dim];
                double c = cos(theta);
                double s = sin(theta);

                if(sy < sx)
                {
                    normal[0] = c;
                    normal[1] = s;
                }
                else
                {
                    normal[0] = -s;
                    normal[1] = c;
                }
                const double radius_i = data(i, m_indexRadius);

                // Setting the new particles
                const int id1 = m_particles->newId();
                const int id2 = m_particles->newId();
                const int col1 = newCol++;
                const int col2 = newCol++;
                copyParticleTo(id1, i, col1);
                copyParticleTo(id2, i, col2);

                for(int d=0; d<m_dim; d++)
                {
                    r(col1, d) = r(i,d) + 0.5*radius_i*normal[d];
                    r(col2, d) = r(i,d) - 0.5*radius_i*normal[d];
                    r0(col1, d) = r0(i,d) + 0.5*radius_i*normal[d];
                    r0(col2, d) = r0(i,d) - 0.5*radius_i*normal[d];
                    v(col1, d) = 0.5*v(col1, d);
                    v(col2, d) = 0.5*v(col1, d);
                }

                data(col1, m_indexRadius) *= 0.5;
                data(col2, m_indexRadius) *= 0.5;
                data(col1, m_indexVolume) *= 0.5;
                data(col2, m_indexVolume) *= 0.5;
                data(i, m_indexNewConnectionId_1) = id1;
                data(i, m_indexNewConnectionId_2) = id2;

                // Setting the new connections for the two particles
                vector<pair<int, vector<double>>> & PDconnections_i = m_particles->pdConnections(id_i);
                vector<pair<int, vector<double>>>  connectionsVector1;
                vector<pair<int, vector<double>>>  connectionsVector2;

                double dr_ij[m_dim];
                for(auto &con:PDconnections_i)
                {
                    const int id_j = con.first;
                    const int j = (*m_idToCol).at(id_j);

                    double dotproduct = 0;
                    for(int d=0; d<m_dim; d++)
                    {
                        dr_ij[d] = r(i, d) - r(j, d);
                        dotproduct += normal[d]*dr_ij[d];
                    }

                    int nId, nCol;

                    if(dotproduct < 0)
                    {
                        nId = id1;
                        nCol = col1;
                    }else
                    {
                        nId = id2;
                        nCol = col2;
                    }

                    double r_len = 0;
                    for(int d=0; d<m_dim; d++)
                    {
                        dr_ij[d] = r(j, d) - r(nCol, d);
                        r_len += dr_ij[d]*dr_ij[d];
                    }
                    r_len = sqrt(r_len);

                    vector<double> newCon = con.second;
                    newCon[m_indexDr0] = r_len;

                    if(nId == id1)
                        connectionsVector1.push_back(pair<int, vector<double>>(id_j, newCon));
                    else
                        connectionsVector2.push_back(pair<int, vector<double>>(id_j, newCon));
                }

                m_particles->setPdConnections(id1, connectionsVector1);
                m_particles->setPdConnections(id2, connectionsVector2);
                data(col1, m_indexBroken) = 0;
                data(col2, m_indexBroken) = 0;
                data(i, m_indexNormal[0]) = normal[0];
                data(i, m_indexNormal[1]) = normal[1];
            }
        }
    }

    m_particles->nParticles(newCol);
    m_particles->totParticles(newCol);
    nParticles = m_particles->nParticles();

    //--------------------------------------------------------------------------
    double n_j[m_dim];
    double dr_ij[m_dim];

    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        if(data(i, m_indexBroken) >= 1)
            continue;

        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = (*m_idToCol)[id_j];

            if(data(j, m_indexUnbreakable) >= 1)
                continue;

            const int broken_j = data(j, m_indexBroken);
            if(broken_j>0)
            {

                const int id1 = data(j, m_indexNewConnectionId_1);
                const int id2 = data(j, m_indexNewConnectionId_2);

                double dotproduct = 0;
                for(int d=0; d<m_dim; d++)
                {
                    dr_ij[d] = r(j, d) - r(i, d);
                    n_j[d] = data(j, m_indexNormal[d]);
                    dotproduct += n_j[d]*dr_ij[d];
                }

                int newId;
                if(dotproduct < 0)
                {
                    newId = id1;
                }else
                {
                    newId = id2;
                }
                const int newCol = (*m_idToCol).at(newId);

                double r_len = 0;
                for(int d=0; d<m_dim; d++)
                {
                    dr_ij[d] = r(i, d) - r(newCol, d);
                    r_len += dr_ij[d]*dr_ij[d];
                }
                r_len = sqrt(r_len);
                con.first = newId;
                con.second[m_indexDr0] = r_len;

                //--------------------------------------------------------------
                // Cheking all the othder bonds
                //--------------------------------------------------------------
                const double radius_j = 1.35*data(j, m_indexRadius);
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

                double center_x = r(j, 0) + 0.5*radius_j*n_j[0];
                double center_y = r(j, 1) + 0.5*radius_j*n_j[1];

                vector<pair<int, vector<double>>> & PDconnections2 = m_particles->pdConnections(id_i);
                for(auto &con2:PDconnections2)
                {
                    if(con2.first == newId)
                        continue;
                    const int id_k = con2.first;
                    const int k = (*m_idToCol).at(id_k);
                    x4 = r(k, 0);
                    y4 = r(k, 1);

                    // If the denominator is close to zero the two lines are parallell of coincident
                    const double threshold = 1e-7;
                    const double det = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);
                    if(fabs(det) > threshold)
                        continue;

                    // Dinding the intersection between the two lines
                    const double px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4))/det;
                    const double py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4))/det;
                    const double dist = (px-center_x)*(px-center_x) + (py-center_y)*(py-center_y);

                    if(dist < r2)
                    {
                        dotproduct = (px - r(i, 0))*(px - r(k, 0));
                        dotproduct += (py - r(i, 1))*(py - r(k, 1));
                        if(dotproduct < 0)
                        {
//                            cout << "#break " << id_i << "(" << id_j << ") - " << id_k << "\t " << px << "\t " << py << "\t (" << x1 << ", " << y1 << ")" << endl;
//                            cout << "px = " << px << "\npy = " << py << endl;
//                            cout << "x1 =" << x1 << "\ny1 = " << y1 << "\nx2 = " << x2 << "\ny2 = " << y2 << "\nx3 = " << x3 << "\ny3 = " << y3 << "\nx4 = " << x4 << "\ny4 = " << y4 << endl;
                            con2.second[m_indexConnected] = 0;
                            vector<pair<int, vector<double>>> & PDconnections_k = m_particles->pdConnections(id_k);
                            for(auto &con_k:PDconnections_k)
                            {
                                if(con_k.first == id_i)
                                {
                                    con_k.second[m_indexConnected] = 0;
                                }
                            }
                        }
                    }
//                    cout << "iwdwe " << id_i << " " << con2.first  << endl;
//                    cout << "break " << id_i << "(" << id_j << ") - " << id_k << "\t " << px << "\t " << py << "\t (" << x1 << ", " << y1 << ")" << endl;
                }
                //--------------------------------------------------------------
            }
        }
    }
    nParticles = m_particles->nParticles();
    //--------------------------------------------------------------------------
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

        if(data(i, m_indexBroken) >= 1)
            continue;

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = (*m_idToCol)[id_j];

            for(int id:m_toBeDeleted)
            {
                if(id == id_j)
                {
                    cout << "id_i:" << id_i << " connected:" << id_j << " (broken)" << endl;
                }
            }
        }
    }

    nParticles = m_particles->nParticles();
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

        if(data(i, m_indexBroken) >= 1)
            continue;

        for(auto &con_i:PDconnections)
        {
            const int id_j = con_i.first;
            const int j = (*m_idToCol).at(id_j);

            vector<pair<int, vector<double>>> & PDconnections_j = m_particles->pdConnections(id_j);

            bool found = false;
            for(auto &con_j:PDconnections_j)
            {
                if(id_i == con_j.first)
                {
                    found = true;
                }
            }
            if(!found)
            {
                cout << "not found id_i:" << id_i << " in " << id_j << " connnection" << endl;
            }
        }
    }



//    if(m_toBeDeleted.size()>0)
//    {
//        cout << m_toBeDeleted.size() << " -"  << endl;
//        for(int id:m_toBeDeleted)
//        {
//            cout << id << " " ;
//        }
//        cout << endl;

//    }

    int np =  m_particles->nParticles();
    for(int id:m_toBeDeleted)
    {
        cout << "del:" << id << endl;
        m_particles->deleteParticleById(id);
        np--;
    }

    m_particles->totParticles(np);
//    int nDel = m_toBeDeleted.size();
//    if( nDel > 0)
//        cout << "Deleted" << endl;
    m_toBeDeleted.clear();
}
//------------------------------------------------------------------------------
void MohrCoulombNodeSplit::evaluateStepOne(const int id_i, const int i)
{
    /*
    mat & data = *m_data;
    const mat & r = m_particles->r();

    if(data(i, m_indexUnbreakable) >= 1)
        return;

    const int broken_i = data(i, m_indexBroken);
    double n_i[m_dim];
    double n_j[m_dim];
    double dr_ij[m_dim];

//    if(broken_i)
//    {
//        for(int d=0; d<m_dim; d++)
//        {
//            n_i[d] = data(i, m_indexNormal[d]);
//        }
//    }

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    for(auto &con:PDconnections)
    {
        const int id_j = con.first;
        const int j = (*m_idToCol)[id_j];

        if(data(j, m_indexUnbreakable) >= 1)
            continue;

//        if(con.second[m_indexConnected] <= 0.5)
//            continue;

        if(broken_i>0)
        {
            con.second[m_indexConnected] = 0;
            continue;
        }

        const int broken_j = data(j, m_indexBroken);
        if(broken_j>0)
        {
            int newId;
            int newCol;
            const int id1 = data(i, m_indexNewConnectionId_1);
            const int id2 = data(i, m_indexNewConnectionId_2);

            double dotproduct = 0;
            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = r(i, d) - r(j, d);
                n_j[d] = data(j, m_indexNormal[d]);
                dotproduct += n_j[d]*dr_ij[d];
            }

            if(dotproduct > 0)
            {
                newId = id1;
            }else
            {
                newId = id2;
            }
            newCol = (*m_idToCol)[newId];

            double r_len;
            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = r(j, d) - r(newCol, d);
                r_len = dr_ij[d]*dr_ij[d];
            }
            r_len = sqrt(r_len);
            con.first = newId;
            con.second[m_indexConnected] = r_len;
            //--------------------------
//            con.second[m_indexConnected] = 0;
        }
    }

    if(m_toBeDeleted.size()>0)
    {
        cout << m_toBeDeleted.size() << " -"  << endl;
        for(int id:m_toBeDeleted)
        {
            cout << id << " " ;
        }
        cout << endl;

    }
    for(int id:m_toBeDeleted)
    {
        cout << "del:" << id << endl;
        m_particles->deleteParticleById(id);
    }
    */
}
//------------------------------------------------------------------------------
void MohrCoulombNodeSplit::copyParticleTo(int id_to, int old_col, int new_col)
{
    ivec & colToId = m_particles->colToId();
    unordered_map<int, int>  & idToCol = m_particles->idToCol();
    mat & r = m_particles->r();
    mat & r0 = m_particles->r0();
    mat & v = m_particles->v();
    mat & F = m_particles->F();
    ivec & isStatic = m_particles->isStatic();
    mat & data = m_particles->data();
    const unordered_map<string, int> &parameters = m_particles->parameters();
    idToCol[id_to] = new_col;
    colToId[new_col] = id_to;

    for(int d=0;d<DIM;d++)
    {
        r(new_col, d) = r(old_col, d);
    }
    for(int d=0;d<DIM;d++)
    {
        r0(new_col, d) = r0(old_col, d);
    }
    for(int d=0;d<DIM;d++)
    {
        v(new_col, d) = v(old_col, d);
    }
    for(int d=0; d<DIM; d++)
    {
        F(new_col, d) = F(old_col, d);
    }
    for(auto &param:parameters)
    {
        const int p = param.second;
        data(new_col, p) = data(old_col, p);
    }
    isStatic(new_col) = isStatic(old_col);
    /*
    const int nPdConnections = recieveData[j++];
    vector<pair<int, vector<double>>> connectionsVector;

    for(int i=0;i<nPdConnections; i++)
    {
        const int con_id = (int) recieveData[j++];
        vector<double> connectionData;

        for(int k=0; k<nPdParameters; k++)
        {
            connectionData.push_back(recieveData[j++]);
        }

        connectionsVector.push_back(pair<int, vector<double>>(con_id, connectionData));
    }
    // Verlet lists
    for(int verletId=0; verletId<nVerletLists; verletId++)
    {
        vector<int> verletlist;
        int nVerletElements = recieveData[j++];
        for(int i=0; i<nVerletElements; i++)
        {
            verletlist.push_back(recieveData[j++]);
        }
        particles.setVerletList(id, verletlist, verletId);
    }
    particles.setPdConnections(id, connectionsVector);
    nParticles++;
    */
}
//------------------------------------------------------------------------------
void MohrCoulombNodeSplit::evaluateStepTwo(const int id_i, const int i)
{
    /*
    (void) id_i;
    // First calculating the total stress on an material point. The stress is averaged
    // all its bonds, and the stress state on a bond is the mean of the
    // stress at each material point.

    if((*m_data)(i, m_indexUnbreakable) >= 1)
        return;
    mat & data = *m_data;
    data(i, m_indexBroken) = 0;

    double cos_theta = cos(M_PI/2. + m_phi);
    double sin_theta = sin(M_PI/2. + m_phi);

    if(m_dim == 2)
    {
        double sx, sy, sxy;
        sx = data(i, m_indexStress[0]);
        sy = data(i, m_indexStress[1]);
        sxy = data(i, m_indexStress[2]);

        const double first = 0.5*(sx + sy);
        const double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

        const double p_2 = first + second;
        const double p_1 = first - second;

        const double shear = fabs(0.5*(p_1 - p_2)*sin_theta);
        const double normal = 0.5*(p_1 + p_2) + 0.5*(p_1 - p_2)*cos_theta;

        if(shear >= fabs(m_C - m_d*normal) && normal < 0)
        {
            data(i, m_indexBroken) = 2;
        }
        else if(p_2 >= m_T)
        {
            data(i, m_indexBroken) = 1;
        }
    }
    */
}
//------------------------------------------------------------------------------
void MohrCoulombNodeSplit::evaluateStepTwo()
{

}

//------------------------------------------------------------------------------
}
