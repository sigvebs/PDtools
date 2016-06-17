#include "meshtopdpartices.h"
#include <PDtools.h>

namespace PDtools
{
//------------------------------------------------------------------------------
PD_Particles convertMshToPdParticles(int dim, int interpolationDegree,
                                     const PdMesh &mesh, Grid &grid)
{
    //--------------------------------------------------------------------------
    // Quads
    //--------------------------------------------------------------------------
    GaussLegendreQuad gl_quad(dim, interpolationDegree);
    const vector<double> & quad_points = gl_quad.x();
    const vector<double> & quad_weights = gl_quad.w();

    const unordered_map<int, int> &idToCol_vert = mesh.getIdToCol_vertices();
    const vector<array<double, 3> > &vertices = mesh.vertices();

    // Setting the nodal shape function at the quadtrature points.
    vector<vector<double>> shapeFunction_quad;

    switch (dim) {
    case 2:
        for(int j=0; j<interpolationDegree; j++) {
            for(int k=0; k<interpolationDegree; k++) {
                const double x = quad_points[j];
                const double y = quad_points[k];
                shapeFunction_quad.push_back(gl_quad.shapeFunction(x, y));
            }
        }
    }

    // For the compuation of the derivaties in the Jakobi determinant, not needed for triangles
    vector<vector<double>> g1;
    vector<vector<double>> g2;

    for(int i=0; i<interpolationDegree; i++) {
        vector<double> a1 = {-(1-quad_points[i]),   1-quad_points[i] , 1 + quad_points[i], -(1 + quad_points[i])};
        vector<double> a2 = {-(1-quad_points[i]), -(1+quad_points[i]), 1 + quad_points[i],  (1 - quad_points[i])};
        for(int j=0; j<4; j++) {
            a1[j] *= 1./4.;
            a2[j] *= 1./4.;
        }
        g1.push_back(a1);
        g2.push_back(a2);
    }
    //--------------------------------------------------------------------------
    // Triangles
    //--------------------------------------------------------------------------
    GaussLegendreTriangle gl_tri(dim, interpolationDegree);
    const vector<double> & tri_points_x = gl_tri.x();
    const vector<double> & tri_points_y = gl_tri.y();
    const vector<double> & tri_weights = gl_tri.w();
    const int tri_elements = gl_tri.nTriPoints();
    vector<vector<double>> shapeFunction_tri;

    for(int i=0; i< tri_elements; i++) {
        const double x = tri_points_x[i];
        const double y = tri_points_y[i];
        shapeFunction_tri.push_back(gl_tri.shapeFunction(x, y));
    }
    //--------------------------------------------------------------------------
    // Initializing particles
    //--------------------------------------------------------------------------
    double h = 1.;
    if(dim = 2) {
        const vector<pair<double, double>> &bondunary = grid.boundary();
        h = bondunary[2].second - bondunary[2].first;
    }
    PD_Particles particles;
    const int nTriangles = mesh.nTriangles();
    const int nQuads = mesh.nQuads();

    const int nodesPerQuad = pow(interpolationDegree, dim);
    const int nodesPerTri = tri_elements;
    int nParticles = nodesPerQuad*nQuads + nodesPerTri*nTriangles;

    particles.maxParticles(nParticles);
    particles.nParticles(nParticles);
    particles.totParticles(nParticles);
    particles.initializeMatrices();
    particles.dim(dim);

    unordered_map<int, int> & idToCol = particles.idToCol();
    arma::ivec & get_id = particles.colToId();
    arma::mat & r0 = particles.r0();
    arma::mat & r = particles.r();
    arma::mat & v = particles.v();
    arma::mat & data = particles.data();
    const int i_volume = particles.registerParameter("volume");
    //--------------------------------------------------------------------------

    int pd_nodeId = 0;
    int pd_col = 0;

    const int myRank = grid.myRank();
    const double L0 = grid.initialPositionScaling();

    switch (dim) {
    case 2:
        for(const vector<int> &element:mesh.elements()) {
            const int elementSize = element.size();

            vector<int> verticeColums;
            verticeColums.reserve(elementSize);
            mat elementVertices(elementSize, dim);

            int i = 0;
            for(int vId:element) {
                int vCol = idToCol_vert.at(vId);
                const array<double, 3> & vert_i = vertices[vCol];

                elementVertices(i, 0) = vert_i[0];
                elementVertices(i, 1) = vert_i[1];
                i++;
            }

            // Computing the Quadrature points
            double xy[dim];

            switch (elementSize) {
                case 3: // TRIANGLE
                {
                    double A = 2.*areaTriangle(elementVertices);
                    for(int i=0; i< tri_elements; i++) {
                        const double w_i = tri_weights[i];
                        const vector<double> &N = shapeFunction_tri[i];
                        xy[0] = 0;
                        xy[1] = 0;
                        for(int j=0; j< 3; j++) {
                            xy[0] += elementVertices(j, 0)*N[j];
                            xy[1] += elementVertices(j, 1)*N[j];
                        }

                        vec3 r_local = {xy[0], xy[1], 0};
                        const int cpuId = grid.particlesBelongsTo(r_local/L0);
                        if(myRank != cpuId)
                        {
                            continue;
                        }

                        const double weight = w_i*A*h;
                        idToCol[pd_nodeId] = pd_col;
                        get_id[pd_col] = pd_col;

                        r(pd_col, 0) = xy[0];
                        r(pd_col, 1) = xy[1];
                        r(pd_col, 2) = 0;

                        r0(pd_col, 0) = xy[0];
                        r0(pd_col, 1) = xy[1];
                        r0(pd_col, 2) = 0;

                        data(pd_col, i_volume) = weight;

                        pd_nodeId++;
                        pd_col++;
                    }
                    break;
                }
                case 4: // QUAD
                {
                    for(int j=0; j<interpolationDegree; j++) {
                        double dxdj = 0;
                        double dydj = 0;
                        for(int l=0; l<4; l++)  {
                            dxdj += elementVertices(l, 0)*g1[j][l];
                            dydj += elementVertices(l, 1)*g1[j][l];
                        }

                        const double x = quad_points[j];
                        const double w_j = quad_weights[j];

                        for(int k=0; k<interpolationDegree; k++) {
                            xy[0] = 0;
                            xy[1] = 0;
                            //                        const vector<double> N = shapeFunction_quad[counter];
                            const double y = quad_points[k];
                            const double w_k = quad_weights[k];

                            const vector<double> N = gl_quad.shapeFunction(x, y);
                            double dxdk = 0;
                            double dydk = 0;

                            for(int l=0; l<4; l++) {
                                xy[0] += elementVertices(l, 0)*N[l];
                                xy[1] += elementVertices(l, 1)*N[l];
                                dxdk += elementVertices(l, 0)*g2[k][l];
                                dydk += elementVertices(l, 1)*g2[k][l];
                            }
                            vec3 r_local = {xy[0], xy[1], 0};
                            const int cpuId = grid.particlesBelongsTo(r_local/L0);
                            if(myRank != cpuId)
                            {
                                continue;
                            }

                            const double J = fabs(dxdj*dydk - dydj*dxdk);
                            const double weight = w_j*w_k*J*h;
                            idToCol[pd_nodeId] = pd_col;
                            get_id[pd_col] = pd_col;

                            r(pd_col, 0) = xy[0];
                            r(pd_col, 1) = xy[1];
                            r(pd_col, 2) = 0;

                            r0(pd_col, 0) = xy[0];
                            r0(pd_col, 1) = xy[1];
                            r0(pd_col, 2) = 0;

                            data(pd_col, i_volume) = weight;

                            pd_nodeId++;
                            pd_col++;
                        }
                    }
                    break;
                } default: {
                    cout << "hello" << endl;
                    break;
                }
            }
        }
        break;
    default:
        break;
    }

    particles.nParticles(pd_col);
    cout << pd_col << " " << nParticles << endl;
    //    particles.totParticles(pd_col);

    return particles;
}
//------------------------------------------------------------------------------
GaussLegendreQuad::GaussLegendreQuad(int _dim, int _n):
    dim(_dim), n(_n)
{
    m_shapeFunction = {0, 0, 0, 0};
    m_x.reserve(_n);
    m_w.reserve(_n);

    // For quads
    switch (_n) {
    case 1:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x1[i]);
            m_w.push_back(w1[i]);
        }
        break;
    case 2:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x2[i]);
            m_w.push_back(w2[i]);
        }
        break;
    case 3:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x3[i]);
            m_w.push_back(w3[i]);
        }
        break;
    case 4:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x4[i]);
            m_w.push_back(w4[i]);
        }
        break;
    case 5:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x5[i]);
            m_w.push_back(w5[i]);
        }
        break;
    default:
        cerr << "Degree of Gauss-Legendre polynomial not implemented: " << _n << endl;
        exit(EXIT_FAILURE);
        break;
    }

}
//------------------------------------------------------------------------------
vector<double> GaussLegendreQuad::shapeFunction(double x, double y)
{
    switch(dim) {
    case 2:
        m_shapeFunction[0] = 1./4.*(1-x)*(1-y);
        m_shapeFunction[1] = 1./4.*(1+x)*(1-y);
        m_shapeFunction[2] = 1./4.*(1+x)*(1+y);
        m_shapeFunction[3] = 1./4.*(1-x)*(1+y);
        break;
    }

    return m_shapeFunction;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreQuad::x() const
{
    return m_x;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreQuad::w() const
{
    return m_w;
}
//------------------------------------------------------------------------------
GaussLegendreTriangle::GaussLegendreTriangle(int _dim, int _n):
    dim(_dim), n(_n)
{
    m_shapeFunction = {0, 0, 0};
    m_x.reserve(_n);
    m_w.reserve(_n);


    // For triangles
    switch (_n) {
    case 1:
        m_nTriPoints = 1;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x1[i]);
            m_y.push_back(y1[i]);
            m_w.push_back(w1[i]);
        }
        break;
    case 2:
        m_nTriPoints = 3;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x2[i]);
            m_y.push_back(y2[i]);
            m_w.push_back(w2[i]);
        }
        break;
    case 3:
        m_nTriPoints = 4;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x3[i]);
            m_y.push_back(y3[i]);
            m_w.push_back(w3[i]);
        }
        break;
    case 4:
        m_nTriPoints = 6;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x4[i]);
            m_y.push_back(y4[i]);
            m_w.push_back(w4[i]);
        }
        break;
    case 5:
        m_nTriPoints = 7;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x5[i]);
            m_y.push_back(y5[i]);
            m_w.push_back(w5[i]);
        }
        break;
    default:
        cerr << "Degree of Gauss-Legendre polynomial for triangles not implemented for degree : " << _n << endl;
        exit(EXIT_FAILURE);
        break;
    }
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::shapeFunction(double x, double y)
{
    switch(dim) {
    case 2:
        m_shapeFunction[0] = 1. - x - y;
        m_shapeFunction[1] = x;
        m_shapeFunction[2] = y;
        break;
    }

    return m_shapeFunction;
}
//------------------------------------------------------------------------------
int GaussLegendreTriangle::nTriPoints() const
{
    return m_nTriPoints;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::x() const
{
    return m_x;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::y() const
{
    return m_y;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::w() const
{
    return m_w;
}
//------------------------------------------------------------------------------
double areaTriangle(const mat &vertices)
{
    const double x1 = vertices(0, 0);
    const double y1 = vertices(0, 1);
    const double x2 = vertices(1, 0);
    const double y2 = vertices(1, 1);
    const double x3 = vertices(2, 0);
    const double y3 = vertices(2, 1);

    return 0.5*fabs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
}
//------------------------------------------------------------------------------
}
