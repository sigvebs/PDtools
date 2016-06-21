#include "boostgeometry_settings.h"

namespace PDtools
{
//------------------------------------------------------------------------------
polygon_2d createCircle(const point_2d &p, const double radius, const int nPoints)
{
    polygon_2d r;
    const double dtheta = 2.*M_PI/nPoints;

    double theta = 0;
    for(int i=0; i<nPoints; i++) {
        r.outer().push_back(point_2d(radius*cos(theta), radius*sin(theta)));
        theta += dtheta;
    }

    bg::correct(r);
    std::cout << std::boolalpha << bg::is_valid(r) << "\n";
    return r;
}

polygon_2d translatePolygon(const polygon_2d &p, const point_2d translate_point)
{
    polygon_2d p_new;

//    trans::translate_transformer<point_2d, point_2d> translate(translate_point.x(),translate_point.y());
//    trans::translate_transformer<point_type, point_type> translate(1.5, 1.5);

    trans::translate_transformer<double, 2, 2> translate(translate_point.x(),translate_point.y());
    boost::geometry::transform(p, p_new, translate);
    return p_new;
}

//------------------------------------------------------------------------------
}
