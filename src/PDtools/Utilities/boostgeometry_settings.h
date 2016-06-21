#ifndef BOOSTGEOMETRY_SETTINGS_H
#define BOOSTGEOMETRY_SETTINGS_H

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
//#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include <boost/geometry/strategies/transform.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>

BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian)

using boost::geometry::dsv;

namespace bg = boost::geometry;
namespace trans = boost::geometry::strategy::transform;
typedef bg::model::d2::point_xy<double> point_2d;
typedef bg::model::polygon<point_2d> polygon_2d;
typedef bg::model::ring<point_2d> ring_2d;
typedef bg::model::box<point_2d> box_2d;
typedef bg::model::segment<point_2d> segment_2d;
typedef bg::model::linestring<point_2d> linestring_2d;

namespace PDtools
{
//------------------------------------------------------------------------------
polygon_2d createCircle(const point_2d &p, const double radius, const int nPoints=32);
polygon_2d translatePolygon(const polygon_2d &p, const point_2d translate_point);
//------------------------------------------------------------------------------
}
#endif // BOOSTGEOMETRY_SETTINGS_H
