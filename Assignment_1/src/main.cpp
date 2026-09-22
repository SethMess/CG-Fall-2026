////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    double result = u(0) * v(1) - v(0) * u(1);
    return result;
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{

    // Find if second line straddles a-b
    double d1 = det(b-a, c-a);
    double d2 = det(b-a, d-a);
    bool l1_straddle = (d1 < 0) != (d2 < 0);

    double d3 = det(d-c, a-c);
    double d4 = det(d-c, b-c);
    bool l2_straddle = (d3 < 0) != (d4 < 0);
    // Find if first line straddles c-d

    return l1_straddle && l2_straddle;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon


    double max_x = poly[0](0), min_x = poly[0](0);
    double max_y = poly[0](1), min_y = poly[0](1);

    for (const auto& point : poly){
        max_x = std::max(max_x, point(0));
        min_x = std::min(min_x, point(0));
        max_y = std::max(max_y, point(1));
        min_y = std::min(min_y, point(1));
    }

    // TODO
    Vector2d outside(max_x + 1, query(1));

    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    int count = 0;
    for (size_t i = 0; i < poly.size(); ++i) {
        size_t next = (i + 1) % poly.size();
        if ( intersect_segment(query, outside, poly[i], poly[next]) ) {
            count++;

        }
    }
    

    return (count % 2) == 1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);

    int num_points;
    in >> num_points;
    double x, y, z;
    while(in >> x >> y >> z){
      Vector2d v;
      v << x, y;
      points.push_back(v);
    }

    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
  std::ofstream out(filename);

  if( ! out ){ return; }

  out << points.size() << "\n";
  for (const Vector2d &point : points) {
    out << point(0) << " " << point(1) << " " << 0 << "\n";
  }

}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    return 0;
}
