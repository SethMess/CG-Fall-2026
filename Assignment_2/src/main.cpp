// C++ include
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;


// This function computes if a ray hits a sphere and return a boolean based on the outcome and t for the value if true.
bool ray_intersect_sphere(const Vector3d &e, const Vector3d &d, const Vector3d &center, double radius, double &t) {
    double a = d.dot(d);
    Vector3d center_to_origin =  e - center;
    double b = 2 * d.dot(center_to_origin);
    double c = center_to_origin.dot(center_to_origin) - (radius * radius);

    double delta = (b * b) - (4 * a * c);

    // discriminant ≥ 0 and a valid t exists
    if (delta >= 0)
    {
        double t1 = (-b - std::sqrt(delta)) / (2 * a);
        double t2 = (-b + std::sqrt(delta)) / (2 * a);

        if (t1 >= 0) {
            t = t1;
            return true;
        } else if (t2 >= 0){
            t = t2;
            return true;
        } else {
            return false;
        }
    }
    return false;
}



void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    const double sphere_radius = 0.9;
    const Vector3d sphere_center(0, 0, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this was a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis

            double a = ray_direction.dot(ray_direction);
            Vector3d center_to_origin =  ray_origin - sphere_center;
            double b = 2 * ray_direction.dot(center_to_origin);
            double c = center_to_origin.dot(center_to_origin) - (sphere_radius * sphere_radius);

            double delta = (b * b) - (4 * a * c);

            // discriminant ≥ 0 and a valid t exists
            if (delta >= 0)
            {
                double t1 = (-b - std::sqrt(delta)) / (2 * a);
                double t2 = (-b + std::sqrt(delta)) / (2 * a);

                double t;
                if (t1 >= 0) {
                    t = t1;
                } else if (t2 >= 0){
                    t = t2;
                } else {
                    continue;
                }
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection = ray_origin + t * ray_direction;

                // Compute normal at the intersection point
                Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
           }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

// This function computes if a ray hits a parallelogram and return a boolean based on the outcome and t for the value if true.
bool ray_intersect_parallelogram(const Vector3d &e, const Vector3d &d, const Vector3d &o, const Vector3d &u, const Vector3d &v, double &t) {
    Matrix3d M;
    M << u, v, -d;
    Vector3d rhs = e - o;

    // check if ray is parallel to the plane
    if (std::abs(M.determinant()) < 1e-9) {
        return false;
    }

    Vector3d solved = M.colPivHouseholderQr().solve(rhs);

    double a = solved(0);
    double b = solved(1);
    t = solved(2);

    if (t >= 0 && (a >= 0 && a <= 1) && (b >= 0 && b <= 1)) {
        return true;
    }

    return false;
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(1, 0.4, 0);
    const Vector3d pgram_v(0, 0.7, -10);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            double t = 0;
            if (ray_intersect_parallelogram(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v, t))
            {
                // The ray hit the parallelogram, compute the exact intersection
                // point
                Vector3d ray_intersection = ray_origin + (t * ray_direction);

                // Compute normal at the intersection point
                // Vector3d ray_normal = ray_intersection.normalized();
                Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(1, 0.4, 0);
    const Vector3d pgram_v(0, 0.7, -10);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            double t = 0;
            if (ray_intersect_parallelogram(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v, t))
            {
                // The ray hit the parallelogram, compute the exact intersection point
                Vector3d ray_intersection = ray_origin + (t * ray_direction);

                // Compute normal at the intersection point
                Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{

  /*
   * Notes about changing the parameters:
   *
   * the specular exponent changes the size of the reflected circle.
   *
   * changing the diffuse color changes the color of the sphere.
   * while changing the specular color changes the color of the highlight.
   *
   * ambient light changes how dark or light the part of the sphere with no light shining on it is.
   */

    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask
    // Store the R G B components of the image
    MatrixXd R = MatrixXd::Zero(800, 800);
    MatrixXd G = MatrixXd::Zero(800, 800);
    MatrixXd B = MatrixXd::Zero(800, 800);

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    const Vector3d light_intesity(1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < A.cols(); ++i)
    {
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // Intersect with the sphere
            // implement the generic ray sphere intersection

            double t;

            if (ray_intersect_sphere(ray_origin, ray_direction, sphere_center, sphere_radius, t))
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection = ray_origin + ray_direction * t;

                // Compute normal at the intersection point
                Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

                // Add shading parameter here
                Vector3d l = (light_position - ray_intersection).normalized();
                Vector3d v = (camera_origin - ray_intersection).normalized();
                Vector3d h = (l + v).normalized();

                const double diffuse = std::max(0., l.dot(ray_normal));
                const double specular = std::pow(std::max(0., h.dot(ray_normal)), specular_exponent);

                // Simple diffuse model
                Vector3d colored_pixel = Vector3d::Constant(ambient) + diffuse_color.cwiseProduct(diffuse * light_intesity) + specular_color.cwiseProduct(specular * light_intesity);

                // Clamp to zero
                colored_pixel = colored_pixel.cwiseMax(0.0).cwiseMin(1.0);

                R(i, j) = colored_pixel(0);
                G(i, j) = colored_pixel(1);
                B(i, j) = colored_pixel(2);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}
