////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <limits>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////

const std::string filename("raytrace.png");

//Camera settings
const double focal_length = 10;
const double field_of_view = 0.7854; //45 degrees
const double image_z = 5;
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 5);
const double camera_aperture = 0.05; // Changes how much blurr happens when using depth of field, bigger number = more out of foucs.
const bool use_dof = false; // Whether to use depth of field or not.
const int dof_samples = 32; // Number of simulated rays per pixel, makes image more consistent but takes more compute time.
// Constants for refraction
const double obj_refraction_index = 1.5;
const bool use_refraction = false;

//Maximum number of recursive calls
const int max_bounce = 5;

// Objects
std::vector<Vector3d> sphere_centers;
std::vector<double> sphere_radii;
std::vector<Matrix3d> parallelograms;

//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.5, 0.1, 0.1, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);
const Vector4d obj_refraction_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector4d> light_colors;
//Ambient light
const Vector4d ambient_light(0.2, 0.2, 0.2, 0);

//Fills the different arrays
void setup_scene()
{
    grid.resize(grid_size + 1);
    for (int i = 0; i < grid_size + 1; ++i)
    {
        grid[i].resize(grid_size + 1);
        for (int j = 0; j < grid_size + 1; ++j)
            grid[i][j] = Vector2d::Random().normalized();
    }

    //Spheres
    sphere_centers.emplace_back(10, 0, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(7, 0.05, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(4, 0.1, 1);
    // sphere_centers.emplace_back(4, 0.1, -5); // having z index of -5 matches the camera focal place so its in focus when using depth of field.
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(1, 0.2, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-2, 0.4, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-5, 0.8, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-8, 1.6, 1);
    sphere_radii.emplace_back(1);

    //parallelograms
    parallelograms.emplace_back();
    parallelograms.back() << -100, -100, 100,
        -1.25, -1.2, 0,
        -100, 100, -100;

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);
}

//We need to make this function visible
Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce);

////////////////////////////////////////////////////////////////////////////////
// Perlin noise code
////////////////////////////////////////////////////////////////////////////////

// Function to linearly interpolate between a0 and a1
// Weight w should be in the range [0.0, 1.0]
double lerp(double a0, double a1, double w)
{
    assert(w >= 0);
    assert(w <= 1);
    //linear and cubic interpolation
    // Linear
    double linear = a0 + w * (a1 - a0);

    // Cubic
    double cubic = (a1 - a0) * (3.0 - w * 2.0) * w * w + a0;

    return cubic;
}

// Computes the dot product of the distance and gradient vectors.
double dotGridGradient(int ix, int iy, double x, double y)
{
    // Compute the distance vector
    // Where x and y are point, and ix, and iy are the ith corner points
    Vector2d dist (x - ix, y - iy);
    // Compute and return the dot-product
    double dot_product = grid[iy][ix].dot(dist);

    return dot_product;
}

// Compute Perlin noise at coordinates x, y
double perlin(double x, double y)
{
    //Determine grid cell coordinates x0, y0
    //computes with floor and ceiling like picture example in the slides
    int x0 = floor(x);
    int x1 = ceil(x);
    int y0 = floor(y);
    int y1 = ceil(y);

    // Determine interpolation weights
    double sx = x - x0;
    double sy = y - y0;

    // Interpolate between grid point gradients
    double n0 = dotGridGradient(x0, y0, x, y);
    double n1 = dotGridGradient(x1, y0, x, y);

    double ix0 = lerp(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);

    double ix1 = lerp(n0, n1, sx);
    double value = lerp(ix0, ix1, sy);

    return value;
}

Vector4d procedural_texture(const double tu, const double tv)
{
    assert(tu >= 0);
    assert(tv >= 0);

    assert(tu <= 1);
    assert(tv <= 1);

    //DONE: uncomment these lines once you implement the perlin noise
    const double color = (perlin(tu * grid_size, tv * grid_size) + 1) / 2;
    return Vector4d(0, color, 0, 0);

    //Example fo checkerboard texture
    // const double color = (int(tu * grid_size) + int(tv * grid_size)) % 2 == 0 ? 0 : 1;
    // return Vector4d(0, color, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

//Compute the intersection between a ray and a sphere, return -1 if no intersection
double ray_sphere_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    // DONE, implement the intersection between the ray and the sphere at index index.
    //return t or -1 if no intersection

    const Vector3d sphere_center = sphere_centers[index];
    const double sphere_radius = sphere_radii[index];


    double a = ray_direction.dot(ray_direction);
    Vector3d center_to_origin =  ray_origin - sphere_center;
    double b = 2 * ray_direction.dot(center_to_origin);
    double c = center_to_origin.dot(center_to_origin) - (sphere_radius * sphere_radius);

    double delta = (b * b) - (4 * a * c);


    double t = -1;

    if (delta < 0)
    {
        return -1;
    }
    else
    {
        double t1 = (-b - std::sqrt(delta)) / (2 * a);
        double t2 = (-b + std::sqrt(delta)) / (2 * a);

        if (t1 >= 0) {
            t = t1;
        } else if (t2 >= 0){
            t = t2;
        } else {
            return -1;
        }

        //DONE set the correct intersection point, update p to the correct value
        p = ray_origin + (ray_direction * t);
        N = (p - sphere_center).normalized();

        return t;
    }

    return -1;
}

//Compute the intersection between a ray and a paralleogram, return -1 if no intersection
double ray_parallelogram_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    // DONE, implement the intersection between the ray and the parallelogram at index index.
    //return t or -1 if no intersection

    const Vector3d pgram_origin = parallelograms[index].col(0);
    const Vector3d A = parallelograms[index].col(1);
    const Vector3d B = parallelograms[index].col(2);
    const Vector3d pgram_u = A - pgram_origin;
    const Vector3d pgram_v = B - pgram_origin;

    Matrix3d M;
    M << pgram_u, pgram_v, -ray_direction;
    Vector3d rhs = ray_origin - pgram_origin;

    bool no_collision = std::abs(M.determinant()) < 1e-9;

    if (no_collision)
    {
        return -1;
    }

    Vector3d solved = M.colPivHouseholderQr().solve(rhs);

    double a = solved(0);
    double b = solved(1);
    double t = solved(2);

    if (t >= 0 && (a >= 0 && a <= 1) && (b >= 0 && b <= 1)) {

        //DONE set the correct intersection point, update p and N to the correct values
        p = ray_origin + (t * ray_direction);
        N = pgram_u.cross(pgram_v).normalized();
        return t;
    }
    return -1;
}

//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
int find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N)
{
    // Find the object in the scene that intersects the ray first
    // we store the index and the 'closest_t' to their expected values
    int closest_index = -1;
    double closest_t = std::numeric_limits<double>::max(); //closest t is "+ infinity"

    Vector3d tmp_p, tmp_N;
    for (int i = 0; i < sphere_centers.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_sphere_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            //The point is before our current closest t
            if (t < closest_t)
            {
                closest_index = i;
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    for (int i = 0; i < parallelograms.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_parallelogram_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            //The point is before our current closest t
            if (t < closest_t)
            {
                closest_index = sphere_centers.size() + i;
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    return closest_index;
}

////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////

//Checks if the light is visible
bool is_light_visible(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &light_position)
{
    // DONE: Determine if the light is visible here
    // Use find_nearest_object
    Vector3d intersection_point;
    Vector3d intersection_normal;

    const int nearest_object = find_nearest_object(ray_origin, ray_direction,intersection_point, intersection_normal);
    if (nearest_object == -1){
        return true;
    } else if ( (light_position - ray_origin).norm() < (intersection_point - ray_origin).norm()){
        return true;
    } else {
        return false;
    }
}

Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const int nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (nearest_object < 0)
    {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    // Ambient light contribution
    const Vector4d ambient_color = obj_ambient_color.array() * ambient_light.array();

    // Punctual lights contribution (direct lighting)
    Vector4d lights_color(0, 0, 0, 0);
    for (int i = 0; i < light_positions.size(); ++i)
    {
        const Vector3d &light_position = light_positions[i];
        const Vector4d &light_color = light_colors[i];

        const Vector3d Li = (light_position - p).normalized();

        // Shoot a shadow ray to determine if the light should affect the intersection point and call is_light_visible
        const bool light_visible = is_light_visible(p + 1.e-4 * Li, Li, light_position);
        if (!light_visible) {
            continue;
        }

        Vector4d diff_color = obj_diffuse_color;

        if (nearest_object == 4)
        {
            //Compute UV coodinates for the point on the sphere
            const double x = p(0) - sphere_centers[nearest_object][0];
            const double y = p(1) - sphere_centers[nearest_object][1];
            const double z = p(2) - sphere_centers[nearest_object][2];
            double tu = acos(z / sphere_radii[nearest_object]) / 3.1415;
            double tv = (3.1415 + atan2(y, x)) / (2 * 3.1415);
            tu = std::min(tu, 1.0);
            tu = std::max(tu, 0.0);

            tv = std::min(tv, 1.0);
            tv = std::max(tv, 0.0);

            diff_color = procedural_texture(tu, tv);
        }

        // Add shading parameters

        // Diffuse contribution
        const Vector4d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution, use obj_specular_color
        const Vector3d Hi = (Li - ray_direction).normalized();
        const Vector4d specular =  obj_specular_color * std::pow(std::max(0.0, Hi.dot(N)), obj_specular_exponent);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    Vector4d refl_color = obj_reflection_color;
    if (nearest_object == 4)
    {
        refl_color = Vector4d(0.5, 0.5, 0.5, 0);
    }
    // DONE: Compute the color of the reflected ray and add its contribution to the current point color.
    // use refl_color

    Vector4d reflection_color(0, 0, 0, 0);

    Vector3d v = -ray_direction;
    // r = 2N(N·v) − v
    Vector3d r = 2 * N * (N.dot(v)) - v;

    if (max_bounce > 0) {
        reflection_color = refl_color.cwiseProduct(shoot_ray(p + 1.e-4 * r, r, max_bounce - 1));
    } // else stop recursion and keep reflection color zero



    // DONE: Compute the color of the refracted ray and add its contribution to the current point color.
    //       Make sure to check for total internal reflection before shooting a new ray.
    Vector4d refraction_color(0, 0, 0, 0);


    // Refraction using snells law: n_1 * sin(theta_1) = n_2 * sin(theta_2), eta = n_1 / n_2
    // cos(theta_1) = -dot(d, N), the angle between the incoming ray and the normal
    // the part of the ray along the surface gets scaled by eta: eta * (d + cos(theta_1) * N)
    // the part along the normal keeps going in: -cos(theta_2) * N, where cos(theta_2) = sqrt(1 - sin^2(theta_2))
    // and sin^2(theta_2) = eta^2 * (1 - cos^2(theta_1))
    // adding them gives: t = eta * d + (eta * cos(theta_1) - cos(theta_2)) * N
    // if 1 - sin^2(theta_2) < 0 there is no refracted ray (total internal reflection), that's k < 0 below
    // only spheres refract since the floor is an infinitely thin plane

    if (use_refraction && max_bounce > 0 && nearest_object < sphere_centers.size()) {
        // Derive snells law variables
        Vector3d d = ray_direction; // ray direction incoming
        Vector3d n = N; // Normal vector from point
        double n1;
        double n2;
        if (d.dot(n) < 0) { // entering the sphere
            n1 = 1;
            n2 = obj_refraction_index;

        } else { // leaving sphere
            n1 = obj_refraction_index;
            n2 = 1;
            n = -n; // inside the sphere, so flip N to face against the ray

        }

        double eta   = n1 / n2;
        double cos_theta1 = -d.dot(n);
        double k = 1 - std::pow(eta, 2) * (1 - std::pow(cos_theta1, 2));
        // k = 1 - sin^2(theta_2) and k < 0, means total internal reflection, and there is no refracted ray
        if (k >= 0) {

            Vector3d t = eta * d + (eta * cos_theta1 - sqrt(k)) * n;
            refraction_color = obj_refraction_color.cwiseProduct(shoot_ray(p + 1.e-4 * t, t, max_bounce - 1));

        }
    }

    // Rendering equation
    Vector4d C = ambient_color + lights_color + reflection_color + refraction_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene()
{
    std::cout << "Simple ray tracer." << std::endl;

    int w = 800;
    int h = 400;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    double image_y = focal_length * std::sin(field_of_view/2)/std::sin(M_PI/2- (field_of_view/2)); //DONE: compute the correct pixels size
    double image_x = image_y * aspect_ratio; //DONE: compute the correct pixels size

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, -image_z);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (unsigned i = 0; i < w; ++i)
    {
        for (unsigned j = 0; j < h; ++j)
        {
            // DONE: Implement depth of field

            Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;
            Vector4d C;

            if (is_perspective && use_dof) {
                Vector4d sum(0,0,0,0);
                for (unsigned s = 0; s < dof_samples; s++) {
                    Vector2d offset (0,0);
                    do {
                        offset = Vector2d::Random(); // gives between -1 and 1 in x and y
                    } while (offset.norm() > 1); // try again if its outside a circle to get a circle camera
                    offset *= camera_aperture;
                    ray_origin = camera_position + Vector3d (offset(0), offset(1), 0);
                    ray_direction = (pixel_center - ray_origin).normalized();
                    sum += shoot_ray(ray_origin, ray_direction, max_bounce);
                }
                 C = sum / dof_samples; // average all samples
            } else {


                if (is_perspective)
                {
                    // Perspective camera

                    ray_origin = camera_position;
                    ray_direction = (pixel_center - camera_position).normalized();
                }
                else
                {
                    // Orthographic camera
                    ray_origin = camera_position + Vector3d(pixel_center[0], pixel_center[1], 0);
                    ray_direction = Vector3d(0, 0, -1);
                }

                C = shoot_ray(ray_origin, ray_direction, max_bounce);

            } // end of dof code

            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    setup_scene();

    raytrace_scene();
    return 0;
}
