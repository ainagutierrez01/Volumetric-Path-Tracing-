#include "coloredvolumetricpathshader.h"
#include "../core/utils.h"
#include "../core/hemisphericalsampler.h"
#include <random>

#define _USE_MATH_DEFINES
#include "math.h"

Colored::Colored() :
    hitColor(Vector3D(1, 0, 0))
{
}

Colored::Colored(Vector3D hitColor_, Vector3D bgColor_) :
    Shader(bgColor_), hitColor(hitColor_)
{
}

// Helper function for uniform random sampling
double randomUniform() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

// Sample distance in participating media using exponential distribution
double sampleDistance(double sigma_t) {
    double xi = randomUniform();
    return -log(1.0 - xi) / sigma_t;
}

// Compute colored transmittance through the medium
Vector3D computeTransmittance(double distance, const Vector3D& extinction) {
    return Vector3D(exp(-extinction.x * distance), exp(-extinction.y * distance), exp(-extinction.z * distance));
}

// Compute inscattering contribution along a ray segment with colored medium
Vector3D Colored::computeInscattering(const Ray& r, double t_max, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, int depth) const {

    // Color-dependent absorption and scattering (creates colored shadows)
    // Red shadows: absorb cyan/blue, let red pass through shadows
    Vector3D mu_a(0.003, 0.025, 0.030);  // Absorption (R, G, B) - absorb cyan/blue heavily
    Vector3D mu_s(0.15, 0.03, 0.02);     // Scattering (R, G, B) - scatter red heavily in light areas
    Vector3D mu_t = mu_a + mu_s;         // Total extinction per channel

    // Use average extinction for distance sampling
    double mu_t_avg = (mu_t.x + mu_t.y + mu_t.z) / 3.0;

    Vector3D L_inscatter(0.0);

    double t_sample = sampleDistance(mu_t_avg);

    if (t_sample < t_max) {                                                 // If scattering event happens before intersection
        Vector3D scatter_pos = r.o + r.d * t_sample;                        // Scattering event position
        Vector3D T_to_scatter = computeTransmittance(t_sample, mu_t);       // Transmittance

        Vector3D L_direct_scatter(0.0);
        for (auto light : lsList) {                                         // Iterate for every light source
            Vector3D light_pos = light->sampleLightPosition();              // Light position
            Vector3D light_dir = (light_pos - scatter_pos).normalized();    // Light direction
            double light_dist = (light_pos - scatter_pos).length();         // Light distance  

            Ray shadow_ray(scatter_pos, light_dir);
            shadow_ray.maxT = light_dist - Epsilon;                         // Cast shadow ray

            if (!Utils::hasIntersection(shadow_ray, objList)) {             // If not occluded, compute contribution --> NEE
                Vector3D light_n = light->getNormal();
                double light_area = light->getArea();
                Vector3D T_to_light = computeTransmittance(light_dist, mu_t);
                double phase = 1.0 / (4.0 * M_PI);
                double geometric = abs(dot(-light_dir, light_n)) / ((light_pos - scatter_pos).lengthSq());
                Vector3D emitted = light->getIntensity() * 2.0;             // Boost light intensity to make volumetric effects more visible

                // Component-wise multiplication for colored transmission
                L_direct_scatter += emitted * T_to_light * phase * geometric * light_area;  // Accumulate contribution
            }
        }

        // Color-dependent albedo
        Vector3D albedo = mu_s / mu_t;
        L_inscatter = L_direct_scatter * albedo * T_to_scatter;                     // Direct scattering contribution

        // Multiple scattering
        if (depth < 5) {
            HemisphericalSampler hs;
            Vector3D random_dir = hs.getSample(Vector3D(0, 1, 0)).normalized();
            Ray scatter_ray(scatter_pos, random_dir, depth + 1);                    // Cast new ray from scattering point
            scatter_ray.maxT = 1000.0;

            Intersection scatter_its;
            if (Utils::getClosestIntersection(scatter_ray, objList, scatter_its)) {
                double dist_to_surface = (scatter_its.itsPoint - scatter_pos).length();
                Vector3D L_indirect = computeInscattering(scatter_ray, dist_to_surface, objList, lsList, depth + 1);        // Recursibe inscattering
                L_inscatter += L_indirect * albedo * T_to_scatter * (1.0 / (4.0 * M_PI));
            }
        }
    }

    return L_inscatter;
}

Vector3D Colored::computeRadiance(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const int MaxDepth) const {
    Intersection its;
    Vector3D Lo(0.0, 0.0, 0.0);

    // Colored extinction coefficients - red shadows
    Vector3D mu_a(0.003, 0.025, 0.030);
    Vector3D mu_s(0.15, 0.03, 0.02);
    Vector3D mu_t = mu_a + mu_s;

    if (Utils::getClosestIntersection(r, objList, its)) {
        double distance = (its.itsPoint - r.o).length();
        Vector3D T = computeTransmittance(distance, mu_t);                                      // Transmittance through heterogeneous media

        Vector3D L_inscatter = computeInscattering(r, distance, objList, lsList, r.depth);      // In-scattering contribution
        Lo += L_inscatter;
        Lo += its.shape->getMaterial().getEmissiveRadiance() * T;                               // Surface emision attenuated

        // Surface reflection
        if (r.depth < MaxDepth) {
            Vector3D wo = -r.d.normalized();
            Vector3D n = its.normal.normalized();
            Vector3D L_reflected = ReflectedRadiance(its, wo, r.depth, objList, lsList);
            Lo += L_reflected * T;
        }
        return Lo;
    }

    // Background with colored fog
    double tau_avg = (mu_t.x + mu_t.y + mu_t.z) / 3.0;
    double exp_r = exp(-mu_t.x);
    double exp_g = exp(-mu_t.y);
    double exp_b = exp(-mu_t.z);
    Vector3D T_light = Vector3D(exp_r > 0.3 ? exp_r : 0.3, exp_g > 0.3 ? exp_g : 0.3, exp_b > 0.3 ? exp_b : 0.3);
    Vector3D ambient = bgColor * 0.2;
    return ambient + bgColor * T_light;
}

Vector3D Colored::ReflectedRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    Vector3D L_dir = DirectRadiance(its, wo, objList, lsList);
    Vector3D L_ind = IndirectRadiance(its, wo, depth, objList, lsList);

    return L_dir + L_ind;
}

Vector3D Colored::DirectRadiance(const Intersection& its, const Vector3D& wo, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    Vector3D L_direct = Vector3D(0.0);

    // Colored extinction
    Vector3D mu_a(0.005, 0.015, 0.025);
    Vector3D mu_s(0.08, 0.04, 0.02);
    Vector3D mu_t = mu_a + mu_s;

    for (int light_idx = 0; light_idx < lsList.size(); light_idx++) {
        double light_area = lsList[light_idx]->getArea();
        Vector3D y = lsList[light_idx]->sampleLightPosition();

        Vector3D n = its.normal.normalized();
        Vector3D wi = (y - its.itsPoint).normalized();

        double p_y = 1.0 / light_area;

        Vector3D reflectance = its.shape->getMaterial().getReflectance(n, wo, wi);

        Ray ray_wi(its.itsPoint + n * Epsilon, wi);
        ray_wi.maxT = (y - its.itsPoint).length() - Epsilon;

        size_t V;

        if (Utils::hasIntersection(ray_wi, objList)) {
            V = 0;
        }
        else {
            V = 1;
        }

        Vector3D light_n = lsList[light_idx]->getNormal();

        double geometric_term = (dot(n, wi) * abs(dot(-wi, light_n))) / (y - its.itsPoint).lengthSq();
        Vector3D emited_light = lsList[light_idx]->getIntensity();

        // Colored volumetric transmittance
        double dist = (y - its.itsPoint).length();
        Vector3D T_light = computeTransmittance(dist, mu_t);

        // Ensure positive dot product
        double cos_term = std::max(0.0, dot(n, wi));

        L_direct += (emited_light * reflectance * T_light * cos_term * abs(dot(-wi, light_n)) * V) / ((y - its.itsPoint).lengthSq() * p_y);
    }
    return L_direct;
}

Vector3D Colored::IndirectRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const size_t maxDepth) const {
    HemisphericalSampler hs;
    Vector3D L_indirect = Vector3D(0.0);

    // Cosine-weighted hemisphere sampling
    double r1 = randomUniform();
    double r2 = randomUniform();

    double phi = 2.0 * M_PI * r1;
    double cosTheta = sqrt(1.0 - r2);
    double sinTheta = sqrt(r2);

    // Local coordinate system
    Vector3D n = its.normal.normalized();
    Vector3D tangent, bitangent;
    if (fabs(n.x) > fabs(n.y)) {
        tangent = Vector3D(-n.z, 0, n.x).normalized();
    }
    else {
        tangent = Vector3D(0, n.z, -n.y).normalized();
    }
    bitangent = cross(n, tangent);

    // Transform to world space
    Vector3D wi = (tangent * cos(phi) * sinTheta + bitangent * sin(phi) * sinTheta + n * cosTheta).normalized();

    double p_wi = cosTheta / M_PI;

    Ray new_r(its.itsPoint + n * Epsilon, wi, depth + 1);
    new_r.maxT = 1000.0;

    // Recursive path tracing
    if (depth < maxDepth) {
        Intersection y;
        if (Utils::getClosestIntersection(new_r, objList, y)) {
            Vector3D reflectance = its.shape->getMaterial().getReflectance(n, wo, wi);
            double cos_term = std::max(0.0, dot(n, wi));

            Vector3D incoming_radiance = computeRadiance(new_r, objList, lsList, maxDepth);

            L_indirect += incoming_radiance * reflectance * cos_term / p_wi;
        }
    }

    return L_indirect;
}

Vector3D Colored::computeColor(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    size_t N = 100;
    Vector3D vector_total(0.0, 0.0, 0.0);

    for (size_t idx = 0; idx < N; idx++)
        vector_total += computeRadiance(r, objList, lsList, 5);

    return vector_total / N;
}