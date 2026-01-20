#include "homogenousvolumetricpathshader.h"
#include "../core/utils.h"
#include "../core/hemisphericalsampler.h"
#include <random>

#define _USE_MATH_DEFINES
#include "math.h"

Homogenous::Homogenous():
    hitColor(Vector3D(1, 0, 0))
{
}

Homogenous::Homogenous(Vector3D hitColor_, Vector3D bgColor_) :
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

// Compute inscattering contribution along a ray segment
Vector3D Homogenous::computeInscattering(const Ray& r, double t_max, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, int depth) const {

    // Absorption and scattering coefficients
    double mu_a = 0.05;
    double mu_s = 0.1;
    double mu_t = mu_a + mu_s;  // 0.05 - 0.5     // 0.1 - 1.0

    Vector3D L_inscatter(0.0);

    double t_sample = sampleDistance(mu_t);

    if (t_sample < t_max) {     // if a scattering event occurs before surface hit
        Vector3D scatter_pos = r.o + r.d * t_sample;
        double T_to_scatter = exp(-mu_t * t_sample);

        Vector3D L_direct_scatter(0.0);
        // -------------------------------------------------------------------------------
        double pScatter = mu_s / mu_t;
        double xi = randomUniform();
    
        if (xi > pScatter) {    // ABSORPTION -> ray dies, contributes no radiance

            return Vector3D(0.0);

        }
        else {                  // SCATTERING EVENT -> continue in-scattering
            for (auto light : lsList) {
                Vector3D light_pos = light->sampleLightPosition();              // Sample light position
                Vector3D light_dir = (light_pos - scatter_pos).normalized();    // Direction to light
                double light_dist = (light_pos - scatter_pos).length();         // Distance to light

                Ray shadow_ray(scatter_pos, light_dir);
                shadow_ray.maxT = light_dist - Epsilon;                         // Shadow ray to light

                if (!Utils::hasIntersection(shadow_ray, objList)) {             // If not occluded, compute contribution
                    Vector3D light_n = light->getNormal();
                    double light_area = light->getArea();
                    double T_to_light = exp(-mu_t * light_dist);                // Transmittance to light
                    double phase = 1.0 / (4.0 * M_PI);                          // Isotropic phase function
                    double geometric = abs(dot(-light_dir, light_n)) / ((light_pos - scatter_pos).lengthSq());
                    Vector3D emitted = light->getIntensity() * 2.0;             // Boost light intensity to make volumetric effects more visible 
                    L_direct_scatter += emitted * phase * geometric * T_to_light * light_area;  // Accumulate contribution
                }
            }
        }
        // -------------------------------------------------------------------------------
        double albedo = mu_s / mu_t;                                            // Single scattering albedo
        L_inscatter = L_direct_scatter * albedo * T_to_scatter;                 // Direct inscattering contribution

        // Multiple scattering
        if (depth < 10) {
            HemisphericalSampler hs;
            Vector3D random_dir = hs.getSample(Vector3D(0, 1, 0)).normalized();
            Ray scatter_ray(scatter_pos, random_dir, depth + 1);                // Cast new ray from scattering point
            scatter_ray.maxT = 1000.0;

            Intersection scatter_its;
            if (Utils::getClosestIntersection(scatter_ray, objList, scatter_its)) {
                double dist_to_surface = (scatter_its.itsPoint - scatter_pos).length();
                Vector3D L_indirect = computeInscattering(scatter_ray, dist_to_surface, objList, lsList, depth + 1);    // Recursive inscattering 
                L_inscatter += L_indirect * albedo * T_to_scatter * (1.0 / (4.0 * M_PI));
            }
        }
    }

    return L_inscatter;

}

Vector3D Homogenous::computeRadiance(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const int MaxDepth) const {
    Intersection its;
    Vector3D Lo(0.0, 0.0, 0.0);

    // Absorption and scattering coefficients
    double mu_a = 0.05;
    double mu_s = 0.1;
    double mu_t = mu_a + mu_s;

    if (Utils::getClosestIntersection(r, objList, its)) {
        double distance = (its.itsPoint - r.o).length();
        double T = exp(-mu_t * distance);                                                       // Transmittance from ray origin to intersection

        Vector3D L_inscatter = computeInscattering(r, distance, objList, lsList, r.depth);      // In-scattering contribution along the ray
        Lo += L_inscatter;
        Lo += its.shape->getMaterial().getEmissiveRadiance() * T;                               // Emissive contribution attenuated by transmittance

        if (r.depth < MaxDepth) {                                                               // Indirect and direct lighting --> recursive path tracing
            Vector3D wo = -r.d.normalized();
            Vector3D n = its.normal.normalized();
            Vector3D L_reflected = ReflectedRadiance(its, wo, r.depth, objList, lsList);
            Lo += L_reflected * T;
        }
        return Lo;
    }

    // If no intersection, return background color with volumetric attenuation
    double tau = mu_t * r.maxT;
    double T_light = std::max(0.3, exp(-tau));  // Increased base light  
    Vector3D ambient = bgColor * 0.2;           // Increased ambient  
    return ambient + bgColor * T_light;

}

Vector3D Homogenous::ReflectedRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    Vector3D L_dir = DirectRadiance(its, wo, objList, lsList);
    Vector3D L_ind = IndirectRadiance(its, wo, depth, objList, lsList);

    return L_dir + L_ind;
}

Vector3D Homogenous::DirectRadiance(const Intersection& its, const Vector3D& wo, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    Vector3D L_direct = Vector3D(0.0);

    // Absorption and scattering coefficients
    double mu_a = 0.05;
    double mu_s = 0.1;
    double mu_t = mu_a + mu_s;

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

        // Volumetric transmittance
        double dist = (y - its.itsPoint).length();
        double tau = mu_t * dist;
        double T_light = exp(-tau);

        // Ensure dot product is positive for proper lighting
        double cos_term = std::max(0.0, dot(n, wi));

        L_direct += (emited_light * reflectance * cos_term * abs(dot(-wi, light_n)) * V * T_light) / ((y - its.itsPoint).lengthSq() * p_y);
    }
    return L_direct;
}

Vector3D Homogenous::IndirectRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const size_t maxDepth) const {
    HemisphericalSampler hs;
    Vector3D L_indirect = Vector3D(0.0);

    // Cosine-weighted hemisphere sampling for better convergence
    double r1 = randomUniform();
    double r2 = randomUniform();

    double phi = 2.0 * M_PI * r1;
    double cosTheta = sqrt(1.0 - r2);
    double sinTheta = sqrt(r2);

    // Build local coordinate system
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
    Vector3D wi = (tangent * cos(phi) * sinTheta +
        bitangent * sin(phi) * sinTheta +
        n * cosTheta).normalized();

    double p_wi = cosTheta / M_PI;  // Cosine-weighted PDF

    Ray new_r(its.itsPoint + n * Epsilon, wi, depth + 1);
    new_r.maxT = 1000.0;

    if (depth < maxDepth) {
        Intersection y;
        if (Utils::getClosestIntersection(new_r, objList, y)) {
            Vector3D reflectance = its.shape->getMaterial().getReflectance(n, wo, wi);
            double cos_term = std::max(0.0, dot(n, wi));

            // Recursively compute radiance
            Vector3D incoming_radiance = computeRadiance(new_r, objList, lsList, maxDepth);
            
            L_indirect += incoming_radiance * reflectance * cos_term / p_wi ;

        }
    }

    return L_indirect;
}

Vector3D Homogenous::computeColor(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    size_t N = 512;
    Vector3D vector_total(0.0, 0.0, 0.0);

    for (size_t idx = 0; idx < N; idx++)
        vector_total += computeRadiance(r, objList, lsList, 5);

    return vector_total / N;

}