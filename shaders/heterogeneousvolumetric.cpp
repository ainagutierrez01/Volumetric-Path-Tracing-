#include "heterogeneousvolumetric.h"
#include "../core/utils.h"
#include "../core/hemisphericalsampler.h"
#include <random>

#define _USE_MATH_DEFINES
#include "math.h"

Heterogeneous::Heterogeneous() :
    hitColor(Vector3D(1, 0, 0))
{
}

Heterogeneous::Heterogeneous(Vector3D hitColor_, Vector3D bgColor_) :
    Shader(bgColor_), hitColor(hitColor_)
{
}

// Helper function for uniform random sampling
static double randomUniform() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

// ==================== NOISE FUNCTIONS ====================

// Permutation function for noise
double permute(double x) {
    return fmod(((x * 34.0) + 1.0) * x, 289.0);
}

// Vector3D permutation
Vector3D permute(const Vector3D& x) {
    return Vector3D(
        fmod(((x.x * 34.0) + 1.0) * x.x, 289.0),
        fmod(((x.y * 34.0) + 1.0) * x.y, 289.0),
        fmod(((x.z * 34.0) + 1.0) * x.z, 289.0)
    );
}

// Taylor inverse square root approximation
double taylorInvSqrt(double r) {
    return 1.79284291400159 - 0.85373472095314 * r;
}

// 3D Simplex Noise (Adapted from GLSL implementation)
double snoise(const Vector3D& v) {
    const double C1 = 1.0 / 6.0;
    const double C2 = 1.0 / 3.0;

    // First corner
    double s = (v.x + v.y + v.z) * C2;
    Vector3D i = Vector3D(floor(v.x + s), floor(v.y + s), floor(v.z + s));

    double t = (i.x + i.y + i.z) * C1;
    Vector3D x0 = v - (i - Vector3D(t, t, t));

    // Other corners
    Vector3D g = Vector3D(
        x0.x >= x0.y ? 1.0 : 0.0,
        x0.y >= x0.z ? 1.0 : 0.0,
        x0.z >= x0.x ? 1.0 : 0.0
    );
    Vector3D l = Vector3D(1.0, 1.0, 1.0) - g;

    Vector3D i1 = Vector3D(
        std::min(g.x, l.z),
        std::min(g.y, l.x),
        std::min(g.z, l.y)
    );
    Vector3D i2 = Vector3D(
        std::max(g.x, l.z),
        std::max(g.y, l.x),
        std::max(g.z, l.y)
    );

    Vector3D x1 = x0 - i1 + Vector3D(C1, C1, C1);
    Vector3D x2 = x0 - i2 + Vector3D(2.0 * C1, 2.0 * C1, 2.0 * C1);
    Vector3D x3 = x0 - Vector3D(1.0, 1.0, 1.0) + Vector3D(3.0 * C1, 3.0 * C1, 3.0 * C1);

    // Permutations
    i = Vector3D(fmod(i.x, 289.0), fmod(i.y, 289.0), fmod(i.z, 289.0));

    double p0 = permute(permute(permute(i.z) + i.y) + i.x);
    double p1 = permute(permute(permute(i.z + i1.z) + i.y + i1.y) + i.x + i1.x);
    double p2 = permute(permute(permute(i.z + i2.z) + i.y + i2.y) + i.x + i2.x);
    double p3 = permute(permute(permute(i.z + 1.0) + i.y + 1.0) + i.x + 1.0);

    // Simplified gradient calculation
    double n_ = 1.0 / 7.0;
    double j0 = p0 - 49.0 * floor(p0 * n_ * n_);
    double j1 = p1 - 49.0 * floor(p1 * n_ * n_);
    double j2 = p2 - 49.0 * floor(p2 * n_ * n_);
    double j3 = p3 - 49.0 * floor(p3 * n_ * n_);

    // Create gradients
    auto makeGradient = [n_](double j) -> Vector3D {
        double x_ = floor(j * n_);
        double y_ = floor(j - 7.0 * x_);
        double x = x_ * n_ + n_;
        double y = y_ * n_ + n_;
        double h = 1.0 - abs(x) - abs(y);
        double sx = x >= 0 ? 1.0 : -1.0;
        double sy = y >= 0 ? 1.0 : -1.0;
        if (h < 0.0) {
            x = x - sx * h;
            y = y - sy * h;
        }
        return Vector3D(x, y, h);
        };

    Vector3D grad0 = makeGradient(j0);
    Vector3D grad1 = makeGradient(j1);
    Vector3D grad2 = makeGradient(j2);
    Vector3D grad3 = makeGradient(j3);

    // Normalize
    double norm0 = taylorInvSqrt(dot(grad0, grad0));
    double norm1 = taylorInvSqrt(dot(grad1, grad1));
    double norm2 = taylorInvSqrt(dot(grad2, grad2));
    double norm3 = taylorInvSqrt(dot(grad3, grad3));

    grad0 = grad0 * norm0;
    grad1 = grad1 * norm1;
    grad2 = grad2 * norm2;
    grad3 = grad3 * norm3;

    // Mix contributions
    double m0 = std::max(0.6 - dot(x0, x0), 0.0);
    double m1 = std::max(0.6 - dot(x1, x1), 0.0);
    double m2 = std::max(0.6 - dot(x2, x2), 0.0);
    double m3 = std::max(0.6 - dot(x3, x3), 0.0);

    m0 *= m0; m0 *= m0;
    m1 *= m1; m1 *= m1;
    m2 *= m2; m2 *= m2;
    m3 *= m3; m3 *= m3;

    return 42.0 * (m0 * dot(grad0, x0) + m1 * dot(grad1, x1) +
        m2 * dot(grad2, x2) + m3 * dot(grad3, x3));
}

// Fractional Brownian Motion
double fbm(const Vector3D& p, double scale = 1.0, int octaves = 3) {
    double value = 0.0;
    double amplitude = 0.5;
    Vector3D position = p / scale;
    Vector3D shift = Vector3D(100.0, 100.0, 100.0);

    for (int i = 0; i < octaves; i++) {
        value += amplitude * snoise(position);
        position = position * 2.0 + shift;
        amplitude *= 0.5;
    }

    return fabs(value);
}

// Get density at a point in space using FBM
double getDensity(const Vector3D& pos, double base_density = 0.8) {
    // Single FBM call with 2 octaves for good balance of detail and speed
    double noise = fbm(pos, 3.0, 2);

    // Remap from [-1, 1] to [0, 1]
    noise = noise * 0.5 + 0.5;

    // Apply power function to create patchier distribution
    // Values close to 0 become very sparse, values close to 1 stay dense
    double density = pow(noise, 2.5);

    // Add threshold to create clear empty regions
    if (density < 0.15) {
        density = 0.0;
    }

    // Scale by base density
    density *= base_density;

    return std::max(0.0, std::min(1.0, density));
}

// ==================== VOLUMETRIC RENDERING ====================

// Compute transmittance between two points through heterogeneous media
double computeTransmittance(const Vector3D& p1, const Vector3D& p2, double mu_s, double mu_a) {
    Vector3D dir = (p2 - p1);
    double dist = dir.length();
    dir = dir.normalized();

    const int STEPS = 16;  // Reduced from 50 for performance
    const double step_size = dist / STEPS;

    double tau = 0.0;
    for (int i = 0; i < STEPS; i++) {
        Vector3D pos = p1 + dir * (i * step_size);
        double density = getDensity(pos);
        tau += (mu_a + mu_s) * density * step_size;
    }

    return exp(-tau);
}

// ------------------ Delta (Woodcock) Tracking ------------------
double sampleDistanceDeltaTracking(const Ray& r, double t_max, double& transmittance, Vector3D& sample_pos, double mu_s, double mu_a)
{
    const double sigma_t_max = (mu_s + mu_a);

    double t = 0.0;

    while (true)
    {
        // 1. Sample exponential distance
        double xi = randomUniform();
        t += -log(1.0 - xi) / sigma_t_max;

        // If we step outside the volume, no collision
        if (t >= t_max)
        {
            sample_pos = r.o + r.d * t_max;
            return t_max;   // "no real collision inside"
        }

        // 2. Compute real extinction at sampled position
        Vector3D pos = r.o + r.d * t;
        double density = getDensity(pos);                    // in [0,1]
        double sigma_t_local = (mu_s + mu_a) * density;      // local sigma_t(x)

        // 3. Acceptance / rejection
        if (randomUniform() < (sigma_t_local / sigma_t_max))
        {
            sample_pos = pos;   // real event
            return t;           // return distance to collision
        }

        // else: null-collision --> keep looping
    }
}

// Compute inscattering contribution along a ray segment
Vector3D Heterogeneous::computeInscattering(const Ray& r, double t_max, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, int depth) const {

    double mu_a = 0.05;  // Absorption - reduced for brighter patches
    double mu_s = 0.1;   // Scattering - increased for more visible fog patches

    Vector3D L_inscatter(0.0);

    Vector3D scatter_pos;
    double T_to_scatter;
    double t_sample = sampleDistanceDeltaTracking(r, t_max, T_to_scatter, scatter_pos, mu_s, mu_a);     // Delta tracking sampling --> distance to next event

    if (t_sample < t_max) {                                 // if a scattering event occurs before intersection hit
        double local_density = getDensity(scatter_pos);     // Get density at scattering point (since it's heterogeneous media)
        double local_mu_a = mu_a * local_density;           // Local absorption coefficient
        double local_mu_s = mu_s * local_density;           // Local scattering coefficient
        double local_mu_t = (mu_a + mu_s) * local_density;  // Local extinction coefficient

        if (local_mu_t > 0.0001) {                          
            // Direct lighting from light sources
            Vector3D L_direct_scatter(0.0);

            double pScatter = local_mu_s / local_mu_t;      // Local single scattering albedo
            double xi = randomUniform();

            if (xi > pScatter) {    // ABSORPTION -> ray dies, contributes no radiance

                return Vector3D(0.0);

            }
            else {      // SCATTERING EVENT -> continue in-scattering
                for (auto light : lsList) {
                    Vector3D light_pos = light->sampleLightPosition();              // Sample light position
                    Vector3D light_dir = (light_pos - scatter_pos).normalized();    // Direction to light
                    double light_dist = (light_pos - scatter_pos).length();         // Distance to light

                    Ray shadow_ray(scatter_pos, light_dir);
                    shadow_ray.maxT = light_dist - Epsilon;                         // Shadow ray to light

                    if (!Utils::hasIntersection(shadow_ray, objList)) {             // If not occluded, compute contribution
                        Vector3D light_n = light->getNormal();
                        double light_area = light->getArea();

                        // Transmittance through heterogeneous media
                        double T_to_light = computeTransmittance(scatter_pos, light_pos, local_mu_s, local_mu_a);

                        // Isotropic phase function
                        double phase = 1.0 / (4.0 * M_PI);
                        double geometric = abs(dot(-light_dir, light_n)) / ((light_pos - scatter_pos).lengthSq());
                        Vector3D emitted = light->getIntensity() * 3.0;             // Boost light intensity to make volumetric effects more visible

                        L_direct_scatter += emitted * phase * geometric * T_to_light * light_area;  // Accumulate contribution
                    }
                }
            }

            double albedo = local_mu_s / local_mu_t;                            // Local single scattering albedo
            L_inscatter = L_direct_scatter * albedo * T_to_scatter;             // Direct inscattering contribution

            // Multiple scattering
            if (depth < 5) {
                HemisphericalSampler hs;
                Vector3D random_dir = hs.getSample(Vector3D(0, 1, 0)).normalized();
                Ray scatter_ray(scatter_pos, random_dir, depth + 1);            // Cast new ray from scattering point
                scatter_ray.maxT = 1000.0;

                Intersection scatter_its;
                if (Utils::getClosestIntersection(scatter_ray, objList, scatter_its)) {
                    double dist_to_surface = (scatter_its.itsPoint - scatter_pos).length();
                    Vector3D L_indirect = computeInscattering(scatter_ray, dist_to_surface,     // Recursive inscattering
                        objList, lsList, depth + 1);
                    L_inscatter += L_indirect * albedo * T_to_scatter * (1.0 / (4.0 * M_PI));
                }
            }
        }
    }

    return L_inscatter;
}

Vector3D Heterogeneous::computeRadiance(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const int MaxDepth) const {
    Intersection its;
    Vector3D Lo(0.0, 0.0, 0.0);

    double mu_a = 0.05;
    double mu_s = 0.1;

    if (Utils::getClosestIntersection(r, objList, its)) {
        double distance = (its.itsPoint - r.o).length();

        // Compute transmittance through heterogeneous media
        double T = computeTransmittance(r.o, its.itsPoint, mu_s, mu_a);

        // In-scattering contribution
        Vector3D L_inscatter = computeInscattering(r, distance, objList, lsList, r.depth);
        Lo += L_inscatter;

        // Surface emission
        Lo += its.shape->getMaterial().getEmissiveRadiance() * T;

        // Surface reflection
        if (r.depth < MaxDepth) {
            Vector3D wo = -r.d.normalized();
            Vector3D n = its.normal.normalized();
            Vector3D L_reflected = ReflectedRadiance(its, wo, r.depth, objList, lsList);
            Lo += L_reflected * T;
        }
        return Lo;
    }

    // Background with atmospheric scattering
    double avg_density = 0.3;
    double tau = (mu_a + mu_s) * avg_density * 10.0;
    double T_light = std::max(0.2, exp(-tau));
    Vector3D ambient = bgColor * 0.15;
    return ambient + bgColor * T_light;
}

Vector3D Heterogeneous::ReflectedRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    Vector3D L_dir = DirectRadiance(its, wo, objList, lsList);
    Vector3D L_ind = IndirectRadiance(its, wo, depth, objList, lsList);

    return L_dir + L_ind;
}

Vector3D Heterogeneous::DirectRadiance(const Intersection& its, const Vector3D& wo, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    Vector3D L_direct = Vector3D(0.0);

    double mu_a = 0.05;
    double mu_s = 0.1;

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

        // Heterogeneous volumetric transmittance
        double T_light = computeTransmittance(its.itsPoint, y, mu_s, mu_a);

        double cos_term = std::max(0.0, dot(n, wi));

        L_direct += (emited_light * reflectance * cos_term * abs(dot(-wi, light_n)) * V * T_light) / ((y - its.itsPoint).lengthSq() * p_y);
    }
    return L_direct;
}

Vector3D Heterogeneous::IndirectRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const size_t maxDepth) const {
    HemisphericalSampler hs;
    Vector3D L_indirect = Vector3D(0.0);

    // Cosine-weighted hemisphere sampling
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

    double p_wi = cosTheta / M_PI;

    Ray new_r(its.itsPoint + n * Epsilon, wi, depth + 1);
    new_r.maxT = 1000.0;

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

Vector3D Heterogeneous::computeColor(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    size_t N = 128;
    Vector3D vector_total(0.0, 0.0, 0.0);

    for (size_t idx = 0; idx < N; idx++)
        vector_total += computeRadiance(r, objList, lsList, 5);

    return vector_total / N;
}