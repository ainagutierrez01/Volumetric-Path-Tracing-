#include "shader.h"
#include "hemisphericalshader.h"
#include "../materials/phong.h"
#include "../lightsources/pointlightsource.h"
#include "../core/hemisphericalsampler.h"
#include "../core/utils.h"

#define M_PI 3.1416

HemisphericalShader::HemisphericalShader() :
    Color(Vector3D(1, 0, 0))
{
}

HemisphericalShader::HemisphericalShader(Vector3D Color_, Vector3D bgColor_) :
    Shader(bgColor_), Color(Color_)
{
}

Vector3D HemisphericalShader::computeColor(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const {
    Intersection its;
    Vector3D Color = Vector3D(0, 0, 0);

    if (!Utils::getClosestIntersection(r, objList, its)) {
        return bgColor;
    }

    // Handle emissive materials - they emit light themselves
    if (its.shape->getMaterial().isEmissive()) {
        Color = its.shape->getMaterial().getEmissiveRadiance();
    }

    Vector3D wo = -r.d; // Viewing direction (from surface to camera)
    Vector3D n = its.normal.normalized();

    // We use hemispherical sampling for diffuse/glossy materials
    if (its.shape->getMaterial().hasDiffuseOrGlossy()) {
        Vector3D ambient = Vector3D(0.25);
        Vector3D diffuse_coefficient = its.shape->getMaterial().getDiffuseReflectance();
        int numSamples = 512;

        Color += ambient * diffuse_coefficient;

        for (int j = 0; j < numSamples; j++) {
            // Sample uniform hemisphere around normal
            Vector3D wi = HemisphericalSampler().getSample(n);
            double pdf = 1.0 / (2.0 * M_PI); // Uniform hemisphere PDF

            // Test visibility in direction wi
            Ray testRay(its.itsPoint, wi);
            Intersection lightIts;

            Vector3D Li_dir = Vector3D(0, 0, 0);
            if (Utils::getClosestIntersection(testRay, objList, lightIts) && lightIts.shape->getMaterial().isEmissive()) {
                Li_dir = lightIts.shape->getMaterial().getEmissiveRadiance();
            }

            // BRDF evaluation
            Vector3D fr = its.shape->getMaterial().getReflectance(n, wo, wi);

            // Cosine term
            double cosTerm = std::max(0.0, dot(n, wi));

            // Monte Carlo estimatation
            Color += (Li_dir * fr * cosTerm) / (pdf * numSamples);
        }

    }
    // We use recursive ray tracing for perfect reflection (specular materials)
    else if (its.shape->getMaterial().hasSpecular()) {
        // We calculate the ideal reflection direction
        Vector3D in_light_direction = -(r.d).normalized();   // Incident light direction (ωi)
        Vector3D normal = its.normal.normalized(); // The normal (n)

        double scale_factor = 2.0 * dot(normal, in_light_direction);
        Vector3D reflection = scale_factor * normal;
        Vector3D ideal_reflection = reflection - in_light_direction; // Final computation of the ideal reflection direction (wr)

        // We build the reflection ray from the intersection point to the ideal reflection direction
        Ray reflection_ray(its.itsPoint, ideal_reflection, r.depth + 1);
        // Recursive Ray Tracing
        Color = computeColor(reflection_ray, objList, lsList);
    }
    // We use recursive ray tracing for refraction
    else if (its.shape->getMaterial().hasTransmission()) {
        //We compute the transmissive refaction wt following the formula:

        Vector3D normal = its.normal.normalized(); // The normal (n)
        Vector3D view_direction = -(r.d).normalized(); // viewing direction (wo)
        double refraction_ratio = 0.75;

        // If the ray is inside the object (we are going from medium µ2 to µ1) we need to flip the normal and invert the refraction ratio
        if (dot(view_direction, normal) < 0) {
            refraction_ratio = 1.0 / refraction_ratio;
            normal = -normal;
        }

        // Cosine of angle between incident ray and surface normal
        double cos_angle = dot(view_direction, normal);

        // Discriminant of the equation
        double discriminant = 1.0 - refraction_ratio * refraction_ratio * (1.0 - (cos_angle * cos_angle));

        Vector3D transmissive_refraction;

        if (discriminant >= 0) {

            double refract_term = sqrt(discriminant);
            Vector3D first_term = (refraction_ratio * view_direction);
            Vector3D second_term = refraction_ratio * cos_angle - refract_term;

            transmissive_refraction = -first_term + normal * second_term;
        }
        else {
            // Happens total internal reflection so the ray cannot exit the material
            // We compute specular reflection instead
            double scale_factor = 2.0 * dot(normal, view_direction);
            Vector3D reflection = scale_factor * normal;
            Vector3D transmissive_refraction = reflection - view_direction;
        }


        Ray refracted_ray(its.itsPoint, transmissive_refraction, r.depth + 1);

        // Recursively trace the refracted ray to get its color contribution
        Color = computeColor(refracted_ray, objList, lsList);
    }
        return Color;
}
