#include "shader.h"
#include "whittedshader.h"
#include "../materials/phong.h"
#include "../lightsources/pointlightsource.h"
#include "../core/utils.h"

WhittedShader::WhittedShader():
    WgColor(Vector3D(1, 0, 0))
{ }

WhittedShader::WhittedShader(Vector3D wgColor_, Vector3D bgColor_):
    Shader(bgColor_),WgColor(wgColor_)
{ }

Vector3D WhittedShader::computeColor(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const
{

    Intersection its;
    
    if (Utils::getClosestIntersection(r, objList, its)) {

        Vector3D Color = Vector3D(0, 0, 0);

        //5.2
        if (its.shape->getMaterial().hasDiffuseOrGlossy()) {
            const Vector3D ambient_constant = Vector3D(0.25, 0.25, 0.25); 
            Vector3D diffuse_coefficient = its.shape->getMaterial().getDiffuseReflectance();

            //We initialize the luminance with the ambient constant and the diffuse coefficient
            Color += ambient_constant * diffuse_coefficient;

            //For the summation part of the equation
            for (const auto& lightSource : lsList) {

                // We calculate the reflectance, to do so, we need:
                Vector3D normal = its.normal.normalized(); // The normal (n)
                Vector3D view_direction = -r.d; // The viewing direction (wo)
                Vector3D in_light_direction = (lightSource->sampleLightPosition() - its.itsPoint).normalized(); //the incident light direction

                Vector3D reflectance = its.shape->getMaterial().getReflectance(normal, view_direction, in_light_direction);

                // To calculate the visibility term Vi(x):
                // We construct a shadow ray
                Ray shadow_ray(its.itsPoint, in_light_direction);
                shadow_ray.maxT = (shadow_ray.o - lightSource->sampleLightPosition()).length();

                // If the shadow ray hits another object before reaching the light, the point is in a shadow.
                bool is_in_shadow = Utils::hasIntersection(shadow_ray, objList);

                int visibility_factor;
                if (is_in_shadow) {
                    visibility_factor = 0;
                }
                else {
                    visibility_factor = 1;
                }

                // We compute the final term:
                Color += lightSource->getIntensity() * reflectance * dot(normal, in_light_direction) * visibility_factor;
            }
        }

        // 5.3
        else if (its.shape->getMaterial().hasSpecular()) {
            // We calculate the ideal reflection direction
            Vector3D in_light_direction = -(r.d).normalized();   // Incident light direction (ωi)
            Vector3D normal = its.normal.normalized(); // The normal (n

            double scale_factor = 2.0 * dot(normal, in_light_direction);
            Vector3D reflection = scale_factor * normal;
            Vector3D ideal_reflection = reflection - in_light_direction; // Final computation of the ideal reflection direction (wr)

            // We build the reflection ray from the intersection point to the ideal reflection direction
            Ray reflection_ray(its.itsPoint, ideal_reflection, r.depth);
            // Recursive Ray Tracing
            Color = computeColor(reflection_ray, objList, lsList);
            
        }

        // 5.4
        else if (its.shape->getMaterial().hasTransmission()) {
            
            //We compute the transmissive refaction wt following the formula:
            
            Vector3D normal = its.normal.normalized(); // The normal (n)
            Vector3D view_direction = - (r.d).normalized(); // viewing direction (wo)
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
                Vector3D second_term = refraction_ratio * cos_angle - refract_term ;
                
                transmissive_refraction = -first_term + normal * second_term;
            }
            else {
                // Happens total internal reflection so the ray cannot exit the material
                // We compute specular reflection instead
                double scale_factor = 2.0 * dot(normal, view_direction);
                Vector3D reflection = scale_factor * normal;
                Vector3D transmissive_refraction = reflection - view_direction;
            }

       
            Ray refracted_ray(its.itsPoint, transmissive_refraction, r.depth );

            // Recursively trace the refracted ray to get its color contribution
            Color += computeColor(refracted_ray, objList, lsList);
        }

        return Color;
    } 
    
    
    else {
        return bgColor;
    }
}
