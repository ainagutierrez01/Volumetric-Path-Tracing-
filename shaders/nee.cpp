#include "nee.h"
#include "../materials/phong.h"
#include "../lightsources/pointlightsource.h"
#include "../lightsources/arealightsource.h"
#include "../core/hemisphericalsampler.h"
#include "../core/utils.h"
#include "../core/scene.h"
#include <math.h>

#define M_PI 3.1416

NEE::NEE() :
    Color(Vector3D(1, 0, 0))
{
}

NEE::NEE(Vector3D Color_, Vector3D bgColor_) :
    Shader(bgColor_), Color(Color_)
{
}

Vector3D NEE::computeRadiance(
    const Ray& r,
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList,
    const Intersection& its,
    const int maxDepth) const
{
    // Termination check
    if (r.depth >= maxDepth)
        return Vector3D(0.0);

    Vector3D Le(0.0), Lr(0.0);

    // Emitted radiance only if emissive surface
    if (its.shape->getMaterial().isEmissive()) {
        Le = its.shape->getMaterial().getEmissiveRadiance();
    }

    // Outgoing direction towards the camera
    Vector3D w_o = -r.d.normalized();

    // Reflected contribution
    Lr = reflectedRadiance(objList, lsList, its, w_o, r.depth, maxDepth);

    return Le + Lr;
}

Vector3D NEE::reflectedRadiance(
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList,
    const Intersection& x,
    const Vector3D& w_o,
    int depth,
    const int maxDepth) const
{
    Vector3D Ldir = directRadiance(objList, lsList, x, w_o);
    Vector3D Lind = indirectRadiance(objList, lsList, x, w_o, depth, maxDepth);
    return Ldir + Lind;
}


Vector3D NEE::directRadiance(
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList,
    const Intersection& x,
    const Vector3D& w_o) const
{
    Vector3D direct = Vector3D(0.0);
    Vector3D n = x.normal.normalized();
    int numLightSamples = 4; // Number of samples per light

    for (int lightId = 0; lightId < lsList.size(); lightId++) {
        double lightArea = lsList.at(lightId)->getArea();

        for (int sample = 0; sample < numLightSamples; sample++) {
            Vector3D p = lsList.at(lightId)->sampleLightPosition();
            double pdf = 1.0 / lightArea;

            Vector3D wi = (p - x.itsPoint).normalized();
            double distance = (p - x.itsPoint).length();
            double dSquared = distance * distance;

            // Shadow test
            Ray shadowRay(x.itsPoint, wi);
            shadowRay.maxT = distance - 1e-4;
            bool is_in_shadow = Utils::hasIntersection(shadowRay, objList);

            if (!is_in_shadow) {
                Vector3D radiance = lsList.at(lightId)->getIntensity();

                // BRDF direction (n, wi, wo)
                Vector3D fr = x.shape->getMaterial().getReflectance(n, wi, w_o);

                Vector3D ny = lsList.at(lightId)->getNormal().normalized();
                double g = (dot(n, wi) * dot(-wi, ny)) / dSquared;

                // Divide by number of samples
                direct += (radiance * fr * g) / (pdf * numLightSamples);
            }
        }
    }
    return direct;
}



Vector3D NEE::indirectRadiance(
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList,
    const Intersection& x,
    const Vector3D& w_o,
    int depth,
    const int maxDepth) const
{
    if (depth >= maxDepth)
        return Vector3D(0.0);

    Vector3D Lind(0.0);
    int numSamples = 4; // hemisphere samples

    Vector3D n = x.normal.normalized();
    HemisphericalSampler hS;

    for (int i = 0; i < numSamples; ++i) {
        // Sample random direction ωi over hemisphere
        Vector3D wi = hS.getSample(n);
        double pdf = 1.0 / (2.0 * M_PI); // uniform hemisphere PDF
        double cosTheta = std::max(0.0, dot(n, wi));

        // Create new ray
        Vector3D origin = x.itsPoint + n * 1e-4;
        Ray newRay(origin, wi, depth + 1);

        Intersection newIts;
        if (Utils::getClosestIntersection(newRay, objList, newIts)) {
            // Recursive reflected radiance at next bounce
            Vector3D Li = reflectedRadiance(objList, lsList, newIts, -wi, depth + 1, maxDepth);
            Vector3D brdf = x.shape->getMaterial().getReflectance(n, wi, w_o);

            // Lind = Li * BRDF * cosθ / pdf
            Lind += (Li * brdf * cosTheta) / (pdf * numSamples);
        }
    }

    return Lind;
}

Vector3D NEE::computeColor(
    const Ray& r,
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList) const
{
    Intersection its;
    if (!Utils::getClosestIntersection(r, objList, its))
        return bgColor;

    return computeRadiance(r, objList, lsList, its, 3);
}

