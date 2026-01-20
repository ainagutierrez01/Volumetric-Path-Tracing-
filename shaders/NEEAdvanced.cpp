#include "NEEAdvanced.h"
#include "../materials/phong.h"
#include "../lightsources/pointlightsource.h"
#include "../lightsources/arealightsource.h"
#include "../core/hemisphericalsampler.h"
#include "../core/utils.h"
#include "../core/scene.h"
#include <math.h>

#define M_PI 3.1416

NEEAdvanced::NEEAdvanced()
    : Color(Vector3D(1, 0, 0))
{
}

NEEAdvanced::NEEAdvanced(Vector3D Color_, Vector3D bgColor_)
    : NEE(Color_, bgColor_)
{
}

Vector3D NEEAdvanced::computeRadiance(
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

    // We handle emissive surfaces
    // Firstly, primary or specular/transmissive rays can see lights
    if (its.shape->getMaterial().isEmissive()) {
        if (r.depth == 0 || r.type == SPECULAR || r.type == TRANSMISSIVE) {
            Le = its.shape->getMaterial().getEmissiveRadiance();
        }
    }

    // Outgoing direction towards the camera
    Vector3D w_o = -r.d.normalized();

    // Reflected/indirect contribution
    Lr = reflectedRadiance(objList, lsList, its, w_o, r.depth, maxDepth);

    return Le + Lr;
}


Vector3D NEEAdvanced::reflectedRadiance(
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList,
    const Intersection& x,
    const Vector3D& w_o,
    int depth,
    const int maxDepth) const
{
    Vector3D Ldir(0.0);

    if (x.shape->getMaterial().hasDiffuseOrGlossy()) {
        Ldir = directRadiance(objList, lsList, x, w_o);
    }

    Vector3D Lind = indirectRadiance(objList, lsList, x, w_o, depth, maxDepth);
    return Ldir + Lind;
}


Vector3D NEEAdvanced::directRadiance(
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList,
    const Intersection& x,
    const Vector3D& w_o) const
{
    Vector3D direct = Vector3D(0.0);
    Vector3D n = x.normal.normalized();
    if (dot(n, w_o) < 0.0) n = -n;
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
            shadowRay.maxT = distance - 1e-3;
            bool is_in_shadow = Utils::hasIntersection(shadowRay, objList);

            if (!is_in_shadow) {
                Vector3D radiance = lsList.at(lightId)->getIntensity();

                // BRDF direction (n, wi, wo)
                Vector3D fr = x.shape->getMaterial().getReflectance(n, wi, w_o);

                Vector3D lightNormal = lsList.at(lightId)->getNormal().normalized();
                double cosThetaX = dot(n, wi);
                double cosThetaY = dot(lightNormal, -wi);
                double G = (cosThetaX * cosThetaY) / dSquared;

                // Divide by number of samples
                direct += (radiance * fr * G) / (pdf * numLightSamples);
            }
        }
    }
    return direct;
}


Vector3D NEEAdvanced::indirectRadiance(
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList,
    const Intersection& x,
    const Vector3D& w_o,
    int depth,
    const int maxDepth) const
{
    if (depth >= maxDepth)
        return Vector3D(0.0);

    Vector3D n = x.normal.normalized();
    if (dot(n, w_o) < 0.0) n = -n;

    // We handle perfect mirror
    if (x.shape->getMaterial().hasSpecular()) {
        Vector3D wi = Utils::computeReflectionDirection(w_o, n);
        Ray reflectedRay(x.itsPoint - n * 1e-5, wi, depth + 1, SPECULAR);
        reflectedRay.minT = 1e-3;

        Intersection newIts;
        if (Utils::getClosestIntersection(reflectedRay, objList, newIts)) {
            Vector3D Li = computeRadiance(reflectedRay, objList, lsList, newIts, maxDepth);
            return Li; // Mirror perfectly reflects one direction
        }
        return Vector3D(0.0);
    }

    // We handle transmissive (refraction)
    if (x.shape->getMaterial().hasTransmission()) {
        double coeficcient = 0.75;
        if (dot(n, w_o) < 0.0) coeficcient = 1/0.75;

        Vector3D wi = Utils::computeRefractionDirection(-w_o, n, coeficcient); // η = n1/n2
        Ray refractedRay(x.itsPoint - n * 1e-3, wi, depth + 1, TRANSMISSIVE);

        Intersection newIts;
        if (Utils::getClosestIntersection(refractedRay, objList, newIts)) {
            Vector3D Li = computeRadiance(refractedRay, objList, lsList, newIts, maxDepth);
            return Li;
        }
        return Vector3D(0.0);
    }

    // We handle diffuse or glossy
    if (x.shape->getMaterial().hasDiffuseOrGlossy()) {
        Vector3D Lind(0.0);
        int numSamples = 4;
        HemisphericalSampler hS;
        for (int i = 0; i < numSamples; ++i) {
            Vector3D wi = hS.getSample(n);
            double pdf = 1.0 / (2.0 * M_PI);
            double cosTheta = std::max(0.0, dot(n, wi));

            Vector3D origin = x.itsPoint + n * 1e-3;
            Ray newRay(origin, wi, depth + 1, DIFFUSE);

            Intersection newIts;
            if (Utils::getClosestIntersection(newRay, objList, newIts)) {
                Vector3D Li = reflectedRadiance(objList, lsList, newIts, -wi, depth + 1, maxDepth);
                Vector3D brdf = x.shape->getMaterial().getReflectance(n, wi, w_o);
                Lind += (Li * brdf * cosTheta) / (pdf * numSamples);
            }
        }

        return Lind;
    }
    return Vector3D(0.0);
}

Vector3D NEEAdvanced::computeColor(
    const Ray& r,
    const std::vector<Shape*>& objList,
    const std::vector<LightSource*>& lsList) const
{
    Intersection its;
    if (!Utils::getClosestIntersection(r, objList, its))
        return bgColor;

    return computeRadiance(r, objList, lsList, its, 3);
}
