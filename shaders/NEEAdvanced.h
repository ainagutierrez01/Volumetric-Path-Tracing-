#ifndef NEEADVANCED_H
#define NEEADVANCED_H

#include "shader.h"
#include "nee.h"
#include "../core/scene.h"
#include "../core/hemisphericalsampler.h"

class NEEAdvanced : public NEE
{
public:
    NEEAdvanced();
    NEEAdvanced(Vector3D Color_, Vector3D bgColor_);

    Vector3D computeColor(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const;
    virtual Vector3D computeRadiance(
        const Ray& r,
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList,
        const Intersection& its,
        const int maxDepth) const;

    virtual Vector3D reflectedRadiance(
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList,
        const Intersection& x,
        const Vector3D& w_o,
        int depth,
        const int maxDepth) const;

    virtual Vector3D directRadiance(
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList,
        const Intersection& x,
        const Vector3D& w_o) const;

    virtual Vector3D indirectRadiance(
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList,
        const Intersection& x,
        const Vector3D& w_o,
        int depth,
        const int maxDepth) const;

    Vector3D bgColor;
private:
    Vector3D Color;
};

#endif // NEEADVANCED_H