#ifndef PATHTRACER_H
#define PATHTRACER_H

#include "shader.h"
#include "../core/scene.h"
#include "../core/hemisphericalsampler.h"

class PathTracer : public Shader
{
public:
    PathTracer();
    PathTracer(Vector3D Color_, Vector3D bgColor_);

    Vector3D computeColor(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const;
    virtual Vector3D computeRadiance(const Ray& r,
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList,
        const int maxDepth) const;

    Vector3D bgColor;
private:
    Vector3D Color;
};

#endif // PATHTRACER_H