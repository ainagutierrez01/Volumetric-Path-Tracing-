#ifndef AREAINTEGRATOR_H
#define AREAINTEGRATOR_H

#include "shader.h"


class AreaIntegrator : public Shader
{
public:
    AreaIntegrator();
    AreaIntegrator(Vector3D Color_, Vector3D bgColor_);

    Vector3D computeColor(const Ray& r,
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList) const;

    Vector3D bgColor;
private:
    Vector3D Color;
};

#endif AREAINTEGRATOR_H