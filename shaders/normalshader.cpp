#include "normalshader.h"

#include "../core/utils.h"

NormalShader::NormalShader() :
    color(Vector3D(1, 0, 0))
{ }

NormalShader::NormalShader(Vector3D hitColor_, Vector3D bgColor_) :
    Shader(bgColor_), color(hitColor_)
{ }

Vector3D NormalShader::computeColor(const Ray &r, const std::vector<Shape*> &objList, const std::vector<LightSource*> &lsList) const
{
    Intersection its;

    if (Utils::getClosestIntersection(r, objList, its)) {

        return (its.normal + Vector3D(1, 1, 1)) * 0.5;
    }

    return bgColor;
}
