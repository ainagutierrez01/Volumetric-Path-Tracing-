#include "depthshader.h"

#include "../core/utils.h"

DepthShader::DepthShader() :
    color(Vector3D(1, 0, 0))
{ }

DepthShader::DepthShader(Vector3D hitColor_, double maxDist_, Vector3D bgColor_) :
    Shader(bgColor_), maxDist(maxDist_), color(hitColor_)
{ }

Vector3D DepthShader::computeColor(const Ray &r, const std::vector<Shape*> &objList, const std::vector<LightSource*> &lsList) const
{
    Intersection its;

    if (Utils::getClosestIntersection(r, objList, its)) {
        Vector3D origin = r.o;
        Vector3D dest = its.itsPoint;

        float dist = (origin - dest).length();
        float dist_n = dist / maxDist;

        return Vector3D(0, 1 - dist_n, 0);
    }


    return bgColor;
}
