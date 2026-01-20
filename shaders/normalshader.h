#ifndef NORMALSHADER_H
#define NORMALSHADER_H

#include "shader.h"


class NormalShader : public Shader
{
public:
    NormalShader();
    NormalShader(Vector3D color_, Vector3D bgColor_);

    Vector3D computeColor(const Ray &r,
                             const std::vector<Shape*> &objList,
                             const std::vector<LightSource*> &lsList) const;

private:
    double maxDist;
    Vector3D color;
};

#endif // NORMALSHADER_H
