#ifndef HEMISPHERICALSHADER_H
#define HEMISPHERICALSHADER_H

#include "shader.h"


class HemisphericalShader : public Shader
{
public:
    HemisphericalShader();
    HemisphericalShader(Vector3D Color_, Vector3D bgColor_);

    Vector3D computeColor(const Ray& r,
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList) const;

    Vector3D bgColor;
private:
    Vector3D Color;
};

#endif HEMISPHERICALSHADER_H