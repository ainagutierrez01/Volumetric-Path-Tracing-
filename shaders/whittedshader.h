#ifndef WHITTEDSHADER_H
#define WHITTEDSHADER_H

#include "shader.h"

//to define the class WhittedShader
class WhittedShader : public Shader
{
public:
    WhittedShader();
    WhittedShader(Vector3D wgColor_, Vector3D bgColor_);

    Vector3D computeColor(const Ray& r,
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList) const;

    Vector3D WgColor;
private:
    Vector3D color;
};


#endif // WhittedShader_H