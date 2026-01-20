#ifndef COLOREDVOLUMETRICPATHSHADER_H
#define COLOREDVOLUMETRICPATHSHADER_H

#include "shader.h"

class Colored : public Shader
{
public:
    Colored();
    Colored(Vector3D hitColor, Vector3D bgColor_);

    Vector3D DirectRadiance(const Intersection& its, const Vector3D& wo, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const;
    Vector3D IndirectRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const size_t maxDepth = 3) const;
    Vector3D ReflectedRadiance(const Intersection& its, const Vector3D& wo, const size_t depth, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList) const;
    Vector3D computeInscattering(const Ray& r, double t_max,
        const std::vector<Shape*>& objList,
        const std::vector<LightSource*>& lsList,
        int depth) const;

    Vector3D computeRadiance(const Ray& r, const std::vector<Shape*>& objList, const std::vector<LightSource*>& lsList, const int MaxDepth = 3) const;

    Vector3D computeColor(const Ray &r,
                             const std::vector<Shape*> &objList,
                             const std::vector<LightSource*> &lsList) const;

    Vector3D hitColor;
};

#endif //Colored
