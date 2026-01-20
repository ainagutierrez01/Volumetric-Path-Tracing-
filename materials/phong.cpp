#include "phong.h"

#include <iostream>
#define M_PI 3.1416

Phong::Phong()
{ }

Phong::Phong(Vector3D Kd_, Vector3D Ks_, float alpha_):
rho_d(Kd_), Ks(Ks_), alpha(alpha_){}


Vector3D Phong::getReflectance(const Vector3D& n, const Vector3D& wo,
    const Vector3D& wi) const {

    //We calculate the ideal reflection direction
    Vector3D wr = 2.0 * dot(n, wi) * n - wi;

    //First term of the reflectance
    Vector3D diffuse = (rho_d / M_PI); //A1

    //Second term of the reflectance
    double scale_factor = ((alpha + 2) / (2 * M_PI));
    Vector3D ref = scale_factor * (Ks * pow(dot(wo, wr), alpha));

    return Vector3D (diffuse + ref);

};

double Phong::getIndexOfRefraction() const
{
    std::cout << "Warning! Calling \"Material::getIndexOfRefraction()\" for a non-transmissive material"
              << std::endl;

    return -1;
}


Vector3D Phong::getEmissiveRadiance() const
{
    return Vector3D(0.0);
}


Vector3D Phong::getDiffuseReflectance() const
{
    return rho_d;
}

