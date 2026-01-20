#include "arealightsource.h"

AreaLightSource::AreaLightSource(Square* areaLightsource_) :
    myAreaLightsource(areaLightsource_)
{ }



Vector3D AreaLightSource::getIntensity() const
{
    return myAreaLightsource->getMaterial().getEmissiveRadiance();
}

Vector3D AreaLightSource::sampleLightPosition() const
{
    // Get two random numbers between 0-1
    double u = (double)std::rand() / RAND_MAX;
    double v = (double)std::rand() / RAND_MAX;

    // Map u,v from [0,1] to [-0.5,0.5] to center around the square's position
    u = u - 0.5;
    v = v - 0.5;

    // Calculate the sampled point on the square
    Vector3D sampledPoint = myAreaLightsource->corner + (u * myAreaLightsource->v1) + (v * myAreaLightsource->v2);

    return sampledPoint;
}

