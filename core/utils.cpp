#include "utils.h"

Utils::Utils()
{ }

double Utils::degreesToRadians(double degrees)
{
    return degrees * M_PI / 180.0;
}



bool Utils::hasIntersection(const Ray& cameraRay, const std::vector<Shape*>& objectsList) //or Shadow Ray
{

    // For each object on the scene...
    for(size_t objIndex = 0; objIndex < objectsList.size(); objIndex ++)
    {
          // Get the current object
          const Shape *obj = objectsList.at(objIndex);
          if (obj->rayIntersectP(cameraRay))
              return true;
    }   

    return false;
}



bool Utils::getClosestIntersection(const Ray& cameraRay, const std::vector<Shape*>& objectsList, Intersection& its) //or Closest Hit Ray
{
    //std::cout << "Need to implement the function Utils::getClosestIntersection() in the file utils.cpp" << std::endl;

    bool hasIntersection = false;

    for (size_t objIndex = 0; objIndex < objectsList.size(); objIndex++)
    {
        // Get the current object
        const Shape* obj = objectsList.at(objIndex);
        if (obj->rayIntersect(cameraRay, its))
            hasIntersection =  true;
    }
    return hasIntersection;
}

double interpolate(double val, double y0, double x0, double y1, double x1 )
{
    return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

double getRed(double value)
{
    if (value > 0.5)
        return interpolate( value, 0.0, 0.5, 1.0, 1.0 );
    else
        return 0.0;
}

double getGreen(double value)
{
    if (value < 0.25)
        return 0.0;
    else if (value < 0.5)
        return interpolate(value, 0.0, 0.25, 1, 0.5);
    else if (value < 0.75)
        return interpolate(value, 1, 0.5, 0, 0.75);
    else
        return 0;
}

double getBlue(double value)
{
    if (value < 0.5)
        return interpolate(value, 1.0, 0.0, 0.0, 0.5);
    else
        return 0;
}


Vector3D Utils::scalarToRGB(double scalar)
{
    return Vector3D( getRed(scalar),
                getGreen(scalar),
                getBlue(scalar) );
}

Vector3D Utils::computeReflectionDirection(const Vector3D &Direction, const Vector3D &normal)
{
    // Compute the perfect reflection direction FILL(...) [OPTIONAL]
    return 2.0 * dot(normal, Direction) * normal - Direction;
}

Vector3D Utils::computeRefractionDirection(const Vector3D& wi_in, const Vector3D& n_in, double eta)
{
    Vector3D n = n_in.normalized();          // make a modifiable copy
    Vector3D wi = wi_in.normalized();

    double cos_i = dot(-wi, n);
    bool entering = cos_i > 0.0;

    if (!entering) {
        // inside → outside: flip normal and invert ratio
        n = -n;
        eta = 1.0 / eta;
        cos_i = dot(-wi, n);
    }

    double k = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
    if (k < 0.0) {
        // total internal reflection
        return Utils::computeReflectionDirection(wi, n);
    }
    else {
        // refraction direction
        return eta * wi + (eta * cos_i - sqrt(k)) * n;
    }
}