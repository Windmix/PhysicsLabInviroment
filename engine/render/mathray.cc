#include <config.h>
#include "mathray.h"


MathRay::MathRay()
{
}

MathRay::MathRay(const glm::vec3& origin, glm::vec3& direction) : m_origin(origin), m_direction(glm::normalize(direction))
{
}

MathRay::~MathRay()
{
}


 