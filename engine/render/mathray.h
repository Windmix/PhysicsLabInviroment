#pragma once
#include <cmath>
#include <limits>

class MathRay
{
public:
    MathRay();
    MathRay(const glm::vec3& origin, glm::vec3& direction);
    ~MathRay();

    // Getters
    glm::vec3 GetOrigin() const { return m_origin; }
    glm::vec3 GetDirection() const { return m_direction; }
    float GetRayLength() const { return m_rayLength; }

    // Setters
    void SetOrigin(const  glm::vec3& origin) { m_origin = origin; }
    void SetDirection(const  glm::vec3& direction) { m_direction = direction; }
    void SetRayLength(const float raylength) { m_rayLength = raylength; }

  
private:
    glm::vec3 m_origin;
    glm::vec3 m_direction;
    float m_rayLength;
};

