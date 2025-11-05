#pragma once
#include <cmath>
#include <glm.hpp>
#include <gtc/epsilon.hpp>




class MathPlane
{
public:
    MathPlane();
    MathPlane(float a, float b, float c, float d);
    MathPlane(const glm::vec3& normal, const glm::vec3& point);
    MathPlane(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3);
    ~MathPlane();

    glm::vec3 GetNormal() const { return m_normal; }
    glm::vec4 GetEquation() const { return glm::vec4(m_normal, m_d); }
    float GetDistance() const { return m_d; }

    float DistanceToPoint(const glm::vec3& point) const;
    bool IsPointOnPlane(const glm::vec3& point, float epsilon = 1e-6f) const;
    bool IsPointInFront(const glm::vec3& point) const;

    void Normalize();
    void SetFromPoints(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3);
    void SetFromNormalAndPoint(const glm::vec3& normal, const glm::vec3& point);

private:
    glm::vec3 m_normal;
    float m_d;
};


