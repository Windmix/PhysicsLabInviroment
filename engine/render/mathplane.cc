#include <config.h>
#include "mathplane.h"

MathPlane::MathPlane() : m_normal(0, 1, 0), m_d(0) {}

MathPlane::MathPlane(float a, float b, float c, float d)
    : m_normal(a, b, c), m_d(d)
{
    Normalize();
}

MathPlane::MathPlane(const glm::vec3& normal, const glm::vec3& point)
{
    SetFromNormalAndPoint(normal, point);
}

MathPlane::MathPlane(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3)
{
    SetFromPoints(p1, p2, p3);
}

MathPlane::~MathPlane() {}

float MathPlane::DistanceToPoint(const glm::vec3& point) const
{
    return glm::dot(m_normal, point) + m_d;
}

bool MathPlane::IsPointOnPlane(const glm::vec3& point, float epsilon) const
{
    return fabs(DistanceToPoint(point)) <= epsilon;
}

bool MathPlane::IsPointInFront(const glm::vec3& point) const
{
    return DistanceToPoint(point) > 0.0f;
}

void MathPlane::Normalize()
{
    float len = glm::length(m_normal);
    if (len > 0.0f)
    {
        m_normal /= len;
        m_d /= len;
    }
}

void MathPlane::SetFromPoints(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3)
{
    glm::vec3 edge1 = p2 - p1;
    glm::vec3 edge2 = p3 - p1;
    m_normal = glm::normalize(glm::cross(edge1, edge2));
    m_d = -glm::dot(m_normal, p1);
}

void MathPlane::SetFromNormalAndPoint(const glm::vec3& normal, const glm::vec3& point)
{
    m_normal = glm::normalize(normal);
    m_d = -glm::dot(m_normal, point);
}





