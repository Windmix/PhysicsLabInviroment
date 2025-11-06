#pragma once
#include "config.h"
#include "quad.h"

Quad::Quad(const glm::vec3& _v0, const glm::vec3& _v1, const glm::vec3& _v2, const glm::vec3& _v3, const glm::vec4& _color, const MathPlane _plane) : v0(_v0), v1(_v1), v2(_v2), v3(_v3), color(_color), plane(_plane)
{
	isHit = false;
	position = glm::vec3(0);
    rotation = glm::vec3(0);
}

Quad::~Quad()
{
}

void Quad::SetRotation(const glm::vec3 direction, float angle)
{
    // Update Euler angles
    auto rotMat = glm::rotate(angle, direction);
    glm::vec3 euler = glm::eulerAngles(glm::quat_cast(rotMat));
    rotation = glm::degrees(euler);

    // Compute quad center
    glm::vec3 center = (v0 + v1 + v2 + v3) * 0.25f;

    // Rotate each vertex around the center
    v0 = center + glm::vec3(rotMat * glm::vec4(v0 - center, 1.0f));
    v1 = center + glm::vec3(rotMat * glm::vec4(v1 - center, 1.0f));
    v2 = center + glm::vec3(rotMat * glm::vec4(v2 - center, 1.0f));
    v3 = center + glm::vec3(rotMat * glm::vec4(v3 - center, 1.0f));

    // Update the plane
    glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
    plane = MathPlane(normal, center);

    transform = glm::translate(glm::mat4(1.0f), position) * rotMat;

}

void Quad::SetPosition(const glm::vec3 pos)
{
    glm::vec3 delta = pos - position;
    position = pos;

    v0 += delta;
    v1 += delta;
    v2 += delta;
    v3 += delta;

     glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

     glm::vec3 center = (v0 + v1 + v2 + v3) * 0.25f;

     plane = MathPlane(normal, center);
    transform = glm::translate(glm::mat4(1.0f), position);
}

glm::vec3 Quad::GetPosition()
{
    return glm::vec3(transform[3]);
}
glm::vec3 Quad::GetRotation()
{
    return rotation;
}









