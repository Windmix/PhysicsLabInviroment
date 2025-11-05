#pragma once
#include "config.h"
#include "quad.h"

Quad::Quad(const glm::vec3& _v0, const glm::vec3& _v1, const glm::vec3& _v2, const glm::vec3& _v3, const glm::vec4& _color) : v0(_v0), v1(_v1), v2(_v2), v3(_v3), color(_color)
{

}

Quad::~Quad()
{
}








