#pragma once
#include <glm.hpp>
#include <vector>
#include <GL/glew.h>
#include "render/mathplane.h"

class Quad
{
public:
    glm::vec3 v0, v1, v2, v3; // quad vertices
    glm::vec4 color;

    Quad(const glm::vec3& _v0, const glm::vec3& _v1,
        const glm::vec3& _v2, const glm::vec3& _v3,
        const glm::vec4& _color = glm::vec4(1.0f));
    ~Quad();



private:


};
