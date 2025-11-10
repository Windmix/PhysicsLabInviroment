#pragma once
#include <glm.hpp>
#include <vector>
#include <GL/glew.h>
#include "render/mathplane.h"
#include "betterphysics.h"

class Quad
{
public:
    glm::vec3 v0, v1, v2, v3; // quad vertices

    std::vector<Physics::ColliderMesh::Triangle> triangles;
    MathPlane plane;
    Physics::AABB aabb;
   
    glm::mat4 transform;



    Quad(const glm::vec3& _v0, const glm::vec3& _v1, const glm::vec3& _v2, const glm::vec3& _v3, const glm::vec4 _color, const MathPlane _plane = MathPlane());

    ~Quad();
    void SetRotation(const glm::vec3 direction, float angle);
    void SetPosition(const glm::vec3 pos);

    glm::vec3 GetPosition();
    glm::vec3 GetRotation();
    void UpdateAndDrawAABB();

   
    bool CheckRayHit(Quad& myQuad, MathRay& ray, Physics::RayProperties& rayproperties);
private:
    void DrawAABB();

    glm::vec3 position; // center position
    glm::vec3 rotation;//  center rotation

};
   
