#pragma once
#include "config.h"
#include "quad.h"
#include "debugrender.h"

Quad::Quad(const glm::vec3& _v0, const glm::vec3& _v1, const glm::vec3& _v2, const glm::vec3& _v3, const glm::vec4 _color, const MathPlane _plane) : v0(_v0), v1(_v1), v2(_v2), v3(_v3), plane(_plane)
{
   

    auto triangle = Physics::ColliderMesh::Triangle();

    //first tirangle
    triangle.verticies[0] = v0;
    triangle.verticies[1] = v1;
    triangle.verticies[2] = v2;
    triangle.normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
    triangle.selectedColor = glm::vec4(1, 1, 0, 1);
    triangle.og_color = _color;
    triangle.color = _color;

    triangles.push_back(triangle);
    

   
    auto triangle1 = Physics::ColliderMesh::Triangle();
    //second triangle
    triangle1.verticies[0] = v0;
    triangle1.verticies[1] = v2;
    triangle1.verticies[2] = v3;
    triangle1.normal = glm::normalize(glm::cross(v2 - v0, v3 - v0));

    triangle1.selectedColor = glm::vec4(1, 1, 0, 1);
    triangle1.og_color = _color;
    triangle1.color = _color;

    triangles.push_back(triangle1);

	position = glm::vec3(0);
    rotation = glm::vec3(0);
    aabb = Physics::AABB();
  


    UpdateAndDrawAABB();
}

Quad::~Quad()
{
}



void Quad::SetRotation(const glm::vec3 direction, float angle)
{
    float angleRad = glm::radians(angle);

    // Compute quad center
    glm::vec3 center = (v0 + v1 + v2 + v3) * 0.25f;

    // Build rotation matrix around center
    glm::mat4 rotMat = glm::translate(glm::mat4(1.0f), center) *
        glm::rotate(glm::mat4(1.0f), angleRad, direction) *
        glm::translate(glm::mat4(1.0f), -center);

    // Rotate vertices
    v0 = glm::vec3(rotMat * glm::vec4(v0, 1.0f));
    v1 = glm::vec3(rotMat * glm::vec4(v1, 1.0f));
    v2 = glm::vec3(rotMat * glm::vec4(v2, 1.0f));
    v3 = glm::vec3(rotMat * glm::vec4(v3, 1.0f));

    // Update two triangles
    triangles[0].verticies[0] = v0;
    triangles[0].verticies[1] = v1;
    triangles[0].verticies[2] = v2;
    triangles[0].normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

    triangles[1].verticies[0] = v0;
    triangles[1].verticies[1] = v2;
    triangles[1].verticies[2] = v3;
    triangles[1].normal = glm::normalize(glm::cross(v2 - v0, v3 - v0));

    // Update rotation Euler angles
    glm::mat3 rot3x3(rotMat);
    glm::vec3 euler = glm::eulerAngles(glm::quat_cast(rot3x3));
    rotation = glm::degrees(euler);

    // Update plane
    glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
    plane = MathPlane(normal, center);

    // Update transform
    transform = rotMat;
}

void Quad::SetPosition(const glm::vec3 pos)
{
    glm::vec3 delta = pos - position;
    position = pos;

    v0 += delta;
    v1 += delta;
    v2 += delta;
    v3 += delta;

    // Update two triangles
    triangles[0].verticies[0] = v0;
    triangles[0].verticies[1] = v1;
    triangles[0].verticies[2] = v2;
    triangles[0].normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

    triangles[1].verticies[0] = v0;
    triangles[1].verticies[1] = v2;
    triangles[1].verticies[2] = v3;
    triangles[1].normal = glm::normalize(glm::cross(v2 - v0, v3 - v0));

    // Update plane
    glm::vec3 center = (v0 + v1 + v2 + v3) * 0.25f;
    glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
    plane = MathPlane(normal, center);

    // Update transform
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
void Quad::UpdateAndDrawAABB()
{

    // Start with extreme values
    aabb.min = glm::vec3(std::numeric_limits<float>::max());
    aabb.max = glm::vec3(std::numeric_limits<float>::lowest());

    // Expand AABB using triangle vertices
    for (int t = 0; t < 2; ++t) // or 2 if you have two triangles
    {
        for (int v = 0; v < 3; ++v)
        {
            aabb.Expand(triangles[t].verticies[v]);
        }
    }

    // Draw the box
    DrawAABB();
}

bool Quad::CheckRayHit(Quad& myQuad, MathRay& ray, Physics::RayProperties& rayproperties)
{
    bool hitAny = false;

    for (int i = 0; i < myQuad.triangles.size(); i++)
    {
        auto& tri = myQuad.triangles[i];

        // Ray-Triangle intersection using Möller–Trumbore
        const glm::vec3& v0 = tri.verticies[0];
        const glm::vec3& v1 = tri.verticies[1];
        const glm::vec3& v2 = tri.verticies[2];

        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;
        glm::vec3 h = glm::cross(ray.GetDirection(), edge2);
        float a = glm::dot(edge1, h);

        if (fabs(a) < 1e-6f) continue; // parallel, skip

        float f = 1.0f / a;
        glm::vec3 s = ray.GetOrigin() - v0;
        float u = f * glm::dot(s, h);
        if (u < 0.0f || u > 1.0f) continue;

        glm::vec3 q = glm::cross(s, edge1);
        float v = f * glm::dot(ray.GetDirection(), q);
        if (v < 0.0f || u + v > 1.0f) continue;

        float t = f * glm::dot(edge2, q);
        if (t > 1e-6f && t <= ray.GetRayLength())
        {
            // HIT
            rayproperties.intersection = ray.GetOrigin() + ray.GetDirection() * t;
            rayproperties.normalEnd = rayproperties.intersection + tri.normal;

            // Make this triangle glow
            tri.SetSelected(true);
            hitAny = true;
        }
    
    }

    return hitAny;

}


void Quad::DrawAABB()
{
    glm::vec3 min = aabb.min;
    glm::vec3 max = aabb.max;

    glm::vec3 corners[8] = 
    {
        {min.x, min.y, min.z},
        {max.x, min.y, min.z},
        {max.x, max.y, min.z},
        {min.x, max.y, min.z},
        {min.x, min.y, max.z},
        {max.x, min.y, max.z},
        {max.x, max.y, max.z},
        {min.x, max.y, max.z}
    };

    // Draw edges (12 lines)
    Debug::DrawLine(corners[0], corners[1], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[1], corners[2], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[2], corners[3], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[3], corners[0], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));

    Debug::DrawLine(corners[4], corners[5], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[5], corners[6], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[6], corners[7], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[7], corners[4], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));

    Debug::DrawLine(corners[0], corners[4], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[1], corners[5], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[2], corners[6], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[3], corners[7], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
}










