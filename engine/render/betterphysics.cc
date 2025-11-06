#include <config.h>
#include "betterphysics.h"



MathRay Physics::ScreenPointToRay(glm::vec2& mousePos, float ScreenWidth, float ScreenHeight)
{
    MathRay ray;
    Render::Camera* camera = Render::CameraManager::GetCamera(CAMERA_MAIN);

    float clipX = (mousePos.x / ScreenWidth) * 2.0f - 1.0f;
    float clipY = 1.0f - (mousePos.y / ScreenHeight) * 2.0f; // flip Y


    glm::vec4 clipCoordNear(clipX, clipY, -1.0f, 1.0f);  // Near plane
    glm::vec4 clipCoordFar(clipX, clipY, 1.0f, 1.0f);    // Far plane

    // Clip → Eye space
    glm::mat4 invProj = glm::inverse(camera->projection);
    glm::vec4 eyeCoordNear = invProj * clipCoordNear;
    glm::vec4 eyeCoordFar = invProj * clipCoordFar;

    // Perspective divide
    eyeCoordNear /= eyeCoordNear.w;
    eyeCoordFar /= eyeCoordFar.w;

    // Eye → World space
    glm::mat4 invView = glm::inverse(camera->view);
    glm::vec4 worldCoordNear = invView * eyeCoordNear;
    glm::vec4 worldCoordFar = invView * eyeCoordFar;

    // Ray origin is the camera position (extracted from inverse view matrix)
    glm::vec3 origin = glm::vec3(invView[3]);

    // Direction is from near point to far point (through the mouse position)
    glm::vec3 direction = glm::normalize(glm::vec3(worldCoordFar) - glm::vec3(worldCoordNear));

    ray.SetOrigin(origin);
    ray.SetDirection(direction);
    ray.SetRayLength(1000.0f);
    return ray;
}

bool Physics::CheckRayHit(Quad& myQuad, MathRay& ray, RayProperties& rayproperties)
{
    //normal on quad
    glm::vec3 QuadNormal = glm::cross((myQuad.v2 - myQuad.v0), (myQuad.v1 - myQuad.v0));


    // Get ray parameters
    glm::vec3 rayOrigin = ray.GetOrigin();
    glm::vec3 rayDirection = ray.GetDirection();
    float rayLength = ray.GetRayLength();

    //MathPlane methods for intersection
    glm::vec3 planeNormal = myQuad.plane.GetNormal();
    float planeD = myQuad.plane.GetDistance();

    float denom = glm::dot(planeNormal, rayDirection);

    // Slightly larger epsilon for better stability
    if (fabs(denom) > 1e-6f) // not parallel to plane
    {
        // Calculate intersection distance 
        float t = -(glm::dot(planeNormal, rayOrigin) + planeD) / denom;

        //a small margin to handle edge cases
        if (t >= -1e-6f && t <= rayLength + 1e-6f)
        {
            rayproperties.intersection = rayOrigin + rayDirection * t;

            // Check if intersection point is inside the quad
            glm::vec3 v0 = myQuad.v1 - myQuad.v0;
            glm::vec3 v1 = myQuad.v3 - myQuad.v0;
            glm::vec3 v2 = rayproperties.intersection - myQuad.v0;

            float dot00 = glm::dot(v0, v0);
            float dot01 = glm::dot(v0, v1);
            float dot02 = glm::dot(v0, v2);
            float dot11 = glm::dot(v1, v1);
            float dot12 = glm::dot(v1, v2);

            float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
            float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

            // Allow small tolerance for edge hits
            const float epsilon = 1e-5f;


            bool isOnPlane = myQuad.plane.IsPointOnPlane(rayproperties.intersection);
            float distance = myQuad.plane.DistanceToPoint(rayproperties.intersection);
            if (u >= -epsilon && v >= -epsilon && u <= 1.0f + epsilon && v <= 1.0f + epsilon)
            {
                // HIT!

                printf("QUAD HIT! at (%.3f, %.3f, %.3f)\n", rayproperties.intersection.x, rayproperties.intersection.y, rayproperties.intersection.z);
                rayproperties.normalEnd = rayproperties.intersection + planeNormal;

                return true;
            }
            else
            {
                // MISS

                printf("MISS");
                printf("Intersection point: (%.6f, %.6f, %.6f)\n", rayproperties.intersection.x, rayproperties.intersection.y, rayproperties.intersection.z);
                printf("Barycentric coords: u=%.6f, v=%.6f\n", u, v);


            }

            printf("Distance to plane: %.6f, On plane: %s\n",
                distance, isOnPlane ? "YES" : "NO");




        }
    }
    return false;
}

inline void Physics::AABB::Expand(const glm::vec3& point)
{
    min = glm::min(min, point);
    max = glm::max(max, point);
}
Physics::AABB::AABB() : min(std::numeric_limits<float>::max()),max(std::numeric_limits<float>::lowest())
{
    model = Render::LoadModel("assets/space/Cube.glb");
}

glm::vec3 Physics::AABB::GetABBCenter() const
{
    return (min + max) * 0.5f;
}

glm::vec3 Physics::AABB::GetABBSize() const
{
    return max - min;
}

void Physics::AABB::drawBox(glm::mat4 transform)
{
    Render::RenderDevice::Draw(model, transform);
}
