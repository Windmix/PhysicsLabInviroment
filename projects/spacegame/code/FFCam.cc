#include "config.h"
#include "FFCam.h"
#include "render/input/inputserver.h"
#include "render/cameramanager.h"
#include "render/debugrender.h"
#include "render/particlesystem.h"

using namespace Input;
using namespace glm;
using namespace Render;

namespace Game
{
    FFCam::FFCam()
    {
    
    }

    void
        FFCam::Update(float dt)
    {
        Mouse* mouse = Input::GetDefaultMouse();
        Keyboard* kbd = Input::GetDefaultKeyboard();

        Camera* cam = CameraManager::GetCamera(CAMERA_MAIN);

        // Handle mouse look when right mouse button is held
        if (mouse->held[Mouse::RightButton])
        {
            // Get mouse movement delta
            float mouseDeltaX = mouse->delta.x;
            float mouseDeltaY = mouse->delta.y;

            // Apply mouse sensitivity
            const float mouseSensitivity = 1.0f;
            float mouseRotX = -mouseDeltaY * mouseSensitivity * dt;
            float mouseRotY = -mouseDeltaX * mouseSensitivity * dt;

            // Use mouse input for rotation instead of keyboard
            rotXSmooth = mix(rotXSmooth, mouseRotY, dt * cameraSmoothFactor);
            rotYSmooth = mix(rotYSmooth, mouseRotX, dt * cameraSmoothFactor);
        }
        else
        {
            //Dampning
            float rotX = 0.0f;
            float rotY = 0.0f;

            const float rotationSpeed = 1.8f * dt;
            rotXSmooth = mix(rotXSmooth, rotX * rotationSpeed, dt * cameraSmoothFactor);
            rotYSmooth = mix(rotYSmooth, rotY * rotationSpeed, dt * cameraSmoothFactor);
        }

        // Handle movement
        vec3 desiredVelocity = vec3(0);
        vec3 copyVelocity = vec3(0);

    

        // Forward/Backward movement
        if (kbd->held[Key::W] || kbd->held[Key::S])
        {
            float moveDirectionFB = kbd->held[Key::W] ? 1.0f : -1.0f;
            float targetSpeedFB;

            if (kbd->held[Key::Shift])
                targetSpeedFB = this->boostSpeed * moveDirectionFB;
            else
                targetSpeedFB = this->normalSpeed * moveDirectionFB;

            this->currentSpeedFB = mix(this->currentSpeedFB, targetSpeedFB, std::min(1.0f, dt * (kbd->held[Key::Shift] ? 30.0f : 90.0f)));
            copyVelocity.z += this->currentSpeedFB;
        }

        // Left/Right movement  
        if (kbd->held[Key::A] || kbd->held[Key::D])
        {
            float moveDirectionLR = kbd->held[Key::A] ? 1.0f : -1.0f;
            float targetSpeedLR;

            if (kbd->held[Key::Shift])
                targetSpeedLR = this->boostSpeed * moveDirectionLR;
            else
                targetSpeedLR = this->normalSpeed * moveDirectionLR;

            this->currentSpeedLR = mix(this->currentSpeedLR, targetSpeedLR, std::min(1.0f, dt * (kbd->held[Key::Shift] ? 30.0f : 90.0f)));
            copyVelocity.x += this->currentSpeedLR;
        }

        // Up/Down movement
        if (kbd->held[Key::E] || kbd->held[Key::Q])
        {
            float moveDirectionTD = kbd->held[Key::E] ? 1.0f : -1.0f;
            float targetSpeedTD;

            if (kbd->held[Key::Shift])
                targetSpeedTD = this->boostSpeed * moveDirectionTD;
            else
                targetSpeedTD = this->normalSpeed * moveDirectionTD;

            this->currentSpeedTD = mix(this->currentSpeedTD, targetSpeedTD, std::min(1.0f, dt * (kbd->held[Key::Shift] ? 30.0f : 90.0f)));
            copyVelocity.y += this->currentSpeedTD;
        }

        // Only set speed to zero if no movement keys are pressed
        if (!(kbd->held[Key::W] || kbd->held[Key::S] || kbd->held[Key::A] || kbd->held[Key::D] || kbd->held[Key::E] || kbd->held[Key::Q]))
        {
            this->currentSpeedFB = 0.0f;
            this->currentSpeedLR = 0.0f;
            this->currentSpeedTD = 0.0f;
        }

        desiredVelocity = copyVelocity + desiredVelocity;
        desiredVelocity = this->transform * vec4(desiredVelocity, 0.0f);
        this->linearVelocity = mix(this->linearVelocity, desiredVelocity, dt * accelerationFactor);

        this->position += this->linearVelocity * dt * 20.0f;
        
     

        // Apply rotation
        quat localOrientation = quat(vec3(-rotYSmooth, rotXSmooth, rotZSmooth));
        this->orientation = this->orientation * localOrientation;
        this->rotationZ -= rotXSmooth;

        mat4 T = translate(this->position) * (mat4)this->orientation;
        this->transform = T * (mat4)quat(vec3(0, 0, rotationZ));
        this->rotationZ = mix(this->rotationZ, 0.0f, dt * cameraSmoothFactor);
 

        // Update camera view transform
        vec3 desiredCamPos = this->position + vec3(this->transform * vec4(0, camOffsetY, -4.0f, 0));
        this->camPos = mix(this->camPos, desiredCamPos, dt * cameraSmoothFactor);
        
        cam->view = lookAt(this->camPos, this->camPos + vec3(this->transform[2]), vec3(this->transform[1]));

        this->forward = position + vec3(transform[1]);
       

    }




}

MathRay Physics::ScreenPointToRay(glm::vec2& mousePos, float ScreenWidth, float ScreenHeight)
{
    MathRay ray;
    Render::Camera* camera = Render::CameraManager::GetCamera(CAMERA_MAIN);

    float clipX = (mousePos.x / ScreenWidth) * 2.0f - 1.0f;
    float clipY = 1.0f - (mousePos.y / ScreenHeight) * 2.0f; // flip Y

    // Create two points: one on near plane and one on far plane
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

bool Physics::CheckRayHit(Quad& myQuad, MathRay& ray, MathPlane& plane, RayProperties& rayproperties)
{
    // Create plane from quad points
             glm::vec3 QuadNormal = glm::cross((myQuad.v2 - myQuad.v0), (myQuad.v1 - myQuad.v0));
           

             // Get ray parameters
             glm::vec3 rayOrigin = ray.GetOrigin();
             glm::vec3 rayDirection = ray.GetDirection();
             float rayLength = ray.GetRayLength();

             // Use your MathPlane methods for intersection
             glm::vec3 planeNormal = plane.GetNormal();
             float planeD = plane.GetDistance();

             float denom = glm::dot(planeNormal, rayDirection);

             // Slightly larger epsilon for better stability
             if (fabs(denom) > 1e-6f) // not parallel to plane
             {
                 // Calculate intersection distance - THIS IS CORRECT!
                 float t = -(glm::dot(planeNormal, rayOrigin) + planeD) / denom;

                 // Add a small margin to handle edge cases
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

                     
                    
                     if (u >= -epsilon && v >= -epsilon && u <= 1.0f + epsilon && v <= 1.0f + epsilon)
                     {
                         // HIT!
                         float distance = plane.DistanceToPoint(rayproperties.intersection);
                         bool isOnPlane = plane.IsPointOnPlane(rayproperties.intersection);
                         bool isInFront = plane.IsPointInFront(rayproperties.intersection);

                         rayproperties.normalEnd = rayproperties.intersection + planeNormal * 2.0f;
                         printf("QUAD HIT! at (%.3f, %.3f, %.3f)\n", rayproperties.intersection.x, rayproperties.intersection.y, rayproperties.intersection.z);
                         printf("Distance to plane: %.6f, On plane: %s, In front: %s\n",
                             distance, isOnPlane ? "YES" : "NO", isInFront ? "YES" : "NO");
                         return true;
                     }
                     else
                     {
                         // MISS
                         float distance = plane.DistanceToPoint(rayproperties.intersection);
                         printf("MISS - Point is on plane but outside quad bounds\n");
                         printf("Intersection point: (%.6f, %.6f, %.6f)\n",
                             rayproperties.intersection.x, rayproperties.intersection.y, rayproperties.intersection.z);
                         printf("Barycentric coords: u=%.6f, v=%.6f\n", u, v);
                         printf("Distance to plane: %.10f\n", distance);
    
                     }
                   

                         
                    
                 }
             }
             return false;
}
