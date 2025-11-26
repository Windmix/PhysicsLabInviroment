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
            const float mouseSensitivity = 0.0020f;
            float targetYaw = -mouseDeltaX * mouseSensitivity;
            float targetPitch = mouseDeltaY * mouseSensitivity;

            // Use mouse input for rotation instead of keyboard
            rotXSmooth = mix(rotXSmooth, targetYaw, dt * cameraSmoothFactor);
            rotYSmooth = mix(rotYSmooth, targetPitch, dt * cameraSmoothFactor);

        }
        else
        {
            rotXSmooth = mix(rotXSmooth, 0.0f, dt * cameraSmoothFactor);
            rotYSmooth = mix(rotYSmooth, 0.0f, dt * cameraSmoothFactor);
        }

        yaw += rotXSmooth;
        pitch += rotYSmooth;

        //vec3 desiredVelocity = vec3(0);
        //vec3 copyVelocity = vec3(0);

    
        // Clamp pitch to avoid flipping
        const float pitchLimit = radians(89.0f);
        pitch = clamp(pitch, -pitchLimit, pitchLimit);

        // Build final orientation = yaw * pitch
        quat qPitch = angleAxis(pitch, vec3(1, 0, 0));
        quat qYaw = angleAxis(yaw, vec3(0, 1, 0));

        orientation = qYaw * qPitch;   // NO roll!

        // Update transform (no Z-roll anymore)
        transform = translate(position) * mat4(orientation);


        // Movement (Unity Style)

        vec3 localInput = vec3(0);

        if (kbd->held[Key::W]) localInput.z += 1.0f;
        if (kbd->held[Key::S]) localInput.z -= 1.0f;

        if (kbd->held[Key::D]) localInput.x -= 1.0f;
        if (kbd->held[Key::A]) localInput.x += 1.0f;

        if (kbd->held[Key::E]) localInput.y += 1.0f;
        if (kbd->held[Key::Q]) localInput.y -= 1.0f;

        float baseSpeed = normalSpeed;
        float boost = boostSpeed;

        float targetSpeed = (kbd->held[Key::Shift]) ? boost : baseSpeed;

        // Smooth acceleration
        linearVelocity = mix(linearVelocity, localInput * targetSpeed, dt * accelerationFactor);

        // Move using camera orientation (Unity style)
        vec3 worldMovement = vec3(orientation * vec4(linearVelocity, 0.0f));
        position += worldMovement * dt;

        transform = translate(position) * mat4(orientation);

        // Camera View Update
        vec3 desiredCamPos = position + vec3(transform * vec4(0, camOffsetY, -4.0f, 0));
        camPos = mix(camPos, desiredCamPos, dt * cameraSmoothFactor);

        cam->view = lookAt
        (
            camPos,
            camPos + vec3(transform[2]),    // forward
            vec3(transform[1])              // up
        );
    }




}


