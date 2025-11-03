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

        if (kbd->held[Key::W] || kbd->held[Key::S])
        {
            float moveDirectionFB = kbd->held[Key::W] ? 1.0f : -1.0f;

            if (kbd->held[Key::Shift])
                this->currentSpeed = mix(this->currentSpeed, this->boostSpeed * moveDirectionFB, std::min(1.0f, dt * 30.0f));
            else
                this->currentSpeed = mix(this->currentSpeed, this->normalSpeed * moveDirectionFB, std::min(1.0f, dt * 90.0f));

            desiredVelocity = vec3(0, 0, this->currentSpeed);
        }
        else if (kbd->held[Key::A] || kbd->held[Key::D])
        {
            float moveDirectionLR = kbd->held[Key::A] ? 1.0f : -1.0f;

            if (kbd->held[Key::Shift])
                this->currentSpeed = mix(this->currentSpeed, this->boostSpeed * moveDirectionLR, std::min(1.0f, dt * 30.0f));
            else
                this->currentSpeed = mix(this->currentSpeed, this->normalSpeed * moveDirectionLR, std::min(1.0f, dt * 90.0f));

            desiredVelocity = vec3(this->currentSpeed, 0, 0);
           
        }
        else if (kbd->held[Key::E] || kbd->held[Key::Q])
        {
            float moveDirectionTD = kbd->held[Key::E] ? 1.0f : -1.0f;

            if (kbd->held[Key::Shift])
                this->currentSpeed = mix(this->currentSpeed, this->boostSpeed * moveDirectionTD, std::min(1.0f, dt * 30.0f));
            else
                this->currentSpeed = mix(this->currentSpeed, this->normalSpeed * moveDirectionTD, std::min(1.0f, dt * 90.0f));

            desiredVelocity = vec3(0, this->currentSpeed, 0);

        }
        else
        {

            this->currentSpeed = 0;

        }
        desiredVelocity = this->transform * vec4(desiredVelocity, 0.0f);
        this->linearVelocity = mix(this->linearVelocity, desiredVelocity, dt * accelerationFactor);

        this->position += this->linearVelocity * dt * 20.0f;

        // Apply rotation
        quat localOrientation = quat(vec3(-rotYSmooth, rotXSmooth, rotZSmooth));
        this->orientation = this->orientation * localOrientation;
        this->rotationZ -= rotXSmooth;
        this->rotationZ = clamp(this->rotationZ, -45.0f, 45.0f);

        mat4 T = translate(this->position) * (mat4)this->orientation;
        this->transform = T * (mat4)quat(vec3(0, 0, rotationZ));
        this->rotationZ = mix(this->rotationZ, 0.0f, dt * cameraSmoothFactor);

        // Update camera view transform
        vec3 desiredCamPos = this->position + vec3(this->transform * vec4(0, camOffsetY, -4.0f, 0));
        this->camPos = mix(this->camPos, desiredCamPos, dt * cameraSmoothFactor);
        cam->view = lookAt(this->camPos, this->camPos + vec3(this->transform[2]), vec3(this->transform[1]));
    }


}