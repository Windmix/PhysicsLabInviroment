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


