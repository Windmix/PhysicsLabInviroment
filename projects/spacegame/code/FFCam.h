#pragma once
#include "render/model.h"
#include "render/mathray.h"
#include "render/mathplane.h"
#include <render/quad.h>


namespace Render
{
    struct ParticleEmitter;
}

namespace Game
{

struct FFCam
{
    FFCam();
    
    glm::vec3 position = glm::vec3(0);
    glm::vec3 forward = glm::vec3(0);
    glm::quat orientation = glm::identity<glm::quat>();
    glm::vec3 camPos = glm::vec3(0, 1.0f, -2.0f);
    glm::mat4 transform = glm::mat4(1);
    glm::vec3 linearVelocity = glm::vec3(0);


    const float normalSpeed = 1.0f;
    const float boostSpeed = normalSpeed * 2.0f;
    const float accelerationFactor = 50.0f;
    const float camOffsetY = 0.0f;
    const float cameraSmoothFactor = 50.0f;

    float currentSpeedFB = 0.0f;
    float currentSpeedLR = 0.0f;
    float currentSpeedTD = 0.0f;


    float rotationZ = 0;
    float rotXSmooth = 0;
    float rotYSmooth = 0;
    float rotZSmooth = 0;

    void Update(float dt);

   


    

};

}
