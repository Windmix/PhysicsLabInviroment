#pragma once
#include "cameramanager.h"
#include "mathray.h"
#include "quad.h"
#include "renderdevice.h"
#include "model.h"
namespace Physics
{
    struct RayProperties
    {
        glm::vec3 intersection, normalEnd;
    };
    class AABB
    {
    public:
        Render::ModelId model;
        glm::vec3 min;
        glm::vec3 max;

        AABB();
        void Expand(const glm::vec3& point);
        glm::vec3 GetABBCenter() const;
        glm::vec3 GetABBSize() const;
        void drawBox(glm::mat4 transform);

    };

    MathRay ScreenPointToRay(glm::vec2& mousePos, float ScreenWidth, float ScreenHeight);
    bool CheckRayHit(Quad& myQuad, MathRay& ray, RayProperties& rayproperties);



};
