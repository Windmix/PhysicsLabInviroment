#pragma once
#include "cameramanager.h"
#include "mathray.h"
#include "renderdevice.h"
#include "model.h"
#include "gltf.h"



namespace Physics
{

    struct RayProperties
    {
        //AABB normals
        glm::vec3 AABBintersection = glm::vec3(0), AABBnormalEnd = glm::vec3(0);
        //original Normals
        glm::vec3 intersection = glm::vec3(0), normalEnd  = glm::vec3(0);

        //is inside
        bool isInsideObject = false;
    };

    struct ColliderMesh
    {
        struct Triangle
        {

            glm::vec4 color;       // current display color
            glm::vec4 selectedColor; // selected color
            glm::vec4 og_color; // og color

            glm::vec3 verticies[3];
            glm::vec3 normal;
            bool selected = false;
            void SetSelected(bool s);
        };
    };

    class AABB
    {
    public:
     
        glm::vec3 min;
        glm::vec3 max;

        AABB();
        void UpdateAndDrawAABB();
        void DrawAABB();

        void Expand(const glm::vec3& point);
        glm::vec3 GetABBCenter() const;
        glm::vec3 GetABBSize() const;
 

    };



    void LoadFromIndexBuffer(fx::gltf::Document  doc, std::vector<Physics::ColliderMesh::Triangle>& refTriangles, Physics::AABB& aabb);

    MathRay ScreenPointToRay(glm::vec2& mousePos, float ScreenWidth, float ScreenHeight);

    bool CheckAABBCollision(const AABB& a, const AABB& b);

    std::vector<std::pair<Physics::AABB, Physics::AABB>> PlaneSweepOverlaps(std::vector<AABB>& aabbs);
    bool CheckRayHitAABB(AABB& aabb, MathRay& ray, RayProperties& rayproperties);
};

namespace SortingAlgorithm
{
    // Merge function
    void Merge(std::vector<Physics::AABB>& arr, int left, int mid, int right);

    // Merge sort recursive function
    void MergeSort(std::vector<Physics::AABB>& arr, int left, int right);
}

