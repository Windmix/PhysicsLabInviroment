#pragma once
#include "cameramanager.h"
#include "mathray.h"
#include "renderdevice.h"
#include "model.h"
#include "gltf.h"

class Object; //forward declaration

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

        Object* owner;
        bool ishit;

        AABB();
        void UpdateAndDrawAABB();
        void DrawAABB();

        void Expand(const glm::vec3& point);
        glm::vec3 GetABBCenter() const;
        glm::vec3 GetABBSize() const;
 

    };
    struct SupportPoint
    {
        glm::vec3 point;   // Minkowski difference
        glm::vec3 pointA;  // vertex on A
        glm::vec3 pointB;  // vertex on B
    };

    struct Simplex
    {
        std::vector<SupportPoint> pts;
        void push_front(const SupportPoint& sp);
        void DrawSimplex(const Simplex& simplex, const glm::vec4& color = glm::vec4(1, 1, 0, 1));
    };

    struct CollisionInfo
    {
        glm::vec3 normal;     // collision normal from A → B
        float penetration;    // depth of penetration
        glm::vec3 contactPoint; // approximate contact point
    };

    struct EPAFace
    {
        SupportPoint a, b, c;
        glm::vec3 normal;
        float distance;

        EPAFace(const SupportPoint& pa, const SupportPoint& pb, const SupportPoint& pc);
    };
    

    //render
    void LoadFromIndexBuffer(fx::gltf::Document  doc, std::vector<Physics::ColliderMesh::Triangle>& refTriangles, Physics::AABB& aabb);
    MathRay ScreenPointToRay(glm::vec2& mousePos, float ScreenWidth, float ScreenHeight);

    //physics
    bool CheckAABBCollision(const AABB& a, const AABB& b);
    std::vector<std::pair<Physics::AABB, Physics::AABB>> PlaneSweepOverlaps(std::vector<AABB>& aabbs);
    bool CheckRayHitAABB(AABB& aabb, MathRay& ray, RayProperties& rayproperties);


    bool EPA(const std::vector<SupportPoint>& simplex, const std::vector<glm::vec3>& vertsA, const std::vector<glm::vec3>& vertsB,
        glm::vec3& outNormal, float& outPenetration, glm::vec3& outPoint);

    //GJK
    bool GJK_Intersect(const std::vector<glm::vec3>& A, const std::vector<glm::vec3>& B, Simplex& simplex, glm::vec3& outCollisionPoint);
    glm::vec3 Support(const std::vector<glm::vec3>& verts, const glm::vec3& dir);
    SupportPoint SupportMinkowski(const std::vector<glm::vec3>& A, const std::vector<glm::vec3>& B, const glm::vec3& dir);
    bool DoSimplex(Simplex& simplex, glm::vec3& dir);

    //draw
    void DrawGJKDirection(const glm::vec3& dir, const glm::vec3& origin = glm::vec3(0), const glm::vec4& color = glm::vec4(0, 1, 0, 1));

};

namespace SortingAlgorithm
{
    // Merge function
    template<typename T>
    void Merge(std::vector<T>& arr, int left, int mid, int right);

    // Merge sort recursive function
    template<typename T>
    void MergeSort(std::vector<T>& arr, int left, int right);
}

