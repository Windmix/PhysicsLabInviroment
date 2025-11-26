#include <config.h>
#include <limits>
#include "betterphysics.h"
#include "debugrender.h"
#include <imgui_impl_opengl3.h>





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



// Simple AABB collision check (can be used as narrow-phase for now)
bool Physics::CheckAABBCollision(const AABB& a, const AABB& b)
{
    bool overlapX = (a.min.x <= b.max.x && a.max.x >= b.min.x);
    bool overlapY = (a.min.y <= b.max.y && a.max.y >= b.min.y);
    bool overlapZ = (a.min.z <= b.max.z && a.max.z >= b.min.z);

    return overlapX && overlapY && overlapZ;
}

std::vector<std::pair<Physics::AABB, Physics::AABB>> Physics::PlaneSweepOverlaps(std::vector<AABB>& aabbs)
{
    std::vector<Physics::AABB> sortedAABBs = aabbs;
    SortingAlgorithm::MergeSort(sortedAABBs, 0, sortedAABBs.size() - 1);

    std::vector<std::pair<Physics::AABB, Physics::AABB>> candidates;

    for (int i = 0; i < sortedAABBs.size(); i++)
    {
        Physics::AABB a = sortedAABBs[i];
        for (int j = i + 1; j < sortedAABBs.size(); j++)
        {
            Physics::AABB b = sortedAABBs[j];

            // Stop early if b is completely to the right of a
            if (b.min.x > a.max.x)
                break;

            // Check overlap in y and z axes
            bool overlapY = (a.min.y <= b.max.y && a.max.y >= b.min.y);
            bool overlapZ = (a.min.z <= b.max.z && a.max.z >= b.min.z);

            if (overlapY && overlapZ)
            {
                // constructs pair(a,b) directly inside vector
                candidates.emplace_back(a, b);
            }
        }
    }

    return candidates;
}

bool Physics::CheckRayHitAABB(AABB& aabb, MathRay& ray, RayProperties& rayproperties)
{
    const float epsilon = 1e-5f;
    glm::vec3 rayOrigin = ray.GetOrigin();
    glm::vec3 rayDir = ray.GetDirection();
    float rayLength = ray.GetRayLength();

    bool hit = false;
    float closestT = rayLength;

    // Iterate over each face of the AABB (6 planes)
    struct Face { glm::vec3 normal; glm::vec3 point; };
    Face faces[6] = {
        { glm::vec3(1, 0, 0), glm::vec3(aabb.max.x, 0, 0) },
        { glm::vec3(-1, 0, 0), glm::vec3(aabb.min.x, 0, 0) },
        { glm::vec3(0,  1, 0), glm::vec3(0, aabb.max.y, 0) },
        { glm::vec3(0, -1, 0), glm::vec3(0, aabb.min.y, 0) },
        { glm::vec3(0, 0,  1), glm::vec3(0, 0, aabb.max.z) },
        { glm::vec3(0, 0, -1), glm::vec3(0, 0, aabb.min.z) }
    };

    for (int i = 0; i < 6; ++i)
    {
        Face& f = faces[i];
        float denom = glm::dot(f.normal, rayDir);

        // if parallel
        if (fabs(denom) < epsilon) continue;

        float t = glm::dot(f.point - rayOrigin, f.normal) / denom;
        // behind or farther than previous hit
        if (t < 0 || t > closestT) continue;

        glm::vec3 intersect = rayOrigin + rayDir * t;

        // Check if intersection point is inside face bounds
        if (intersect.x + epsilon >= aabb.min.x && intersect.x - epsilon <= aabb.max.x &&
            intersect.y + epsilon >= aabb.min.y && intersect.y - epsilon <= aabb.max.y &&
            intersect.z + epsilon >= aabb.min.z && intersect.z - epsilon <= aabb.max.z)
        {
            closestT = t;
            rayproperties.AABBintersection = intersect;
            rayproperties.AABBnormalEnd = intersect + f.normal;
            hit = true;
        }
    }

    if (hit)
    {/*
        printf("AABB HIT!\n");
        printf("Intersection point: (%.3f, %.3f, %.3f)\n",
            rayproperties.AABBintersection.x,
            rayproperties.AABBintersection.y,
            rayproperties.AABBintersection.z);
        printf("Face normal end: (%.3f, %.3f, %.3f)\n",
            rayproperties.AABBnormalEnd.x,
            rayproperties.AABBnormalEnd.y,
            rayproperties.AABBnormalEnd.z);*/
    }


    return hit;
}
void Physics::LoadFromIndexBuffer(fx::gltf::Document doc, std::vector<Physics::ColliderMesh::Triangle>& refTriangles, Physics::AABB& aabb)
{
    fx::gltf::Primitive const& primitive = doc.meshes[0].primitives[0];

    fx::gltf::Accessor const& ibAccessor = doc.accessors[primitive.indices];
    fx::gltf::BufferView const& ibView = doc.bufferViews[ibAccessor.bufferView];
    fx::gltf::Buffer const& ib = doc.buffers[ibView.buffer];

    fx::gltf::Accessor const& vbAccessor = doc.accessors[primitive.attributes.find("POSITION")->second];
    fx::gltf::BufferView const& vbView = doc.bufferViews[vbAccessor.bufferView];
    fx::gltf::Buffer const& vb = doc.buffers[vbView.buffer];

    size_t numIndices = ibAccessor.count;
    uint16_t const* indexBuffer = (uint16_t const*)&ib.data[ibAccessor.byteOffset + ibView.byteOffset];
    float const* vertexBuffer = (float const*)&vb.data[vbAccessor.byteOffset + vbView.byteOffset];

    size_t vSize = (vbAccessor.type == fx::gltf::Accessor::Type::Vec3) ? 3 : 4;

    for (size_t i = 0; i < numIndices; i += 3)
    {
        Physics::ColliderMesh::Triangle tri;
        tri.verticies[0] = glm::vec3(vertexBuffer[vSize * indexBuffer[i] + 0],
            vertexBuffer[vSize * indexBuffer[i] + 1],
            vertexBuffer[vSize * indexBuffer[i] + 2]);
        tri.verticies[1] = glm::vec3(vertexBuffer[vSize * indexBuffer[i + 1] + 0],
            vertexBuffer[vSize * indexBuffer[i + 1] + 1],
            vertexBuffer[vSize * indexBuffer[i + 1] + 2]);
        tri.verticies[2] = glm::vec3(vertexBuffer[vSize * indexBuffer[i + 2] + 0],
            vertexBuffer[vSize * indexBuffer[i + 2] + 1],
            vertexBuffer[vSize * indexBuffer[i + 2] + 2]);

        // Expand the AABB 
        aabb.Expand(tri.verticies[0]);
        aabb.Expand(tri.verticies[1]);
        aabb.Expand(tri.verticies[2]);

        // Compute normal
        glm::vec3 AB = tri.verticies[1] - tri.verticies[0];
        glm::vec3 AC = tri.verticies[2] - tri.verticies[0];
        tri.normal = glm::cross(AC, AB);

        refTriangles.push_back(std::move(tri));
    }

}



void Physics::AABB::Expand(const glm::vec3& point)
{
    min = glm::min(min, point);
    max = glm::max(max, point);
}
Physics::AABB::AABB() : min(std::numeric_limits<float>::max()),max(std::numeric_limits<float>::lowest())
{
}

glm::vec3 Physics::AABB::GetABBCenter() const
{
    return (min + max) * 0.5f;
}

glm::vec3 Physics::AABB::GetABBSize() const
{
    return max - min;
}

void Physics::ColliderMesh::Triangle::SetSelected(bool s)
{
    selected = s;

    color = selected ? selectedColor : og_color;
}

void SortingAlgorithm::Merge(std::vector<Physics::AABB>& arr, int left, int mid, int right)
{
    int n1 = mid - left + 1;
    int n2 = right - mid;

    std::vector<Physics::AABB> L(n1);
    std::vector<Physics::AABB> R(n2);

    for (int i = 0; i < n1; i++)
        L[i] = arr[left + i];

    for (int j = 0; j < n2; j++)
        R[j] = arr[mid + 1 + j];

    int i = 0;
    int j = 0;
    int k = left;

    while (i < n1 && j < n2)
    {
        if (L[i].min.x <= R[j].min.x)  // sorting by min.x
        { 
            arr[k] = L[i];
            i++;
        }
        else 
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1)
    {
        arr[k] = L[i];
        i++; 
        k++;
    }
    while (j < n2) 
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void SortingAlgorithm::MergeSort(std::vector<Physics::AABB>& arr, int left, int right)
{
    if (left < right)
    {
        int mid = left + (right - left) / 2;

        SortingAlgorithm::MergeSort(arr, left, mid);
        SortingAlgorithm::MergeSort(arr, mid + 1, right);
        SortingAlgorithm::Merge(arr, left, mid, right);
    }
}
