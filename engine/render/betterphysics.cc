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

bool Physics::GJK_Intersect(const std::vector<glm::vec3>& A, const std::vector<glm::vec3>& B)
{
    //an initial direction
    glm::vec3 dir(1, 0, 0);

    //Build first simplex point
    Simplex simplex;
    simplex.pts.push_back(SupportMinkowski(A, B, dir));

    // New direction = toward origin
    dir = -simplex.pts[0];

    // 3. Iterate
    for (int iter = 0; iter < 50; iter++)
    {
        // New support point
        glm::vec3 newPt = SupportMinkowski(A, B, dir);

        // If new point does not pass origin in direction → no collision
        if (glm::dot(newPt, dir) < 0)
            return false;

        // Add point to simplex
        simplex.push_front(newPt);
        glm::mat4 transform = glm::mat4(0.5);
        transform[3] = glm::vec4(newPt, 1);
        Debug::DrawBox(transform, glm::vec4(1, 1, 0, 1));
        // Draw current simplex
        simplex.DrawSimplex(simplex);
        //Draw search direction
        DrawGJKDirection(dir);

        // Process simplex; dir will be updated
        if (DoSimplex(simplex, dir))
        {
            return true; // Origin inside → collision
        }
    }

    return false; // rare fallback
}

glm::vec3 Physics::Support(const std::vector<glm::vec3>& verts, const glm::vec3& dir)
{
    float best = -FLT_MAX;
    glm::vec3 bestV = glm::vec3(0);

    for (auto& v : verts)
    {
        float d = glm::dot(v, dir);
        if (d > best) 
        {
            best = d;
            bestV = v;
        }
    }
    return bestV;
}

glm::vec3 Physics::SupportMinkowski(const std::vector<glm::vec3>& A, const std::vector<glm::vec3>& B, const glm::vec3& dir)
{
    return Support(A, dir) - Support(B, -dir);
}

bool Physics::DoSimplex(Simplex& simplex, glm::vec3& dir)
{
    // Number of points in simplex
    switch (simplex.pts.size())
    {
        // --------------------------
        //       LINE (2 points)
        // --------------------------
    case 2:
    {
        glm::vec3 a = simplex.pts[0];
        glm::vec3 b = simplex.pts[1];

        glm::vec3 ab = b - a;
        glm::vec3 ao = -a;

        // If AB still points toward origin → keep both points
        if (glm::dot(ab, ao) > 0)
        {
            // New direction is perpendicular to AB toward origin
            dir = glm::cross(glm::cross(ab, ao), ab);
        }
        else
        {
            // Remove B, keep only A
            simplex.pts = { a };
            dir = ao;
        }

        return false;
    }

    // --------------------------
    //     TRIANGLE (3 points)
    // --------------------------
    case 3:
    {
        glm::vec3 a = simplex.pts[0];
        glm::vec3 b = simplex.pts[1];
        glm::vec3 c = simplex.pts[2];

        glm::vec3 ab = b - a;
        glm::vec3 ac = c - a;
        glm::vec3 ao = -a;

        glm::vec3 abc = glm::cross(ab, ac);

        // Check if origin is in region AB
        if (glm::dot(glm::cross(abc, ac), ao) > 0)
        {
            simplex.pts = { a, c };
            dir = glm::cross(glm::cross(ac, ao), ac);
        }
        // Check if origin is in region AC
        else if (glm::dot(glm::cross(ab, abc), ao) > 0)
        {
            simplex.pts = { a, b };
            dir = glm::cross(glm::cross(ab, ao), ab);
        }
        else
        {
            // Origin is either above or below triangle
            if (glm::dot(abc, ao) > 0)
            {
                dir = abc;
            }
            else
            {
                // Flip winding to maintain correct orientation
                std::swap(simplex.pts[1], simplex.pts[2]);
                dir = -abc;
            }
        }
        return false;
    }

    // --------------------------
    //    TETRAHEDRON (4 points)
    // --------------------------
    case 4:
    {
        glm::vec3 a = simplex.pts[0];
        glm::vec3 b = simplex.pts[1];
        glm::vec3 c = simplex.pts[2];
        glm::vec3 d = simplex.pts[3];

        glm::vec3 ao = -a;

        glm::vec3 ab = b - a;
        glm::vec3 ac = c - a;
        glm::vec3 ad = d - a;

        glm::vec3 abc = glm::cross(ab, ac);
        glm::vec3 acd = glm::cross(ac, ad);
        glm::vec3 adb = glm::cross(ad, ab);

        // Check each face if origin lies outside
        if (glm::dot(abc, ao) > 0)
        {
            simplex.pts = { a, b, c };
            dir = abc;
            return false;
        }

        if (glm::dot(acd, ao) > 0)
        {
            simplex.pts = { a, c, d };
            dir = acd;
            return false;
        }

        if (glm::dot(adb, ao) > 0)
        {
            simplex.pts = { a, d, b };
            dir = adb;
            return false;
        }

        // Origin is inside the tetrahedron → COLLISION
        return true;
    }
    }

    return false;
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
void Physics::DrawGJKDirection(const glm::vec3& dir, const glm::vec3& origin, const glm::vec4& color)
{
    Debug::DrawLine(origin, origin + dir, 2.0f, color, color);
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

void Physics::Simplex::push_front(const glm::vec3& p)
{
    pts.insert(pts.begin(), p);
}
void Physics::Simplex::DrawSimplex(const Simplex& simplex, const glm::vec4& color)
{
    // Draw edges between all points in the simplex
    for (int i = 0; i < simplex.pts.size(); i++)
    {
        for (int j = i + 1; j < simplex.pts.size(); j++)
        {
            Debug::DrawLine(simplex.pts[i], simplex.pts[j], 2.0f, color, color);
        }
    }
}