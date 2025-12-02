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

bool Physics::GJK_Intersect(const std::vector<glm::vec3>& A, const std::vector<glm::vec3>& B, Simplex& simplex, glm::vec3& outCollisionPoint)
{
    glm::vec3 dir(1, 0, 0);


    simplex.push_front(SupportMinkowski(A, B, dir));

    dir = -simplex.pts[0].point;

    for (int iter = 0; iter < 50; iter++)
    {
        SupportPoint sp = SupportMinkowski(A, B, dir);

        // If support point does not pass origin in search direction, no collision
        if (glm::dot(sp.point, dir) < 0)
            return false;

        simplex.push_front(sp);

        // Optional debug drawing
        // Debug::DrawBox(glm::translate(glm::mat4(1.0f), sp.point) * glm::scale(glm::mat4(1.0f), glm::vec3(0.5f)), glm::vec4(1,1,0,1));
        // simplex.DrawSimplex(glm::vec4(1,0,1,1));

        // Update simplex and direction
        if (DoSimplex(simplex, dir))
        {
            glm::vec3 sum(0.0f);
            for (auto& pt : simplex.pts)
                sum += pt.pointA;  // average the contributing points on A
            outCollisionPoint = sum / float(simplex.pts.size());
            // Collision found, simplex contains last tetrahedron
            return true;
        }
    }

    return false;
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

Physics::SupportPoint Physics::SupportMinkowski(const std::vector<glm::vec3>& A, const std::vector<glm::vec3>& B, const glm::vec3& dir)
{
    glm::vec3 pA = Support(A, dir);   // support
    glm::vec3 pB = Support(B, -dir);

    SupportPoint sp;
    sp.point = pA - pB;
    sp.pointA = pA;
    sp.pointB = pB;
    return sp;
}

bool Physics::DoSimplex(Simplex& simplex, glm::vec3& dir)
{
    switch (simplex.pts.size())
    {
    case 2: // Line
    {
        glm::vec3 a = simplex.pts[0].point;
        glm::vec3 b = simplex.pts[1].point;

        glm::vec3 ab = b - a;
        glm::vec3 ao = -a;

        if (glm::dot(ab, ao) > 0)
            dir = glm::cross(glm::cross(ab, ao), ab);
        else
        {
            simplex.pts = { simplex.pts[0] }; // keep only A
            dir = ao;
        }

        return false;
    }

    case 3: // Triangle
    {
        glm::vec3 a = simplex.pts[0].point;
        glm::vec3 b = simplex.pts[1].point;
        glm::vec3 c = simplex.pts[2].point;

        glm::vec3 ab = b - a;
        glm::vec3 ac = c - a;
        glm::vec3 ao = -a;

        glm::vec3 abc = glm::cross(ab, ac);

        if (glm::dot(glm::cross(abc, ac), ao) > 0)
        {
            simplex.pts = { simplex.pts[0], simplex.pts[2] };
            dir = glm::cross(glm::cross(ac, ao), ac);
        }
        else if (glm::dot(glm::cross(ab, abc), ao) > 0)
        {
            simplex.pts = { simplex.pts[0], simplex.pts[1] };
            dir = glm::cross(glm::cross(ab, ao), ab);
        }
        else
        {
            if (glm::dot(abc, ao) > 0)
                dir = abc;
            else
            {
                std::swap(simplex.pts[1], simplex.pts[2]);
                dir = -abc;
            }
        }
        return false;
    }

    case 4: // Tetrahedron
    {
        glm::vec3 a = simplex.pts[0].point;
        glm::vec3 b = simplex.pts[1].point;
        glm::vec3 c = simplex.pts[2].point;
        glm::vec3 d = simplex.pts[3].point;

        glm::vec3 ao = -a;

        glm::vec3 ab = b - a;
        glm::vec3 ac = c - a;
        glm::vec3 ad = d - a;

        glm::vec3 abc = glm::cross(ab, ac);
        glm::vec3 acd = glm::cross(ac, ad);
        glm::vec3 adb = glm::cross(ad, ab);

        if (glm::dot(abc, ao) > 0)
        {
            simplex.pts = { simplex.pts[0], simplex.pts[1], simplex.pts[2] };
            dir = abc;
            return false;
        }

        if (glm::dot(acd, ao) > 0)
        {
            simplex.pts = { simplex.pts[0], simplex.pts[2], simplex.pts[3] };
            dir = acd;
            return false;
        }

        if (glm::dot(adb, ao) > 0)
        {
            simplex.pts = { simplex.pts[0], simplex.pts[3], simplex.pts[1] };
            dir = adb;
            return false;
        }

        return true; // origin inside tetrahedron
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
    //this bool is for change color
    aabb.ishit = hit;

    return hit;
}

bool Physics::EPA(const std::vector<SupportPoint>& simplex, const std::vector<glm::vec3>& vertsA, const std::vector<glm::vec3>& vertsB, glm::vec3& outNormal, float& outPenetration, glm::vec3& outPoint)
{
    if (simplex.size() != 4)
        return false;

    std::vector<EPAFace> faces;
    faces.emplace_back(simplex[0], simplex[1], simplex[2]);
    faces.emplace_back(simplex[0], simplex[3], simplex[1]);
    faces.emplace_back(simplex[0], simplex[2], simplex[3]);
    faces.emplace_back(simplex[1], simplex[3], simplex[2]);

    const float tolerance = 0.0001f;
    const int maxIterations = 50;
    int iterations = 0;

    while (iterations < maxIterations && !faces.empty())
    {
        // 1. Find face closest to origin
        auto closestIt = std::min_element(faces.begin(), faces.end(),
            [](const EPAFace& a, const EPAFace& b) { return a.distance < b.distance; });
        EPAFace& closestFace = *closestIt;

        // 2. Get new support point in face normal direction
        SupportPoint newPoint = SupportMinkowski(vertsA, vertsB, closestFace.normal);
        float newDist = glm::dot(newPoint.point, closestFace.normal);

        // Check termination
        if (newDist - closestFace.distance < tolerance)
        {
            outNormal = closestFace.normal;
            outPenetration = newDist;

            // Approximate contact point on object A
            outPoint = (closestFace.a.pointA + closestFace.b.pointA + closestFace.c.pointA) / 3.0f;
            return true;
        }

        // 3. Find all faces visible from new point
        std::vector<int> visibleFaces;
        for (int i = 0; i < faces.size(); ++i)
        {
            if (glm::dot(faces[i].normal, newPoint.point - faces[i].a.point) > 0.0f)
                visibleFaces.push_back(i);
        }

        if (visibleFaces.empty())
        {
            // No visible faces, cannot expand further
            break;
        }

        // 4. Build horizon edges (edges on boundary between visible and non-visible faces)
        std::vector<std::pair<SupportPoint, SupportPoint>> horizonEdges;

        auto edgeMatches = [](const std::pair<SupportPoint, SupportPoint>& e1,
            const std::pair<SupportPoint, SupportPoint>& e2) -> bool
            {
                return (e1.first.point == e2.second.point && e1.second.point == e2.first.point);
            };

        for (int idx : visibleFaces)
        {
            EPAFace& f = faces[idx];
            std::pair<SupportPoint, SupportPoint> edges[3] = { {f.a, f.b}, {f.b, f.c}, {f.c, f.a} };

            for (auto& e : edges)
            {
                bool isShared = false;
                for (int vidx : visibleFaces)
                {
                    if (vidx == idx) continue;
                    EPAFace& other = faces[vidx];
                    std::pair<SupportPoint, SupportPoint> otherEdges[3] = { {other.a, other.b}, {other.b, other.c}, {other.c, other.a} };
                    for (auto& oe : otherEdges)
                    {
                        if (edgeMatches(e, oe))
                        {
                            isShared = true;
                            break;
                        }
                    }
                    if (isShared) break;
                }
                if (!isShared)
                    horizonEdges.push_back(e);
            }
        }

        // 5. Remove visible faces
        std::sort(visibleFaces.rbegin(), visibleFaces.rend());
        for (int idx : visibleFaces)
            faces.erase(faces.begin() + idx);

        // 6. Create new faces from horizon edges and new point
        for (auto& edge : horizonEdges)
        {
            EPAFace newFace(edge.first, edge.second, newPoint);

            // Ensure normal points away from origin
            if (glm::dot(newFace.normal, newFace.a.point) < 0.0f)
            {
                std::swap(newFace.b, newFace.c);
                newFace.normal = -newFace.normal;
                newFace.distance = -newFace.distance;
            }

            faces.push_back(newFace);
        }

        iterations++;
    }

    return false; // Failed to converge
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
    owner = nullptr;
    ishit = false;
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
template<typename T>
void SortingAlgorithm::Merge(std::vector<T>& arr, int left, int mid, int right)
{
    int n1 = mid - left + 1;
    int n2 = right - mid;

    std::vector<T> L(n1);
    std::vector<T> R(n2);

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
template<typename T>
void SortingAlgorithm::MergeSort(std::vector<T>& arr, int left, int right)
{
    if (left < right)
    {
        int mid = left + (right - left) / 2;

        SortingAlgorithm::MergeSort(arr, left, mid);
        SortingAlgorithm::MergeSort(arr, mid + 1, right);
        SortingAlgorithm::Merge(arr, left, mid, right);
    }
}

void Physics::Simplex::push_front(const SupportPoint& sp)
{
    pts.insert(pts.begin(), sp);
}
void Physics::Simplex::DrawSimplex(const Simplex& simplex, const glm::vec4& color)
{
    // Draw edges between all points in the simplex
    for (int i = 0; i < simplex.pts.size(); i++)
    {
        for (int j = i + 1; j < simplex.pts.size(); j++)
        {
            Debug::DrawLine(simplex.pts[i].point, simplex.pts[j].point, 2.0f, color, color);
        }
    }
}

inline Physics::EPAFace::EPAFace(const SupportPoint& pa, const SupportPoint& pb, const SupportPoint& pc) : a(pa), b(pb), c(pc)
{
    normal = glm::normalize(glm::cross(b.point - a.point, c.point - a.point));
    distance = glm::dot(normal, a.point);
    if (distance < 0.0f)
    {
        normal = -normal;
        distance = -distance;
        std::swap(b, c);
    }
}
