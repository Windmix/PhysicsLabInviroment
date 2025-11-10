#include <config.h>
#include "object.h"
#include "gltf.h"
#include "debugrender.h"

Object::Object()
{

	modelID = 0;
	rotation = glm::vec3(0);
	position = glm::vec3(0);
	transform = glm::mat4(1);
	scale = glm::vec3(1.0f);
    plane = MathPlane();

}

Object::~Object()
{
}
void Object::createObject(std::string filePath)

{
    modelID = Render::LoadModel(filePath);

    // Initialize transform/position/scale
    transform = glm::mat4(1.0f);
    position = glm::vec3(0.0f);
    rotation = glm::vec3(0.0f);
    scale = glm::vec3(1.0f);
}


void Object::DrawAABBOnObject()
{
    glm::vec3 min = aabb.min;
    glm::vec3 max = aabb.max;

    glm::vec3 corners[8] =
    {
        {min.x, min.y, min.z},
        {max.x, min.y, min.z},
        {max.x, max.y, min.z},
        {min.x, max.y, min.z},
        {min.x, min.y, max.z},
        {max.x, min.y, max.z},
        {max.x, max.y, max.z},
        {min.x, max.y, max.z}
    };

    // Draw edges (12 lines)
    Debug::DrawLine(corners[0], corners[1], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[1], corners[2], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[2], corners[3], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[3], corners[0], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));

    Debug::DrawLine(corners[4], corners[5], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[5], corners[6], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[6], corners[7], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[7], corners[4], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));

    Debug::DrawLine(corners[0], corners[4], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[1], corners[5], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[2], corners[6], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
    Debug::DrawLine(corners[3], corners[7], 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
}

bool Object::CheckRayHit(Object& myObj, MathRay& ray, Physics::RayProperties& rayproperties)
{

    bool hitAny = false;
    float closestT = std::numeric_limits<float>::max(); // store closest hit
    Physics::ColliderMesh::Triangle* closestTri = nullptr;

    // Build normal matrix to rotate normals correctly (ignores translation)
    glm::mat3 normalMatrix = glm::transpose(glm::inverse(glm::mat3(myObj.transform)));

    for (auto& tri : myObj.triangles)
    {
        // Transform triangle vertices from local space to world space
        glm::vec3 v0 = glm::vec3(myObj.transform * glm::vec4(tri.verticies[0], 1.0f));
        glm::vec3 v1 = glm::vec3(myObj.transform * glm::vec4(tri.verticies[1], 1.0f));
        glm::vec3 v2 = glm::vec3(myObj.transform * glm::vec4(tri.verticies[2], 1.0f));

        // Compute the edges of the triangle
        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;

        // Start Möller–Trumbore algorithm to detect ray-triangle intersection
        glm::vec3 h = glm::cross(ray.GetDirection(), edge2);
        float a = glm::dot(edge1, h);

        // Skip triangle if ray is parallel to it
        if (fabs(a) < 1e-6f) continue;

        float f = 1.0f / a;
        glm::vec3 s = ray.GetOrigin() - v0;
        float u = f * glm::dot(s, h);

        // Reject intersection if it is outside the triangle (u parameter)
        if (u < 0.0f || u > 1.0f) continue;

        glm::vec3 q = glm::cross(s, edge1);
        float v = f * glm::dot(ray.GetDirection(), q);

        // Reject intersection if it is outside the triangle (v parameter)
        if (v < 0.0f || u + v > 1.0f) continue;

        // Compute the distance along the ray
        float t = f * glm::dot(edge2, q);

        // Reject intersections behind the ray or beyond its length
        if (t < 1e-6f || t > ray.GetRayLength()) continue;

        // Keep the closest hit only
        if (t < closestT)
        {
            closestT = t;
            closestTri = &tri; // Save pointer to closest triangle

            // Compute intersection point in world space
            glm::vec3 intersection = ray.GetOrigin() + ray.GetDirection() * t;

            // Rotate normal to world space
            glm::vec3 normal = glm::normalize(normalMatrix * tri.normal);

            // Make sure the normal always faces against the ray
            if (glm::dot(ray.GetDirection(), normal) > 0.0f)
                normal = -normal;

            // Store intersection info for later use
            rayproperties.intersection = intersection;
            rayproperties.normalEnd = intersection + normal;

            hitAny = true;
        }

    }
    // Only mark the closest triangle
    if (closestTri)
    {
        closestTri->SetSelected(true);
      
    }

    return hitAny;
}

void Object::UpdateAndDrawAABBObject()
{
   // Reset AABB
    aabb.min = glm::vec3(std::numeric_limits<float>::max());
    aabb.max = glm::vec3(std::numeric_limits<float>::lowest());

    // Expand AABB by transformed triangle vertices
    for (auto& tri : triangles)
    {
       
        aabb.Expand(glm::vec3(transform * glm::vec4(tri.verticies[0], 1.0f)));
        aabb.Expand(glm::vec3(transform * glm::vec4(tri.verticies[1], 1.0f)));
        aabb.Expand(glm::vec3(transform * glm::vec4(tri.verticies[2], 1.0f)));
    }



    DrawAABBOnObject();
}
void Object::SetOBjectRotation(glm::vec3 direction, float angle)
{
    // Normalize the axis
    glm::vec3 axis = glm::normalize(direction);

    // Build rotation matrix
    glm::mat4 r = glm::rotate(glm::mat4(1.0f), angle, axis);

    // Rebuild the transform matrix: T * R * S
    glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
    glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);

    transform = t * r * s;
}

void Object::SetOBjectPosition(glm::vec3 newPosition)
{


    position = newPosition; // update stored position

    // keep the same rotation and scale
    glm::mat4 r = glm::mat4(1.0f);
    glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
    glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);

    transform = t * r * s;
}


void Object::drawObject() const
{
	Render::RenderDevice::Draw(modelID, transform);
}

void Object::SetScale(glm::vec3 _scale)
{
	scale = _scale;
	glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
	glm::mat4 r = glm::mat4(1.0f);
	glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);
	transform = t * r * s;
}



