#include <config.h>
#include "object.h"
#include "gltf.h"
#include "debugrender.h"
#include <core/cvar.h>
#include <gtx/matrix_decompose.hpp>

Object::Object()
{

	modelID = 0;
	position = glm::vec3(0);
	transform = glm::mat4(1);
	scale = glm::vec3(1.0f);
    plane = MathPlane();
    velocity = glm::vec3(0.0f);
    acceleration = glm::vec3(0.0f);
    force = glm::vec3(0.0f);
    torque = glm::vec3(0.0f);
    accumulatedForce = glm::vec3(0.0f);
    angularVelocity = glm::vec3(0.0f);
    accumulatedTorque = glm::vec3(0.0f);
    orientation = glm::quat(1, 0, 0, 0);
    forceMagnitude = 10000.0f;
    storedHitindex = -1;
    mass = 20.0f;
    totalMass = 0.0f;
    centerOfMass = glm::vec3(0);
    inertiaTensor = glm::mat3(0);
    totalMeshArea = 0.0f;
    inertiaTensorInv = glm::mat3(0);
    previousOrientation = glm::quat();
    previousPosition = glm::vec3(0);



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
    if (aabb.ishit)
    {
        //make it green
        glm::vec4 green = glm::vec4(0, 1, 0, 1);

        Debug::DrawLine(corners[0], corners[1], 3.0f, green, green);
        Debug::DrawLine(corners[1], corners[2], 3.0f, green, green);
        Debug::DrawLine(corners[2], corners[3], 3.0f, green, green);
        Debug::DrawLine(corners[3], corners[0], 3.0f, green, green);

        Debug::DrawLine(corners[4], corners[5], 3.0f, green, green);
        Debug::DrawLine(corners[5], corners[6], 3.0f, green, green);
        Debug::DrawLine(corners[6], corners[7], 3.0f, green, green);
        Debug::DrawLine(corners[7], corners[4], 3.0f, green, green);

        Debug::DrawLine(corners[0], corners[4], 3.0f, green, green);
        Debug::DrawLine(corners[1], corners[5], 3.0f, green, green);
        Debug::DrawLine(corners[2], corners[6], 3.0f, green, green);
        Debug::DrawLine(corners[3], corners[7], 3.0f, green, green);
    }
    else
    {
        glm::vec4 c1 = glm::vec4(0, 1, 1, 1); // cyan
        glm::vec4 c2 = glm::vec4(1, 0, 1, 1); // magenta

        Debug::DrawLine(corners[0], corners[1], 3.0f, c1, c2);
        Debug::DrawLine(corners[1], corners[2], 3.0f, c1, c2);
        Debug::DrawLine(corners[2], corners[3], 3.0f, c1, c2);
        Debug::DrawLine(corners[3], corners[0], 3.0f, c1, c2);

        Debug::DrawLine(corners[4], corners[5], 3.0f, c1, c2);
        Debug::DrawLine(corners[5], corners[6], 3.0f, c1, c2);
        Debug::DrawLine(corners[6], corners[7], 3.0f, c1, c2);
        Debug::DrawLine(corners[7], corners[4], 3.0f, c1, c2);

        Debug::DrawLine(corners[0], corners[4], 3.0f, c1, c2);
        Debug::DrawLine(corners[1], corners[5], 3.0f, c1, c2);
        Debug::DrawLine(corners[2], corners[6], 3.0f, c1, c2);
        Debug::DrawLine(corners[3], corners[7], 3.0f, c1, c2);
    }
}
 
std::vector<glm::vec3> Object::GetWorldVertices() const
{
    std::vector<glm::vec3> worldVerts;
    worldVerts.reserve(triangles.size() * 3); // pre-allocate memory

    for (const auto& tri : triangles)
    {
        worldVerts.push_back(glm::vec3(transform * glm::vec4(tri.verticies[0], 1.0f)));
        worldVerts.push_back(glm::vec3(transform * glm::vec4(tri.verticies[1], 1.0f)));
        worldVerts.push_back(glm::vec3(transform * glm::vec4(tri.verticies[2], 1.0f)));
    }

    return worldVerts;
}

std::vector<glm::vec3> Object::GetWorldInterpolationVertices() const
{
    std::vector<glm::vec3> worldVerts;
    worldVerts.reserve(triangles.size() * 3); // pre-allocate memory

    for (const auto& tri : triangles)
    {
        worldVerts.push_back(glm::vec3(renderTransform * glm::vec4(tri.verticies[0], 1.0f)));
        worldVerts.push_back(glm::vec3(renderTransform * glm::vec4(tri.verticies[1], 1.0f)));
        worldVerts.push_back(glm::vec3(renderTransform * glm::vec4(tri.verticies[2], 1.0f)));
    }

    return worldVerts;
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

            rayproperties.isInsideObject = myObj.IsPointInsideMesh(ray.GetOrigin());
            if (rayproperties.isInsideObject)
            {
                rayproperties.normalEnd = intersection - normal;
            }
            else
            {
                rayproperties.normalEnd = intersection + normal;
            }
            rayproperties.intersection = intersection;

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

bool Object::IsPointInsideMesh(const glm::vec3& point)
{
    glm::vec3 dir = glm::vec3(1.0f, 0.0f, 0.0f); // arbitrary direction
    MathRay testRay(point, dir);

    int hits = 0;

    for (auto& tri : triangles)
    {
        // Transform triangle vertices to world space
        glm::vec3 v0 = glm::vec3(transform * glm::vec4(tri.verticies[0], 1.0f));
        glm::vec3 v1 = glm::vec3(transform * glm::vec4(tri.verticies[1], 1.0f));
        glm::vec3 v2 = glm::vec3(transform * glm::vec4(tri.verticies[2], 1.0f));

        // --- Same Möller–Trumbore intersection logic ---
        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;
        glm::vec3 h = glm::cross(testRay.GetDirection(), edge2);
        float a = glm::dot(edge1, h);

        if (fabs(a) < 1e-6f) continue; // parallel

        float f = 1.0f / a;
        glm::vec3 s = point - v0;
        float u = f * glm::dot(s, h);
        if (u < 0.0f || u > 1.0f) continue;

        glm::vec3 q = glm::cross(s, edge1);
        float v = f * glm::dot(testRay.GetDirection(), q);
        if (v < 0.0f || u + v > 1.0f) continue;

        float t = f * glm::dot(edge2, q);
        if (t > 1e-6f) // ignore self-hit
            hits++;
    }

    return (hits % 2 == 1); // odd = inside
}

void Object::UpdateAABBObject()
{
   // Reset AABB
    aabb.min = glm::vec3(std::numeric_limits<float>::max());
    aabb.max = glm::vec3(std::numeric_limits<float>::lowest());

    // Expand AABB by transformed triangle vertices
    for (auto& tri : triangles)
    {
       
        aabb.Expand(glm::vec3(transform * glm::vec4(tri.verticies[0], 1.0f)  ));
        aabb.Expand(glm::vec3(transform * glm::vec4(tri.verticies[1], 1.0f)));
        aabb.Expand(glm::vec3(transform * glm::vec4(tri.verticies[2], 1.0f)));
    }
    // save pointer to abb as this object
    aabb.owner = this;
    
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

void Object::SetOBjectPosition(glm::vec3 newPosition) // sends a copy of position
{


    position = newPosition; // update stored position

    // keep the same rotation and scale
    glm::mat4 r = glm::mat4(1.0f);
    glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
    glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);

    transform = t * r * s;
}

void Object::ApplyGravityForce(glm::vec3 gravitationalVector)
{


    accumulatedForce += gravitationalVector * mass * 9.82f;
    this->force = gravitationalVector * mass * 9.82f;
}

void Object::ApplyGravityFrom(Object& other)
{
    if (&other == this) return; // do not apply gravity to itself

    // F(i,j) = G * (M(i) * M(j)) / r^2 * r(vector)
    glm::vec3 dir = other.position - this->position;
    float distSq = glm::dot(dir, dir);

    // avoid infinite force when objects overlap
    float minDistSq = 0.01f;
    distSq = glm::max(distSq, minDistSq);
    float dist = sqrt(distSq);

    //r(vector) direction
    glm::vec3 dirNorm = dir / dist;

    //G * (M(i) * M(j)) / r^2
    float forceMag = G * (this->mass * other.mass) / distSq;
    glm::vec3 force = dirNorm * forceMag;

    this->force += force;
}

void Object::ApplyForce(const glm::vec3& force, glm::vec3& forcehitPoint)
{
 
    accumulatedForce += force;
    this->force = force;

    glm::vec3 r = forcehitPoint - centerOfMass; // vector from COM to application point

    // torque = r × F
    torque = glm::cross(r, force);
    accumulatedTorque += glm::cross(r, force);

}   

void Object::Integrate(float dt)
{
    Core::CVar* r_freeze_pos = Core::CVarCreate(Core::CVarType::CVar_Int, "r_freeze_pos", "0");
    int freezePosBool = Core::CVarReadInt(r_freeze_pos);

    Core::CVar* r_freeze_rot = Core::CVarCreate(Core::CVarType::CVar_Int, "r_freeze_rot", "0");
    int freezeRotBool = Core::CVarReadInt(r_freeze_rot);

    // --- LINEAR MOTION ---
    // Semi-implicit Euler
    if (freezePosBool == 0)
    {
        acceleration = this->force / mass; // F = kg/s^2 --> a = F/kg
        velocity += acceleration * dt;
        position += velocity * dt;
    }
   
  
   
 
    if (freezeRotBool == 0)
    {

        // --- ANGULAR MOTION ---
        glm::vec3 angularAcceleration = inertiaTensorInv * torque;
        angularVelocity += angularAcceleration * dt;

        glm::quat deltaRot = glm::quat(0.0f, angularVelocity * dt);
        orientation += 0.5f * deltaRot * orientation;
        orientation = glm::normalize(orientation);

        // Update transform for rendering
        glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
        glm::mat4 r = glm::mat4_cast(orientation);
        glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);
        transform = t * r * s;

    }
   
    //clear force and torque
    force = glm::vec3(0);
    torque = glm::vec3(0);
   
    
}
void Object::Interpolate(float alpha)
{
    // Decompose previousTransform
    glm::vec3 prevPos, prevScale;
    glm::quat prevRot;
    glm::vec3 skew; glm::vec4 perspective;
    glm::decompose(previousTransform, prevScale, prevRot, prevPos, skew, perspective);

    // Decompose current transform
    glm::vec3 currPos, currScale;
    glm::quat currRot;
    glm::decompose(transform, currScale, currRot, currPos, skew, perspective);

    // Interpolate
    glm::vec3 interpPos = glm::mix(prevPos, currPos, alpha);
    glm::quat interpRot = glm::slerp(prevRot, currRot, alpha);
    glm::vec3 interpScale = glm::mix(prevScale, currScale, alpha);

    // Build visual transform
    renderTransform = glm::translate(glm::mat4(1.0f), interpPos) * glm::mat4_cast(interpRot) * glm::scale(glm::mat4(1.0f), interpScale);
}
void Object::drawForceDirection(glm::vec3 intersect, glm::vec3 dir)
{
    Debug::DrawLine(intersect, dir, 10.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 0, 0, 1));
}

void Object::drawAnglularAxis(glm::vec3 angularDir)
{

   

    glm::vec3 dir = glm::normalize(angularDir); 
    float displayLength = 30.0f;                 // arbitrary scale for visibility

    glm::vec3 start = centerOfMass + dir * displayLength;
    glm::vec3 end = centerOfMass - dir * displayLength;

    Debug::DrawLine(start, end, 5.0f, glm::vec4(1, 0, 0, 1), glm::vec4(0, 0, 1, 1));
}

void Object::RotateAxisAngle(const glm::vec3& axis, float angle)
{
    glm::vec3 ax = glm::normalize(axis);
    glm::quat dq = glm::angleAxis(angle, ax);

    // update orientation state (accumulative)
    orientation = dq * orientation;
    orientation = glm::normalize(orientation);

    // rebuild transform
    glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
    glm::mat4 r = glm::mat4_cast(orientation);
    glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);
    transform = t * r * s;
}

void Object::calculateCenterOfMass()
{
    // Reset
    totalMeshArea = 0.0f;
    centerOfMass = glm::vec3(0.0f);
    totalMass = 0.0f;

    // Build object transform (without rotation for CoM)
    glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
    glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);
    glm::mat4 transformMatrix = t * s;

    // First, compute total mesh area
    for (auto& tri : triangles)
    {
        glm::vec3 v0 = glm::vec3(transformMatrix * glm::vec4(tri.verticies[0], 1.0f));
        glm::vec3 v1 = glm::vec3(transformMatrix * glm::vec4(tri.verticies[1], 1.0f));
        glm::vec3 v2 = glm::vec3(transformMatrix * glm::vec4(tri.verticies[2], 1.0f));

        totalMeshArea += 0.5f * glm::length(glm::cross(v1 - v0, v2 - v0));
    }

    // Now compute weighted center of mass
    for (auto& tri : triangles)
    {
        glm::vec3 v0 = glm::vec3(transformMatrix * glm::vec4(tri.verticies[0], 1.0f));
        glm::vec3 v1 = glm::vec3(transformMatrix * glm::vec4(tri.verticies[1], 1.0f));
        glm::vec3 v2 = glm::vec3(transformMatrix * glm::vec4(tri.verticies[2], 1.0f));

        float area = 0.5f * glm::length(glm::cross(v1 - v0, v2 - v0));
        float triMass = mass * (area / totalMeshArea);

        glm::vec3 centroid = (v0 + v1 + v2) / 3.0f;

        centerOfMass += centroid * triMass;
        totalMass += triMass;
    }

    centerOfMass /= totalMass;

}

void Object::calculateInertiaTensor()
{
   
    //      |   ∑mi​(yi^2​ + zi^2​)  −∑mi​xi​yi​       −∑mi​xi​zi​​          |
    // I =  |  −∑mi​xi​yi       ​∑mi​(xi^2​ + zi^2​)   −∑mi​yi​zi​​          |   
    //      |  −∑mi​xi​zi           ​−∑mi​yi​zi​     ∑mi​(xi^2​ + yi^2​)​    |
    // where i is for index

    for (auto& tri : triangles)
    {
        glm::vec3 r0 = tri.verticies[0] - centerOfMass;
        glm::vec3 r1 = tri.verticies[1] - centerOfMass;
        glm::vec3 r2 = tri.verticies[2] - centerOfMass;

        float area = 0.5f * glm::length(glm::cross(r1 - r0, r2 - r0));
        float triMass = mass * (area / totalMeshArea);

        auto addTri = [&](const glm::vec3& r)
            {
            inertiaTensor += triMass / 3.0f * glm::mat3
            (r.y * r.y + r.z * r.z, -r.x * r.y, -r.x * r.z,
             -r.x * r.y, r.x * r.x + r.z * r.z, -r.y * r.z,
             -r.x * r.z, -r.y * r.z, r.x * r.x + r.y * r.y
            );
            };

        addTri(r0); 
        addTri(r1); 
        addTri(r2);
    }

    inertiaTensorInv = glm::inverse(inertiaTensor);
}

void Object::drawObject() const
{
	/*Render::RenderDevice::Draw(modelID, transform);
    Render::RenderDevice::Draw(modelID, renderTransform);*/
}

void Object::SetScale(glm::vec3 _scale)
{
	scale = _scale;
	glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
	glm::mat4 r = glm::mat4(1.0f);
	glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);
	transform = t * r * s;
}


ObjectGlobalData* ObjectGlobalData::instance = nullptr;

ObjectGlobalData::ObjectGlobalData()
{
    ammountOfObjects = 0;
}

ObjectGlobalData& ObjectGlobalData::GetInstance()
{
    if (!instance)
    {
        instance = new ObjectGlobalData();
    }
    return *instance;
}

void ObjectGlobalData::ClearInstance()
{
    delete instance;
    instance = nullptr;
}
