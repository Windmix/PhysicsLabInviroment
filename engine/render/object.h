#pragma once
#include "renderdevice.h"
#include "model.h"
#include "betterphysics.h"
#include "mathplane.h"


class Object
{
public:
	Physics::AABB aabb;
	Render::Model Model;
	Render::ModelId modelID;
	glm::vec3 position;
	glm::mat4 transform;
	glm::vec3 scale;

	MathPlane plane;

	// Store triangles
	std::vector<Physics::ColliderMesh::Triangle> triangles;

	// Physics
	// 
	// --- Linear ---
	glm::vec3 velocity;
	glm::vec3 accumulatedForce;
	glm::vec3 force;

	// --- Angular ---
	glm::vec3 angularVelocity;
	glm::vec3 torque;
	glm::vec3 accumulatedTorque;


	glm::vec3 forceDirection;
	glm::vec3 acceleration;

	// Orientation
	glm::quat orientation;

	// Store the applied force
	float forceMagnitude;
	int storedHitindex;

	//G from the gravity formula
	float G = 0.001f;

	// Inertia
	float mass;
	float totalMass;
	glm::vec3 centerOfMass;

	glm::mat3 inertiaTensor;
	glm::mat3 inertiaTensorInv;
	float totalMeshArea;

	Object();
	~Object();

	//render
	void createObject(std::string filePath);
	void UpdateAABBObject();
	void drawObject() const;
	void DrawAABBOnObject();

	//ray checks
	bool CheckRayHit(Object& myObj, MathRay& ray, Physics::RayProperties& rayproperties);
	bool IsPointInsideMesh(const glm::vec3& point);

	//pos and rot
	void SetOBjectRotation(glm::vec3 direction, float angle);
	void SetOBjectPosition(glm::vec3 newPosition);
	void ApplyGravityForce(glm::vec3 force);
	void ApplyGravityFrom(Object& other);

	// Physics
	void ApplyForce(const glm::vec3& force, glm::vec3& forcehitPoint); // apply force
	void Integrate(float dt);
	void drawForceDirection(glm::vec3 intersect, glm::vec3 dir);
	void drawAnglularAxis(glm::vec3 angularDir);
	void RotateAxisAngle(const glm::vec3& axis, float angle);
	void calculateCenterOfMass(); // calculate center of mass with weighted vertices
	void calculateInertiaTensor();



	void SetScale(glm::vec3 s);
private:

};
class ObjectGlobalData
{
	// Private constructor - prevents external instantiation
	ObjectGlobalData();




	static ObjectGlobalData* instance;  // Static instance pointer
public:


	static ObjectGlobalData& GetInstance();
	void ClearInstance();

	std::vector<Object> phyiscsObjects;
	std::vector<Object> Objects;

	int ammountOfObjects;
};


	
	

