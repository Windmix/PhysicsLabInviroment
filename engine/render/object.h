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
	glm::vec3 rotation;
	glm::mat4 transform;
	glm::vec3 scale;
	MathPlane plane;

	// Store triangles similar to Quad
	std::vector<Physics::ColliderMesh::Triangle> triangles;
	Object();
	~Object();

	//render
	void createObject(std::string filePath);
	void UpdateAndDrawAABBObject();
	void drawObject() const;
	void DrawAABBOnObject();
	bool CheckRayHit(Object& myObj, MathRay& ray, Physics::RayProperties& rayproperties);
	bool IsPointInsideMesh(const glm::vec3& point);

	//pos and rot
	void SetOBjectRotation(glm::vec3 direction, float angle);
	void SetOBjectPosition(glm::vec3 newPosition);


	void SetScale(glm::vec3 s);
private:

};

	
	

