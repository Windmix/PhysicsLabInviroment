#pragma once
#include "renderdevice.h"
#include "model.h"
#include "betterphysics.h"



class Object
{
public:
	Physics::AABB aabb;
	std::vector<glm::vec3> vertices;
	Render::ModelId model;
	glm::vec3 position;
	glm::vec3 rotation;
	glm::mat4 transform;
	glm::vec3 scale;

	Object();
	~Object();

	void createObject(std::string filePath);
	void drawObject() const;
	void SetScale(glm::vec3 s);

};

