#include <config.h>
#include "object.h"
#include "gltf.h"

Object::Object()
{
	aabb = Physics::AABB();
	model = 0;
	rotation = glm::vec3(0);
	position = glm::vec3(0);
	transform = glm::mat4(1);
	scale = glm::vec3(1.0f);

}

Object::~Object()
{
}
void Object::createObject(std::string filePath)

{    // Load the model for rendering
    model = Render::LoadModel(filePath);


    

}
void Object::drawObject() const
{
	Render::RenderDevice::Draw(model, transform);
}

void Object::SetScale(glm::vec3 _scale)
{
	scale = _scale;
	glm::mat4 t = glm::translate(glm::mat4(1.0f), position);
	glm::mat4 r = glm::mat4(1.0f);
	glm::mat4 s = glm::scale(glm::mat4(1.0f), scale);
	transform = t * r * s;
}


