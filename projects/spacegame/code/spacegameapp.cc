//------------------------------------------------------------------------------
// spacegameapp.cc
// (C) 2022 Individual contributors, see AUTHORS file
//------------------------------------------------------------------------------
#include "config.h"
#include "spacegameapp.h"
#include <cstring>
#include "imgui.h"
#include "render/renderdevice.h"
#include "render/shaderresource.h"
#include <vector>
#include "render/textureresource.h"
#include "render/model.h"
#include "render/cameramanager.h"
#include "render/lightserver.h"
#include "render/debugrender.h"
#include "core/random.h"
#include "render/input/inputserver.h"
#include "core/cvar.h"
#include <chrono>
#include "FFCam.h"
#include "render/mathray.h"
#include <render/quad.h>
#include "render/betterphysics.h"
#include "string.h"
#include "render/object.h"

using namespace Display;
using namespace Render;

namespace Game
{

//------------------------------------------------------------------------------
/**
*/
SpaceGameApp::SpaceGameApp()
{
    // empty
}

//------------------------------------------------------------------------------
/**
*/
SpaceGameApp::~SpaceGameApp()
{
	// empty
}

//------------------------------------------------------------------------------
/**
*/
bool
SpaceGameApp::Open()
{
	App::Open();
	this->window = new Display::Window;
    this->window->SetSize(2560, 1440);

    if (this->window->Open())
	{
		// set clear color to gray
		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

        RenderDevice::Init();

		// set ui rendering function
		this->window->SetUiRender([this]()
		{
			this->RenderUI();
		});
        
        return true;
	}
	return false;
}

//------------------------------------------------------------------------------
/**
*/
void
SpaceGameApp::Run()
{
    int w;
    int h;
    this->window->GetSize(w, h);
    glm::mat4 projection = glm::perspective(glm::radians(90.0f), float(w) / float(h), 0.01f, 1000.f);
    Camera* cam = CameraManager::GetCamera(CAMERA_MAIN);
    cam->projection = projection;
    
    // Setup skybox
    std::vector<const char*> skybox
    {
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png"
    };
    TextureResourceId skyboxId = TextureResource::LoadCubemap("skybox", skybox, true);
    RenderDevice::SetSkybox(skyboxId);
    
    Input::Keyboard* kbd = Input::GetDefaultKeyboard();

    const int numLights = 40;
    Render::PointLightId lights[numLights];
    // Setup lights
    for (int i = 0; i < numLights; i++)
    {
        glm::vec3 translation = glm::vec3(
            Core::RandomFloatNTP() * 20.0f,
            Core::RandomFloatNTP() * 20.0f,
            Core::RandomFloatNTP() * 20.0f
        );
        glm::vec3 color = glm::vec3(
            Core::RandomFloat(),
            Core::RandomFloat(),
            Core::RandomFloat()
        );
        lights[i] = Render::LightServer::CreatePointLight(translation, color, Core::RandomFloat() * 4.0f, 1.0f + (15 + Core::RandomFloat() * 10.0f));
    }

    FFCam ffCam;
    MathRay ray;
    std::vector<Quad> Quads;
    Object obj;
    obj.createObject("assets/space/Cube.glb");
    obj.SetScale(glm::vec3(0.1f) );

    glm::vec3 v0(-1.0f, 0.0f, -1.0f);
    glm::vec3 v1(1.0f, 0.0f , -1.0f);
    glm::vec3 v2(1.0f, 0.0f, 1.0f);
    glm::vec3 v3(-1.0f, 0.0f, 1.0f);

    for (int i = 0; i < 2; i++)
    {


        glm::vec4 color(
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX,
            1.0f
        );

      
        auto quad = Quad(v0, v1, v2, v3, color);
       
       
      
       auto newPos = glm::vec3(0.0f + i * 10.0f, 0.0f + i * 10.0f, 0.0f);
       quad.SetPosition(newPos);

       Quads.push_back(quad);
  

       
    }
    bool drawRay = false;
    float angle = 1.0f;
  

    glm::vec3 SavedOrigin, SavedEnd;
    Physics::RayProperties rayProperties;

    std::clock_t c_start = std::clock();
    double dt = 0.01667f;

    // game loop
    while (this->window->IsOpen())
	{
        auto timeStart = std::chrono::steady_clock::now();
		glClear(GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
        
        this->window->Update();

        // Handle cursor show/hide on right mouse button
        Input::Mouse* mouse = Input::GetDefaultMouse();
        if (mouse->pressed[mouse->RightButton])
        {
            glfwSetInputMode(this->window->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        }
        if (mouse->released[mouse->RightButton])
        {
            glfwSetInputMode(this->window->window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        }

        if (kbd->pressed[Input::Key::Code::End])
        {
            ShaderResource::ReloadShaders();
        }

        //Free flight cam
        ffCam.Update(dt);

        //text
        auto text = std::to_string(Quads[0].GetRotation().x) + " , " + std::to_string(Quads[0].GetRotation().y) + " , " + std::to_string(Quads[0].GetRotation().z);
        Debug::DrawDebugText(text.c_str(), glm::vec3(0), {1,0,0,1});

        for (int i = 0; i < Quads.size(); i++)
        {
            Debug::DrawQuad(Quads[i].v0, Quads[i].v1, Quads[i].v2, Quads[i].v3, Quads[i].color, Debug::DoubleSided, 1.0f);
            Quads[i].SetRotation(glm::vec3(1.0f, 1.0f, 0.0f), angle = dt);
            glm::vec3 start = Quads[i].plane.getOrigin();
            glm::vec3 end = start + Quads[i].plane.GetNormal();
            Debug::DrawLine(start, end, 1.0f, glm::vec4(1, 0, 0, 1), glm::vec4(1, 0, 0, 1));
             
        }
        //who get hit?
        for (int i = 0; i < Quads.size(); i++)
        {
           
            if (mouse->held[mouse->LeftButton])
            {
                ray = Physics::ScreenPointToRay(mouse->position, w, h);
                drawRay = true;
                Quads[i].isHit = Physics::CheckRayHit(Quads[i], ray, rayProperties);
            }
            
            if (Quads[i].isHit)
            {
                Debug::DrawLine(rayProperties.intersection, rayProperties.normalEnd, 3.0f, glm::vec4(0, 1, 1, 1), glm::vec4(1, 0, 1, 1));
                break;
            }

        }
        if (mouse->released[mouse->LeftButton])
        {

            drawRay = false;

        }
        
        //drawLazer
        if (drawRay)
        {
            Debug::DrawLine(SavedOrigin = ray.GetOrigin(), SavedEnd = ray.GetOrigin() + ray.GetDirection() * ray.GetRayLength(), 1.0f, glm::vec4(0.0f, 1.0f, 0.0f, 1.0f), glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
           
        }
        else
        {
            Debug::DrawLine(SavedOrigin, SavedEnd, 1.0f, glm::vec4(0.0f, 1.0f, 0.0f, 1.0f), glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
        }
        obj.aabb.drawBox(glm::mat4(1));
        obj.drawObject();
       
    

     

  
           
        // Draw some debug text

      
        // Execute the entire rendering pipeline
        RenderDevice::Render(this->window, dt);

		// transfer new frame to window
		this->window->SwapBuffers();

        auto timeEnd = std::chrono::steady_clock::now();
        dt = std::min(0.04, std::chrono::duration<double>(timeEnd - timeStart).count());

        if (kbd->pressed[Input::Key::Code::Escape])
            this->Exit();
	}
}

//------------------------------------------------------------------------------
/**
*/
void
SpaceGameApp::Exit()
{
    this->window->Close();
}

//------------------------------------------------------------------------------
/**
*/
void
SpaceGameApp::RenderUI()
{
	if (this->window->IsOpen())
	{
        ImGui::Begin("Debug");
        Core::CVar* r_draw_light_spheres = Core::CVarGet("r_draw_light_spheres");
        int drawLightSpheres = Core::CVarReadInt(r_draw_light_spheres);
        if (ImGui::Checkbox("Draw Light Spheres", (bool*)&drawLightSpheres))
            Core::CVarWriteInt(r_draw_light_spheres, drawLightSpheres);
        
        Core::CVar* r_draw_light_sphere_id = Core::CVarGet("r_draw_light_sphere_id");
        int lightSphereId = Core::CVarReadInt(r_draw_light_sphere_id);
        if (ImGui::InputInt("LightSphereId", (int*)&lightSphereId))
            Core::CVarWriteInt(r_draw_light_sphere_id, lightSphereId);
        
        ImGui::End();

        Debug::DispatchDebugTextDrawing();
	}
}

} // namespace Game