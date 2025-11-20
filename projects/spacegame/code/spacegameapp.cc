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

    const int numLights = 3;
    Render::PointLightId lights[numLights];
    // Setup lights
    for (int i = 0; i < numLights; i++)
    {
        glm::vec3 translation = glm::vec3(0.0f + i * 10.0f, 0.0f + i * 10.0f, 0.0f);
        glm::vec3 color = glm::vec3(Core::RandomFloat(), Core::RandomFloat(), Core::RandomFloat());
        lights[i] = Render::LightServer::CreatePointLight(translation, color, Core::RandomFloat() * 4.0f, 1.0f + (15 + Core::RandomFloat() * 10.0f));
    }

    FFCam ffCam;
    MathRay ray;
    bool drawRay = false;
    static float angle = 0.0f; // keep angle across frames
    glm::vec3 SavedOrigin, SavedEnd;
    Physics::RayProperties rayProperties;
     
    glm::vec3 v0(-1.0f, 0.0f, -1.0f);
    glm::vec3 v1(1.0f, 0.0f , -1.0f);
    glm::vec3 v2(1.0f, 0.0f, 1.0f);
    glm::vec3 v3(-1.0f, 0.0f, 1.0f);
    std::vector<Quad> Quads;
    for (int i = 0; i < 10; i++)
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
       quad.UpdateAndDrawAABB();
       Quads.push_back(quad);
      
  

       
    }

    //Create obj physics
    std::vector<Object> objs;
    //create obj2 model
    std::vector<Object> objs2;
    for (int i = 0; i < 3; i++)
    {
        glm::vec4 color(
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX,
            1.0f);

        Object obj;
       
        std::string filePath = "assets/space/spaceship_physics.glb";
        obj.createObject(filePath);

        fx::gltf::Document doc;
        if (filePath.substr(filePath.find_last_of(".") + 1) == "glb")
            doc = fx::gltf::LoadFromBinary(filePath);
        else
            doc = fx::gltf::LoadFromText(filePath);

        Physics::LoadFromIndexBuffer(doc, obj.triangles, obj.aabb);
        obj.SetScale(glm::vec3(1.0f));
        auto newPos = glm::vec3(0.0f + i * 10.0f, 0.0f + i * 10.0f, 0.0f);
    
        for (auto& tri : obj.triangles)
        {
            tri.color = glm::vec4(0.0f, 1.0f, 1.0f, 1.0f);
            tri.og_color = glm::vec4(0.0f, 1.0f, 1.0f, 1.0f);
            tri.selectedColor = glm::vec4(1, 1, 0, 1);

        }
        obj.SetOBjectPosition(newPos);
        obj.UpdateAndDrawAABBObject();
        objs.push_back(obj);

       
        Object obj2;

        std::string filePath2 = "assets/space/spaceship.glb";
        obj2.createObject(filePath2);

        obj2.SetScale(glm::vec3(3.0f));

        obj2.SetOBjectPosition(newPos);
       
        objs2.push_back(obj2);




    }

    std::clock_t c_start = std::clock();
    double dt = 0.01667f;

    float timeSinceAxisChange = 0.0f;

    //RandomAxis
    glm::vec3 randomAxis;
    randomAxis = glm::normalize(glm::vec3(
        static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f,
        static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f,
        static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f));
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
        /*auto text = std::to_string(Quads[0].GetRotation().x) + " , " + std::to_string(Quads[0].GetRotation().y) + " , " + std::to_string(Quads[0].GetRotation().z);
        Debug::DrawDebugText(text.c_str(), glm::vec3(0), {1,0,0,1});*/

        //render draw
        //for (int i  = 0; i < Quads.size(); i++)
        //{
        //    
        //       
        //    for (auto& tri : Quads[i].triangles)
        //    {
        //        Debug::DrawTriangle(tri.verticies[0], tri.verticies[1], tri.verticies[2], tri.color, Debug::DoubleSided, 1.0f);
        //        tri.SetSelected(false);
        //    }
        //        

        //   Quads[i].SetRotation(glm::vec3(1.0f, 1.0f, 0.0f), angle = dt);
        //    Quads[i].UpdateAndDrawAABB();
        //    glm::vec3 start = Quads[i].plane.getOrigin();
        //    glm::vec3 end = start + Quads[i].plane.GetNormal();
        //    Debug::DrawLine(start, end, 1.0f, glm::vec4(1, 0, 0, 1), glm::vec4(1, 0, 0, 1));
        //     
        //}

      
        //int hitIndex = -1;
        //float closestDistance = std::numeric_limits<float>::max();

        ////closes one who get hit?
        //for (int i = 0; i < Quads.size(); i++)
        //{
        //   
        //    if (mouse->pressed[mouse->LeftButton])
        //    {
        //        ray = Physics::ScreenPointToRay(mouse->position, w, h);

        //        drawRay = true;
        //    }
        //    
        //    if (Physics::CheckRayHitAABB(Quads[i].aabb, ray, rayProperties))
        //    {
        //   
        //        if (rayProperties.AABBintersection.length() < closestDistance) // distance from ray origin to intersection
        //        {
        //            closestDistance = rayProperties.AABBintersection.length();
        //            hitIndex = i;
        //        }
        //    }
        //}
        //// After loop, only one hit
        //if (hitIndex != -1)
        //{
        //    Quads[hitIndex].CheckRayHit(Quads[hitIndex], ray, rayProperties);
        //    Debug::DrawLine(rayProperties.intersection, rayProperties.normalEnd, 3.0f, glm::vec4(0, 0, 1, 1), glm::vec4(0, 1, 1, 1));
        //}
        //if (mouse->released[mouse->LeftButton])
        //{

        //    drawRay = false;

        //}
      
        //render
        angle += dt;  
        timeSinceAxisChange += dt;
        
        // random -1..1 axis
        static glm::vec3 targetAxis;

        float blendSpeed =  2.0 * dt; // adjust speed (0.5–2.0 for slower/faster blend)
        if (timeSinceAxisChange >= 2.0f)
        {
            targetAxis += glm::normalize(glm::vec3(
                static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f,
                static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f,
                static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f));
            timeSinceAxisChange = 0.0f;
        }
      
        randomAxis = glm::normalize(glm::mix(randomAxis, targetAxis, blendSpeed));

        for (int i = 0; i < objs.size(); i++)
        {
           
            // Update AABB and render
            objs[i].UpdateAndDrawAABBObject();
            objs[i].drawObject();
            //objs[i].SetOBjectRotation(randomAxis, angle);
            
            

            // Draw triangles with transform applied
            for (auto& tri : objs[i].triangles)
            {
               
                glm::vec3 v0 = glm::vec3(objs[i].transform * glm::vec4(tri.verticies[0], 1.0f));
                glm::vec3 v1 = glm::vec3(objs[i].transform * glm::vec4(tri.verticies[1], 1.0f));
                glm::vec3 v2 = glm::vec3(objs[i].transform * glm::vec4(tri.verticies[2], 1.0f));

                // If this triangle was hit -> draw filled, else wireframe
                if (tri.selected)
                {
                    Debug::DrawTriangle(v0, v1, v2, glm::vec4(0, 1, 1, 1), Debug::Filled, 2.0f); 
                }
                else
                {
                    Debug::DrawTriangle(v0, v1, v2, tri.color, Debug::WireFrame, 1.0f);
                }
              
                tri.SetSelected(false);

            }  
        }
        for (int i = 0; i < objs2.size(); i++)
        {

           
           
           //drawOBJ
            objs2[i].drawObject();
            //objs2[i].SetOBjectRotation(randomAxis, angle);

            
        }


        int hitIndex = -1;
       float closestDistance = std::numeric_limits<float>::max(); // give largest value as possible

       //closes one who get hit?
       for (int i = 0; i < objs.size(); i++)
       {
           if (mouse->held[mouse->LeftButton])
           {
               ray = Physics::ScreenPointToRay(mouse->position, w, h);
               drawRay = true;
           }
          
           if (Physics::CheckRayHitAABB(objs[i].aabb, ray, rayProperties))
           {
               // distance from ray origin to intersection
               if (rayProperties.AABBintersection.length() < closestDistance) 
               {
                   closestDistance = rayProperties.AABBintersection.length();
                   hitIndex = i;
               }
           }
       }
       
       // After loop, only one hit
       if (hitIndex != -1)
       {
           if(objs[hitIndex].CheckRayHit(objs[hitIndex], ray, rayProperties))
               Debug::DrawLine(rayProperties.intersection, rayProperties.normalEnd, 3.0f, glm::vec4(0, 0, 1, 1), glm::vec4(0, 1, 1, 1));
       }
       if (mouse->released[mouse->LeftButton])
       {
           drawRay = false;
       }
        
        //draw Lazer
        if (drawRay)
        {
            Debug::DrawLine(SavedOrigin = ray.GetOrigin(), SavedEnd = ray.GetOrigin() + ray.GetDirection() * ray.GetRayLength(), 1.0f, glm::vec4(0.0f, 1.0f, 0.0f, 1.0f), glm::vec4(1.0f, 0.0f, 0.0f, 1.0f)); 
        }
        else
        {
            Debug::DrawLine(SavedOrigin, SavedEnd, 1.0f, glm::vec4(0.0f, 1.0f, 0.0f, 1.0f), glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
        }
      
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