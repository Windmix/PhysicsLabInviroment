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
    this->window->SetSize(1920, 1080);

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


    ObjectGlobalData& globalData = ObjectGlobalData::GetInstance();
    globalData.ammountOfObjects = 4;

    for (int i = 0; i < globalData.ammountOfObjects; i++)
    {
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

        glm::vec3 startPos = glm::vec3(1 + i * 10.0f, 1 + i * 10.0f, 1 + i * 10.0f);
        for (auto& tri : obj.triangles)
        {
            tri.color = glm::vec4(0.0f, 1.0f, 1.0f, 1.0f);
            tri.og_color = glm::vec4(0.0f, 1.0f, 1.0f, 1.0f);
            tri.selectedColor = glm::vec4(1, 1, 0, 1);

        }
        if (i == 0)
        {
            obj.mass = 50000.0f;
        }
    
        obj.SetOBjectPosition(startPos); 
        obj.UpdateAABBObject();

        obj.calculateCenterOfMass();
        obj.calculateInertiaTensor();
        globalData.phyiscsObjects.push_back(obj);

       
     /*   Object obj2;
        std::string filePath2 = "assets/space/Cube.glb";
        obj2.createObject(filePath2);
        obj2.SetScale(glm::vec3(1.0f));
        obj2.SetOBjectPosition(newPos);
        obj2.calculateCenterOfMass();
        obj2.calculateInertiaTensor();
        globalData.phyiscsObjects2.push_back(obj2);*/




    }

    std::clock_t c_start = std::clock();
    double dt = 0.01667f;
    //------------------------------------------------------------------------------------------------------------
    //cvar create
    //physics
    auto CvarGravity = Core::CVarCreate(Core::CVarType::CVar_Int, "r_apply_gravity", "0");
    auto Cvar_gravity_direction_x = Core::CVarCreate(Core::CVarType::CVar_Float, "r_gravity_direction_x", "0");
    auto Cvar_gravity_direction_y = Core::CVarCreate(Core::CVarType::CVar_Float, "r_gravity_direction_y", "0");
    auto Cvar_gravity_direction_z = Core::CVarCreate(Core::CVarType::CVar_Float, "r_gravity_direction_z", "0");
    Core::CVar* r_freeze_rot = Core::CVarCreate(Core::CVarType::CVar_Int, "r_freeze_rot", "0");
    int freezeRotBool = Core::CVarReadInt(r_freeze_rot);
    Core::CVar* r_freeze_pos = Core::CVarCreate(Core::CVarType::CVar_Int, "r_freeze_pos", "0");
    int freezePosBool = Core::CVarReadInt(r_freeze_pos);

    //Draws
    auto CvarDrawAABBid = Core::CVarCreate(Core::CVarType::CVar_Int, "r_draw_AABB_id", "0");
    auto CvarDrawAABB = Core::CVarCreate(Core::CVarType::CVar_Int, "r_draw_AABB", "0");
    //------------------------------------------------------------------------------------------------------------
    while (this->window->IsOpen())// game loop
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

        // random -1..1 axis
        static glm::vec3 targetAxis;

        //------------------------------------------------------------------------------------------------------------
        //cvar create
        //physics
        int applyGravity = Core::CVarReadInt(CvarGravity);
        float gravity_direction_x = Core::CVarReadFloat(Cvar_gravity_direction_x);
        float gravity_direction_y = Core::CVarReadFloat(Cvar_gravity_direction_y);
        float gravity_direction_z = Core::CVarReadFloat(Cvar_gravity_direction_z);
        glm::vec3 gravityDirection(gravity_direction_x, gravity_direction_y, gravity_direction_z);

        int IdDrawABB = Core::CVarReadInt(CvarDrawAABBid);
        bool cvarDrawAABB = Core::CVarReadInt(CvarDrawAABB);
        //-------------------------------------------------------------------------------------------------------------

        

        for (int i = 0; i < globalData.phyiscsObjects.size(); i++) // renderloop
        {

            // Physics first
            globalData.phyiscsObjects[i].calculateCenterOfMass();

            if(freezeRotBool == 0)
                globalData.phyiscsObjects[i].drawAnglularAxis(globalData.phyiscsObjects[i].angularVelocity);

            if (applyGravity == 1)
            {
                globalData.phyiscsObjects[i].ApplyGravityForce(gravityDirection);
            }
            
            if (freezeRotBool == 0 || freezeRotBool == 0)
            {
                globalData.phyiscsObjects[i].Integrate(dt);
            }
            // Update AABB after moving
            globalData.phyiscsObjects[i].UpdateAABBObject();


           
            //also collect all objects
            std::vector<Physics::AABB> sweepAABBEntries;


            for (auto& obj : globalData.phyiscsObjects)
            {
                sweepAABBEntries.push_back(obj.aabb);
            }
            //contains all pairs of AABBs that overlap in 3D space
            //Broad-phase collision detection 
            std::vector<std::pair<Physics::AABB, Physics::AABB>> candidates;
            candidates = Physics::PlaneSweepOverlaps(sweepAABBEntries);

            //Narrow-phase collision detection and response
            for (auto& pair : candidates)
            {
                Object* objA = pair.first.owner;
                Object* objB = pair.second.owner;

                // Broad-phase
                if (Physics::CheckAABBCollision(pair.first, pair.second))
                {
                    // Narrow-phase
                    const std::vector<glm::vec3>& vertsA = objA->GetWorldVertices();
                    const std::vector<glm::vec3>& vertsB = objB->GetWorldVertices();


                    if (Physics::GJK_Intersect(vertsA, vertsB))
                    {
                        // Collision detected
                        // TODO: call EPA to get penetration vector
                        // For now, apply simple separation impulse

                     /*   glm::vec3 collisionDir = glm::normalize(objB->centerOfMass - objA->centerOfMass);
                        float impulseStrength = 30.0f;

                        objA->ApplyForce(-collisionDir * impulseStrength, objA->centerOfMass);
                        objB->ApplyForce(collisionDir * impulseStrength, objB->centerOfMass);*/
                    }

                
                }
            }

            // Draw AABB
            if (cvarDrawAABB)
            {
                globalData.phyiscsObjects[IdDrawABB].DrawAABBOnObject();
            }

            // Draw the object at new position
            globalData.phyiscsObjects[i].drawObject();

            // Draw triangles with applied transform
            for (auto& tri : globalData.phyiscsObjects[i].triangles)
            {
                glm::vec3 v0 = glm::vec3(globalData.phyiscsObjects[i].transform * glm::vec4(tri.verticies[0], 1.0f));
                glm::vec3 v1 = glm::vec3(globalData.phyiscsObjects[i].transform * glm::vec4(tri.verticies[1], 1.0f));
                glm::vec3 v2 = glm::vec3(globalData.phyiscsObjects[i].transform * glm::vec4(tri.verticies[2], 1.0f));

                //if triangle is hit by ray, it will recolor it
                if (tri.selected)
                    Debug::DrawTriangle(v0, v1, v2, glm::vec4(0, 1, 1, 1), Debug::AlwaysOnTop, 2.0f);
                else
                    Debug::DrawTriangle(v0, v1, v2, tri.color, Debug::WireFrame, 1.0f);

                tri.SetSelected(false);
            }
        }
        //for (int i = 0; i < globalData.phyiscsObjects2.size(); i++)
        //{


        //    globalData.phyiscsObjects2[i].Integrate(dt);
        //    //Draw model
        //    globalData.phyiscsObjects2[i].drawObject();


        //}

        //Newton’s Law of Gravity, using first cube as a planet with big mass.
        for (int i = 0; i < globalData.phyiscsObjects.size(); i++)
        {
            for (int j = 0; j < globalData.phyiscsObjects.size(); j++)
            {
                if (i == j) continue;
                globalData.phyiscsObjects[i].ApplyGravityFrom(globalData.phyiscsObjects[j]);
                globalData.phyiscsObjects[j].Integrate(dt);

            }
            
        }

        int hitIndex = -1;
        float closestDistance = std::numeric_limits<float>::max(); // give largest value as possible

        
        if (mouse->held[mouse->LeftButton])
        {
            ray = Physics::ScreenPointToRay(mouse->position, w, h);
            drawRay = true;
        }
        if (mouse->released[mouse->LeftButton])
        {
            drawRay = false;

            //reset the ray's values
            ray = MathRay();
        }

        for (int i = 0; i < globalData.phyiscsObjects.size(); i++)
        {
            //closes ABB who get hit?
            if (Physics::CheckRayHitAABB(globalData.phyiscsObjects[i].aabb, ray, rayProperties))
            {
                // distance from ray origin to intersection
                if (rayProperties.AABBintersection.length() < closestDistance)
                {
                    //need this check for if its inside aabb and want to intersect another object
                    if (globalData.phyiscsObjects[i].CheckRayHit(globalData.phyiscsObjects[i], ray, rayProperties))
                    {
                        closestDistance = rayProperties.AABBintersection.length();
                        hitIndex = i;
                    }
                }
            }
        }
      
        // After loop, only one hit
        if (hitIndex != -1 && drawRay != false)
        {
          

            //draw normal
            if (globalData.phyiscsObjects[hitIndex].CheckRayHit(globalData.phyiscsObjects[hitIndex], ray, rayProperties))
            {
                Debug::DrawLine(rayProperties.intersection, rayProperties.normalEnd, 3.0f, glm::vec4(0, 0, 1, 1), glm::vec4(0, 1, 1, 1));

              

                //drawing the direction later
                glm::vec3 forceDirection = glm::normalize(ray.GetDirection());
                globalData.phyiscsObjects[hitIndex].forceDirection = forceDirection;
                globalData.phyiscsObjects[hitIndex].storedHitindex = hitIndex;
                globalData.phyiscsObjects[hitIndex].ApplyForce(forceDirection* globalData.phyiscsObjects[hitIndex].forceMagnitude, rayProperties.intersection);
                
                //center of mass drawing
                Debug::DrawLine(rayProperties.intersection, globalData.phyiscsObjects[hitIndex].centerOfMass, 3.0f, glm::vec4(1, 0, 1, 1), glm::vec4(1, 1, 1, 1));

               
              
            }
           


            // globalData.phyiscsObjects2[hitIndex].forceDirection = forceDirection;
           // globalData.phyiscsObjects2[hitIndex].ApplyForce(forceDirection * globalData.phyiscsObjects2[hitIndex].forceMagnitude, rayProperties.intersection);
        }

        //for (int i = 0; i < globalData.phyiscsObjects.size(); i++)
        //{
        //    if (globalData.phyiscsObjects[i].storedHitindex != -1)
        //    {
        //        
        //       // globalData.phyiscsObjects2[i].drawForceDirection(globalData.phyiscsObjects2[i].position, (globalData.phyiscsObjects2[i].position + globalData.phyiscsObjects2[i].forceDirection * glm::length(globalData.phyiscsObjects2[i].accumulatedForce * 0.1f)));
        //        //debug draw How much force applied on
        //        glm::vec3 currentForce = globalData.phyiscsObjects[i].forceDirection * globalData.phyiscsObjects[i].forceMagnitude;
        //      /*  std::ostringstream ss;
        //        ss << std::fixed << std::setprecision(2);
        //        ss << "Current applied Force: " << glm::length(currentForce) << " N \n"
        //           << "Accumulated Force: " << glm::length(globalData.phyiscsObjects[i].accumulatedForce) << " N";
        //        Debug::DrawDebugText(ss.str().c_str(),globalData.phyiscsObjects[i].position,glm::vec4(1, 1, 1, 1)
        //        );*/
        //    }
        //}
        
        //draw Lazer
        if (drawRay)
        {
            Debug::DrawLine(SavedOrigin = ray.GetOrigin(), SavedEnd = ray.GetOrigin() + ray.GetDirection() * ray.GetRayLength(), 1.0f, glm::vec4(0.0f, 1.0f, 0.0f, 1.0f), glm::vec4(1.0f, 0.0f, 0.0f, 1.0f)); 
        }
        else
        {
           // Debug::DrawLine(SavedOrigin, SavedEnd, 1.0f, glm::vec4(0.0f, 1.0f, 0.0f, 1.0f), glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
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
    ObjectGlobalData& globalData = ObjectGlobalData::GetInstance();
    globalData.ClearInstance();
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
        ObjectGlobalData& globalData = ObjectGlobalData::GetInstance();
      

        ImGui::Begin("Debug");
        Core::CVar* r_draw_light_spheres = Core::CVarGet("r_draw_light_spheres");
        int drawLightSpheres = Core::CVarReadInt(r_draw_light_spheres);
        if (ImGui::Checkbox("Draw Light Spheres", (bool*)&drawLightSpheres))
            Core::CVarWriteInt(r_draw_light_spheres, drawLightSpheres);
        
        Core::CVar* r_draw_light_sphere_id = Core::CVarGet("r_draw_light_sphere_id");
        int lightSphereId = Core::CVarReadInt(r_draw_light_sphere_id);
        if (ImGui::InputInt("LightSphereId", (int*)&lightSphereId))
            Core::CVarWriteInt(r_draw_light_sphere_id, lightSphereId);


        //freeze pos
        Core::CVar* r_freeze_pos = Core::CVarGet("r_freeze_pos");
        int freezePos = Core::CVarReadInt(r_freeze_pos);
        if (ImGui::Checkbox("freezePos", (bool*)&freezePos))
        {
            Core::CVarWriteInt(r_freeze_pos, freezePos);

        }
        Core::CVar* r_freeze_rot = Core::CVarGet("r_freeze_rot");
        int freezeRot = Core::CVarReadInt(r_freeze_rot);
        if (ImGui::Checkbox("freezeRot", (bool*)&freezeRot))
        {
            Core::CVarWriteInt(r_freeze_rot, freezeRot);

        }
      

        //gravity
        Core::CVar * r_apply_gravity = Core::CVarGet("r_apply_gravity");
        int applyGravity = Core::CVarReadInt(r_apply_gravity);
        if (ImGui::Checkbox("Apply Gravity", (bool*)&applyGravity))
        {

            Core::CVarWriteInt(r_apply_gravity, applyGravity);
        }

        // --- Gravity direction ---
        float gravityDir[3];
        gravityDir[0] = Core::CVarReadFloat(Core::CVarGet("r_gravity_direction_x"));
        gravityDir[1] = Core::CVarReadFloat(Core::CVarGet("r_gravity_direction_y"));
        gravityDir[2] = Core::CVarReadFloat(Core::CVarGet("r_gravity_direction_z"));
        if (ImGui::InputFloat3("Gravity Direction", gravityDir))
        {
            Core::CVarWriteFloat(Core::CVarGet("r_gravity_direction_x"), gravityDir[0]);
            Core::CVarWriteFloat(Core::CVarGet("r_gravity_direction_y"), gravityDir[1]);
            Core::CVarWriteFloat(Core::CVarGet("r_gravity_direction_z"), gravityDir[2]);
        }

        Core::CVar* r_draw_AABB = Core::CVarGet("r_draw_AABB");
        int DrawAABB = Core::CVarReadInt(r_draw_AABB);
        if (ImGui::Checkbox("allow draw AABB", (bool*)&DrawAABB))
        {
            Core::CVarWriteInt(r_draw_AABB, DrawAABB);

        }

        Core::CVar* r_draw_AABB_id = Core::CVarGet("r_draw_AABB_id");
        int drawAABBId = Core::CVarReadInt(r_draw_AABB_id);
        if (ImGui::InputInt("drawAABBId", (int*)&drawAABBId))
        {

            //draw AABBgiv
            Core::CVarWriteInt(r_draw_AABB_id, drawAABBId);
            if (drawAABBId > (globalData.ammountOfObjects - 1))
            {

                Core::CVarWriteInt(r_draw_AABB_id, 0);
            }
            else if (drawAABBId < 0)
            {
                Core::CVarWriteInt(r_draw_AABB_id, (globalData.ammountOfObjects - 1));
            }
          
         
           
        }
          


        ImGui::End();
        Debug::DispatchDebugTextDrawing();
	}
}

} // namespace Game