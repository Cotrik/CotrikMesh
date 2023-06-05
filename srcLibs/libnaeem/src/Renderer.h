#ifndef RENDERER_H_
#define RENDERER_H_

#include "Mesh.h"
#include <thread>

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkPolygon.h>
#include <vtkQuad.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCallbackCommand.h>
#include <vtkRendererCollection.h>

enum CustomEventIds {
    RENDER_EVENT = vtkCommand::UserEvent + 1
};

class vtkRenderCallback_ : public vtkCommand {
    public:
        static vtkRenderCallback_* New() {return new vtkRenderCallback_;}

        virtual void Execute(vtkObject* caller, unsigned long eventId, void* callData) {
            if (mesh_ == nullptr) return;
            vtkRenderWindowInteractor* interactor_ = static_cast<vtkRenderWindowInteractor*>(caller);
            vtkSmartPointer<vtkRenderer> renderer_ = interactor_->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

            for (auto& v: mesh_->V) {
                points->InsertNextPoint(v.x, v.y, v.z);
            }

            for (auto& f: mesh_->F) {
                if (f.Vids.size() != 4 || f.Vids.empty()) continue;
                vtkSmartPointer<vtkQuad> quad = vtkQuad::New();
                for (int id = 0; id < f.Vids.size(); id++) quad->GetPointIds()->SetId(id, f.Vids[id]);
                cells->InsertNextCell(quad);
            }

            polyData->SetPoints(points);
            polyData->SetPolys(cells);

            vtkPolyDataMapper* mapper_ = vtkPolyDataMapper::New();
            mapper_->SetInputData(polyData);
            
            vtkActor* actor_ = vtkActor::New();
            actor_->SetMapper(mapper_);
            {
                vtkProperty* property_ = vtkProperty::New();
                property_->SetEdgeVisibility(true);
                actor_->SetProperty(property_);
            }
    
            renderer_->RemoveAllViewProps();
            renderer_->AddActor(actor_);
            // if (!render_) {
                // interactor_->Start();
                // render_ = true;
            // } else {
                interactor_->Render();
            // }
        }

        Mesh* mesh_;
        bool render_ = false;
};

class Renderer {
    public:
        Renderer() {}
        Renderer(Mesh& mesh) {
            mesh_ = &mesh;
            Begin();
        }
        ~Renderer() {}

        void SetMesh(Mesh& mesh) {
            mesh_ = &mesh;
            Begin();
        }

        void Render() {
            render_ = true;
        }
    private:
        Mesh* mesh_;
        vtkSmartPointer<vtkRenderer> renderer_;
        vtkSmartPointer<vtkRenderWindow> renderWindow_;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor_;
        bool render_ = false;
        bool isRendered_ = false;

        void Begin() {
            
            std::cout << "Starting thread for rendering..." << std::endl;
            std::thread t(&Renderer::Update, this);
            // std::thread t([&](){
            //     renderer_ = vtkRenderer::New();
            //     renderWindow_ = vtkRenderWindow::New();
            //     renderWindowInteractor_ = vtkRenderWindowInteractor::New();
            //     renderer_->SetBackground(0.8, 0.8, 0.8);
            //     renderWindow_->AddRenderer(renderer_);
            //     renderWindow_->SetSize(640, 480);
            //     renderWindowInteractor_->SetRenderWindow(renderWindow_);
            //     renderWindowInteractor_->Initialize();
            //     renderWindowInteractor_->SetInteractorStyle(vtkInteractorStyleTrackballCamera::New());
            //     vtkNew<vtkRenderCallback_> renderCallback_;
            //     renderCallback_->mesh_ = mesh_;
            //     renderWindowInteractor_->AddObserver(CustomEventIds::RENDER_EVENT, renderCallback_);
            //     renderWindowInteractor_->ProcessEvents();
            //     while (true) {
            //         if (render_) {
            //             render_ = false;
            //             renderWindowInteractor_->InvokeEvent(CustomEventIds::RENDER_EVENT);
            //         }
            //     }
            // });
            std::cout << "After starting thread for rendering..." << std::endl;
            t.detach();
        }

        void Update() {
            renderer_ = vtkRenderer::New();
            renderWindow_ = vtkRenderWindow::New();
            renderWindowInteractor_ = vtkRenderWindowInteractor::New();
            renderer_->SetBackground(0.8, 0.8, 0.8);
            renderWindow_->AddRenderer(renderer_);
            renderWindow_->SetSize(640, 480);
            renderWindowInteractor_->SetRenderWindow(renderWindow_);
            renderWindowInteractor_->Initialize();
            renderWindowInteractor_->SetInteractorStyle(vtkInteractorStyleTrackballCamera::New());
            vtkNew<vtkRenderCallback_> renderCallback_;
            renderCallback_->mesh_ = mesh_;
            // renderWindowInteractor_->AddObserver(vtkCommand::TimerEvent, timerCallback_);
            // renderWindowInteractor_->CreateRepeatingTimer(100);
            renderWindowInteractor_->AddObserver(CustomEventIds::RENDER_EVENT, renderCallback_);
            // renderWindowInteractor_->Start();
            // std::thread interactorThread(&vtkRenderWindowInteractor::Start, renderWindowInteractor_);
            // interactorThread.detach();
            /*std::cout << "Thread started" << std::endl;
            while (true) {
                // std::cout << "checking render_..." << std::endl;
                if (render_) {
                    render_ = false;
                    renderWindowInteractor_->InvokeEvent(CustomEventIds::RENDER_EVENT);
                    renderWindowInteractor_->ProcessEvents();
                }
                if (render_) {
                    // std::cout << "rendering is set to true" << std::endl;
                    render_ = false;
                    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
                    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
                    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

                    for (auto& v: mesh_->V) {
                        points->InsertNextPoint(v.x, v.y, v.z);
                    }

                    for (auto& f: mesh_->F) {
                        if (f.Vids.size() != 4 || f.Vids.empty()) continue;
                        vtkSmartPointer<vtkQuad> quad = vtkQuad::New();
                        for (int id = 0; id < f.Vids.size(); id++) quad->GetPointIds()->SetId(id, f.Vids[id]);
                        cells->InsertNextCell(quad);
                    }

                    polyData->SetPoints(points);
                    polyData->SetPolys(cells);

                    vtkPolyDataMapper* mapper_ = vtkPolyDataMapper::New();
                    mapper_->SetInputData(polyData);
                    
                    vtkActor* actor_ = vtkActor::New();
                    actor_->SetMapper(mapper_);
                    {
                        vtkProperty* property_ = vtkProperty::New();
                        property_->SetEdgeVisibility(true);
                        actor_->SetProperty(property_);
                    }
            
                    renderer_->RemoveAllViewProps();
                    renderer_->AddActor(actor_);

                    // if (!isRendered_) {
                    //     renderWindowInteractor_->Initialize();
                    //     renderWindowInteractor_->Start();
                    //     isRendered_ = true;
                    // }
                    // renderWindow_->Render();
                    // std::this_thread::sleep_for(std::chrono::milliseconds(500));
                }
            }*/
        }
};

#endif