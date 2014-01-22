#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <algorithm>

#include <boost/thread/thread.hpp>

#include <pangolin/pangolin.h>

#include "bspline.h"

using namespace cascade;
using namespace pangolin;
using namespace std;

uint8_t bg_colour = 255;
size_t LOD = 200;

void DrawBsplineCtrlPt(Image<uint8_t>& img, Bspline<Cubic,double,2> const& bspline)
{

    int const r = 2;

    for(int i = 0; i < bspline.GetNumCtrlPoints(); ++i)
    {

        boost::array<double,2> ctrl_pt = bspline.GetCtrlPoint(i);
        size_t const x = ctrl_pt[0];
        size_t const y = ctrl_pt[1];

        for(int dx = -r; dx <= r; ++dx)
            for(int dy = -r; dy <= r; ++dy)
            {
                if(dx*dx+dy*dy <= (r+0.5)*(r+0.5))
                {
                    int w_x = x+dx; int w_y = y+dy;
                    if(w_x >= 0 && w_x < img.w && w_y >= 0 && w_y < img.h)
                    {
                        img.ptr[(w_y*img.w+w_x)*4] = 255;
                        img.ptr[(w_y*img.w+w_x)*4+1] = 0;
                        img.ptr[(w_y*img.w+w_x)*4+2] = 0;
                        img.ptr[(w_y*img.w+w_x)*4+3] = 255;
                    }
                }
            }
    }

}

void DrawBspline(Image<uint8_t>& img, Bspline<Cubic,double,2> const& bspline, double const LOD)
{

    for(int pt_idx = -1, seg_idx = 0; seg_idx != bspline.GetNumCtrlPoints()+1; ++pt_idx, ++seg_idx)
    {
        for(int d = 0; d < LOD; ++d)
        {
            double t = d/double(LOD-1);

            boost::array<double,2> pt = bspline.Interpolate(pt_idx, t);
            int const x = round(pt[0]);
            int const y = round(pt[1]);

            if(x >= 0 && x < img.w && y >= 0 && y < img.h)
            {
                img.ptr[(y*img.w+x)*4] = 0;
                img.ptr[(y*img.w+x)*4+1] = 0;
                img.ptr[(y*img.w+x)*4+2] = 255;
                img.ptr[(y*img.w+x)*4+3] = 255;
            }
        }
    }

}

int main( int /*argc*/, char* argv[] )
{

    size_t const img_w = 640;
    size_t const img_h = 480;

    const int UI_WIDTH = 180;

    // Create OpenGL window in single line thanks to GLUT
    pangolin::CreateWindowAndBind("Main",UI_WIDTH+img_w, img_h);

    // 3D Mouse handler requires depth testing to be enabled
    glEnable(GL_DEPTH_TEST);

    // Issue specific OpenGl we might need
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Add named OpenGL viewport to window and provide 3D Handler
    View& canvas_view = pangolin::CreateDisplay()
                        .SetBounds(0.0, 1.0, Attach::Pix(UI_WIDTH), 1.0, -640.0f/480.0f)
                        .SetHandler(new Handler2D);

    // Add named Panel and bind to variables beginning 'ui'
    // A Panel is just a View with a default layout and input handling
    pangolin::CreatePanel("ui")
            .SetBounds(0.0, 1.0, 0.0, Attach::Pix(UI_WIDTH));

    Var<bool> open_bspline_button("ui.Open B-spline", false, false);
    Var<bool> closed_bspline_button("ui.Closed B-spline", false, false);

    Var<bool> reset("ui.Reset", false, false);
    Var<bool> save_canvas("ui.Save Canvas", false, false);

    Image<uint8_t> img;
    img.Alloc(img_w, img_h, 4*sizeof(unsigned char)*img_w);
    std::fill(img.ptr, img.ptr+img.h*img.pitch, bg_colour);

    Bspline<Cubic,double,2> bspline;    

    GlTexture* canvas = new GlTexture(img_w,img_h);
    canvas->Upload(img.ptr, GL_RGBA, GL_UNSIGNED_BYTE);

    // Default hooks for exiting (Esc) and fullscreen (tab).
    while( !pangolin::ShouldQuit() )
    {
        // Clear entire screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        if(Pushed(open_bspline_button))
            bspline.SetBsplineType(OPEN);

        if(Pushed(closed_bspline_button))
            bspline.SetBsplineType(CLOSED);

        if(Pushed(reset))
        {
            bspline.ClearCtrlPoints();
            std::fill(img.ptr, img.ptr+img.h*img.pitch, bg_colour);
            canvas->Upload(img.ptr, GL_RGBA, GL_UNSIGNED_BYTE);
        }

        if(((Handler2D*)canvas_view.handler)->IsLeftButtonClicked())
        {
            std::fill(img.ptr, img.ptr+img.h*img.pitch, bg_colour);

            size_t x = ((Handler2D*)canvas_view.handler)->GetLastPos()[0] - canvas_view.vp.l;
            size_t y = ((Handler2D*)canvas_view.handler)->GetLastPos()[1] - canvas_view.vp.b;

            double pt[2] = {x, y};
            bspline.AddCtrlPoint(pt);

            cout << "Add control points: (" << x << "," << y << ")" << endl;
            cout << "Number of control points: " << bspline.GetNumCtrlPts() << endl;

            /* Draw B-Spline */
            DrawBsplineCtrlPt(img, bspline);

            if(bspline.GetNumCtrlPts() >= 4)
                DrawBspline(img, bspline, LOD);

            canvas->Upload(img.ptr, GL_RGBA, GL_UNSIGNED_BYTE);

        }

        // Activate efficiently by object
        canvas_view.ActivateAndScissor();
        canvas->RenderToViewport();

        if( Pushed(save_canvas) )
            canvas_view.SaveOnRender("canvas");

        // Swap frames and Process Events
        pangolin::FinishFrame();
    }

    img.Dealloc();

    delete canvas;

    return 0;
}
