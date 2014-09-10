#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <algorithm>

#include <pangolin/pangolin.h>

#include "bspline.h"

using namespace pangolin;
using namespace std;

////////////////////////////////////////////////////////////////////////////
//  Constant colours
////////////////////////////////////////////////////////////////////////////
float colour_roi[3] = {1.0, 0.0, 0.0};
float colour_knot_pt[3] = {1.0, 0.0, 0.0};
float colour_ctrl_pt[3] = {0.0, 1.0, 1.0};
float colour_spline[3] = {1.0, 1.0, 1.0};

////////////////////////////////////////////////////////////////////////////
//  Pangolin UI
////////////////////////////////////////////////////////////////////////////
Var<int>* LOD; // Level of details, i.e., sampling numbers for B-spline
Var<bool>* check_knot_mode;
Var<bool>* button_open_bspline;
Var<bool>* button_closed_bspline;
Var<bool>* button_reset;
Var<bool>* button_save_canvas;

////////////////////////////////////////////////////////////////////////////
//  Global functions
////////////////////////////////////////////////////////////////////////////
void draw_spline(Bspline<float,2> const& bspline)
{

    for(int pt_idx = -1, seg_idx = 0; seg_idx < bspline.get_num_ctrl_pts()+1; ++pt_idx, ++seg_idx)
    {
        for(int d = 0; d < *LOD; ++d)
        {
            array<float,2> pt1 = bspline.cubic_intplt(pt_idx, d/float(*LOD));
            array<float,2> pt2 = bspline.cubic_intplt(pt_idx, (d+1)/float(*LOD));

            glColor3fv(colour_spline);
            glBegin(GL_LINES);
            glVertex2f(pt1[0], pt1[1]);
            glVertex2f(pt2[0], pt2[1]);
            glEnd();
        }
    }
}

void draw_knot_pts(Bspline<float,2> const& bspline)
{
    for(int k = 0; k < bspline.get_num_knot_pts(); ++k)
    {
        array<float,2> const& pt = bspline.get_knot_pt(k);

        glColor3fv(colour_knot_pt);
        glPointSize(5);
        glBegin(GL_POINTS);
        glVertex2f(pt[0], pt[1]);
        glEnd();
    }
}

void draw_ctrl_pts(Bspline<float,2> const& bspline)
{
    for(int k = 0; k < bspline.get_num_ctrl_pts(); ++k)
    {
        array<float,2> const& pt = bspline.get_ctrl_pt(k);

        glColor3fv(colour_ctrl_pt);
        glPointSize(5);
        glBegin(GL_POINTS);
        glVertex2f(pt[0], pt[1]);
        glEnd();
    }
}

////////////////////////////////////////////////////////////////////////////
//  Main function
////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

    size_t const img_w = 640;
    size_t const img_h = 480;

    const int ui_width = 180;

    // Create OpenGL window in single line thanks to GLUT
    pangolin::CreateWindowAndBind("B-spline Canvas", ui_width+img_w, img_h);

    // 3D Mouse handler requires depth testing to be enabled
    glEnable(GL_DEPTH_TEST);

    // Issue specific OpenGl we might need
    glEnable(GL_POINT_SMOOTH);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0.0, 0.0, 0.0, 1.0);

    // Add named OpenGL viewport to window and provide 3D Handler
    View& canvas_view = pangolin::CreateDisplay()
                        .SetBounds(0.0, 1.0, Attach::Pix(ui_width), 1.0, -640.0f/480.0f)
                        .SetHandler(new Handler2D);

    // Add named Panel and bind to variables beginning 'ui'
    // A Panel is just a View with a default layout and input handling
    pangolin::CreatePanel("ui")
            .SetBounds(0.0, 1.0, 0.0, Attach::Pix(ui_width));

    LOD = new Var<int>("ui.Level of Details", 10, 1, 50);
    check_knot_mode = new Var<bool>("ui.Assign Knots", true, true, false);
    button_open_bspline = new Var<bool>("ui.Open B-spline", false, false);
    button_closed_bspline = new Var<bool>("ui.Closed B-spline", false, false);
    button_reset = new Var<bool>("ui.Reset", false, false);
    button_save_canvas = new Var<bool>("ui.Save Canvas", false, false);

    Bspline<float,2> bspline;

    while(!pangolin::ShouldQuit())
    {
        // Clear entire screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if(Pushed(*button_open_bspline))
            bspline.set_bspline_type(Bspline<float,2>::OPEN);

        if(Pushed(*button_closed_bspline))
            bspline.set_bspline_type(Bspline<float,2>::CLOSED);

        if(Pushed(*button_reset))
            bspline.clear();

        if(((Handler2D*)canvas_view.handler)->IsLeftButtonClicked())
        {
            size_t x = ((Handler2D*)canvas_view.handler)->GetLastPos()[0] - canvas_view.vp.l;
            size_t y = ((Handler2D*)canvas_view.handler)->GetLastPos()[1] - canvas_view.vp.b;

            array<float,2> pt = {(float)x, (float)y};
            if(*check_knot_mode)
            {
                bspline.add_knot_pt(pt);

                cout << "Add node points: (" << x << "," << y << ")" << endl;
                cout << "Number of node points: " << bspline.get_num_knot_pts() << endl;

            }
            else
            {
                bspline.add_ctrl_pt(pt);

                cout << "Add control points: (" << x << "," << y << ")" << endl;
                cout << "Number of control points: " << bspline.get_num_ctrl_pts() << endl;
            }
        }

        // Activate efficiently by object
        canvas_view.ActivatePixelOrthographic();

        if( Pushed(*button_save_canvas) )
            canvas_view.SaveOnRender("canvas");

        draw_knot_pts(bspline);
        draw_ctrl_pts(bspline);
        draw_spline(bspline);

        // Swap frames and Process Events
        pangolin::FinishFrame();

    }

    return 0;
}
