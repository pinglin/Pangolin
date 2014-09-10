/* This file is part of the Pangolin Project.
 * http://github.com/stevenlovegrove/Pangolin
 *
 * Copyright (c) 2014 Steven Lovegrove
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef PLOTTER_H
#define PLOTTER_H

#include <pangolin/gl.h>
#include <pangolin/view.h>
#include <pangolin/handler.h>
#include <pangolin/glsl.h>
#include <pangolin/datalog.h>
#include <pangolin/colour.h>
#include <pangolin/glfont.h>
#include <pangolin/range.h>

#include <set>

namespace pangolin
{

struct Marker
{
    enum Direction
    {
        Horizontal,
        Vertical
    };

    enum Equality
    {
        LessThan = -1,
        Equal = 0,
        GreaterThan = 1
    };

    Marker(Direction d, float value, Equality leg = Equal, Colour c = Colour() )
        : direction(d), value(value), leg(leg), colour(c)
    {
    }

    Direction direction;
    float value;
    Equality leg;
    Colour colour;
};

class PANGOLIN_EXPORT Plotter : public View, Handler
{
public:
    Plotter(
        DataLog* log,
        float left=0, float right=600, float bottom=-1, float top=1,
        float tickx=30, float ticky=0.5,
        Plotter* linked_plotter_x = 0,
        Plotter* linked_plotter_y = 0
    );

    virtual ~Plotter();

    void Render();

    XYRange& GetSelection();

    XYRange& GetDefaultView();
    void SetDefaultView(const XYRange& range);

    XYRange& GetView();
    void SetView(const XYRange& range);
    void SetViewSmooth(const XYRange& range);

    void ScrollView(float x, float y);
    void ScrollViewSmooth(float x, float y);

    void ScaleView(float x, float y, float cx, float cy);
    void ScaleViewSmooth(float x, float y, float cx, float cy);

    void ResetView();

    void SetTicks(float tickx, float ticky);

    void Track(const std::string& x="$i", const std::string& y = "");
    void ToggleTracking();

    void Trigger(const std::string& x="$0", int edge = -1, float value = 0.0f);
    void ToggleTrigger();

    void SetBackgroundColour(const Colour& col);
    void SetAxisColour(const Colour& col);
    void SetTickColour(const Colour& col);

    void ScreenToPlot(int xpix, int ypix, float &xplot, float &yplot);
    void Keyboard(View&, unsigned char key, int x, int y, bool pressed);
    void Mouse(View&, MouseButton button, int x, int y, bool pressed, int mouse_state);
    void MouseMotion(View&, int x, int y, int mouse_state);
    void PassiveMouseMotion(View&, int x, int y, int button_state);
    void Special(View&, InputSpecial inType, float x, float y, float p1, float p2, float p3, float p4, int button_state);

    Marker& AddMarker(Marker::Direction d, float value, Marker::Equality leg = Marker::Equal, Colour c = Colour() );
    void ClearMarkers();

protected:
    struct PANGOLIN_EXPORT Tick
    {
        float val;
        float factor;
        std::string symbol;
    };

    struct PANGOLIN_EXPORT PlotAttrib
    {
        PlotAttrib(std::string name, int plot_id, int location = -1)
            : name(name), plot_id(plot_id), location(location) { }

        std::string name;
        int plot_id;
        int location;
    };

    struct PANGOLIN_EXPORT PlotSeries
    {
        PlotSeries();
        void CreatePlot(const std::string& x, const std::string& y, Colour c, std::string title);

        GlSlProgram prog;
        GlText title;
        bool contains_id;
        std::vector<PlotAttrib> attribs;
        GLenum drawing_mode;
        Colour colour;
        bool used;
    };

    struct PANGOLIN_EXPORT PlotImplicit
    {
        // Assign to gl_FragColor
        void CreatePlot(const std::string& code);

        // Expression uses x,y and assignes colours [0,1] to r,g,b,a
        void CreateColouredPlot(const std::string& code);

        // Expression uses x,y and evaluates to true/false;
        void CreateInequality(const std::string& ie, Colour c);

        // Expression uses x,y and evaluates to a number
        void CreateDistancePlot(const std::string& dist);

        GlSlProgram prog;
    };

    void FixSelection();
    void UpdateView();
    Tick FindTickFactor(float tick);

    DataLog* log;

    ColourWheel colour_wheel;
    Colour colour_bg;
    Colour colour_tk;
    Colour colour_ax;

    GlSlProgram prog_lines;
    GlSlProgram prog_text;

    std::vector<PlotSeries> plotseries;
    std::vector<Marker> plotmarkers;
    std::vector<PlotImplicit> plotimplicits;

    Tick tick[2];
    XYRange rview_default;
    XYRange rview;
    XYRange target;
    XYRange selection;

    void ComputeTrackValue( float track_val[2] );

    bool track;
    std::string track_x;
    std::string track_y;
    float last_track_val[2];

    // -1: falling, -0:disable, 1: rising
    int trigger_edge;
    float trigger_value;
    std::string trigger;

    float hover[2];
    int last_mouse_pos[2];

    Plotter* linked_plotter_x;
    Plotter* linked_plotter_y;
};

} // namespace pangolin

#endif // PLOTTER_H
