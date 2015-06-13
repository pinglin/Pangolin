#include <pangolin/pangolin.h>
#include <pangolin/video/video_record_repeat.h>
#include <pangolin/gl/gltexturecache.h>

struct GlFormat
{
    GlFormat() {}

    GlFormat(const pangolin::VideoPixelFormat& fmt)
    {
        switch( fmt.channels) {
        case 1: glformat = GL_LUMINANCE; break;
        case 3: glformat = (fmt.format == "BGR24") ? GL_BGR : GL_RGB; break;
        case 4: glformat = (fmt.format == "BGRA24") ? GL_BGRA : GL_RGBA; break;
        default: throw std::runtime_error("Unable to display video format");
        }

        switch (fmt.channel_bits[0]) {
        case 8: gltype = GL_UNSIGNED_BYTE; break;
        case 16: gltype = GL_UNSIGNED_SHORT; break;
        case 32: gltype = GL_FLOAT; break;
        default: throw std::runtime_error("Unknown channel format");
        }
    }

    GLint glformat;
    GLenum gltype;
};

void RenderToViewport(
    pangolin::Image<unsigned char>& image,
    const GlFormat& fmt, bool flipx=false, bool flipy=false, bool linear_sampling = true
) {
    pangolin::GlTexture& tex = pangolin::TextureCache::I().GlTex(image.w, image.h, GL_RGBA, GL_RGBA, fmt.gltype);
    tex.Bind();
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, linear_sampling ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, linear_sampling ? GL_LINEAR : GL_NEAREST);
    tex.Upload(image.ptr,0,0, image.w, image.h, fmt.glformat, fmt.gltype);
    tex.RenderToViewport(pangolin::Viewport(0,0,image.w, image.h), flipx, flipy);
}

void VideoViewer(const std::string& input_uri, const std::string& output_uri)
{
    // Open Video by URI
    pangolin::VideoRecordRepeat video(input_uri, output_uri);
    int total_frames = std::numeric_limits<int>::max();

    if(video.Streams().size() == 0) {
        pango_print_error("No video streams from device.\n");
        return;
    }

    // Check if video supports VideoPlaybackInterface
    pangolin::VideoPlaybackInterface* video_playback = video.Cast<pangolin::VideoPlaybackInterface>();
    if( video_playback ) {
        total_frames = video_playback->GetTotalFrames();
        std::cout << "Video length: " << total_frames << " frames" << std::endl;
    }

    std::vector<unsigned char> buffer;
    buffer.resize(video.SizeBytes()+1);

    // Create OpenGL window - guess sensible dimensions
    pangolin::CreateWindowAndBind( "VideoViewer",
        video.Width() * video.Streams().size(), video.Height()
    );

    // Setup resizable views for video streams
    std::vector<GlFormat> glfmt;
    pangolin::DisplayBase().SetLayout(pangolin::LayoutEqual);
    for(unsigned int d=0; d < video.Streams().size(); ++d) {
        pangolin::View& view = pangolin::CreateDisplay().SetAspect(video.Streams()[d].Aspect());
        pangolin::DisplayBase().AddDisplay(view);
        glfmt.push_back(GlFormat(video.Streams()[d].PixFormat()));
    }

    const int FRAME_SKIP = 30;
    int frame = 0;
    pangolin::Var<int>  max_frame("max_frame", total_frames );
    pangolin::Var<bool> linear_sampling("linear_sampling", true );
    pangolin::Var<float> int16_scale("int16.scale", 20.0 );
    pangolin::Var<float> int16_bias("int16.bias", 0.0 );

#ifdef CALLEE_HAS_CPP11
    // Show/hide streams
    for(size_t v=0; v < pangolin::DisplayBase().NumChildren() && v < 9; v++) {
        pangolin::RegisterKeyPressCallback('1'+v, [v](){
            pangolin::DisplayBase()[v].ToggleShow();
        } );
    }

    pangolin::RegisterKeyPressCallback('r', [&](){
        if(!video.IsRecording()) {
            video.Record();
            pango_print_info("Started Recording.\n");
        }else{
            video.Stop();
            pango_print_info("Finished recording.\n");
        }
        fflush(stdout);
    });
    pangolin::RegisterKeyPressCallback('p', [&](){
        video.Play();
        max_frame = std::numeric_limits<int>::max();
        pango_print_info("Playing from file log.\n");
        fflush(stdout);
    });
    pangolin::RegisterKeyPressCallback('s', [&](){
        video.Source();
        max_frame = std::numeric_limits<int>::max();
        pango_print_info("Playing from source input.\n");
        fflush(stdout);
    });
    pangolin::RegisterKeyPressCallback(' ', [&](){
        max_frame = (frame < max_frame) ? frame : std::numeric_limits<int>::max();
    });
    pangolin::RegisterKeyPressCallback(pangolin::PANGO_SPECIAL + pangolin::PANGO_KEY_LEFT, [&](){
        if(video_playback) {
            const int frame = std::min(video_playback->GetCurrentFrameId()-FRAME_SKIP, video_playback->GetTotalFrames()-1);
            video_playback->Seek(frame);
        }else{
            // We can't go backwards
        }
    });
    pangolin::RegisterKeyPressCallback(pangolin::PANGO_SPECIAL + pangolin::PANGO_KEY_RIGHT, [&](){
        if(video_playback) {
            const int frame = std::max(video_playback->GetCurrentFrameId()+FRAME_SKIP, 0);
            video_playback->Seek(frame);
        }else{
            // Pause at this frame
            max_frame = frame+1;
        }
    });
    pangolin::RegisterKeyPressCallback('l', [&](){ linear_sampling = true; });
    pangolin::RegisterKeyPressCallback('n', [&](){ linear_sampling = false; });
#endif

    std::vector<pangolin::Image<unsigned char> > images;

    // Stream and display video
    while(!pangolin::ShouldQuit())
    {
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 1.0f, 1.0f);

        if (frame == 0 || frame < max_frame) {
            images.clear();
            if (video.Grab(&buffer[0], images) ){
                ++frame;
            }
        }

        for(unsigned int i=0; i<images.size(); ++i)
        {
            if(pangolin::DisplayBase()[i].IsShown()) {
                pangolin::DisplayBase()[i].Activate();
                if(glfmt[i].gltype == GL_UNSIGNED_SHORT) {
                    pangolin::GlSlUtilities::Scale(int16_scale, int16_bias);
                    RenderToViewport(images[i], glfmt[i], false, true, linear_sampling);
                    pangolin::GlSlUtilities::UseNone();
                }else{
                    RenderToViewport(images[i], glfmt[i], false, true, linear_sampling);
                }
            }
        }

        pangolin::FinishFrame();
    }
}


int main( int argc, char* argv[] )
{
    const std::string dflt_output_uri = "pango://video.pango";

    if( argc > 1 ) {
        const std::string input_uri = std::string(argv[1]);
        const std::string output_uri = (argc > 2) ? std::string(argv[2]) : dflt_output_uri;
        try{
            VideoViewer(input_uri, output_uri);
        } catch (pangolin::VideoException e) {
            std::cout << e.what() << std::endl;
        }
    }else{
        const std::string input_uris[] = {
            "dc1394:[fps=30,dma=10,size=640x480,iso=400]//0",
            "convert:[fmt=RGB24]//v4l:///dev/video0",
            "convert:[fmt=RGB24]//v4l:///dev/video1",
            "openni:[img1=rgb]//",
            "test:[size=160x120,n=1,fmt=RGB24]//"
            ""
        };

        std::cout << "Usage  : VideoViewer [video-uri]" << std::endl << std::endl;
        std::cout << "Where video-uri describes a stream or file resource, e.g." << std::endl;
        std::cout << "\tfile:[realtime=1]///home/user/video/movie.pvn" << std::endl;
        std::cout << "\tfile:///home/user/video/movie.avi" << std::endl;
        std::cout << "\tfiles:///home/user/seqiemce/foo%03d.jpeg" << std::endl;
        std::cout << "\tdc1394:[fmt=RGB24,size=640x480,fps=30,iso=400,dma=10]//0" << std::endl;
        std::cout << "\tdc1394:[fmt=FORMAT7_1,size=640x480,pos=2+2,iso=400,dma=10]//0" << std::endl;
        std::cout << "\tv4l:///dev/video0" << std::endl;
        std::cout << "\tconvert:[fmt=RGB24]//v4l:///dev/video0" << std::endl;
        std::cout << "\tmjpeg://http://127.0.0.1/?action=stream" << std::endl;
        std::cout << "\topenni:[img1=rgb]//" << std::endl;
        std::cout << std::endl;

        // Try to open some video device
        for(int i=0; !input_uris[i].empty(); ++i )
        {
            try{
                pango_print_info("Trying: %s\n", input_uris[i].c_str());
                VideoViewer(input_uris[i], dflt_output_uri);
                return 0;
            }catch(pangolin::VideoException) { }
        }
    }

    return 0;
}
