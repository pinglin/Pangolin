#include <pangolin/video/uvc.h>

namespace pangolin
{

UvcVideo::UvcVideo()
    : ctx_(NULL),
      dev_(NULL),
      devh_(NULL)
{
    res_ = uvc_init(&ctx_, NULL);
    if(res_ != UVC_SUCCESS) {
        uvc_perror(res_, "uvc_init");
        throw VideoException("Unable to open UVC Context");
    }

    InitDevice(0, 0, NULL, 640, 480, 30);

    Start();
}

UvcVideo::~UvcVideo()
{
    DeinitDevice();

//    if (ctx_) {
//        // Work out how to kill this properly
//        uvc_exit(ctx_);
//        ctx_ = 0;
//    }   
}

void UvcVideo::InitDevice(int vid, int pid, const char* sn, int width, int height, int fps)
{    

    res_ = uvc_find_device(ctx_, &dev_, vid, pid, sn );
    if(res_ != UVC_SUCCESS) {
        uvc_perror(res_, "uvc_find_device");
        throw VideoException("Unable to find UVC Device");
    }

    res_ = uvc_open(dev_, &devh_);
    if(res_ != UVC_SUCCESS) {
        uvc_perror(res_, "uvc_open");
        throw VideoException("Unable to open UVC Device");
    }

    // Print out all avaliable configuration
    uvc_print_diag(devh_, stderr);

    res_ = uvc_get_stream_ctrl_format_size(devh_, &ctrl_, UVC_FRAME_FORMAT_UNCOMPRESSED, width, height, fps);
    if(res_ != UVC_SUCCESS) {
        uvc_perror(res_, "uvc_get_stream_ctrl_format_size");
        uvc_close(devh_);
        uvc_unref_device(dev_);
        throw VideoException("Unable to make the device mode.");
    }
    
    uvc_print_stream_ctrl(&ctrl_, stderr);

    // assume bgr format
    const VideoPixelFormat pfmt = VideoFormatFromString("RGB24");
    const StreamInfo stream_info(pfmt, width, height, (width*pfmt.bpp)/8, 0);
    streams.push_back(stream_info);

    size_bytes = width*height*3;

}

void UvcVideo::DeinitDevice()
{
    Stop();    
}

void UvcVideo::Start()
{

    res_ = uvc_stream_open_ctrl(devh_, &strmh_, &ctrl_);
    if(res_ != UVC_SUCCESS) {
        uvc_perror(res_, "uvc_stream_open_ctrl");
        uvc_close(devh_);
        uvc_unref_device(dev_);
        throw VideoException("Unable to open a new stream.");
    }

    res_ = uvc_stream_start(strmh_, NULL, NULL, 0);
    if (res_ != UVC_SUCCESS) {
        uvc_perror(res_, "uvc_stream_start");
        uvc_close(devh_);
        uvc_unref_device(dev_);
        throw VideoException("Unable to start streaming.");
    }

}

void UvcVideo::Stop()
{
    if(devh_) {
        uvc_stop_streaming(devh_);
    }
}

size_t UvcVideo::SizeBytes() const
{
    return size_bytes;
}

const std::vector<StreamInfo>& UvcVideo::Streams() const
{
    return streams;
}

bool UvcVideo::GrabNext(unsigned char* image, bool wait)
{

    uvc_frame_t* frame_ = NULL;
    res_ = uvc_stream_get_frame(strmh_, &frame_, 0);
    if (res_ != UVC_SUCCESS) {
        uvc_perror(res_, "uvc_stream_get_frame");
        uvc_close(devh_);
        uvc_unref_device(dev_);
        throw VideoException("Unable to get frame.");
    }

    if(frame_) {

        uvc_frame_t frame_rgb;

        frame_rgb.data = image;
        frame_rgb.data_bytes = streams[0].Width()*streams[0].Height()*3;

        res_ = uvc_any2rgb(frame_, &frame_rgb);
        if (res_ != UVC_SUCCESS) {
            uvc_perror(res_, "uvc_any2rgb");
            uvc_close(devh_);
            uvc_unref_device(dev_);
            throw VideoException("Unable to convert yuv422 to rgb.");
        }

        return true;

    } else {
        std::cerr << "No data..." << std::endl;
        return false;
    }

}

bool UvcVideo::GrabNewest( unsigned char* image, bool wait )
{
    return GrabNext(image, wait);
}

}
