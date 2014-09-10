/* This file is part of the Pangolin Project.
 * http://github.com/stevenlovegrove/Pangolin
 *
 * Copyright (c) 2014 Richard Newcombe
 *               2014 Steven Lovegrove
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

#include <pangolin/video/openni2.h>

#include <PS1080.h>
#include <OniVersion.h>

namespace pangolin
{

VideoPixelFormat VideoFormatFromOpenNI2(openni::PixelFormat fmt)
{
    std::string pvfmt;

    switch (fmt) {
    case openni::PIXEL_FORMAT_DEPTH_1_MM:   pvfmt = "GRAY16LE"; break;
    case openni::PIXEL_FORMAT_DEPTH_100_UM: pvfmt = "GRAY16LE"; break;
    case openni::PIXEL_FORMAT_SHIFT_9_2:    pvfmt = "GRAY16LE"; break; // ?
    case openni::PIXEL_FORMAT_SHIFT_9_3:    pvfmt = "GRAY16LE"; break; // ?
    case openni::PIXEL_FORMAT_RGB888:       pvfmt = "RGB24"; break;
    case openni::PIXEL_FORMAT_GRAY8:        pvfmt = "GRAY8"; break;
    case openni::PIXEL_FORMAT_GRAY16:       pvfmt = "GRAY16LE"; break;
    case openni::PIXEL_FORMAT_YUV422:       pvfmt = "YUYV422"; break;
#if ONI_VERSION_MAJOR >= 2 && ONI_VERSION_MINOR >= 2
    case openni::PIXEL_FORMAT_YUYV:         pvfmt = "Y400A"; break;
#endif
    default:
        throw VideoException("Unknown OpenNI pixel format");
        break;
    }

    return VideoFormatFromString(pvfmt);
}

void OpenNiVideo2::PrintOpenNI2Modes(openni::SensorType sensorType)
{
    // Query supported modes for device
    const openni::Array<openni::VideoMode>& modes =
            device.getSensorInfo(sensorType)->getSupportedVideoModes();

    switch (sensorType) {
    case openni::SENSOR_COLOR: pango_print_info("OpenNI Colour Modes:\n"); break;
    case openni::SENSOR_DEPTH: pango_print_info("OpenNI Depth Modes:\n"); break;
    case openni::SENSOR_IR:    pango_print_info("OpenNI IR Modes:\n"); break;
    }

    for(int i = 0; i < modes.getSize(); i++) {
        std::string sfmt = "PangolinUnknown";
        try{
            sfmt = VideoFormatFromOpenNI2(modes[i].getPixelFormat()).format;
        }catch(VideoException){}
        pango_print_info( "  %dx%d, %d fps, %s\n",
            modes[i].getResolutionX(), modes[i].getResolutionY(),
            modes[i].getFps(), sfmt.c_str()
        );
    }
}

openni::VideoMode OpenNiVideo2::FindOpenNI2Mode(
    openni::SensorType sensorType,
    int width, int height,
    int fps, openni::PixelFormat fmt
) {
    // Query supported modes for device
    const openni::Array<openni::VideoMode>& modes =
            device.getSensorInfo(sensorType)->getSupportedVideoModes();

    // Select last listed mode which matches parameters
    int best_mode = -1;
    for(int i = 0; i < modes.getSize(); i++) {
        if( (!width || modes[i].getResolutionX() == width) &&
            (!height || modes[i].getResolutionY() == height) &&
            (!fps || modes[i].getFps() == fps) &&
            (!fmt || modes[i].getPixelFormat() == fmt)
        ) {
            best_mode = i;
        }
    }

    if(best_mode >= 0) {
        return modes[best_mode];
    }

    throw pangolin::VideoException("Video mode not supported");
}

OpenNiVideo2::OpenNiVideo2(OpenNiSensorType s1, OpenNiSensorType s2, ImageDim dim, int fps)
{
    sensor_type[0] = s1;
    sensor_type[1] = s2;

    use_depth = false;
    use_ir = false;
    use_rgb = false;
    depth_to_color = false;
    use_ir_and_rgb = false;

    const char* deviceURI = openni::ANY_DEVICE;
    fromFile = (deviceURI!=NULL);

    openni::Status rc = openni::OpenNI::initialize();
    if (rc != openni::STATUS_OK) {
        throw VideoException( "Unable to initialise OpenNI library", openni::OpenNI::getExtendedError() );
    }

    rc = device.open(deviceURI);
    if (rc != openni::STATUS_OK) {
        throw VideoException("Failed to open device", openni::OpenNI::getExtendedError());
    }

//    PrintOpenNI2Modes(openni::SENSOR_COLOR);
//    PrintOpenNI2Modes(openni::SENSOR_DEPTH);
//    PrintOpenNI2Modes(openni::SENSOR_IR);

    for(int i=0; i<2; ++i) {
        openni::SensorType sensortype;
        openni::PixelFormat pixelfmt;

        switch( sensor_type[i] ) {
        case OpenNiDepthRegistered:
            depth_to_color = true;
        case OpenNiDepth:
            sensortype = openni::SENSOR_DEPTH;
            pixelfmt = openni::PIXEL_FORMAT_DEPTH_1_MM;
            use_depth = true;
            break;
        case OpenNiIrProj:
        case OpenNiIr:
            sensortype = openni::SENSOR_IR;
            pixelfmt = openni::PIXEL_FORMAT_GRAY16;
            use_ir = true;
            break;
        case OpenNiIr24bit:
            sensortype = openni::SENSOR_IR;
            pixelfmt = openni::PIXEL_FORMAT_RGB888;
            use_ir = true;
            break;
        case OpenNiIr8bitProj:
        case OpenNiIr8bit:
            sensortype = openni::SENSOR_IR;
            pixelfmt = openni::PIXEL_FORMAT_GRAY8;
            use_ir = true;
            break;
        case OpenNiRgb:
            sensortype = openni::SENSOR_COLOR;
            pixelfmt = openni::PIXEL_FORMAT_RGB888;
            use_rgb = true;
            break;
        case OpenNiGrey:
            sensortype = openni::SENSOR_COLOR;
            pixelfmt = openni::PIXEL_FORMAT_GRAY8;
            use_rgb = true;
            break;
        case OpenNiUnassigned:
        default:
            continue;
        }

        openni::VideoMode onivmode;
        try {
            onivmode = FindOpenNI2Mode(sensortype, dim.x, dim.y, fps, pixelfmt);
        }catch(VideoException e) {
            pango_print_error("Unable to find compatible OpenNI Video Mode. Please choose from:\n");
            PrintOpenNI2Modes(sensortype);
            fflush(stdout);
            throw e;
        }

        if(fromFile) {
            // do something with mode?
        }

        rc = video_stream[i].create(device, sensortype);
        if(rc != openni::STATUS_OK)
            throw VideoException("Couldn't create sensor", openni::OpenNI::getExtendedError());

        rc = video_stream[i].setVideoMode(onivmode);
        if(rc != openni::STATUS_OK)
            throw VideoException("Couldn't set OpenNI VideoMode", openni::OpenNI::getExtendedError());

        video_stream[i].setMirroringEnabled(false);

        if(sensortype == openni::SENSOR_COLOR) {
            video_stream[i].getCameraSettings()->setAutoExposureEnabled(true);
            video_stream[i].getCameraSettings()->setAutoWhiteBalanceEnabled(true);
        }

        const VideoPixelFormat fmt = VideoFormatFromOpenNI2(pixelfmt);
        const StreamInfo stream(
            fmt, onivmode.getResolutionX(), onivmode.getResolutionY(),
            (onivmode.getResolutionX() * fmt.bpp) / 8,
            (unsigned char*)0 + sizeBytes
        );

        sizeBytes += stream.SizeBytes();
        streams.push_back(stream);
    }
    use_ir_and_rgb = use_rgb && use_ir;

    if(fromFile) {
        // Go as fast as we can.
        device.getPlaybackControl()->setSpeed(-1);
    }

    if(depth_to_color) {
        device.setImageRegistrationMode(openni::IMAGE_REGISTRATION_DEPTH_TO_COLOR);
    }else{
        device.setImageRegistrationMode(openni::IMAGE_REGISTRATION_OFF);
    }

    Start();
}

void OpenNiVideo2::SetDepthCloseRange(bool enable)
{
    // Set this property on all devices. It doesn't matter if it fails.
    for(int i=0; i<2; ++i) {
        video_stream[i].setProperty(XN_STREAM_PROPERTY_CLOSE_RANGE, enable);
    }
}

void OpenNiVideo2::SetDepthHoleFilter(bool enable)
{
    // Set this property on all devices. It doesn't matter if it fails.
    for(int i=0; i<2; ++i) {
        video_stream[i].setProperty(XN_STREAM_PROPERTY_HOLE_FILTER, enable);
        video_stream[i].setProperty(XN_STREAM_PROPERTY_GAIN,50);
    }
}

void OpenNiVideo2::SetDepthColorSyncEnabled(bool enable)
{
    device.setDepthColorSyncEnabled(enable);
}

void OpenNiVideo2::SetRegisterDepthToImage(bool enable)
{
    if(enable) {
        device.setImageRegistrationMode(openni::IMAGE_REGISTRATION_DEPTH_TO_COLOR);
    }else{
        device.setImageRegistrationMode(openni::IMAGE_REGISTRATION_OFF);
    }
}

OpenNiVideo2::~OpenNiVideo2()
{
    Stop();

    for(int i=0; i<2; ++i) {
        if( video_stream[i].isValid()) {
            video_stream[i].destroy();
        }
    }

    openni::OpenNI::shutdown();
}

size_t OpenNiVideo2::SizeBytes() const
{
    return sizeBytes;
}

const std::vector<StreamInfo>& OpenNiVideo2::Streams() const
{
    return streams;
}

void OpenNiVideo2::Start()
{
    for(int i=0; i<2; ++i) {
        video_stream[i].start();
    }
}

void OpenNiVideo2::Stop()
{
    for(int i=0; i<2; ++i) {
        video_stream[i].stop();
    }
}

bool OpenNiVideo2::GrabNext( unsigned char* image, bool wait )
{
    unsigned char* out_img = image;

    openni::Status rc;

    for(int i=0; i<2; ++i) {
        if(!video_stream[i].isValid()) {
            rc = openni::STATUS_NO_DEVICE;
            continue;
        }

        if(use_ir_and_rgb) video_stream[i].start();

        rc = video_stream[i].readFrame(&video_frame[i]);
        if(rc != openni::STATUS_OK) {
            pango_print_error("Error reading frame:\n%s", openni::OpenNI::getExtendedError() );
        }

        const bool toGreyscale = false;
        if(toGreyscale) {
            const int w = streams[i].Width();
            const int h = streams[i].Height();

            openni::RGB888Pixel* pColour = (openni::RGB888Pixel*)video_frame[i].getData();
            for(int i = 0 ; i  < w*h;i++){
                openni::RGB888Pixel rgb = pColour[i];
                int grey = ((int)(rgb.r&0xFF) +  (int)(rgb.g&0xFF) + (int)(rgb.b&0xFF))/3;
                grey = std::min(255,std::max(0,grey));
                out_img[i] = grey;
            }
        }else{
            memcpy(out_img, video_frame[i].getData(), streams[i].SizeBytes());
        }

        if(use_ir_and_rgb) video_stream[i].stop();

        out_img += streams[i].SizeBytes();
    }

    return rc == openni::STATUS_OK;
}

bool OpenNiVideo2::GrabNewest( unsigned char* image, bool wait )
{
    return GrabNext(image,wait);
}

}
