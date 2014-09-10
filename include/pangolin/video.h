/* This file is part of the Pangolin Project.
 * http://github.com/stevenlovegrove/Pangolin
 *
 * Copyright (c) 2011 Steven Lovegrove
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

#ifndef PANGOLIN_VIDEO_H
#define PANGOLIN_VIDEO_H

// Pangolin video supports various cameras and file formats through
// different 3rd party libraries.
//
// Video URI's take the following form:
//  scheme:[param1=value1,param2=value2,...]//device
//
// scheme = file | dc1394 | v4l | openni | convert | mjpeg
//
// files - read one or more streams from image files
// e.g.  "files://~/data/dataset/img_*.jpg"
// e.g.  "files://~/data/dataset/img_[left,right]_*.pgm"
//
// file/files - read PVN file format (pangolin video) or other formats using ffmpeg
//  e.g. "file:[realtime=1]///home/user/video/movie.pvn"
//  e.g. "file:[stream=1]///home/user/video/movie.avi"
//  e.g. "files:///home/user/sequence/foo%03d.jpeg"
//
// dc1394 - capture video through a firewire camera
//  e.g. "dc1394:[fmt=RGB24,size=640x480,fps=30,iso=400,dma=10]//0"
//  e.g. "dc1394:[fmt=FORMAT7_1,size=640x480,pos=2+2,iso=400,dma=10]//0"
//  e.g. "dc1394:[fmt=FORMAT7_3,deinterlace=1]//0"
//
// v4l - capture video from a Video4Linux (USB) camera (normally YUVY422 format)
//           method=mmap|read|userptr
//  e.g. "v4l:///dev/video0"
//  e.g. "v4l[method=mmap]:///dev/video0"
//
// openni2 - capture video / depth from OpenNI2 SDK  (Kinect / Xtrion etc)
//           imgN=grey|rgb|ir|ir8|ir24|depth|reg_depth
//  e.g. "openni2://'
//  e.g. "openni2:[img1=rgb,img2=depth,coloursync=true]//"
//  e.g. "openni2:[img1=depth,close=closerange,holefilter=true]//"
//  e.g. "openni2:[size=320x240,fps=60,img1=ir]//"
//
// openni - capture video / depth from OpenNI 1.0 SDK (Kinect / Xtrion etc)
//           Sensor modes containing '8' will truncate to 8-bits.
//           Sensor modes containing '+' explicitly enable IR illuminator
//           imgN=rgb|ir|ir8|ir+|ir8+|depth|reg_depth
//  e.g. "openni://'
//  e.g. "openni:[img1=rgb,img2=depth]//"
//  e.g. "openni:[size=320x240,fps=60,img1=ir]//"
//
// depthsense - capture video / depth from DepthSense SDK.
//              DepthSenseViewer can be used to alter capture settings.
//  e.g. "depthsense://"
//
// convert - use FFMPEG to convert between video pixel formats
//  e.g. "convert:[fmt=RGB24]//v4l:///dev/video0"
//  e.g. "convert:[fmt=GRAY8]//v4l:///dev/video0"
//
// mjpeg - capture from (possibly networked) motion jpeg stream using FFMPEG
//  e.g. "mjpeg://http://127.0.0.1/?action=stream"
//
// split - split a single stream video into a multi stream video based on Region of Interest / memory specification
//           roiN=X+Y+WxH
//           memN=Offset:WxH:PitchBytes:Format
//  e.g. "split:[roi1=0+0+640x480,roi2=640+0+640x480]//files:///home/user/sequence/foo%03d.jpeg"
//  e.g. "split:[mem1=307200:640x480:1280:GRAY8,roi2=640+0+640x480]//files:///home/user/sequence/foo%03d.jpeg"
//
// test - output test video sequence
//  e.g. "test://"
//  e.g. "test:[size=640x480,fmt=RGB24]//"

#include <pangolin/video_common.h>
#include <pangolin/image.h>

#include <vector>

namespace pangolin
{

class PANGOLIN_EXPORT StreamInfo
{
public:
    inline StreamInfo()
        : fmt(VideoFormatFromString("GRAY8")) {}

    inline StreamInfo(VideoPixelFormat fmt, const Image<unsigned char> img_offset )
        : fmt(fmt), img_offset(img_offset) {}

    inline StreamInfo(VideoPixelFormat fmt, size_t w, size_t h, size_t pitch, unsigned char* offset)
        : fmt(fmt), img_offset(w,h,pitch,offset) {}
    
    //! Format representing how image is layed out in memory
    inline VideoPixelFormat PixFormat() const { return fmt; }
    
    //! Image width in pixels
    inline size_t Width() const { return img_offset.w; }
    
    //! Image height in pixels
    inline size_t Height() const { return img_offset.h; }
    
    //! Pitch: Number of bytes between one image row and the next
    inline size_t Pitch() const { return img_offset.pitch; }
    
    //! Number of contiguous bytes in memory that the image occupies
    inline size_t SizeBytes() const { return img_offset.h * img_offset.pitch; }
    
    //! Offset in bytes relative to start of frame buffer
    inline unsigned char* Offset() const { return img_offset.ptr; }
    
    //! Return Image wrapper around raw base pointer
    Image<unsigned char> StreamImage(unsigned char* base_ptr) const {
        Image<unsigned char> img = img_offset;
        img.ptr += (size_t)base_ptr;
        return img;
    }
    
protected:
    VideoPixelFormat fmt;        
    Image<unsigned char> img_offset;
};

//! Interface to video capture sources
struct PANGOLIN_EXPORT VideoInterface
{
    virtual ~VideoInterface() {}

    //! Required buffer size to store all frames
    virtual size_t SizeBytes() const = 0;
    
    //! Get format and dimensions of all video streams
    virtual const std::vector<StreamInfo>& Streams() const = 0;
    
    //! Start Video device
    virtual void Start() = 0;
    
    //! Stop Video device
    virtual void Stop() = 0;
    
    //! Copy the next frame from the camera to image.
    //! Optionally wait for a frame if one isn't ready
    //! Returns true iff image was copied
    virtual bool GrabNext( unsigned char* image, bool wait = true ) = 0;
    
    //! Copy the newest frame from the camera to image
    //! discarding all older frames.
    //! Optionally wait for a frame if one isn't ready
    //! Returns true iff image was copied
    virtual bool GrabNewest( unsigned char* image, bool wait = true ) = 0;
};

//! Generic wrapper class for different video sources
struct PANGOLIN_EXPORT VideoInput : public VideoInterface
{
    VideoInput();
    VideoInput(const std::string& uri);
    ~VideoInput();
    
    void Open(const std::string& uri);
    void Reset();
    
    size_t SizeBytes() const;
    const std::vector<StreamInfo>& Streams() const;
    
    // Return details of first stream
    unsigned Width() const;
    unsigned Height() const;
    VideoPixelFormat PixFormat() const;
    const Uri& VideoUri() const;
    
    void Start();
    void Stop();
    bool GrabNext( unsigned char* image, bool wait = true );
    bool GrabNewest( unsigned char* image, bool wait = true );

    // Return pointer to inner video class as VideoType
    template<typename VideoType>
    VideoType* Cast() {
        return dynamic_cast<VideoType*>(video);
    }

    // experimental - not stable
    bool Grab( unsigned char* buffer, std::vector<Image<unsigned char> >& images, bool wait = true, bool newest = false);
    
protected:
    Uri uri;
    VideoInterface* video;
};

//! Open Video Interface from string specification (as described in this files header)
PANGOLIN_EXPORT
VideoInterface* OpenVideo(const std::string& uri);

//! Open Video Interface from Uri specification
PANGOLIN_EXPORT
VideoInterface* OpenVideo(const Uri& uri);

}

#endif // PANGOLIN_VIDEO_H
