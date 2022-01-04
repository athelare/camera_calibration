//
// Created by Jiyu_ on 2021/8/2.
//

#ifndef OPENCV_HELLO_DVPFRAMECAPTURE_H
#define OPENCV_HELLO_DVPFRAMECAPTURE_H
#include "include/DVPCamera.h"
#include <opencv2/opencv.hpp>

class DVPFrameCapture{
    char *friendlyName{};
    dvpStatus status;
    dvpHandle h{};
    dvpFrame frame{};
    void *p{};
public:
    ~DVPFrameCapture(){
        dvpStop(h);
        dvpClose(h);
    }
    explicit DVPFrameCapture(int cameraIndex = 0){
        dvpCameraInfo info[8];
        dvpUint32 count = 0;
        /* 枚举设备 */
        dvpRefresh(&count);
        if (count > 8)
            count = 8;
        for (int i = 0; i < count; i++)
        {
            if (dvpEnum(i, &info[i]) == DVP_STATUS_OK)
            {
                printf("[%d]-Camera FriendlyName : %s\r\n", i, info[i].FriendlyName);
            }
        }
        friendlyName = info[cameraIndex].FriendlyName;
        dvpRegion region;
        double exp;
        float gain;
        bool trigMode = false;

        status = dvpOpenByName(friendlyName, OPEN_NORMAL, &h);
        if(status != DVP_STATUS_OK){
            printf("dvpOpenByName failed with err:%d\r\n", status);
            return;
        }
        status = dvpSetTriggerState(h, false);
        if (status == DVP_STATUS_OK)
        {
            trigMode = false;
        }
        else
        {
            printf("dvpSetTriggerState failed with err:%d\r\n", status);
            return;
        }

        /* 打印ROI信息 */
        status = dvpGetRoi(h, &region);
        if (status != DVP_STATUS_OK)
        {
            printf("dvpGetRoi failed with err:%d\r\n", status);
            return;
        }
        printf("%s, region: x:%d, y:%d, w:%d, h:%d\r\n", friendlyName, region.X, region.Y, region.W, region.H);

        /* 打印曝光增益信息 */
        status = dvpGetExposure(h, &exp);
        if (status != DVP_STATUS_OK)
        {
            printf("dvpGetExposure failed with err:%d\r\n", status);
            return;
        }

        status = dvpGetAnalogGain(h, &gain);
        if (status != DVP_STATUS_OK)
        {
            printf("dvpGetAnalogGain failed with err:%d\r\n", status);
            return;
        }

        printf("%s, exposure: %lf, gain: %f\r\n", friendlyName, exp, gain);


        /* 开始视频流 */
        status = dvpStart(h);
        if (status != DVP_STATUS_OK)
        {
            return;
        }
    }

    bool nextFrame(cv::Mat& img){
        status = dvpGetFrame(h, &frame, &p, 3000);
        img = cv::Mat(frame.iHeight, frame.iWidth, CV_8UC3, p);
        return status == DVP_STATUS_OK;
    }

    bool isOpened(){
        bool isValid;
        status = dvpIsValid(h, &isValid);
        return status == DVP_STATUS_OK && isValid;
    }

    DVPFrameCapture& operator >> (cv::Mat& img){
        nextFrame(img);
        return *this;
    }

    bool savePictureToFile(const char* fileName) {
        return DVP_STATUS_OK == dvpSavePicture(&frame, p, fileName, 100);
    }
};







#endif //OPENCV_HELLO_DVPFRAMECAPTURE_H
