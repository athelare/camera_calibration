echo %1
echo off
scp -r config CMakeLists.txt main.cpp %1:~/calibrate
ssh %1 "cd ~/calibrate; mkdir build; cd build; cmake ..; make; mv camera_calibration ..; cd ..; rm -rf build CMakeLists.txt main.cpp;"