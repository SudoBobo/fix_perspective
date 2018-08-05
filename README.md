## Overview

Simple program which fix photo perspective (from "trapeze" perspective to "scan-like")
Tested with `.jpg` images with dark background. Made for educational purposes using OpenCV.
![alt text](https://github.com/SudoBobo/fix_perspective/blob/master/2.jpg)
![alt text](https://github.com/SudoBobo/fix_perspective/blob/master/2_fixed.jpg)

![alt text](https://github.com/SudoBobo/fix_perspective/blob/master/v3.jpg)
![alt text](https://github.com/SudoBobo/fix_perspective/blob/master/v3_fixed.jpg)

![alt text](https://github.com/SudoBobo/fix_perspective/blob/master/z.jpg)
![alt text](https://github.com/SudoBobo/fix_perspective/blob/master/z_fixed.jpg)


### Requirements
```
CMake
gcc or g++ or Clang
OpenCV => 3.4
```
Instructions on how to install OpenCV can be found here:
[official OpenCV site (instructions for linux)](https://docs.opencv.org/3.4/d7/d9f/tutorial_linux_install.html)

### Usage
```
cd fix_perspective
cmake CMakeLists.txt
make
./fix_perspective image.jpg

```
Result will be saved in fix_perspective/image_fixed.jpg
