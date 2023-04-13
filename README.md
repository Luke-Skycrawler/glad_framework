
## Windows building
Make sure you have [`cmake`](https://cmake.org/download/) and `vcpkg` installed. Detailed instructions can be found [here](https://vcpkg.io/en/getting-started.html).

Navigate to the vcpkg directory, and installed the following packages via vcpkg,
`vcpkg install glfw3 glad glm assimp`

Navigate to the root directory `glad_framework`, and type in  
`cmake -B "build" -S . -DCMAKE_TOOLCHAIN_FILE=[path to vcpkg]/scripts/buildsystems/vcpkg.cmake`
cmake -B "build" -S . -DCMAKE_TOOLCHAIN_FILE=D:/vcpkg/scripts/buildsystems/vcpkg.cmake
`cmake --build build`

## Linux building

Installing on Linux is roughly the same, with a few dependencies to be installed via the system package manager. On Debian and Ubuntu, those packages can be installed with `apt-get install  libxinerama-dev libxcursor-dev libgl1-mesa-dev xorg-dev`. 
Also, you must first make sure you have CMake, Git, and GCC by typing as root (`sudo`) `apt-get install g++ cmake git`. 

With those packages installed, follow the instructions on https://vcpkg.io/en/getting-started.html to install `vcpkg`, then use `vcpkg` to install the packages. Best to do it with `sudo` to suppress possible errors.
`vcpkg install glfw3 glad glm assimp`

Next, navigate root directory `glad_framework`, and config the project with
`cmake -B "build" -S . -DCMAKE_TOOLCHAIN_FILE=[path to vcpkg]/scripts/buildsystems/vcpkg.cmake`.
and build it with `cmake --build build`.

## Mac OS X building
Building on Mac OS X is fairly simple:
```
brew install cmake assimp glm glfw freetype
cmake -S . -B build
cmake --build build -j$(sysctl -n hw.logicalcpu)
```

## Run
Depending on your build system, the executable directory can differ. Make sure you have all the binaries copied directly under `glad_framework\src` (with `assets\` and `shaders\` in it).  


## Code Structure

The code is hugely based on the awesome [online opengl tutorial](https://learnopengl.com/Getting-started/OpenGL/), and you can refer to their website as documentation of specific code. 

### Load and manipulate models
We've provided an example in `main.cpp`.
First, we load a bunny with 
```C++
    Model bunny("assets/bunny/bunny.obj");
```
Then within the main `while (!glfwWindowShouldClose(window))` loop, we set the model matrix in shaders and call `Model.draw(Shader)` to do the rendering. 
```C++
    depthShader.setMat4("model",glm::translate(glm::mat4(1.0f),glm::vec3(0.0f,0.0f,0.0f)));
    bunny.Draw(depthShader);
```
To modify the scene, you can add your code to the lambda expression `renderScene`; 

### Control
The keyboard control function is defined in `void processInput(GLFWwindow *window)` in `main.cpp`. Currently, the keys `w`, `a`, `s`, `d` are occupied for moving and `l`, `c` are occupied by switching light/camera movement and switching on/off of the upper-right corner display.


## Reference
- [LearnOpenGL](https://github.com/JoeyDeVries/LearnOpenGL)