// ----------------------------------------------------------------------------
// main.cpp
//
//  Created on: 24 Jul 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: IGR201 Practical; OpenGL and Shaders (DO NOT distribute!)
//
// Copyright 2020 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Window parameters
GLFWwindow *g_window = nullptr;

// GPU objects
GLuint g_program = 0; // A GPU program contains at least a vertex shader and a fragment shader

// OpenGL identifiers
GLuint g_vao = 0;
GLuint g_posVbo = 0;
GLuint g_colorVbo = 0;
GLuint g_ibo = 0;

// All vertex positions packed in one array [x0, y0, z0, x1, y1, z1, ...]
std::vector<float> g_vertexPositions;
// All triangle indices packed in one array [v00, v01, v02, v10, v11, v12, ...] with vij the index of j-th vertex of the i-th triangle
std::vector<unsigned int> g_triangleIndices;
std::vector<float> g_vertexColors;




// Basic camera model
class Camera {
public:
  inline float getFov() const { return m_fov; }
  inline void setFoV(const float f) { m_fov = f; }
  inline float getAspectRatio() const { return m_aspectRatio; }
  inline void setAspectRatio(const float a) { m_aspectRatio = a; }
  inline float getNear() const { return m_near; }
  inline void setNear(const float n) { m_near = n; }
  inline float getFar() const { return m_far; }
  inline void setFar(const float n) { m_far = n; }
  inline void setPosition(const glm::vec3 &p) { m_pos = p; }
  inline glm::vec3 getPosition() { return m_pos; }

  inline glm::mat4 computeViewMatrix() const {
    return glm::lookAt(m_pos, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
  }

  // Returns the projection matrix stemming from the camera intrinsic parameter.
  inline glm::mat4 computeProjectionMatrix() const {
    return glm::perspective(glm::radians(m_fov), m_aspectRatio, m_near, m_far);
  }

private:
  glm::vec3 m_pos = glm::vec3(0, 0, 0);
  float m_fov = 45.f;        // Field of view, in degrees
  float m_aspectRatio = 1.f; // Ratio between the width and the height of the image
  float m_near = 0.1f; // Distance before which geometry is excluded fromt he rasterization process
  float m_far = 10.f; // Distance after which the geometry is excluded fromt he rasterization process
};
Camera g_camera;

class Mesh {
public:
    void init() { // should properly set up the geometry buffer
        glCreateVertexArrays(1, &m_vao);
        glBindVertexArray(m_vao);

        size_t vertexBufferSize = sizeof(float) * m_vertexPositions.size(); // Gather the size of the buffer from the CPU-side vector
        glCreateBuffers(1, &m_posVbo);
        glNamedBufferStorage(m_posVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
        glNamedBufferSubData(m_posVbo, 0, vertexBufferSize, m_vertexPositions.data()); // Fill the data storage from a CPU array
        glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

        size_t colorBufferSize = sizeof(float) * m_vertexColors.size(); // Gather the size of the buffer from the CPU-side vector
        glCreateBuffers(1, &m_colVbo);
        glNamedBufferStorage(m_colVbo, colorBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
        glNamedBufferSubData(m_colVbo, 0, colorBufferSize, m_vertexColors.data()); // Fill the data storage from a CPU array
        glBindBuffer(GL_ARRAY_BUFFER, m_colVbo);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

        size_t normalBufferSize = sizeof(float) * m_vertexNormals.size(); // Gather the size of the buffer from the CPU-side vector
        glCreateBuffers(1, &m_normalVbo);
        glNamedBufferStorage(m_normalVbo, normalBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
        glNamedBufferSubData(m_normalVbo, 0, normalBufferSize, m_vertexNormals.data()); // Fill the data storage from a CPU array
        glBindBuffer(GL_ARRAY_BUFFER, m_normalVbo);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

        size_t indexBufferSize = sizeof(unsigned int) * m_triangleIndices.size();
        glCreateBuffers(1, &m_ibo);
        glNamedBufferStorage(m_ibo, indexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
        glNamedBufferSubData(m_ibo, 0, indexBufferSize, m_triangleIndices.data());

    }
    
    static std::shared_ptr<Mesh> genSphere(const size_t resolution = 360) {
        std::shared_ptr<Mesh> mesh(new Mesh);

        const size_t numVertice = resolution * resolution;
        const float stepPhi   = 2.f * M_PI / (resolution);
        const float stepTheta = 2.f * M_PI / (resolution);

        
        // mesh->m_triangleIndices.clear();
        for (unsigned int i = 0; i < resolution; i++) {
            for (unsigned int j = 0; j < resolution; j++) {
                // x = r*sin(i)*cos(j)
                // y = r*sin(i)*sin(j)
                // z = r*cos(i)
                mesh->m_vertexPositions.push_back(sin(stepPhi* i) * cos(stepTheta*j));
                mesh->m_vertexPositions.push_back(sin(stepPhi* i) * sin(stepTheta*j));
                mesh->m_vertexPositions.push_back(cos(stepPhi* i));

               
                mesh->m_vertexColors.push_back(sin(stepPhi * i) * cos(stepTheta * j));
                mesh->m_vertexColors.push_back(sin(stepPhi * i) * sin(stepTheta * j));
                mesh->m_vertexColors.push_back(cos(stepPhi * i));
               

                mesh->m_vertexNormals.push_back(sin(stepPhi * i) * cos(stepTheta * j));
                mesh->m_vertexNormals.push_back(sin(stepPhi * i) * sin(stepTheta * j));
                mesh->m_vertexNormals.push_back(cos(stepPhi * i));

            }
        }
        
        for (int i = 0; i < resolution; i++) {
            for (int j = 0; j < resolution; j++) {
                mesh->m_triangleIndices.push_back(resolution * i + j);
                mesh->m_triangleIndices.push_back((resolution * i + (j + 1) % resolution));
                mesh->m_triangleIndices.push_back((resolution * (i + 1) + (j + 1) % resolution));

                mesh->m_triangleIndices.push_back((resolution * (i + 1) + (j + 1) % resolution));
                mesh->m_triangleIndices.push_back((resolution * (i + 1) + j));
                mesh->m_triangleIndices.push_back(resolution * i + j);
            }
        }
        return mesh;
    }
    // ...
private:
    std::vector<float> m_vertexPositions;
    std::vector<float> m_vertexNormals;
    std::vector<unsigned int> m_triangleIndices;
    std::vector<float> m_vertexColors;

    GLuint m_vao = 0;
    GLuint m_posVbo = 0;
    GLuint m_normalVbo = 0;
    GLuint m_ibo = 0;
    GLuint m_colVbo;

    // Constants
    const static float kSizeSun = 1;
    const static float kSizeEarth = 0.5;
    const static float kSizeMoon = 0.25;
    const static float kRadOrbitEarth = 10;
    const static float kRadOrbitMoon = 2;
    // Model transformation matrices
    glm::mat4 g_sun, g_earth, g_moon;

    // ...

};
std::shared_ptr<Mesh> g_mesh;

GLuint loadTextureFromFileToGPU(const std::string &filename) {
  int width, height, numComponents;
  // Loading the image in CPU memory using stbd_image
  unsigned char *data = stbi_load(
    filename.c_str(),
    &width, &height,
    &numComponents, // 1 for a 8 bit greyscale image, 3 for 24bits RGB image, 4 for 32bits RGBA image
    0);

  GLuint texID;
  // TODO:

  // Freeing the now useless CPU memory
  stbi_image_free(data);
  glBindTexture(GL_TEXTURE_2D, 0); // unbind the texture

  return texID;
}

// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window.
void windowSizeCallback(GLFWwindow* window, int width, int height) {
  g_camera.setAspectRatio(static_cast<float>(width)/static_cast<float>(height));
  glViewport(0, 0, (GLint)width, (GLint)height); // Dimension of the rendering region in the window
}

// Executed each time a key is entered.
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
  if(action == GLFW_PRESS && key == GLFW_KEY_W) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  } else if(action == GLFW_PRESS && key == GLFW_KEY_F) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  } else if(action == GLFW_PRESS && (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q)) {
    glfwSetWindowShouldClose(window, true); // Closes the application if the escape key is pressed
  }
}

void errorCallback(int error, const char *desc) {
  std::cout <<  "Error " << error << ": " << desc << std::endl;
}

void initGLFW() {
  glfwSetErrorCallback(errorCallback);

  // Initialize GLFW, the library responsible for window management
  if(!glfwInit()) {
    std::cerr << "ERROR: Failed to init GLFW" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Before creating the window, set some option flags
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create the window
  g_window = glfwCreateWindow(
    1024, 768,
    "Interactive 3D Applications (OpenGL) - Simple Solar System",
    nullptr, nullptr);
  if(!g_window) {
    std::cerr << "ERROR: Failed to open window" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  // Load the OpenGL context in the GLFW window using GLAD OpenGL wrangler
  glfwMakeContextCurrent(g_window);
  glfwSetWindowSizeCallback(g_window, windowSizeCallback);
  glfwSetKeyCallback(g_window, keyCallback);
}

void initOpenGL() {
  // Load extensions for modern OpenGL
  if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cerr << "ERROR: Failed to initialize OpenGL context" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  glCullFace(GL_BACK); // Specifies the faces to cull (here the ones pointing away from the camera)
  glEnable(GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).
  glDepthFunc(GL_LESS);   // Specify the depth test for the z-buffer
  glEnable(GL_DEPTH_TEST);      // Enable the z-buffer test in the rasterization
  glClearColor(0.7f, 0.7f, 0.7f, 1.0f); // specify the background color, used any time the framebuffer is cleared
}

// Loads the content of an ASCII file in a standard C++ string
std::string file2String(const std::string &filename) {
  std::ifstream t(filename.c_str());
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

// Loads and compile a shader, before attaching it to a program
void loadShader(GLuint program, GLenum type, const std::string &shaderFilename) {
  GLuint shader = glCreateShader(type); // Create the shader, e.g., a vertex shader to be applied to every single vertex of a mesh
  std::string shaderSourceString = file2String(shaderFilename); // Loads the shader source from a file to a C++ string
  const GLchar *shaderSource = (const GLchar *)shaderSourceString.c_str(); // Interface the C++ string through a C pointer
  glShaderSource(shader, 1, &shaderSource, NULL); // load the vertex shader code
  glCompileShader(shader);
  glAttachShader(program, shader);
  glDeleteShader(shader);
}

void initGPUprogram() {
  g_program = glCreateProgram(); // Create a GPU program, i.e., two central shaders of the graphics pipeline
  loadShader(g_program, GL_VERTEX_SHADER, "vertexShader.glsl");
  loadShader(g_program, GL_FRAGMENT_SHADER, "fragmentShader.glsl");
  glLinkProgram(g_program); // The main GPU program is ready to be handle streams of polygons

  glUseProgram(g_program);
}

// Define your mesh(es) in the CPU memory
//std::vector<float> g_vertexPositions;
//std::vector<unsigned int> g_triangleIndices;
//std::vector<float> g_vertexColors;
void initCPUgeometry() {
  
        g_vertexPositions = { // the array of vertex positions [x0, y0, z0, x1, y1, z1, ...]
        0.f, 0.f, 0.f,
        1.f, 0.f, 0.f,
        0.f, 1.f, 0.f
        };
        g_triangleIndices = { 0, 1, 2 }; // indices just for one triangle

        g_vertexColors = { // the array of vertex colors [r0, g0, b0, r1, g1, b1, ...]
        1.f, 0.f, 0.f,
        0.f, 1.f, 0.f,
        0.f, 0.f, 1.f
        };
        **/
    //g_mesh = Mesh::genSphere();
   // g_sun = glm::scale(glm::mat4(1.0), glm::vec3(kSizeSun));
}

void initGPUgeometry() {
    
  // Create a single handle that joins together attributes (vertex positions,
  // normals) and connectivity (triangles indices)
  glCreateVertexArrays(1, &g_vao);
  glBindVertexArray(g_vao);

  // Generate a GPU buffer to store the positions of the vertices
  size_t vertexBufferSize = sizeof(float)*g_vertexPositions.size(); // Gather the size of the buffer from the CPU-side vector
  glCreateBuffers(1, &g_posVbo);
  glNamedBufferStorage(g_posVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
  glNamedBufferSubData(g_posVbo, 0, vertexBufferSize, g_vertexPositions.data()); // Fill the data storage from a CPU array
  glBindBuffer(GL_ARRAY_BUFFER, g_posVbo);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), 0);

  // Generate a GPU buffer to store the colors
  size_t colorBufferSize = sizeof(float) * g_vertexColors.size(); // Gather the size of the buffer from the CPU-side vector
  glCreateBuffers(1, &g_colorVbo);
  glNamedBufferStorage(g_colorVbo, colorBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
  glNamedBufferSubData(g_colorVbo, 0, colorBufferSize, g_vertexColors.data()); // Fill the data storage from a CPU array
  glBindBuffer(GL_ARRAY_BUFFER, g_colorVbo);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  glBindVertexArray(0); // deactivate the VAO for now, will be activated at rendering time

  // Same for the index buffer that stores the list of indices of the
  // triangles forming the mesh
  size_t indexBufferSize = sizeof(unsigned int)*g_triangleIndices.size();
  glCreateBuffers(1, &g_ibo);
  glNamedBufferStorage(g_ibo, indexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
  glNamedBufferSubData(g_ibo, 0, indexBufferSize, g_triangleIndices.data());
  
 }

void initCamera() {
  int width, height;
  glfwGetWindowSize(g_window, &width, &height);
  g_camera.setAspectRatio(static_cast<float>(width)/static_cast<float>(height));

  g_camera.setPosition(glm::vec3(0.0, 0.0, 3.0));
  g_camera.setNear(0.1);
  g_camera.setFar(80.1);
}

void init() {
  initGLFW();
  initOpenGL();
//  initCPUgeometry();
  initGPUprogram();
  g_mesh = Mesh::genSphere();
 // initGPUgeometry();
  g_mesh->init();
  initCamera();
}

void clear() {
  glDeleteProgram(g_program);
  g_mesh.reset();

  glfwDestroyWindow(g_window);
  glfwTerminate();
}

// The main rendering call
void render() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Erase the color and z buffers.

  const glm::mat4 viewMatrix = g_camera.computeViewMatrix();
  const glm::mat4 projMatrix = g_camera.computeProjectionMatrix();

  glUniformMatrix4fv(glGetUniformLocation(g_program, "viewMat"), 1, GL_FALSE, glm::value_ptr(viewMatrix)); // compute the view matrix of the camera and pass it to the GPU program
  glUniformMatrix4fv(glGetUniformLocation(g_program, "projMat"), 1, GL_FALSE, glm::value_ptr(projMatrix)); // compute the projection matrix of the camera and pass it to the GPU program

  
  const glm::vec3 camPosition = g_camera.getPosition();
  glUniform3f(glGetUniformLocation(g_program, "camPos"), camPosition[0], camPosition[1], camPosition[2]);


  glBindVertexArray(g_vao);     // bind the VAO storing geometry data
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ibo); // bind the IBO storing geometry data
  glDrawElements(GL_TRIANGLES, g_triangleIndices.size(), GL_UNSIGNED_INT, 0); // Call for rendering: stream the current GPU geometry through the current GPU program
  g_mesh->render();


  

}

// Update any accessible variable based on the current time
void update(const float currentTimeInSec) {
  // std::cout << currentTimeInSec << std::endl;

}

int main(int argc, char ** argv) {
  init(); // Your initialization code (user interface, OpenGL states, scene with geometry, material, lights, etc)
  while(!glfwWindowShouldClose(g_window)) {
    update(static_cast<float>(glfwGetTime()));
    render();
    glfwSwapBuffers(g_window);
    glfwPollEvents();
  }
  clear();
  return EXIT_SUCCESS;
}
