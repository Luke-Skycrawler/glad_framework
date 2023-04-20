#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <spdlog/spdlog.h>
#include "shader.h"
#include "camera.h"
#include "mesh.h"
#include "model.h"
#include "utils.h"
#include "light.h"
#include "dynamics.h"
#include "globals.h"
#include <igl/boundary_facets.h>
#include <igl/edges.h>
//#include <tetgen.h>


void render_arrow(float len);
Globals globals;
unsigned int depthMapFBO,depthMap;
static const int SHADOW_WIDTH=800,SHADOW_HEIGHT=600;
bool model_draw = true, display_corner = true, move_light = false;
// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;
// camera
static const float z_camera = 2.0f;
Camera camera(glm::vec3(0.0f, 0.0f, z_camera));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;
static const float boundf = 1.0f;
// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

bool concentric = true, implicit = false;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

MatrixXi readTetFile(string filename);
MatrixXf readNodeFile(string filename);

vector<unsigned> bar_geometry(vector<vec3> &xcs);
vector<Eigen::Vector3f> read_obj(const string &filename,
                                 vector<unsigned> &indices)
{
    vector<Eigen::Vector3f> vertices;
    ifstream infile(filename);
    indices.resize(0);
    string line;
    while (getline(infile, line))
    {
        istringstream iss(line);
        string prefix;
        iss >> prefix;
        if (prefix == "v")
        {
            double x, y, z;
            iss >> x >> y >> z;
            Eigen::Vector3f vertex(x, y, z);
            // string key = to_string(x) + "," + to_string(y) + "," + to_string(z);
            // auto it = vertexMap.find(key);
            // if (it == vertexMap.end())
            // {
            //     vertexMap[key] = vertices.size();
            //     vertices.push_back(vertex);
            // }
            vertices.push_back(vertex);
        }
        else if (prefix == "f")
        {
            int v1, v2, v3;
            iss >> v1 >> v2 >> v3;
            // Eigen::Vector3i index;

            // index << v1 - 1, v2 - 1, v3 - 1;
            // indices.push_back(index);
            indices.push_back(v1 - 1);
            indices.push_back(v2 - 1);
            indices.push_back(v3 - 1);
        }
    }
    return vertices;
}
shayMesh::shayMesh(const string &filename)
{
    vertices = read_obj(filename, indices);
    nv = vertices.size();
    nf = indices.size() / 3;
    setupMesh();
}

shayMesh::shayMesh(const vector<vec3> &xcs, const vector<unsigned> &indices)
{
    nv = xcs.size();
    nf = indices.size() / 3;
    vertices.resize(nv);
    for (int i = 0; i < nv; i++)
    {
        vertices[i] = Vector3f(xcs[i][0], xcs[i][1], xcs[i][2]);
    }
    this->indices = indices;
    setupMesh();
}
void extract_edges(vector<Edge>& edges, const MatrixXi& T);

shayMesh::shayMesh(string node_file, string tet_file, std::vector<Edge> &edges)
{
    auto tets = readTetFile(tet_file);
    auto nodes = readNodeFile(node_file);
    nv = nodes.rows();
    MatrixXi F, J, K;

    igl::boundary_facets(tets, F, J, K);
    //igl::edges(tets, this -> E);
    extract_edges(edges,tets);
    nf = F.rows();
    
    vertices.resize(nv);
    //_indices.resize(nf * 3);
    indices.resize(nf * 3);
    auto max_coef = nodes.maxCoeff();
    
    for (int i = 0; i < nv; i++)
    {
        vertices[i] = nodes.row(i).transpose() / max_coef;
    }
    for (int i = 0; i < nf; i++)
    {
        indices[3 * i + 0] = F(i,0);
        indices[3 * i + 1] = F(i,1);
        indices[3 * i + 2] = F(i,2);
    }
    setupMesh();
}
void shayMesh::setupMesh()
{
    // create buffers/arrays
    nvn = nf * 3;
    vs.reserve(nvn);
    _indices.resize(nvn);
    for (int i = 0; i < nf * 3; i++)
    {
        vs.push_back({vertices[indices[i]], Vector3f(0, 0, 0)});
        _indices[i] = i;
    }
    compute_normals();

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);

    glBindVertexArray(vao);
    // load data into vertex buffers
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    // A great thing about structs is that their memory layout is sequential for all its items.
    // The effect is that we can simply pass a pointer to the struct and it translates perfectly to a glm::vec3/2 array which
    // again translates to 3/2 floats which translates to a byte array.
    glBufferData(GL_ARRAY_BUFFER, nvn * sizeof(vRender), vs.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, nf * sizeof(int) * 3, _indices.data(), GL_DYNAMIC_DRAW);

    // set the vertex attribute pointers
    // vertex Positions
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vRender), (void *)0);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vRender), (void *)offsetof(vRender, normals));
}

void shayMesh::draw()
{
    glBindVertexArray(vao);
    glDrawElements(GL_TRIANGLES, nf * 3, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}

void shayMesh::update_positions()
{
    for (int i = 0; i < nf * 3; i++)
    {
        vs[i].pos = vertices[indices[i]];
    }
    compute_normals();

    glBindVertexArray(vao);
    // load data into vertex buffers
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    // A great thing about structs is that their memory layout is sequential for all its items.
    // The effect is that we can simply pass a pointer to the struct and it translates perfectly to a glm::vec3/2 array which
    // again translates to 3/2 floats which translates to a byte array.
    glBufferData(GL_ARRAY_BUFFER, nvn * sizeof(vRender), vs.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, nf * sizeof(int) * 3, _indices.data(), GL_DYNAMIC_DRAW);

    // set the vertex attribute pointers
    // vertex Positions
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vRender), (void *)0);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vRender), (void *)offsetof(vRender, normals));
}

// lighting
glm::vec3 LightPositions[]={
    glm::vec3(0.0f, 0.0f, 2.5f),
    glm::vec3(-2.0f, -2.0f, 5.0f),
    glm::vec3(2.0f, 2.0f, 5.0f),
    glm::vec3(2.0f, -2.0f, 5.0f)   
};
glm::vec3 &lightPos(LightPositions[0]);

glm::mat4 from_eigen(Eigen::Matrix3d &eig_matrix){
    glm::mat3 a = glm::make_mat3(eig_matrix.data());
    glm::mat4 ret(a);
    return ret;
}

glm::vec3 from_eigen(Eigen::Vector3d &eig_vec) {
    glm::vec3 v = glm::make_vec3(eig_vec.data());
    return v;
}

void processInput(GLFWwindow* window)
{
    if (move_light) {
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            lightPos += 2.5f * deltaTime * camera.Front;
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            lightPos -= 2.5f * deltaTime * camera.Front;
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            lightPos -= 2.5f * deltaTime * camera.Right;
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            lightPos += 2.5f * deltaTime * camera.Right;
    }
    else {
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            camera.ProcessKeyboard(FORWARD, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            camera.ProcessKeyboard(BACKWARD, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            camera.ProcessKeyboard(LEFT, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            camera.ProcessKeyboard(RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    // if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    // {
    //     auto proj_t = -(camera.Front / camera.Front[2]) * z_camera;
    //     globals.body.force(proj_t[0], proj_t[1]);
    // }
    // else globals.body.f = glm::vec2(0.0);
}

shayMesh *triangle_highlight;
int indices_highlight = -1;
void reset_globals()
{
    std::ifstream f("../src/config.json");
    globals.config = json::parse(f);
    globals.selected = -1;
#define CUBE_CASE
#ifdef CUBE_CASE
    // shayMesh rendered_mesh("../src/assets/cube.obj");
    // auto body2_ptr = new shayMesh("../src/assets/dragon_8kface.obj");
    vector<Edge> edges;
    auto body2_ptr = new shayMesh("../src/assets/dragon/dragon.node", "../src/assets/dragon/dragon.ele", edges);
    // auto body2_ptr = new shayMesh("../src/assets/cube.obj");
    shayMesh &rendered_mesh = *body2_ptr;
    auto &vertices{rendered_mesh.vertices};

    // extract_edges(edges, body.meshes[0].indices);
    #define NODE
    #ifdef NODE
    auto &E {body2_ptr -> E};
    for (int i = 0; i < E.size(); i++)
    {
        edges.push_back({0.0, E(i, 0), E(i,1)});
    }
    #else
    extract_edges(edges, rendered_mesh.indices);
    #endif
    vector<vec3> velocity, position;
    vector<bool> is_static;

    velocity.resize(vertices.size());
    position.resize(vertices.size());
    is_static.resize(vertices.size());

    for (int i = 0; i < vertices.size(); i++)
    {
        auto &p{vertices[i]};
        position[i] = vec3(p[0], p[1], p[2]);
        velocity[i] = vec3(0.0, 0.0, 0.0);
        if (globals.config["static"])
        is_static[i] = position[i][0] == 0.0;
        else is_static[i] = false;
    }
    globals.mesh = new MassSpringMesh{velocity, position, edges, is_static};

    // tetgenio tet;
    // tet.load_node("a.node");
    // tet.load_elem("a.ele");
    // tet.facetlist;
    // tet.pointlist;
    // tet.edgelist;
    triangle_highlight = new shayMesh("../src/assets/tri.obj");

#else
    vector<vec3> xcs, vcs;
    auto bar_indices = bar_geometry(xcs);
    vector<Edge> bar_edges;
    extract_edges(bar_edges, bar_indices);
    shayMesh rendered_mesh(xcs, bar_indices);
    for (int i = 0; i < xcs.size(); i++)
        vcs[i] = vec3(0.0, 0.0, 0.0);
    globals.mesh = new MassSpringMesh{vcs, xcs, bar_edges};
    auto &vertices{rendered_mesh.vertices};
#endif
    init();
    globals.rendered_mesh = body2_ptr;
}

int main()
{
    // glfw: initialize and configure
    // ------------------------------------------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // ------------------------------------------------------------------
    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "glad framework", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetCharCallback(window, char_callback);
    glfwSetMouseButtonCallback(window, click_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    // glad: load all OpenGL function pointers
    // ------------------------------------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // ------------------------------------------------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader programs
    // ------------------------------------------------------------------
    // Shader lightingShader("../src/shaders/cursor/cursor.vert", "../src/shaders/cursor/cursor.frag", "../src/shaders/cursor/cursor.geom");
    Shader lightingShader("../src/shaders/shadow/shadow.vert", "../src/shaders/shadow/shadow.frag");
    // Shader lightingShader("../src/shaders/light/multi_light.vert", "../src/shaders/light/multi_light.frag");
    Shader depthShader("../src/shaders/depth/depth.vert", "../src/shaders/depth/depth.frag");
    Shader cornerShader("../src/shaders/corner/corner.vert", "../src/shaders/corner/corner.frag");
    // select buffers setup
    // ------------------------------------------------------------------
    unsigned int tex, buf;
    // Generate a name for the buffer object, bind it to the
    // GL_TEXTURE_BINDING, and allocate 4K for the buffer
    glGenBuffers(1, &buf);
    glBindBuffer(GL_TEXTURE_BUFFER, buf);
    glBufferData(GL_TEXTURE_BUFFER, sizeof(int), NULL, GL_DYNAMIC_READ);
    // Generate a new name for our texture
    glGenTextures(1, &tex);
    // Bind it to the buffer texture target to create it
    glBindTexture(GL_TEXTURE_BUFFER, tex);
    // Attach the buffer object to the texture and specify format as
    // single channel floating point
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32I, buf);
    // Now bind it for read-write to one of the image units
    glBindImageTexture(0, tex, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32I);
    // ------------------------------------------------------------------
    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    unsigned int quadVBO, quadVAO;
    float corner[] = {
        0.5f, 1.0f, 0.0f, 1.0f,
        0.5f, 0.5f, 0.0f, 0.0f,
        1.0f, 0.5f, 1.0f, 0.0f,
        0.5f, 1.0f, 0.0f, 1.0f,
        1.0f, 0.5f, 1.0f, 0.0f,
        1.0f, 1.0f, 1.0f, 1.0f};
    float quad[] = {
        -1.0f, 1.0f, 0.0f, 1.0f,
        -1.0f, -1.0f, 0.0f, 0.0f,
        1.0f, -1.0f, 1.0f, 0.0f,
        -1.0f, 1.0f, 0.0f, 1.0f,
        1.0f, -1.0f, 1.0f, 0.0f,
        1.0f, 1.0f, 1.0f, 1.0f};
    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &quadVBO);

    glBindVertexArray(quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad), &quad, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    unsigned int cornerVAO, cornerVBO;
    glGenVertexArrays(1, &cornerVAO);
    glGenBuffers(1, &cornerVBO);
    glBindVertexArray(cornerVAO);
    glBindBuffer(GL_ARRAY_BUFFER, cornerVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(corner), &corner, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    reset_globals();

    // TODO:  move this to initialization
    map<array<int, 2>, int> lut;
    int n_mass = globals.mesh->mass_x.size();
    SparseMatrix<double> sparse_matrix(n_mass * 3, n_mass * 3);
    gen_empty_sm(n_mass, globals.mesh->edges, sparse_matrix, lut);

    Light lights(LightPositions, 4);
    // load textures (we now use a utility function to keep the code more organized)
    // ------------------------------------------------------------------
    unsigned int diffuseMap = loadTexture("../src/assets/wood.png");
    unsigned int specularMap = loadTexture("../src/assets/wood_specular.png");

    cornerShader.setInt("screenTexture", 0);

    unsigned int texColorBuffer;
    glGenTextures(1, &texColorBuffer);
    glBindTexture(GL_TEXTURE_2D, texColorBuffer);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texColorBuffer, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, texColorBuffer, 0);

    unsigned int rbo;
    glGenRenderbuffers(1, &rbo);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, SCR_WIDTH, SCR_HEIGHT);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        std::cout << "error: framebuffer\n";
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    gen_preview_framebuffer();
    // ------------------------------------------------------------------
    // shader configuration
    // ------------------------------------------------------------------
    lightingShader.use();
    lightingShader.setInt("material.diffuse", 0);
    lightingShader.setInt("material.specular", 1);
    lightingShader.setInt("shadowMap", 2);

    lightingShader.setFloat("material.shininess", 64);

    // ------------------------------------------------------------------
    // render loop
    // ------------------------------------------------------------------
    int ts = 0;
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = glfwGetTime();
        // lightingShader.setFloat("time",currentFrame);
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        implicit_euler(sparse_matrix, lut);
        auto &vertices{globals.rendered_mesh->vertices};
        for (int i = 0; i < vertices.size(); i++)
        {
        auto &p{globals.mesh->mass_x[i]};
        vertices[i] = Vector3f(p[0], p[1], p[2]);
        }
        globals.rendered_mesh->update_positions();
        // input
        processInput(window);
        // render setup
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_STENCIL_TEST);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        glStencilFunc(GL_ALWAYS, 1, 0XFF);
        glStencilMask(0XFF);
        lightingShader.use();

        // ------------------------------------------------------------------
        // render
        // ------------------------------------------------------------------
        // be sure to activate shader when setting uniforms/drawing objects
        float scale = 1.02;
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 1.0f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();
        glm::mat4 model = glm::mat4(1.0f);
        glm::mat4 tmpmodel = glm::scale(model, glm::vec3(scale, scale, scale));
        glm::vec3 box2Pos(0.3, 0.0, 1.2);
        glm::mat4 lightSpaceTrans = glm::lookAt(lightPos, glm::vec3(0.0f), camera.WorldUp);

        const auto renderScene = [&](Shader &shader)
        {
            // renderPlane();

            // translate and render another
            model = glm::translate(model, glm::vec3(0.0f, 0.0f, -1.0f));
            shader.setMat4("model", model);
            renderPlane();

            model = glm::mat4(1.0f);
            model[1][1] = 0.0f;
            model[2][2] = 0.0f;
            model[1][2] = 1.0f;
            model[2][1] = 1.0f;
            shader.setMat4("model", glm::translate(model, glm::vec3(0.0f, 0.0f, boundf)));
            renderPlane();
            shader.setMat4("model", glm::translate(model, glm::vec3(0.0f, 0.0f, -boundf)));
            renderPlane();

            model = glm::mat4(1.0f);
            model[0][0] = 0.0f;
            model[2][2] = 0.0f;
            model[0][2] = 1.0f;
            model[2][0] = 1.0f;
            shader.setMat4("model", glm::translate(model, glm::vec3(0.0f, 0.0f, boundf)));
            renderPlane();
            shader.setMat4("model", glm::translate(model, glm::vec3(0.0f, 0.0f, -boundf)));
            renderPlane();

            shader.setMat4("model", glm::mat4(1.0f));
            // model = glm::translate(model, from_eigen(globals.body ->S[0].x));
            // shader.setMat4("model", model);
            shader.setVec3("objectColor", 0.8f, 0.8f, 0.8f);

            // body.Draw(shader);
            if (globals.config["wireframe"])
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

            globals.rendered_mesh->draw();
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

            auto &vs = triangle_highlight->vertices;
            if (indices_highlight != -1)
            {
                for (int i = 0; i < 3; i++)
                {
                    auto vid = globals.rendered_mesh->indices[indices_highlight * 3 + i];
                    vs[i] = globals.rendered_mesh->vertices[vid];
                }
                triangle_highlight->update_positions();
                shader.setVec3("objectColor", 1.0f, 0.0f, 0.0f);
                shader.setMat4("model", glm::mat4(1.0f));
                triangle_highlight->draw();
            }
        };
        if (display_corner)
        {
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        depthShader.use();

        // view/projection transformations
        model = glm::mat4(1.0f);

        // depthShader.setMat4("projection",projection);
        depthShader.setMat4("projection", glm::perspective(glm::radians(89.0f), (float)SHADOW_WIDTH / SHADOW_HEIGHT, 0.1f, 10.0f));
        depthShader.setMat4("view", lightSpaceTrans);
        depthShader.setMat4("model", model);
        depthShader.setVec3("viewPos", lightPos);
        renderScene(depthShader);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        model = glm::mat4(1.0f);
        }

        int viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        lightingShader.use();
        lightingShader.setVec2("pickPosition", glm::vec2(lastX / viewport[2] * 2 - 1.0f, (1 - lastY / viewport[3]) * 2 - 1.0f));
        lightingShader.setMat4("lightView", glm::perspective(glm::radians(89.0f), (float)SHADOW_WIDTH / SHADOW_HEIGHT, 0.1f, 10.0f) * lightSpaceTrans);
        view = camera.GetViewMatrix();
        lightingShader.setVec3("objectColor", 1.0f, 0.5f, 0.31f);
        lightingShader.setVec3("lightColor", 1.0f, 1.0f, 1.0f);
        lightingShader.setVec3("lightPos", lightPos);
        lightingShader.setVec3("viewPos", camera.Position);
        // view/projection transformations
        lightingShader.setMat4("projection", projection);
        lightingShader.setMat4("view", view);

        // world transformation
        lightingShader.setMat4("model", model);

        // bind diffuse map
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, diffuseMap);
        // bind specular map
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, specularMap);
        if (display_corner)
        {
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, depthMap);
        }
        renderScene(lightingShader);
        // also draw the lamp object
        // lights.Draw(camera);

        if (display_corner)
        {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glDisable(GL_DEPTH_TEST);
        cornerShader.use();
        glBindVertexArray(cornerVAO);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, depthMap);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        }
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // ------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------
    glDeleteVertexArrays(1, &quadVAO);
    glDeleteBuffers(1, &quadVBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    delete globals.mesh;
    return 0;
}
void gen_preview_framebuffer()
{
    glGenFramebuffers(1, &depthMapFBO);
    glGenTextures(1, &depthMap);
    glBindTexture(GL_TEXTURE_2D, depthMap);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,
                 SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0);
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);
}
// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow *window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;
    Vector2d xm(lastX / SCR_WIDTH, 1 - lastY / SCR_HEIGHT);
    xm = xm.array() * 2.0 - 1.0;
    camera.ProcessMouseMovement(xoffset, yoffset);
    if (globals.selected != -1)
    {
        globals.xm = xm;
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 1.0f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();
        auto P_glm = projection * view;
        Matrix4d P;
        P << P_glm[0][0], P_glm[1][0], P_glm[2][0], P_glm[3][0],
            P_glm[0][1], P_glm[1][1], P_glm[2][1], P_glm[3][1],
            P_glm[0][2], P_glm[1][2], P_glm[2][2], P_glm[3][2],
            P_glm[0][3], P_glm[1][3], P_glm[2][3], P_glm[3][3];
        globals.P = P;
        globals.P_inv = P.inverse();
        cout << xm.transpose() << "\n";
    }
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}

void char_callback(GLFWwindow* window, unsigned int codepoint)
{
    if (codepoint == 'c')
        display_corner=!display_corner;
    if (codepoint == 'l')
        move_light=!move_light;

    if (codepoint == 'r')
    {
        delete globals.mesh;
        delete globals.rendered_mesh;
        reset_globals();
        // init();
    }
}
void renderPlane(){

    static unsigned int planeVBO,planeVAO=0;
    static float planeVertices[] = {
        // positions            // normals         // texcoords
         25.0f,  25.0f, 0.0f,  0.0f, 1.0f, 0.0f,  25.0f,  0.0f,
        -25.0f,  25.0f, 0.0f,  0.0f, 1.0f, 0.0f,   0.0f,  0.0f,
        -25.0f, -25.0f, 0.0f,  0.0f, 1.0f, 0.0f,   0.0f, 25.0f,
         25.0f,  25.0f, 0.0f,  0.0f, 1.0f, 0.0f,  25.0f,  0.0f,
        -25.0f, -25.0f, 0.0f,  0.0f, 1.0f, 0.0f,   0.0f, 25.0f,
         25.0f, -25.0f, 0.0f,  0.0f, 1.0f, 0.0f,  25.0f, 25.0f
    };
    if(planeVAO == 0){
        // plane VAO
        glGenVertexArrays(1, &planeVAO);
        glGenBuffers(1, &planeVBO);
        glBindVertexArray(planeVAO);
        glBindBuffer(GL_ARRAY_BUFFER, planeVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(planeVertices), planeVertices, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
        glBindVertexArray(0);
    }
    glBindVertexArray(planeVAO);
    glDrawArrays(GL_TRIANGLES,0,6);
}
void renderCube(int light){
    static float vertices[] = {
        // positions          // normals           // texture coords
        -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  0.0f,
         0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  0.0f,
         0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
         0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
        -0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  1.0f,
        -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  0.0f,

        -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,
         0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  0.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
        -0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  1.0f,
        -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,

        -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
        -0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
        -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
        -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
        -0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
        -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

         0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
         0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
         0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
         0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
         0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
         0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

        -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,
         0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  1.0f,
         0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
         0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
        -0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  0.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,

        -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f,
         0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  1.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
        -0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,
        -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f
    };
    // first, configure the cube's VAO (and VBO)
    static unsigned int VBO=-1, cubeVAO=-1,lightCubeVAO;
    if(cubeVAO==-1){
        glGenVertexArrays(1, &cubeVAO);
        glGenBuffers(1, &VBO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

        glBindVertexArray(cubeVAO);

        // position attribute
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,8*sizeof(float),(void*)(3*sizeof(float)));
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,8*sizeof(float),(void*)(6*sizeof(float)));
        glEnableVertexAttribArray(2);

        // second, configure the light's VAO (VBO stays the same; the vertices are the same for the light object which is also a 3D cube)
        glGenVertexArrays(1, &lightCubeVAO);
        glBindVertexArray(lightCubeVAO);

        // we only need to bind to the VBO (to link it with glVertexAttribPointer), no need to fill it; the VBO's data already contains all we need (it's already bound, but we do it again for educational purposes)
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }
    glBindVertexArray(light?lightCubeVAO:cubeVAO);
    glDrawArrays(GL_TRIANGLES, 0, 36);
}
void render_arrow(float len){
    float max_len = 0.5f;
    float sqrt3by2 = sqrt(3.0f) / 2, width = 0.5, height = 2.0f, 
        bottom_width = 0.2f, bottom_height = 1.0f, 
        mid_height = 1.0f;
    height = fmax(len, max_len);
    mid_height = bottom_height = len < max_len ? height / 2.0f: len - max_len / 2.0f;    
    width = 0.25f * max_len;
    bottom_width = width * 0.4f;
    // cout << len << " " << " " << height << " " << width << " " << mid_height <<"\n";

    float vertices[] = {
        0.0, height, 0.0,
        width, mid_height, 0.0,
        -width / 2.0,mid_height, sqrt3by2 * width, 

        0.0, height, 0.0,
        width, mid_height, 0.0,
        -width / 2.0,mid_height, -sqrt3by2 * width,

        
        0.0, height, 0.0,
        -width / 2.0,mid_height, sqrt3by2 * width, 
        -width / 2.0,mid_height, -sqrt3by2 * width,

        width  , mid_height, 0.0,
        -width / 2.0,mid_height, sqrt3by2 * width, 
        -width / 2.0,mid_height, -sqrt3by2 * width,

        bottom_width  , mid_height, 0.0,
        -bottom_width / 2.0,mid_height, sqrt3by2 * bottom_width, 
        -bottom_width / 2.0,mid_height, -sqrt3by2 * bottom_width,
        
        bottom_width  , mid_height - bottom_height, 0.0,
        -bottom_width / 2.0, mid_height - bottom_height, sqrt3by2 * bottom_width, 
        -bottom_width / 2.0, mid_height - bottom_height, -sqrt3by2 * bottom_width,


        bottom_width  , mid_height, 0.0,
        bottom_width  , mid_height - bottom_height, 0.0,
        -bottom_width / 2.0, mid_height - bottom_height, sqrt3by2 * bottom_width, 

        bottom_width  , mid_height, 0.0,
        -bottom_width / 2.0, mid_height, sqrt3by2 * bottom_width, 
        -bottom_width / 2.0, mid_height - bottom_height, sqrt3by2 * bottom_width, 


        -bottom_width / 2.0, mid_height, sqrt3by2 * bottom_width, 
        -bottom_width / 2.0, mid_height, -sqrt3by2 * bottom_width,        
        -bottom_width / 2.0, mid_height - bottom_height, -sqrt3by2 * bottom_width,

        -bottom_width / 2.0, mid_height, sqrt3by2 * bottom_width,        
        -bottom_width / 2.0, mid_height - bottom_height, sqrt3by2 * bottom_width, 
        -bottom_width / 2.0, mid_height - bottom_height, -sqrt3by2 * bottom_width,


        bottom_width  , mid_height, 0.0,
        -bottom_width / 2.0, mid_height, -sqrt3by2 * bottom_width,
        bottom_width  , mid_height - bottom_height, 0.0,
        

        -bottom_width / 2.0, mid_height, -sqrt3by2 * bottom_width,
        bottom_width  , mid_height - bottom_height, 0.0,
        -bottom_width / 2.0, mid_height - bottom_height, -sqrt3by2 * bottom_width,

    };
    // first, configure the cube's VAO (and VBO)
    unsigned int vbo=-1, vao=-1;
    if(vao==-1){
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, 36);
}
void click_callback(GLFWwindow* window,int button,int action,int mods){
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        //    auto proj_t = camera.Front / camera.Front[2] * -(z_camera + 0.5f);
        //    globals.body.force(proj_t[0], proj_t[1]);
        Vector2d xm(lastX / SCR_WIDTH, 1 - lastY / SCR_HEIGHT);
        xm = xm.array() * 2.0 - 1.0;
        // normalize xm into range (-1, 1)
        // cout << xm.transpose();
        assert(xm(0) >= -1.0 && xm(0) <= 1.0);
        assert(xm(1) >= -1.0 && xm(1) <= 1.0);
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 1.0f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();
        auto P_glm = projection * view;
        Matrix4d P;
        P << P_glm[0][0], P_glm[1][0], P_glm[2][0], P_glm[3][0],
            P_glm[0][1], P_glm[1][1], P_glm[2][1], P_glm[3][1],
            P_glm[0][2], P_glm[1][2], P_glm[2][2], P_glm[3][2],
            P_glm[0][3], P_glm[1][3], P_glm[2][3], P_glm[3][3];
        auto [ti, z, abc] = raytrace_triangle(xm, P, globals.mesh->mass_x, globals.rendered_mesh->indices);
        if (ti == -1)
        {
            globals.selected = -1;
            return;
        }
        indices_highlight = ti;
        if (globals.config["verbose"])
            spdlog::info("indice = {}", indices_highlight);
        auto l = P_glm * glm::vec4(-0.5f, -0.5f, -0.5f, 1.0f);
        auto u = P_glm * glm::vec4(0.5f, 0.5f, 0.5f, 1.0f);
        if (globals.config["verbose"])
            spdlog::info("l, u = ({}, {}), ({}, {}), xm = ({}, {}) \n\n", l[0] / l[3], l[1] / l[3], u[0] / u[3], u[1] / u[3], xm[0], xm[1]);
        const auto argmax = [](const vec3 &v) -> int
        {
            // returns the max index of v
            if (v[0] > v[1] && v[0] > v[2])
                return 0;
            else if (v[1] > v[2])
                return 1;
            else
                return 2;
        };
        auto t = argmax(abc);
        auto selected_id = globals.rendered_mesh->indices[ti * 3 + t];
        Vector4d xs4d;
        xs4d << globals.mesh->mass_x[selected_id], 1.0;
        auto xs4d_ = (P * xs4d);
        auto xproj = xs4d_.head(2) / xs4d_[3];
        auto dx = xm - xproj;
        if (dx.norm() <globals.config["select_radius"])
        {
            globals.xm = xm;
            globals.selected = selected_id;
            globals.P_inv = P.inverse();
            globals.P = P;
            globals.dxm = dx;
            spdlog::info("selected = {}", globals.selected);
        }
        else
            globals.selected = -1;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        spdlog::info("released");
        globals.selected = -1;
        indices_highlight = -1;
    }
}

vector<unsigned> bar_geometry(vector<vec3> &xcs)
{
    double L = 1.0, w = 0.2;
    int n_x = 20, n_yz = 4;
    double dx = L / n_x;
    static const int faces[][4] = {
        {0, 1, 3, 2},
        {4, 5, 1, 0},
        {2, 3, 7, 6},
        {4, 0, 2, 6},
        {1, 5, 7, 3},
        {5, 4, 6, 7}};
    int n_elements = n_yz * n_yz * n_x;
    int n_nodes = (n_yz + 1) * (n_yz + 1) * (n_x + 1);
    int n_boundary_elements = n_yz * n_x * 4;
    const auto _trans = [=](int I)
    {
        int i = I / (n_yz * n_yz);
        int i_yz = I % (n_yz * n_yz);
        int j = i_yz / n_yz;
        int k = i_yz % n_yz;
        return (n_yz + 1) * (n_yz + 1) * i + (n_yz + 1) * j + k;
    };
    const auto trans = [=](int I) -> unsigned
    {
        int i = I / 4;
        int j = I % 4 / 2;
        int k = I % 2;
        return (n_yz + 1) * (n_yz + 1) * i + (n_yz + 1) * j + k;
    };
    auto nodes = new int[n_elements][8];
    vector<unsigned> indices;
    indices.resize(n_elements * 36);
    for (int e = 0; e < n_elements; e++)
    {
        e = _trans(e);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                {
                    int I = e + i * (n_yz + 1) * (n_yz + 1) + j * (n_yz + 1) + k;
                    int J = i * 4 + j * 2 + k;
                    nodes[e][J] = I;
                }
        for (int i = 0; i < 6; i++)
        {
            indices[e * 36 + i * 6 + 0] = trans(faces[i][0]) + e;
            indices[e * 36 + i * 6 + 1] = trans(faces[i][1]) + e;
            indices[e * 36 + i * 6 + 2] = trans(faces[i][2]) + e;
            indices[e * 36 + i * 6 + 3] = trans(faces[i][3]) + e;
            indices[e * 36 + i * 6 + 4] = trans(faces[i][4]) + e;
            indices[e * 36 + i * 6 + 5] = trans(faces[i][5]) + e;
        }
    }
    xcs.resize(n_nodes);
    const auto xc = [=](int I) -> vec3
    {
        int i = I / ((n_yz + 1) * (n_yz + 1));
        int i_yz = I % ((n_yz + 1) * (n_yz + 1));
        int j = i_yz / (n_yz + 1);
        int k = i_yz % (n_yz + 1);

        return vec3{i * 1.0 , j * 1.0, k * 1.0} * dx;
    };
    for (int i = 0; i < n_nodes; i++)
    {
        xcs[i] = xc(i);
    }
    return indices;
}

void shayMesh::compute_normals(bool clockwise)
{
    for (unsigned int i = 0; i < nvn; i++)
        vs[i].normals = Vector3f(0, 0, 0); // initialize all normals to zero
    auto &f{indices};
    auto &v{vertices};
    for (unsigned int i = 0; i < nf; i++)
    {
        Vector3f N = (vs[i * 3 + 1].pos - vs[i * 3 + 0].pos).cross(vs[i * 3 + 2].pos - vs[i * 3 + 0].pos); // face normal (not normalized)
        if (clockwise)
            N = -N;
        vs[i * 3 + 0].normals += N;
        vs[i * 3 + 1].normals += N;
        vs[i * 3 + 2].normals += N;
    }
    for (unsigned int i = 0; i < nvn; i++)
        vs[i].normals = vs[i].normals.normalized();
}

MatrixXf readNodeFile(string filename) {
    ifstream file(filename);

    if (!file) {
        cerr << "Error: Could not open file " << filename << endl;
        return MatrixXf::Zero(0, 0);
    }

    int num_vertices, dim, num_attrs, marker;
    file >> num_vertices >> dim >> num_attrs >> marker;

    MatrixXf V(num_vertices, 3);

    for (int i = 0; i < num_vertices; i++) {
        float x, y, z;
        int index, is_surface;

        file >> index >> x >> y >> z >> is_surface;

        V.row(i) << x, y, z;
    }

    file.close();
    return V;
}

MatrixXi readTetFile(string filename) {
    vector<Vector4i> tets;
    ifstream file(filename);

    if (!file) {
        cerr << "Error: Could not open file " << filename << endl;
        return MatrixXi::Zero(0, 0);
    }

    int num_tets, num_corners, num_attrs;
    file >> num_tets >> num_corners >> num_attrs;

    MatrixXi T(num_tets, 4);

    for (int i = 0; i < num_tets; i++) {
        int _, v1, v2, v3, v4;

        file >> _ >>  v1 >> v2 >> v3 >> v4;

        // Subtract 1 from vertex indices to convert from 1-indexing to 0-indexing
        T.row(i) << v1 - 1, v2 - 1, v3 - 1, v4 - 1;
    }

    file.close();
    return T;
}
