#version 330 core
layout (location = 0) in vec3 in_Position;    // Entry position.
layout (location = 1) in vec3 in_Normal;      // Normal vector of the vertex.
layout (location = 2) in vec4 in_InstanceColor;  // Color of the instance.
layout (location = 3) in mat4 in_InstanceModel;  // Instance model matrix.

uniform bool instanced;     // Boolean for instanced drawing (if true instanced drawing).
uniform mat4 projection;    // Projection matrix.
uniform mat4 view;          // View matrix.
uniform mat4 model;         // Model matrix.

out float vertIsInstanced;      // Instanced drawing if value == 1.0;
out vec4 vertInstanceColor;  // Instanced color.
out vec3 vertFragmentPosition;
out vec3 vertNormal;

void main() {
    vertNormal = mat3(transpose(inverse(model))) * in_Normal;
    vertFragmentPosition = vec3(model * vec4(in_Position, 1.0));
    
    // Instanced draw of the object.
    if (instanced) {
        vertIsInstanced = 1.0;
        vertInstanceColor = in_InstanceColor;
        gl_Position = projection * view * in_InstanceModel * vec4(in_Position, 1.0);
    }
    // Non-instanced draw of the object.
    else {
        vertIsInstanced = 0.0;
        gl_Position = projection * view * model * vec4(in_Position, 1.0);
    }
}
