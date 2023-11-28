#version 330 core

in float vertIsInstanced;
in vec4 vertInstanceColor;
in vec3 vertFragmentPosition;
in vec3 vertNormal;

uniform vec3 cameraPosition;  // Položaj kamere.
uniform vec3 lightPosition;  // Položaj luèi.
uniform vec3 lightIntensity;  // Intenziteta luèi.
uniform vec3 ambientIntensity;  // Intenziteta ambientne osvetlitve.
uniform float ambientReflectivity;  // Lastnosti materiala.
uniform float diffuseReflectivity;  // Difuzna odbojnost materiala.
uniform float specularReflectivity;  // Spekularna odbojnost materiala.
uniform int exponent;  // Eksponent spekularne osvetlitve.
uniform float shading;  // Spekularna odbojnost materiala.

out vec4 out_Color;  // Out fragment color.
uniform vec4 Color;  // Color vector in RGBA.


void main(void) {
    if (shading > 10.5) {
        vec3 normal = normalize(vertNormal.xyz);
        vec3 v = normalize(cameraPosition - vertFragmentPosition);  // Izraèun vektorja gledišèa.
        vec3 l = normalize(lightPosition - vertFragmentPosition);  // Izraèun vpadnega vektorja luèi.
        vec3 h = normalize(v + l);  // Izraèun normaliziranega vektorja h.

        // Izraèun Phongovega osvetlitvenega modela.
        vec3 ambientComponent = ambientIntensity * ambientReflectivity;
        vec3 diffuseComponent = abs(lightIntensity * diffuseReflectivity * dot(normal, l));
        vec3 specularComponent = lightIntensity * specularReflectivity * pow(abs(dot(normal, l)), exponent);

        vec3 RGB = (ambientComponent + diffuseComponent + specularComponent) / 3;
        out_Color = vec4(RGB.x, RGB.y, RGB.z, 1);
    }
    else {
        // Instanced color if instanced drawing.
        if (vertIsInstanced > 0.5) {
            out_Color = vertInstanceColor;
        }
        else {
            out_Color = Color;
        }
    }       
}
