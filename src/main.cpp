
#include <cmath>
#include <iostream>

#include <Maths/constants.h>
#include <Maths/vec2.h>
#include <Maths/vec3.h>
#include <Maths/vec4.h>
#include <Maths/mat2.h>

//TODO: mat3, mat4

int main(int argc, const char **argv) {

    std::cout << "EPSILON(float) = " << Maths_EPSILON(float) << std::endl;

    std::cout << "dot(12, 74) = " << m::ivec2::dot(m::ivec2(1, 2), m::ivec2(7, 4)) << std::endl;

    std::cout << "cross(123, 712) = " << m::ivec3::cross(m::ivec3(1, 2, 3), m::ivec3(7, 1, 2)) << std::endl;

    std::cout << "magn(1234) = " << m::ivec4(1, 2, 3, 4).magn() << std::endl;

    std::cout << "inverse(1234) = " << m::mat2{1, 2, 3, 4}.inverse() << std::endl;

    return 0;
}
