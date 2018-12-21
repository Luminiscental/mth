
#include <cmath>
#include <iostream>

#include <Maths/constants.h>
#include <Maths/vec2.h>
#include <Maths/vec3.h>
#include <Maths/vec4.h>
#include <Maths/mat2.h>
#include <Maths/mat3.h>
#include <Maths/mat4.h>

// TODO: quaternions

int main(int argc, const char **argv) {

    std::cout << "EPSILON(float) = " << Maths_EPSILON(float) << std::endl;

    std::cout << "dot(12, 74) = " << m::ivec2::dot(m::ivec2(1, 2), m::ivec2(7, 4)) << std::endl;

    std::cout << "cross(123, 712) = " << m::ivec3::cross(m::ivec3(1, 2, 3), m::ivec3(7, 1, 2)) << std::endl;

    std::cout << "magn(1234) = " << m::ivec4(1, 2, 3, 4).magn() << std::endl;

    std::cout << "inverse(1234) = " << m::mat2{1, 2, 3, 4}.inverse() << std::endl;

    std::cout << std::endl << std::endl;

    m::mat3 test{4, 2, 3, 4, 5, 6, 7, 8, 9};

    std::cout << "test = " << test << std::endl;

    std::cout << "minors(0, 1) = " << test.getMinor(0, 1) << std::endl;

    std::cout << "det(test) = " << test.det() << std::endl;

    std::cout << "cofactors(test) = " << test.cofactors() << std::endl;

    std::cout << "inverse(test) = " << test.inverse() << std::endl;

    std::cout << "test * test = " << test * test << std::endl;

    std::cout << "test * inverse(test) = " << test * test.inverse() << std::endl;

    m::mat3 scaleTest = m::mat3::scale(1, 3, 2);

    std::cout << "scaleTest = " << scaleTest << std::endl;

    std::cout << "123 scaled = " << scaleTest * m::vec3(1, 2, 3) << std::endl;

    m::mat4 test4{2, 3, 4, 5, 1, 6, -3, 4, 1, 5, 5, 1, 1, 1, 0, 1};

    std::cout << "test4 = " << test4 << std::endl;

    std::cout << "minors(2, 2) = " << test4.getMinor(2, 2) << std::endl;

    std::cout << "det(test4) = " << test4.det() << std::endl;

    std::cout << "inverse(test4) = " << test4.inverse() << std::endl;

    std::cout << "test4^2 = " << test4 * test4 << std::endl;

    std::cout << "test4 * inverse(test4) = " << test4 * test4.inverse() << std::endl;

    m::mat4 projTest = m::mat4::perspectiveProjection(60.0f, 60.0f, 1.0f, 10.0f);

    std::cout << "projTest = " << projTest << std::endl;

    std::cout << "513 projected = " << projTest * m::vec4(m::vec3(5, 1, 3), 1) << std::endl;

    return 0;
}
