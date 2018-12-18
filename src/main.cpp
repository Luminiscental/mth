
#include <cmath>
#include <iostream>

#include <Maths/constants.h>
#include <Maths/vec2.h>
#include <Maths/vec3.h>
#include <Maths/vec4.h>

// Don't judge my crappy ghetto testing system

template <typename T>
bool test(T a, T b) {

    if (a == b) {

        return true;

    } else {

        std::cerr << a << " != " << b << std::endl;
        return false;
    }
}

bool testVec2() {

    return test(m::ivec2::dot(m::ivec2(2, 10) / 2, -m::ivec2(-3, -4) * 6), 6 * 23) &&
           test(m::ivec2(3, 3).magnSqr(), m::ivec2::dot(m::ivec2(3, 3), m::ivec2(3, 3))) &&
           test(m::ivec2::dot(m::ivec2(3, 3), m::ivec2(3, 3)), 18) &&
           test(std::abs(m::ivec2(1, -1).arg() + 0.7853981634) < Maths_EPSILON(float), true); 
}

bool testVec3() {

    return test(3 * m::ivec3(1, 2, -4), m::ivec3(3, 6, -12)) &&
           test(m::ivec2::dot(m::ivec3(1, -4, 2).xy(), m::ivec2(2, -1)) * m::ivec3(1, m::ivec2(1, -1)), m::ivec3(6, 6, -6)) &&
           test(m::ivec3::cross(m::ivec3(2, 1, 0), m::ivec3(1, -1, -3)).magnSqr(), 54);
}

bool testVec4() {
    return test(2 * m::ivec4(m::ivec3(6, 0, 3), -9) / 3, m::ivec4(4, 0, 2, -6)) &&
           test(m::ivec4(m::ivec2(1, 3), m::ivec2(-2, 1)).magnSqr(), 15) &&
           test(m::ivec4::dot(m::ivec4(2, m::ivec2(1, -1), 3), m::ivec4(1, 0, 1, 0)), 1);
}

int main(int argc, const char **argv) {

    if (!testVec2()) std::cerr << "vec2 test failed!" << std::endl;
    else std::cout << "vec2: passed" << std::endl;
    if (!testVec3()) std::cerr << "vec3 test failed!" << std::endl;
    else std::cout << "vec3: passed" << std::endl;
    if (!testVec4()) std::cerr << "vec4 test failed!" << std::endl;
    else std::cout << "vec4: passed" << std::endl;

    return 0;
}
