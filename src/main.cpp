
#include <iostream>

#include <Maths/vec.h>
#include <Maths/mat.h>
#include <Maths/quat.h>

#include <Maths/constants.h>
#include <Maths/equations.h>

int main(int argc, const char **argv) {

    std::cout << std::endl;

    // x + 3y - 2z + w = 12
    // 2x - 6y + 13z - 2w = 1
    // x + y + z + 2w = 0
    // x - 2y + 3z - 4w = -3

    m::mat4 systemCoeffs{ 1,  2,  1,  1, // x
                          3, -6,  1, -2, // y
                         -2, 13,  1,  3, // z
                          1, -2,  2, -4};// w

    m::vec4 systemOutput{12, 1, 0, -3};

    m::eq::LinearSystem<float, 4> system(systemCoeffs, systemOutput);

    std::cout << "x + 3y - 2z + w = 12" << std::endl
              << "2x - 6y + 13z - 2w = 1" << std::endl
              << "x + y + z + 2w = 0" << std::endl
              << "x - 2y + 3z - 4w = -3" << std::endl << std::endl;

    m::vec4 solution = system.solve();

    std::cout << " -> x = " << solution.x() << ", y = " << solution.y() << ", z = " << solution.z() << ", w = " << solution.w() << std::endl << std::endl;

    // (1 + i) + (2 - i)x + (3 + i)x^2 = 0
    
    m::cvec3 quadraticCoeffs{std::complex<double>(1, 1), std::complex<double>(2, -1), std::complex<double>(3, 1)};

    m::eq::Polynomial<2> quadratic(quadraticCoeffs);

    std::cout << "(1 + i) + (2 - i)x + (3 + i)x^2 = 0" << std::endl << std::endl
              << " -> x = ";

    auto quadraticSolutions = quadratic.solve();
    
    for (std::complex<double> root : quadraticSolutions) {

        std::cout << root.real() << " + " << root.imag() << "i"; // TODO: Replace std::complex

        if (root != *quadraticSolutions.end()) std::cout << " or ";
    }

    std::cout << std::endl << std::endl;

    m::quat xRotation = m::quat::rotation(Maths_PI(float) / 4, m::x_axis<float>);
    m::quat yRotation = m::quat::rotation(Maths_PI(float) / 6, m::y_axis<float>);

    m::quat netRotation = yRotation * xRotation;

    m::vec3 rotatedVector = netRotation.rotate(m::vec3{1, 1, 2});

    std::cout << "(1, 1, 2) rotated 45 degrees around the X-axis and then 30 degrees around the Y-axis = " << rotatedVector << std::endl << std::endl;
    
    return 0;
}
