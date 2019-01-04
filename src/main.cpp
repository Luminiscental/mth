
#include <iostream>

#include <m/m.h>

#include <m/vec.h>
#include <m/mat.h>
#include <m/quat.h>
#include <m/polynomial.h>

int main() {

    m::mat4 transformation = m::mat::rotation(m::PI<float> / 3, m::X_AXIS<float>)
                           * m::mat::translate(m::vec3(1.0f, 13.0f, -2.0f));

    m::vec3 initialPosition(1, 2, 1);
    m::vec4 positionHandle(initialPosition.x(), initialPosition.y(), initialPosition.z(), 1);
    
    m::vec3 transformedPosition = (transformation * positionHandle).xyz();

    m::quat additionalRotation = m::quat::rotation(m::TAU<float> / 12, m::vec3(0.5f, 0.5f, -0.1f));

    additionalRotation *= m::quat(2, 3, 1, 5);
    additionalRotation.real() -= 4;

    transformedPosition = additionalRotation.rotate(transformedPosition);

    std::cout << std::endl << "Transformation matrix:" << std::endl;
    std::cout << std::endl << m::mat::rotation(additionalRotation) * transformation << std::endl;

    std::cout << std::endl << "Went from " << initialPosition << " to " << transformedPosition << std::endl;

    m::tmat<double, 8, 8> bigMatrix(0.0, 3.0, -1.0, 5.0, 7.0, 6.0, 2.0, -3.0,
                                    1.0, 5.0, 2.0, -3.0, 2.0, -2.0, -2.0, 5.0,
                                    7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,
                                    3.0, 2.0, -1.0, 3.0, -5.0, -6.0, -1.0, 1.0,
                                    5.0, 3.0, 4.0, 3.0, 5.0, -2.0, -3.0, 1.0,
                                    6.0, 4.0, 3.0, 2.0, 7.0, -5.0, -7.0, 11.0,
                                    0.0, 2.0, 0.0, 1.5, 2.5, -3.0, 5.0, 4.0,
                                    2.0, 4.0, 8.0, 12.0, -11.0, 2.0, 5.0, -6.0);


    std::cout << std::endl << "A : " << std::endl;
    std::cout << std::endl << bigMatrix << std::endl;

    std::cout << std::endl << "A^{-1} : " << std::endl;
    std::cout << std::endl << bigMatrix.inverse() << std::endl;

    m::tmat<float, 2, 3> rectMat(1, 2,
                                 3, 4,
                                 5, 6);

    auto rectProd = rectMat * m::vec2(-1, -2);

    std::cout << std::endl << "B : " << std::endl;
    std::cout << std::endl << rectMat << std::endl;

    std::cout << std::endl << "B * (-1, -2) : " << std::endl;
    std::cout << std::endl << rectProd << std::endl;

    std::cout << std::endl;

    return 0;
}
