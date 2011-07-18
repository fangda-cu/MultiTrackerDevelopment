/*
 * wStrandTest.cc
 *
 *  Created on: 18/07/2011
 *      Author: jaubry
 */

#include "wStrandTest.hh"
#include "ElasticStrand.hh"
#include "ElasticStrandStaticStepper.hh"

using namespace strandsim;

int main()
{
    static const int nverts = 5;
    static const Scalar totalLength = 20.0;
    static const Scalar radius = 0.1;
    static const Scalar YoungsModulus = 1000.0;
    static const Scalar shearModulus = 100.0;
    static const Scalar density = 1.0;

    ElasticStrandParameters params(radius, YoungsModulus, shearModulus, density);
    VecXd dofs(nverts * 4 - 1);
    for (int i = 0; i < dofs.size(); i += 4)
        dofs[i] = i * 0.25 * totalLength / (nverts - 1);
    ElasticStrand strand(dofs, params);

    ElasticStrandStaticStepper stepper;

    std::cout << strand << '\n';

    for (int iteration = 0; iteration < 5; iteration++)
    {
        std::cout << "Iteration " << iteration << '\n';
        std::cout << "Forces: " << strand.getTotalForces() << '\n';
        stepper.execute(strand);
        std::cout << "New position: " << strand << '\n';
    }

    return 0;
}
