#pragma once
#ifndef BASIC_MATH_HPP
#define BASIC_MATH_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include "gauss_legendre.hpp"

// some fast math algorithms used in the code.
namespace basic_math
{

    // calculates gauss-legendre nodes in [-1,+1].
    std::vector<double> gauss_legendre_nodes(size_t degree)
    {
        std::vector<double> nodes;

        GaussLegendreRule gaussLegendre(degree, -1.0, 1.0);
        for (size_t i = 0; i < degree; i = i + 1)
        {
            double node_i = gaussLegendre.getNode(i);
            nodes.push_back(node_i);
        }
        return nodes;
    }

    // calculates gauss-legendre weights in [-1,+1].
    std::vector<double> gauss_legendre_weights(size_t degree)
    {
        std::vector<double> weights;

        GaussLegendreRule gaussLegendre(degree, -1.0, 1.0);
        for (size_t i = 0; i < degree; i = i + 1)
        {
            double weight_i = gaussLegendre.getWeight(i);
            weights.push_back(weight_i);
        }
        return weights;
    }

    // calculates gauss-legendre nodes in [a,b].
    std::vector<double> gauss_legendre_nodes_interval(double a, double b, size_t degree)
    {
        std::vector<double> nodes_interval;
        auto nodes = gauss_legendre_nodes(degree);
        for (size_t idx = 0; idx < nodes.size(); idx++)
        {
            double xi = nodes[idx];
            nodes_interval.push_back((b - a) / 2.0 * xi + (b + a) / 2.0);
        }

        return nodes_interval;
    }

    // calculates gauss-legendre weights in [a,b].
    std::vector<double> gauss_legendre_weights_interval(double a, double b, size_t degree)
    {
        std::vector<double> weights_interval;
        auto weights = gauss_legendre_weights(degree);
        for (size_t idx = 0; idx < weights.size(); idx++)
        {
            double wi = weights[idx];
            weights_interval.push_back((b - a) / 2.0 * wi);
        }
        return weights_interval;
    }

} // namespace basic_math

#endif // BASIC_MATH_HPP