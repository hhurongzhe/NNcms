#ifndef GAUSSLEGENDRERULE_HPP
#define GAUSSLEGENDRERULE_HPP

namespace GaussLegendre
{
    /**
     * Purpose:
     *	This class defines a Gauss-Legendre "point". Namely, a Gauss-Legendre point is defined to contain an array of nodes and weights and an order (the size of the nodes/weights array's).
     */
    class GaussLegendrePoint
    {
    public:
        GaussLegendrePoint(int order, double *node, double *weight);
        GaussLegendrePoint(int order);
        GaussLegendrePoint();

        double getNode(int i);
        double getWeight(int i);
        int getOrder();

        void setNode(int i, double value);
        void setWeight(int i, double value);

        ~GaussLegendrePoint();

    private:
        double *node;
        double *weight;
        int order;

        bool isNodeNUll();
        bool isWeightNUll();
        void init(int order, double *node, double *weight);
        void validateNodeAndWeight(int i);
        void setNodePtr(double *node);
        void setWeightPtr(double *weight);
        void setOrder(int order);
    };
}

/**
 * Purpose:
 *	Calculates an n-point Gauss-Legendre rule over the interval (a,b).
 */
class GaussLegendreRule
{
public:
    GaussLegendreRule(int n, double a, double b);

    double getNode(int i);
    double getWeight(int i);
    double getOrder();
    double getLowerBound();
    double getUpperBound();

private:
    int n;
    double a;
    double b;
    GaussLegendre::GaussLegendrePoint gaussLegendrePoint;

    void setGaussLegendrePoint();
    void mapGaussLegendrePoint(GaussLegendre::GaussLegendrePoint &negativeZeroOne);
};

#endif