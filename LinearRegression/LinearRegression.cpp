#include "LinearRegression.hpp"


void LinearRegression::calculateCoefficient() {
    const double N = sizeOfData();
    const double numerator = (N * m_SumXY - m_SumX * m_SumY);
    const double denominator = (N * m_SumSquareX - m_SumX * m_SumX);

    m_Coefficient = numerator / denominator;
}

void LinearRegression::calculateConstantTerm() {
    const double N = sizeOfData();
    const double numerator = (m_SumY * m_SumSquareX - m_SumX * m_SumXY);
    const double denominator = (N * m_SumSquareX - m_SumX * m_SumX);

    m_ConstTerm = numerator / denominator;
}

int LinearRegression::sizeOfData() {
    return m_X.size();
}



LinearRegression::LinearRegression(const std::vector<Symmetry::Point>& points) {
    using namespace Symmetry;

    for (const Point& point : points) {
        m_SumXY += point.x * point.y;
        m_SumX += point.x;
        m_SumY += point.y;
        m_SumSquareX += point.x * point.x;
        m_SumSquareY += point.y * point.y;
        m_X.push_back(point.x);
        m_Y.push_back(point.y);
    }
}

std::pair<double, double> LinearRegression::fit(const bool printOutput) {
    // Fitting the line to the points.
    if (m_Coefficient == 0 && m_ConstTerm == 0) {
        calculateCoefficient();
        calculateConstantTerm();
    }

    // Dumping the best fitting line.
    if (printOutput) {
        std::cout << "The best fitting line: y = " << m_Coefficient << "x + " << m_ConstTerm << std::endl;
    }

    return { m_Coefficient, m_ConstTerm };
}

void LinearRegression::showData() {
    for (int i = 0; i < 62; i++) {
        printf("_");
    }
    std::cout << "LINEAR REGRESSION INPUT DATA:" << std::endl;
    printf("\n\n");
    printf("|%15s%5s %15s%5s%20s\n", "X", "", "Y", "", "|");

    for (int i = 0; i < m_X.size(); i++) {
        printf("|%20f %20f%20s\n", m_X[i], m_Y[i], "|");
    }

    for (int i = 0; i < 62; i++) {
        printf("_");
    }

    printf("\n");
}

double LinearRegression::predict(const double x) {
    return m_Coefficient * x + m_ConstTerm;
}

double LinearRegression::errorSquare() {
    double ans = 0;
    for (int i = 0; i < m_X.size(); i++) {
        ans += ((predict(m_X[i]) - m_Y[i]) * (predict(m_X[i]) - m_Y[i]));
    }

    return ans;
}

double LinearRegression::errorIn(const double x) {
    for (int i = 0; i < m_X.size(); i++) {
        if (x == m_X[i]) {
            return (m_Y[i] - predict(m_X[i]));
        }
    }

    return 0;
}