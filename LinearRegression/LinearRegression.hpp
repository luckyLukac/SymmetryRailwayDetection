#include <iostream>
#include <stdio.h>
#include <vector>

#include "../LocalSymmetryDetector/Structs/Point.tpp"


class LinearRegression {
private:
    std::vector<double> m_X;     // Vector of X coordinates.
    std::vector<double> m_Y;     // Vector of Y coordinates.
    double m_Coefficient = 0.0;  // Predicted coefficient k in y=kx+n.
    double m_ConstTerm = 0.0;    // Predicted constant term n in y=kx+n.
    double m_SumXY = 0.0;        // Cummulative sum of X*Y values.
    double m_SumX = 0.0;         // Sum of X values.
    double m_SumY = 0.0;         // Sum of Y values.
    double m_SumSquareX = 0.0;   // Sum of squared X values.
    double m_SumSquareY = 0.0;   // Sum of squared Y values.


    /// <summary>
    /// Calculation of the best coefficent (slope) k in y=kx+n;
    /// </summary>
    void calculateCoefficient();

    /// <summary>
    /// Calculation of the best constant term of the fitting line.
    /// </summary>
    void calculateConstantTerm();

    /// <summary>
    /// Retrievement of point count, used for fitting the line.
    /// </summary>
    /// <returns>Number of points</returns>
    int sizeOfData();


public:
    /// <summary>
    /// Main constructor of LinearRegression class.
    /// </summary>
    /// <param name="points"></param>
    LinearRegression(const std::vector<Symmetry::Point>& points);


    /// <summary>
    /// Calculation of best fitting line to input points.
    /// </summary>
    /// <param name="printOutput">: printing the output to console</param>
    /// <returns>Pair of line coefficient and constant term</returns>
    std::pair<double, double> fit(const bool printOutput = false);

    /// <summary>
    /// Dump of input point data.
    /// </summary>
    void showData();

    /// <summary>
    /// Prediction of the value according to the linear function.
    /// </summary>
    /// <param name="x">: desired position by X</param>
    /// <returns>Y value at the given X</returns>
    double predict(const double x);

    /// <summary>
    /// Calculation of sum of squares of errors.
    /// </summary>
    /// <returns>Sum of squares of errors</returns>
    double errorSquare();

    /// <summary>
    /// Calculation of the difference between the actual value and model-predicted value.
    /// </summary>
    /// <param name="x">: the desired X coordinate</param>
    /// <returns>Error</returns>
    double errorIn(const double x);
};
