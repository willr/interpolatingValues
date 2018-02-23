#include <iostream>
#include <vector>

template<typename T>
uint64_t
getSensorAverage(std::vector<T> list)
{
    int counter = 0;
    uint64_t last = 0;
    uint64_t accum = 0;
    for (auto sensor : list)
    {
        auto time = sensor.first;

        if (last != 0)
        {
            accum += time - last;
        }
        last = time;

        counter++;
    }

    return accum / (counter - 1);
}

double
step(double val, double lastVal, uint64_t timestamp, uint64_t lastTime, double increment, int incrementPos)
{
    double result = (((val - lastVal) / (timestamp - lastTime)) * ((lastTime + (increment * incrementPos)) - lastTime)) + val;

    return result;
}

std::vector<double>
getImuValues(std::vector<std::pair<uint64_t, double>> imus, const int *position, uint64_t endTime, int numParts)
{
    std::vector<double> results{};

    double increment = 1.0 / (double)numParts;
    bool first = true;
    double lastVal = 0;
    uint64_t lastTime = 0;
    for (int i = *position; i < imus.size(); i++)
    {
        auto imuVal = imus[i];
        uint64_t timestamp = imuVal.first;
        double val = imuVal.second;

        if (timestamp > endTime)
        {
            break;
        }

        for (int k = 0; k < numParts; k++)
        {
            double a = step(lastVal, val, lastTime, timestamp, increment, k);
            double result = a + val;
            results.push_back(result);
        }
    }

    return results;
}

int main()
{
    std::cout << "Hello, World!" << std::endl;

    std::vector<std::pair<uint64_t, double>> imuVals{{7090906, 0.244929984211922},
                                                {7190913, 0.191808849573135},
                                                {7290907, 0.00410881638526917},
                                                {7391908, -0.8545823097229},
                                                {7490907, -0.537616312503815},
                                                {7590907, 0.0135817481204867},
                                                {7690907, 0.0205639861524105},
                                                {7790908, -0.248326748609543},
                                                {7890913, -0.284558653831482},
                                                {7990907, -0.187815010547638},
                                                {8091938, -0.240342035889626}};

    std::vector<std::pair<uint64_t, double>> gpsValues{
            {7009129, 0},
            {7042462, 0},
            {7075795, 0},
            {7109128, 0},
            {7142461, 0},
            {7175794, 0},
            {7206967, 0},
            {7240300, 0},
            {7273633, 0},
            {7306966, 0},
            {7340299, 0},
            {7373632, 0},
            {7402809, 0},
            {7436142, 0},
            {7469475, 0},
            {7502808, 0},
            {7536141, 0},
            {7569474, 0},
            {7599573, 0},
            {7632906, 0},
            {7666239, 0},
            {7699572, 0},
            {7732905, 0},
            {7766238, 0},
            {7798438, 0},
            {7831771, 0},
            {7865104, 0},
            {7898437, 0},
            {7931770, 0},
            {7965103, 0},
            {7988182, 0}
    };

    uint64_t imuAvg = getSensorAverage<std::pair<uint64_t, double>>(imuVals);
    uint64_t gpsAvg = getSensorAverage<std::pair<uint64_t, double>>(gpsValues);;

    double gps2imuRate = (double)gpsAvg / (double)imuAvg;
    int numParts = 0;
    if (gps2imuRate < 4)
    {
        numParts = 4 / gps2imuRate;
    }

    int imuPos = 0;
    std::vector<double> results;
    for(int j = 0; j < gpsValues.size(); j++)
    {
        auto gpsTime = gpsValues[j].first;
        results = getImuValues(imuVals, &imuPos, gpsTime, numParts);
    }

    int i = 0;

    return 0;
}
