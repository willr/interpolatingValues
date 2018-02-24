#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>


std::vector<std::tuple<uint64_t, double, double, double>>
imuStorage()
{
    std::vector<std::tuple<uint64_t, double, double, double>> imuVals {
            std::make_tuple(7090906, 0.24492998421192200, 0.19367285072803500, -10.58722114562990),
            std::make_tuple(7190913, 0.19180884957313500, -0.48705899715423600, -9.60891246795654),
            std::make_tuple(7290907, 0.00410881638526917, -0.24062308669090300, -8.07599163055420),
            std::make_tuple(7391908, -0.85458230972290000, -0.51548421382904100, -8.20866394042969),
            std::make_tuple(7490907, -0.53761631250381500, 0.05440217256546020, -9.01467323303223),
            std::make_tuple(7590907, 0.01358174812048670, 0.20375525951385500, -9.27913379669189),
            std::make_tuple(7690907, 0.02056398615241050, -0.35768926143646200, -9.14729595184326),
            std::make_tuple(7790908, -0.24832674860954300, 0.03896053135395050, -9.25193214416504),
            std::make_tuple(7890913, -0.28455865383148200, -0.13461312651634200, -9.39185905456543),
            std::make_tuple(7990907, -0.18781501054763800, -0.40911358594894400, -9.21934127807617),
            std::make_tuple(8091938, -0.24034203588962600, -0.29952114820480300, -9.37930870056152),
            std::make_tuple(8190908, 0.03914649039506910, -0.21077662706375100, -9.17698478698731),
            std::make_tuple(8290907, -0.17362956702709200, -0.13741925358772300, -9.20533084869385),
            std::make_tuple(8390907, -0.13457714021205900, -0.39932504296302800, -9.21733283996582)

    };
    return imuVals;
}

std::vector<std::tuple<uint64_t, double, double, double>>
gpsStorage()
{
    std::vector<std::tuple<uint64_t, double, double, double>> gpsValues{
            std::make_tuple(7009129, 0, 0, 0),
            std::make_tuple(7042462, 0, 0, 0),
            std::make_tuple(7075795, 0, 0, 0),
            std::make_tuple(7109128, 0, 0, 0),
            std::make_tuple(7142461, 0, 0, 0),
            std::make_tuple(7175794, 0, 0, 0),
            std::make_tuple(7206967, 0, 0, 0),
            std::make_tuple(7240300, 0, 0, 0),
            std::make_tuple(7273633, 0, 0, 0),
            std::make_tuple(7306966, 0, 0, 0),
            std::make_tuple(7340299, 0, 0, 0),
            std::make_tuple(7373632, 0, 0, 0),
            std::make_tuple(7402809, 0, 0, 0),
            std::make_tuple(7436142, 0, 0, 0),
            std::make_tuple(7469475, 0, 0, 0),
            std::make_tuple(7502808, 0, 0, 0),
            std::make_tuple(7536141, 0, 0, 0),
            std::make_tuple(7569474, 0, 0, 0),
            std::make_tuple(7599573, 0, 0, 0),
            std::make_tuple(7632906, 0, 0, 0),
            std::make_tuple(7666239, 0, 0, 0),
            std::make_tuple(7699572, 0, 0, 0),
            std::make_tuple(7732905, 0, 0, 0),
            std::make_tuple(7766238, 0, 0, 0),
            std::make_tuple(7798438, 0, 0, 0),
            std::make_tuple(7831771, 0, 0, 0),
            std::make_tuple(7865104, 0, 0, 0),
            std::make_tuple(7898437, 0, 0, 0),
            std::make_tuple(7931770, 0, 0, 0),
            std::make_tuple(7965103, 0, 0, 0),
            std::make_tuple(7988182, 0, 0, 0),
    };

    return gpsValues;
}

template<typename T>
uint64_t
getSensorAverage(std::vector<T> list)
{
    const int AverageLimit = 10;
    int counter = 0;
    uint64_t last = 0;
    uint64_t accum = 0;
    for (auto sensor : list)
    {
        auto time = std::get<0>(sensor);

        if (last != 0)
        {
            accum += time - last;
        }
        last = time;

        counter++;
        if (counter > AverageLimit)
        {
            break;
        }
    }

    return accum / (counter - 1);
}

double
slope(double beginVal, double endVal, uint64_t beginTime, uint64_t endTime)
{
    // m = (y2 - y1) / (x2 - x1)
    // y = m x + b
    // y = ((y2 - y1) / (x2 - x1)) * x + b
    // b = (-1 * ((y2 - y1) / (x2 - x1)) * x) + y

    double x1 = beginTime;
    double y1 = beginVal;
    double x2 = endTime;
    double y2 = endVal;

//    double b = (-1 * (((y2 - y1) / (x2 - x1)) * x1)) + y1;
    double m = (y2 - y1) / (x2 - x1);
    double b = (-1 * (m * x1)) + y1;
//    double y = ((y2 - y1) / (x2 - x1)) * xInc + b;
//    std::cout << "imuStep: " << std::fixed << xInc << ", " << std::fixed << y << std::endl;

    return b;
};

double
step(double beginVal, double endVal, uint64_t beginTime, uint64_t endTime, double b, int xInc)
{
    // m = (y2 - y1) / (x2 - x1)
    // y = m x + b
    // y = ((y2 - y1) / (x2 - x1)) * x + b
    // b = (-1 * ((y2 - y1) / (x2 - x1)) * x) + y

    double x1 = beginTime;
    double y1 = beginVal;
    double x2 = endTime;
    double y2 = endVal;

//    double b = (-1 * (((y2 - y1) / (x2 - x1)) * x1)) + y1;
    double y = ((y2 - y1) / (x2 - x1)) * xInc + b;
//    std::cout << "imuStep: " << std::fixed << xInc << ", " << std::fixed << y << std::endl;

    return y;
};

double
step(double beginVal, double endVal, uint64_t beginTime, uint64_t endTime, int xInc)
{
    // m = (y2 - y1) / (x2 - x1)
    // y = m x + b
    // y = ((y2 - y1) / (x2 - x1)) * x + b
    // b = (-1 * ((y2 - y1) / (x2 - x1)) * x) + y

    double x1 = beginTime;
    double y1 = beginVal;
    double x2 = endTime;
    double y2 = endVal;

//    double b = (-1 * (((y2 - y1) / (x2 - x1)) * x1)) + y1;
    double m = (y2 - y1) / (x2 - x1);
    double b = (-1 * (m * x1)) + y1;

//    double y = ((y2 - y1) / (x2 - x1)) * xInc + b;
    double y = (m) * xInc + b;

//    std::cout << "imuStep: " << std::fixed << xInc << ", " << std::fixed << y << std::endl;

    return y;
};

std::vector<std::tuple<uint64_t, double, double, double>>
getImuValuesWithTimePrint(const std::vector<std::tuple<uint64_t, double, double, double>> &imus, int *position, uint64_t endTime)
{
    std::vector<std::tuple<uint64_t, double, double, double>> results{};

    int pos = *position;
    auto imuTime = std::get<0>(imus[pos]);
    while (imuTime <= endTime)
    {
//        std::cout << "imu: " << std::fixed << imus[pos].first << " , " << imus[pos].second << std::endl;
        results.push_back(std::make_tuple(std::get<0>(imus[pos]), std::get<1>(imus[pos]), std::get<2>(imus[pos]), std::get<3>(imus[pos])));

        pos++;
        *position = pos;
        imuTime = std::get<0>(imus[pos]);
    }

    return results;
}

void imuFromGps1(const std::vector<std::tuple<uint64_t, double, double, double>> &imuVals,
                 const std::vector<std::tuple<uint64_t, double, double, double>> &gpsValues)
{
    int imuPos = 0;

    for (auto &gpsValue : gpsValues)
    {
        auto gpsTime = std::get<0>(gpsValue);
        auto imuEntries = getImuValuesWithTimePrint(imuVals, &imuPos, gpsTime);
        for (auto &imuEntry : imuEntries)
        {
            std::cout << "imu: " << std::fixed << std::get<0>(imuEntry) << " , " << std::get<1>(imuEntry) << " , " << std::get<2>(imuEntry) << " , " << std::get<3>(imuEntry) << std::endl;
        }
        std::cout << "gps time: " << std::fixed << gpsTime << std::endl;

        int i = 0;
    }
    int i = 0;
}

std::vector<std::tuple<uint64_t, double, double, double>>
getImuValuesBetweenTimeBand(const std::vector<std::tuple<uint64_t, double, double, double>> &imus,
                            int *position,
                            uint64_t startTime,
                            uint64_t endTime,
                            uint64_t timeIncrement)
{
    std::vector<std::tuple<uint64_t, double, double, double>> results{};

    int startPos = *position;
    int nextPos = startPos + 1;

    auto startImu = imus[startPos];
    auto nextImu = imus[nextPos];

    // find the start of the generated imu within the timeband
    auto currentImu = startImu;
    auto currentTimestamp = std::get<0>(startImu);
    auto nextTimestamp = std::get<0>(nextImu);
    while (currentTimestamp < startTime)
    {
        currentTimestamp += timeIncrement;

        if (currentTimestamp > nextTimestamp)
        {
            *position += startPos;

            startPos = *position;
            nextPos = startPos + 1;

            startImu = imus[startPos];
            nextImu = imus[nextPos];

            currentTimestamp = std::get<0>(startImu);
            nextTimestamp = std::get<0>(nextImu);
        }
    }

    // now we are at least at the same or greater position of the startTime
    auto x_x1 = std::get<0>(startImu);
    auto x_y1 = std::get<1>(startImu);
    auto x_x2 = std::get<0>(nextImu);
    auto x_y2 = std::get<1>(nextImu);
//    auto x_b = slope(x_x1, x_x2, x_y1, x_y2);

    auto y_x1 = std::get<0>(startImu);
    auto y_y1 = std::get<2>(startImu);
    auto y_x2 = std::get<0>(nextImu);
    auto y_y2 = std::get<2>(nextImu);
//    auto y_b = slope(x_x1, x_x2, x_y1, x_y2);

    auto z_x1 = std::get<0>(startImu);
    auto z_y1 = std::get<3>(startImu);
    auto z_x2 = std::get<0>(nextImu);
    auto z_y2 = std::get<3>(nextImu);
//    auto z_b = slope(x_x1, x_x2, x_y1, x_y2);

    while (currentTimestamp < endTime)
    {
        // now do x pos
        auto currentXval = step(x_y1, x_y2, x_x1, x_x2, currentTimestamp);
//        auto currentXval = step(x_y1, x_y2, x_x1, x_x2, x_b, currentTimestamp);

        // now do y pos
        auto currentYval = step(y_y1, y_y2, y_x1, y_x2, currentTimestamp);
//        auto currentYval = step(y_y1, y_y2, y_x1, y_x2, y_b, currentTimestamp);

        // now do z pos
        auto currentZval = step(z_y1, z_y2, z_x1, z_x2, currentTimestamp);
//        auto currentZval = step(z_y1, z_y2, z_x1, z_x2, z_b, currentTimestamp);

        results.push_back(std::make_tuple(currentTimestamp,currentXval, currentYval, currentZval));
        currentTimestamp += timeIncrement;

        if (currentTimestamp > nextTimestamp)
        {
            *position = startPos + 1;

            startPos = *position;
            nextPos = startPos + 1;

            startImu = imus[startPos];
            nextImu = imus[nextPos];

            currentTimestamp = std::get<0>(startImu);
            nextTimestamp = std::get<0>(nextImu);
        }
    }

    return results;
}

std::vector<std::tuple<uint64_t, double, double, double>>
imuFromGps2(const std::vector<std::tuple<uint64_t, double, double, double>> &imuVals,
                 const std::vector<std::tuple<uint64_t, double, double, double>> &gpsValues,
                 uint64_t timeIncrement)
{
    int imuPos = 0;

    std::vector<std::tuple<uint64_t, double, double, double>> results{};

    std::tuple<uint64_t, double, double, double> lastGps;
    bool firstSeen = false;
    for (auto &gpsValue : gpsValues)
    {
        if (! firstSeen)
        {
            lastGps = gpsValue;
            firstSeen = true;
            continue;
        }
        auto beginGpsTime = std::get<0>(lastGps);
        auto endGpsTime = std::get<0>(gpsValue);
        auto imuEntries = getImuValuesBetweenTimeBand(imuVals, &imuPos, beginGpsTime, endGpsTime, timeIncrement);

        lastGps = gpsValue;

        std::cout << "gps time >>> : " << std::fixed << beginGpsTime << std::endl;
        for (auto &imuEntry : imuEntries)
        {
            std::cout << "imu: " << std::fixed << std::get<0>(imuEntry) << " , " << std::get<1>(imuEntry) << " , " << std::get<2>(imuEntry) << " , " << std::get<3>(imuEntry) << std::endl;
            results.push_back(imuEntry);
        }
        std::cout << "gps time <<< : " << std::fixed << endGpsTime << std::endl;

        int i = 0;
    }
    int i = 0;

    return results;
}

int main()
{
    std::cout << "Hello, World!" << std::endl;

    double x1 = -3;
    double y1 = 3;
    double x2 = 3;
    double y2 = -1;

//    step(y1, y2, x1, x2, x2);

    auto imuVals = imuStorage();
    auto gpsValues = gpsStorage();

    std::cout.precision(17);

//    uint64_t imuAvg = getSensorAverage<std::tuple<uint64_t, double, double, double>>(imuVals);
//    uint64_t gpsAvg = getSensorAverage<std::tuple<uint64_t, double, double, double>>(gpsValues);
//
//    double gps2imuRate = (double)gpsAvg / (double)imuAvg;
//    int numParts = 0;
//    if (gps2imuRate < 4)
//    {
//        numParts = 4 / gps2imuRate;
//    }

    int numParts = 3;
    uint64_t gpsAvg = getSensorAverage<std::tuple<uint64_t, double, double, double>>(gpsValues);

    std::vector<std::tuple<uint64_t, double, double, double>> results;
//    imuFromGps1(imuVals, gpsValues);

    uint64_t timeIncrement = gpsAvg / numParts;
    results = imuFromGps2(imuVals, gpsValues, timeIncrement);

    std::cout << "results: " << results.size() << std::endl;
    std::ofstream outfile ("imu.csv", std::ofstream::binary);
    outfile.precision(17);

    outfile << "timestamp, " <<
            "x, " <<
            "y, " <<
            "z "  << std::endl;
    for (auto r : results)
    {
        auto imuResult = r;
        outfile << std::fixed <<
                   std::get<0>(imuResult) << ", " <<
                   std::get<1>(imuResult) << ", " <<
                   std::get<2>(imuResult) << ", " <<
                   std::get<3>(imuResult) << std::endl;
    }

    return 0;
}
