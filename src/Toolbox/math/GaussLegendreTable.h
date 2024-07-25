#pragma once
const double GaussLegendreKnot[8][8] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {-0.5773502691896257, 0.5773502691896257, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {-0.7745966692414834, 0.0, 0.7745966692414834, 0.0, 0.0, 0.0, 0.0, 0.0},
    {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, 0.0, 0.0, 0.0, 0.0},
    {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640, 0.0, 0.0, 0.0},
    {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521, 0.0, 0.0},
    {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585, 0.0},
    {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363}
};

const double GaussLegendreCoff[8][8] = {
    {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.5555555555555556, 0.8888888888888888, 0.5555555555555556, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539, 0.0, 0.0, 0.0, 0.0},
    {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891, 0.0, 0.0, 0.0},
    {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704, 0.0, 0.0},
    {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697, 0.0},
    {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763}
};
