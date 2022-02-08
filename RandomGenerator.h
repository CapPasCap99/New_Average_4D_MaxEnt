//
// Created by Grzegorz Gradziuk on 19.08.21.
//

#ifndef INC_4DCHROM_LUCAS_RANDOMGENERATOR_H
#define INC_4DCHROM_LUCAS_RANDOMGENERATOR_H

#include <random>

class RandomGenerator{
public:
    RandomGenerator(int seed, int len);
    int unidir();
    int unidir_loop();
    int unimove();
    int unimove2();
    int unipol();
    int unisite(); //now can be easily modified for a better choice of random sites (to have less uncounted moves)
    double disReal();

private:
    int len; //polymer length
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> unireal {0,1};
    std::uniform_int_distribution<int> uniint01 {0,1};
    std::uniform_int_distribution<int> uniint02 {0,2};
    std::uniform_int_distribution<int> uniint15 {1,5};
};

RandomGenerator::RandomGenerator(int seed, int len) {
    this -> len = len;
    gen.seed(seed);
}
int RandomGenerator::unidir() {
    return uniint02(gen);
}
int RandomGenerator::unidir_loop() {
    return uniint15(gen);
}
int RandomGenerator::unimove() {
    return uniint02(gen);
}
int RandomGenerator::unimove2() {
    return uniint01(gen);
}
int RandomGenerator::unipol() {
    return uniint01(gen);
}
int RandomGenerator::unisite() {
    return std::uniform_int_distribution<int>{0, len-1}(gen);
}
double RandomGenerator::disReal() {
    return unireal(gen);
}

#endif //INC_4DCHROM_LUCAS_RANDOMGENERATOR_H
