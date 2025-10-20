//
// Created by Grzegorz Gradziuk on 19.08.21.
//

#ifndef INC_4DCHROM_LUCAS_RANDOMGENERATOR_H
#define INC_4DCHROM_LUCAS_RANDOMGENERATOR_H

#include <boost/random.hpp>

class RandomGenerator{
public:
    RandomGenerator(int seed, int len, int lin_len, int OriC);
    int unidir();
    int unidir_loop();
    int unimove();
    int unimove2();
    int weightedpol();
    int unisitering(); //now can be easily modified for a better choice of random sites (to have less uncounted moves)
    int unisitelin();
    double disReal();

private:
    int len; //polymer length
    int linlen; //linear polymer is length 2*linlen+1
    int OriC;
    boost::random::mt19937_64 gen;
    boost::random::uniform_01<double> unireal;
    boost::uniform_int<int> uniint01{0,1};
    boost::uniform_int<int> uniint02{0,2};
    boost::uniform_int<int> uniint15{1,5};
};

RandomGenerator::RandomGenerator(int seed, int len, int lin_len, int OriC) {
    this -> len = len;
    this -> linlen=lin_len;
    this -> OriC = OriC;
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
int RandomGenerator::weightedpol() {
    return boost::random::discrete_distribution<int,int>{len,2*linlen}(gen);
}
int RandomGenerator::unisitering() {
    return std::uniform_int_distribution<int>{0, len-1}(gen);
}
// int RandomGenerator::unisitelin() {
//     int raw_site = std::uniform_int_distribution<int>{OriC - linlen, OriC + linlen - 1}(gen);
//     return (raw_site + len) % len;  // ensures wrap-around into valid [0, len-1], if you change OriC from 1220
// }

int RandomGenerator::unisitelin() {
    return (std::uniform_int_distribution<int>{OriC-linlen,OriC+linlen-1}(gen))%len;
}

double RandomGenerator::disReal() {
    return unireal(gen);
}

#endif //INC_4DCHROM_LUCAS_RANDOMGENERATOR_H


