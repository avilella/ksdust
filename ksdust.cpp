#include <string>
#include <iostream>
#include <math.h>
#include <map>

// Taken from SGA preprocess dust calculation requested by Albert Vilella circa 2011
// Dust scoring scheme as given by:
// Morgulis A. "A fast and symmetric DUST implementation to Mask
// Low-Complexity DNA Sequences". J Comp Bio.
double calculateDustScore(const std::string& seq)
{
    std::map<std::string, int> scoreMap;
    
    // Cannot calculate dust scores on very short reads
    if(seq.size() < 3)
        return 0.0f;

    // Slide a 3-mer window over the sequence and insert the sequences into the map
    for(size_t i = 0; i < seq.size() - 3; ++i)
    {
        std::string triMer = seq.substr(i, 3);
        scoreMap[triMer]++;
    }

    // Calculate the score by summing the square of every element in the map
    double sum = 0;
    std::map<std::string, int>::iterator iter = scoreMap.begin();
    for(; iter != scoreMap.end(); ++iter)
    {
        int tc = iter->second;
        double score = (double)(tc * (tc - 1)) / 2.0f;
        sum += score;
    }
    return sum / (seq.size() - 2);
}

int main(int argc, char** argv)
{

    for (std::string line; std::getline(std::cin, line);) {
        double dustScore = calculateDustScore(line);
        std::cout << dustScore << std::endl;
    }

    return 0;
}
