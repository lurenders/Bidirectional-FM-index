
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "customtypedefs.h"

void printHelp() {
    
        std::cout << "Usage: ./BWT_build <base filename>\n\n";
        std::cout << "Following files are required:\n";
        std::cout << "\t<base filename>.txt: input text T\n";
        std::cout << "\t<base filename>.sa: suffix array of T\n";
        std::cout << "\t<base filename>.rev.sa: suffix array of reverse of T\n\n";
 
}

void createBWT(const std::string& text, const std::vector<length_t>& sa,
               std::string& bwt, std::vector<length_t>& counts,
               std::vector<unsigned char>& indexToChar,
               std::vector<int>& charToIndex) {

    std::vector<length_t> charToCount(256, 0);
    bwt.reserve(sa.size());

    for (auto saIt = sa.cbegin(); saIt != sa.cend(); saIt++) {
        // get the previous character
        char c;
        if ((*saIt) != 0) {
            c = text[*saIt - 1];
        } else {
            // if the entry in the suffix array is zero, then the final
            // character is the special end character
            c = '$';
        }

        // append this char to the bwt
        bwt += c;

        // increase count of c
        charToCount[c]++;
    }
    // create the alphabet
    std::vector<int> allChars(256, -1);
    length_t charNumber = 0;
    for (int c = 0; c < 256; c++) {
        if (charToCount[(unsigned char)c] > 0) {
            allChars[(unsigned char)c] = charNumber++;
            indexToChar.push_back(c);
        }
    }

    charToIndex.swap(allChars);
    counts.swap(charToCount);
}

void createIndex(const std::string& name, const std::string& output) {
    std::cout << "Read in the textfile" << std::endl;
    std::string text = readString(name + ".txt");
    if (text[text.size() - 1] == '\n') {
        text.erase(text.end() - 1);
    }
    std::cout << "Textsize: " << text.size() << std::endl;
    std::vector<length_t> sa;
    bool dollarPresent = text[text.size() - 1] == '$';
    if (!dollarPresent) {
        throw std::runtime_error("NO dollar at end of file");
    }

    std::cout << "Read in suffix array" << std::endl;
    try {
        readArray(name + ".sa", text.size(), sa);
    } catch (const std::exception& e) {
        // retry with .1 extensions
        readArray(name + ".sa.1", text.size(), sa);
    }

    std::cout << "Size of SA: " << sa.size() << std::endl;
    if (sa.size() != text.size()) {
        throw std::runtime_error("files do not have same size!");
    }

    std::cout << "Creating alphabet and BWT" << std::endl;
    std::string bwt;
    std::vector<length_t> counts;
    std::vector<int> charToIndex;
    std::vector<unsigned char> indexToChar;
    createBWT(text, sa, bwt, counts, indexToChar, charToIndex);

    std::cout << "Alphabet contains " << indexToChar.size() << " characters"
              << std::endl;

    std::cout << "Writing counts to disk" << std::endl;
    std::ofstream ofs(output + ".cct", std::ios::binary);

    ofs.write((char*)&counts[0], counts.size() * sizeof(length_t));
    ofs.close();

    std::cout << "Writing BWT to disk" << std::endl;
    ofs = std::ofstream(output + ".bwt");
    ofs.write((char*)&bwt[0], bwt.size());
    ofs.close();

    // create sparse suffix arryas
    for (int saSF = 1; saSF <= 256; saSF *= 2) {
        std::vector<length_t> spSA((sa.size() + saSF - 1) / saSF);

        for (size_t i = 0; i < spSA.size(); i++) {
            spSA[i] = sa[i * saSF];
        }

        ofs = std::ofstream(output + ".sa." + std::to_string(saSF));
        ofs.write((char*)&spSA[0], spSA.size() * sizeof(length_t));

        ofs.close();
        std::cout << "Wrote sparse file for factor: " << saSF << std::endl;
    }

    // no need for SA anymore
    sa.clear();

    // create prefix occ table
    std::cout << "Writing prefix occurences" << std::endl;
    PrefixOccurences prefixOccurences(indexToChar.size(), charToIndex, bwt);
    prefixOccurences.write(output + ".poc");

    // do the same for the interleaveds
    PrefixOccurrences4 prefixOccurrences4(indexToChar.size(), charToIndex, bwt);
    prefixOccurrences4.write(output + ".poc.dna");
    // no need for bwt
    // bwt.clear();

    std::cout << "Read reverse sa" << std::endl;
    std::vector<length_t> revSA;
    readArray(name + ".rev.sa", text.size(), revSA);
    if (revSA.size() != text.size()) {
        throw std::runtime_error("reverse SA has different length");
    }
    std::cout << "build reverse bwt" << std::endl;
    std::string revBWT;
    revBWT.reserve(text.size());

    for (auto saIt = revSA.begin(); saIt != revSA.end(); saIt++) {
        // get the previous character
        char c;
        if ((*saIt) != 0) {
            c = text[text.size() - *saIt];
        } else {
            // if the entry in the suffix array is zero, take the first
            // character of the text (= final of rev text)
            c = text[0];
        }

        // append this char to the reversed bwt
        revBWT += c;
    }

    std::cout << "Writing reversed prefix occurences" << std::endl;
    PrefixOccurences revPrefixOccurences(indexToChar.size(), charToIndex,
                                         revBWT);

    revPrefixOccurences.write(output + ".rev.poc");
    PrefixOccurrences4 revPrefixOccurrences4(indexToChar.size(), charToIndex,
                                             revBWT);
    revPrefixOccurrences4.write(output + ".rev.poc.dna");

    // sanity check

    for (size_t c = 0; c < indexToChar.size(); c++) {
        for (size_t i = 0; i < bwt.size() + 1; i++) {
            if (prefixOccurrences4.getNumberOfPrefixOccurences(c, i) !=
                prefixOccurences.getNumberOfPrefixOccurences(c, i)) {
                std::cout << "Error at position " << i << ": "
                          << prefixOccurrences4.getNumberOfOccurences(c, i)
                          << " vs "
                          << prefixOccurences.getNumberOfPrefixOccurences(c, i)
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    for (size_t c = 0; c < indexToChar.size(); c++)
        for (size_t i = 0; i < bwt.size() + 1; i++)
            if (revPrefixOccurrences4.getNumberOfPrefixOccurences(c, i) !=
                revPrefixOccurences.getNumberOfPrefixOccurences(c, i)) {
                std::cout << "Error at position " << i << ": "
                          << revPrefixOccurrences4.getNumberOfOccurences(c, i)
                          << " vs "
                          << revPrefixOccurences.getNumberOfPrefixOccurences(c,
                                                                             i)
                          << std::endl;
                exit(EXIT_FAILURE);
            }
}

int main(int argc, char* argv[]) {

    if (argc != 2 && argc != 3) {
        std::cout << argc;
        printHelp();
        return EXIT_FAILURE;
    }

    std::string name = argv[1];
    std::string output = name;
    if (argc == 3) {
        output = argv[2];
    }

    try {
        createIndex(name, output);
    } catch (const std::exception& e) {
        std::cerr << "Error while creating index: " << e.what() << std::endl;
        printHelp();
        return EXIT_FAILURE;
    }
    std::cout << "successs" << std::endl;
    return EXIT_SUCCESS;
}
