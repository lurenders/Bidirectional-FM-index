

#include "customtypedefs.h"

// ============================================================================
// HELPERS FOR BWT
// ============================================================================

Range makeRange(length_t b, length_t e) {
    Range r;
    r.begin = b;
    r.end = e;
    return r;
}

AppMatchSA makeAppMatchSA(Range range, int ED, length_t depth) {
    AppMatch match;
    match.editDist = ED;
    match.range = range;
    AppMatchSA msa;
    msa.match = match;
    msa.depth = depth;

    return msa;
}

// ============================================================================
// HELPERS FOR BIDIRECTIONAL BWT
// ============================================================================

const BiAppMatchSA makeBiAppMatchSA(SARangePair ranges, int ED,
                                    length_t depth) {
    BiAppMatchSA m;
    m.ranges = ranges;
    m.editDist = ED;
    m.depth = depth;

    return m;
}

std::string tostring(Range range) {
    return "[" + std::to_string(range.begin) + "," + std::to_string(range.end) +
           ")";
}

AppMatchSA discardReverseRange(const BiAppMatchSA& biappmatch) {
    return makeAppMatchSA(biappmatch.ranges.rangeSA, biappmatch.editDist,
                          biappmatch.depth);
}

std::vector<AppMatchSA>
discardReverseRanges(const std::vector<BiAppMatchSA>& bivector) {
    std::vector<AppMatchSA> vector;
    vector.reserve(bivector.size());

    std::transform(bivector.begin(), bivector.end(), std::back_inserter(vector),
                   discardReverseRange);
    return vector;
}

void setPieces(Search& s, const std::vector<Substring>& pieces) {
    s.pieces = pieces;

    // set the directions for the pieces
    for (length_t i = 1; i < s.pieceOrder.size(); i++) {
        int pieceNumber = s.pieceOrder[i];
        Direction d = s.directions[i - 1];
        Substring& piece = s.pieces[pieceNumber];
        piece.setDirection(d);
    }
}

Search makeSearch(std::vector<int> pieceOrder, std::vector<int> lowerBounds,
                  std::vector<int> upperBounds) {
    if (pieceOrder.size() != lowerBounds.size() ||
        pieceOrder.size() != upperBounds.size()) {
        throw std::runtime_error(
            "Could not create search, the sizes of all vectors are not equal");
    }
    Search s;
    s.upperBounds.swap(upperBounds);
    s.lowerBounds.swap(lowerBounds);
    s.pieceOrder.swap(pieceOrder);

    for (length_t i = 1; i < s.pieceOrder.size(); i++) {
        Direction d =
            (s.pieceOrder[i] > s.pieceOrder[i - 1]) ? FORWARD : BACKWARD;
        s.directions.push_back(d);
    }
    return s;
}

std::ifstream getStream(const std::string& file) {
    std::ifstream ifs(file);
    if (!ifs) {
        throw std::runtime_error("Cannot open file: " + file);
    }

    return ifs;
}

std::string readString(const std::string& s) {
    std::ifstream ifs = getStream(s);
    std::string ret;
    ifs.seekg(0, std::ios::end);
    ret.resize(ifs.tellg());
    ifs.seekg(0, std::ios::beg);
    ifs.read((char*)&ret[0], ret.size());
    ifs.close();
    return ret;
}

std::vector<length_t> readArray(size_t length, std::ifstream& ifs) {

    std::vector<length_t> sa;
    sa.reserve(length);

    length_t element;
    while (ifs >> element) {
        sa.push_back(element);
    }
    ifs.close();
    return sa;
}

void readArray(const std::string& s, size_t length, std::vector<length_t>& sa) {

    std::ifstream ifs = getStream(s);
    ifs.seekg(0, std::ios::end);
    size_t elementsInBinary = ifs.tellg() / sizeof(length_t);
    ifs.seekg(0, std::ios::beg);
    if (elementsInBinary == length) {
        // truely a binary file
        sa.resize(elementsInBinary);

        ifs.read((char*)&sa[0], sa.size() * sizeof(length_t));
        ifs.close();

    } else {
        // text file
        sa = readArray(length, ifs);
    }
}