
#include "BWT.h"

using namespace std;

// ============================================================================
// CLASS READ MAPPER
// ============================================================================

// ----------------------------------------------------------------------------
// PRIVATE ROUTINES FOR READING IN THE FILES
// ----------------------------------------------------------------------------

void BWT::fromFiles(const string& prefix) {
    cout << "Reading in files with prefix " << prefix << endl;
    // read the text
    cout << "Reading " << prefix << ".txt" << endl;
    string text = readString(prefix + ".txt");
    textLength =
        (text[text.size() - 1] == '\n') ? text.size() - 1 : text.size();

    cout << "Done reading text (size: " << text.size() << ")" << endl;

    // read the counts table
    counts.assign(256, 0);
    cout << "Reading " << prefix << ".cct" << endl;
    readArray(prefix + ".cct", 256, counts);

    vector<length_t> cumCounts;
    charToIndex = vector<int>(256, -1);
    length_t cumCount = 0; // cumulative character counts
    for (size_t i = 0; i < counts.size(); i++) {
        if (counts[i] == 0)
            continue;
        cumCounts.push_back(cumCount);
        cumCount += counts[i];
        charToIndex[i] = indexToChar.size();
        indexToChar.push_back(char(i));
    }

    counts.swap(cumCounts);
    cout << "Text contains " << indexToChar.size() << " unique characters"
         << endl;

    for (auto i : indexToChar) {
        cout << (unsigned char)i << "\t";
    }
    cout << endl;

    // read the SA
    cout << "Reading " << prefix << ".sa." << sparseFactorSA << endl;
    string saFilename = prefix + ".sa." + to_string(sparseFactorSA);
    length_t saLength = (textLength + sparseFactorSA - 1) / sparseFactorSA;
    readArray(saFilename, saLength, sa);

    cout << "Done reading suffix array (sparseness factor: " << sparseFactorSA
         << ")" << endl;

    // read the BWT
    cout << "Reading " << prefix << ".bwt" << endl;
    bwt = readString(prefix + ".bwt");

    cout << "Done reading BWT (size: " << bwt.size() << ")" << endl;

    // scan the text for word delimiting characters
    cout << "Scanning for delimiting characters" << endl;
    for (auto it = text.cbegin(); it != text.cend(); it++) {
        if (*it == DOLLAR) {
            dollarPos = it - text.cbegin();
        } else if (*it == '%' || *it == '#') {
            delimitingPos.push_back(it - text.cbegin());
        }
    }
    cout << "Finished scanning" << endl;

    // read the prefix occurrence table
    cout << "Reading " << prefix << ".poc.dna" << endl;

    if (!prefixOccurrences.read(prefix + ".poc.dna"))
        throw runtime_error("Cannot open file: " + prefix + ".poc.dna");
    cout << "Done reading prefix occurrence table" << endl;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR PATTERN MATCHING
// ----------------------------------------------------------------------------
length_t BWT::getNumberOfPrefOcc(char symbol, length_t index) const {
    return prefixOccurrences.getNumberOfPrefixOccurences(charToIndex[symbol],
                                                         index);
}

length_t BWT::getNumberOfOcc(char symbol, length_t index) const {
    return prefixOccurrences.getNumberOfOccurences(charToIndex[symbol], index);
}

length_t BWT::findLF(length_t k) const {
    return counts[charToIndex[(unsigned char)bwt[k]]] +
           getNumberOfOcc(bwt[k], k);
}

length_t BWT::findSA(length_t index) const {
    // if index modulo the sparsefactor is zero then this row is in the sparse
    // suffix array
    if ((index & (sparseFactorSA - 1)) == 0) {
        return sa[index >> logSparseFactorSA];
    }

    // else iterate over LF mappings untill an index is found that is in the
    // sparse suffix array.
    length_t j = 0;
    while ((index & (sparseFactorSA - 1)) != 0) {
        j++;
        index = findLF(index);
    }

    // return the entry of the new found index plus the number of iterations
    // modulo the lenght of the bwt
    length_t ret = sa[index >> logSparseFactorSA] + j;
    if (ret >= bwt.size()) {
        return ret - bwt.size();
    }
    return ret;
}

vector<length_t> BWT::exactMatches(const string& s) {
    // find the range in the suffix array that matches the string
    Range range = matchString(s);

    // declare the return vector
    vector<length_t> positions;
    positions.reserve(range.width());

    // fill in the vector with all values in this range in the suffix array
    for (length_t i = range.begin; i < range.end; i++) {
        positions.push_back(findSA(i));
    }

    // sort the vector and return
    sort(positions.begin(), positions.end());
    return positions;
}

Range BWT::matchString(const string& s) {
    // start at the end
    auto it = s.crbegin();

    // find the range for this intitial character in the BWT string
    length_t positionInAlphabet = charToIndex[(unsigned char)*it];

    length_t start = counts[positionInAlphabet];
    length_t end;
    if (positionInAlphabet != indexToChar.size() - 1) {
        end = counts[positionInAlphabet + 1];
    } else {
        end = bwt.size();
    }
    nodeCounter++;

    // iterate starting from the second character over the string
    for (++it; it != s.crend(); it++) {
        // find number of occurences of this char before and after and so the
        // new range is found
        length_t startOfChar = counts[charToIndex[(unsigned char)*it]];
        start = getNumberOfOcc(*it, start) + startOfChar;
        end = getNumberOfOcc(*it, end) + startOfChar;
        nodeCounter++;
        if (start == end) {
            // no matches found
            return Range();
        }
    }

    // return this range as a pair
    return Range(start, end);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<AppMatch> BWT::approxMatches(const string& pattern, int maxED) {
    resetCounters();
    vector<AppMatchSA> occurences;

    // Intialize a band-diagonal matrix. The length of the diagonal should
    // be the (P.size + maxED + 1) because (P.size + maxED) is the largest
    // possible sequence we will align. The +1 is for the extra 0th row/col.
    EditMatrix matrix(pattern.size() + maxED + 1, maxED, 0, false, false);

    // set the intial range this is the entire bwt
    Range initialRange = Range(0, bwt.size());

    recApproxMatchesNaive(initialRange, pattern, matrix, maxED, occurences, 0);

    // in the occurence vector are ranges in the suffix array that have their
    // correspondent edit distance these indexes correspond to values that are
    // the positions in the original string
    return mapOccurencesInSAToOccurencesInText(occurences, maxED, pattern);
}

vector<pair<Range, char>> BWT::getCharExtensions(const Range& range) {
    vector<pair<Range, char>> nextChars;
    nextChars.reserve(range.width());

    // iterate over the entire alphabet
    for (const auto& c : indexToChar) {

        length_t occBef = getNumberOfOcc(c, range.begin);
        length_t occAft = getNumberOfOcc(c, range.end);

        // check if this character occurs in the specified range
        if (occAft > occBef) {
            length_t s = occBef + counts[charToIndex[(unsigned char)c]];
            length_t e = occAft + counts[charToIndex[(unsigned char)c]];
            Range newRange = Range(s, e);

            // push this range and character for the next iteration
            nextChars.push_back(make_pair(newRange, c));
        }
    }

    return nextChars;
}

void BWT::recApproxMatchesNaive(const Range& range, const string& P,
                                EditMatrix& M, const int& maxED,
                                vector<AppMatchSA>& occ, length_t depth) {

    // store all character extensions ((characters before this character)
    vector<pair<Range, char>> nextChar = getCharExtensions(range);
    nodeCounter += nextChar.size();

    // report if fully alligned and also check the children
    // check the last entry in the matrix at collumn corresponding to first
    // character of P
    EditDistance ED =
        (depth + maxED >= P.size()) ? M(depth, P.size()) : maxED + 1;

    if (ED.getED() <= maxED) {
        // push if this is better than th allowed edit distance
        occ.push_back(makeAppMatchSA(range, ED.getED(), depth));
    }

    // checking the children so going one row down
    depth++;

    // Matrix range we're filling in: M[i, startj:endj]
    // index = vertical (index i) ; pattern P = horizontal (index j)
    length_t firstCol = max<int>(1, depth - maxED);
    length_t lastCol = min<int>(P.size(), depth + maxED);

    for (length_t ci = 0; ci < nextChar.size(); ci++) {
        char c = nextChar[ci].second;

        // compute the next edit distance matrix column
        unsigned int minChdED = maxED + 1;
        for (length_t j = firstCol; j <= lastCol; j++) {
            // search is in reverse order so the at operator needs to be
            // inversed
            M.updateMatrix(P[P.size() - j] != c, depth, j);
            matrixElementCounter++;
            minChdED = min(minChdED, M(depth, j).getED());
        }

        // edit distance threshold exceeded: backtracking
        if (minChdED > maxED) {
            continue;
        }
        // continue the search on this child
        recApproxMatchesNaive(nextChar[ci].first, P, M, maxED, occ, depth);
    }
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMAtE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<AppMatch>
BWT::mapOccurencesInSAToOccurencesInText(const vector<AppMatchSA>& occ,
                                         length_t ED, const Substring& p) {
    vector<AppMatch> occurencesInText;
    occurencesInText.reserve(bwt.size());
    // map the startposition to the best match (lowest edit distance)
    map<length_t, AppMatch> posToBestMatch;

    if (occ.size() == 0) {
        return {};
    }
    if (occ.size() == 1) {
        // all occ are distinct
        positionsInPostProcessingCounter = occ[0].match.range.width();
        return convertToMatchesInText(occ[0], p, ED);
    }

    // more than 1 reported occurence, could be redundant
    for (const auto& it : occ) {

        Range range = it.match.range;
        positionsInPostProcessingCounter += range.width();

        auto matchesInTextToCheck = convertToMatchesInText(it, p, ED);
        for (const auto& match : matchesInTextToCheck) {

            length_t startPos = match.range.begin;

            // check starthas already been seen (within a range of ED + 1)
            const auto& old = posToBestMatch.find(startPos);
            if (old != posToBestMatch.end()) {
                // this particular position has already been seen, check if
                // update is needed
                if (match.editDist < old->second.editDist) {
                    // better ED so update
                    posToBestMatch[startPos] = match;

                } else if (match.editDist == old->second.editDist &&
                           match.range.width() < old->second.range.width()) {
                    // same ed but better depth -> update
                    posToBestMatch[startPos] = match;
                }
                continue;
            }
            posToBestMatch[startPos] = match;
        }
    }

    // neighbours could be redundant

    // not perfect yet, example pattern is CAAAG and stored neighbours have TAAA
    // and AAAG, AAAG is better but will not push out neighbour because
    // different ends
    for (auto it = posToBestMatch.begin(); it != posToBestMatch.end(); it++) {

        const auto& entry = it->second;
        length_t i = it->first;

        // check the lower neighbours if they have same end
        length_t lowestN = (i > ED) ? i - ED : 0;

        for (length_t n = lowestN; n < i; n++) {
            auto nEntry = posToBestMatch.find(n);
            if (nEntry != posToBestMatch.end()) {
                if (entry.range.end == nEntry->second.range.end) {
                    // they have the same end so redundant, only keep the
                    // one with best ed

                    if (entry.editDist <= nEntry->second.editDist) {
                        // current ED better -> remove neighbour
                        posToBestMatch.erase(n);
                    } else {
                        // neighbour has better ED + shorter -> remove current
                        posToBestMatch.erase(i);
                        break;
                    }
                }
            }
        }
    }

    // add all the best matches to the return vector
    for (auto it = posToBestMatch.begin(); it != posToBestMatch.end(); it++) {
        occurencesInText.push_back(it->second);
    }

    return occurencesInText;
}

AppMatch BWT::verifyInText(length_t s, length_t e, int remED,
                           const Substring& pat) const {

    string text = getText();
    // the part of the text that needs to be checked
    Substring textAt = Substring(text, s, min(e, s + remED + pat.size()));

    // need to allign startOfText to string
    // make matrix
    length_t m = textAt.size();
    length_t n = pat.size();
    EditMatrix matrix(n + remED + 1, remED, 0, s != 0, true);

    length_t currentBestRow = 0;
    length_t currentBestED = remED + 1;
    // fill in the matrix
    for (length_t i = 1; i <= m; i++) {
        int bestEDRow = remED + 1; // collumns of the matrix for this row
        length_t firstCol = max<int>(1, i - remED);
        length_t lastCol = min<int>(n, i + remED);

        for (length_t j = firstCol; j <= lastCol; j++) {
            matrix.updateMatrix(textAt[i - 1] != pat[j - 1], i, j);
            bestEDRow = min<length_t>(bestEDRow, matrix(i, j).getED());
        }

        if (bestEDRow > remED) {
            // no allignment possible -> return empty match
            AppMatch m;
            m.range = Range();
            m.editDist = remED + 1;
            return m;
        }

        if (lastCol == n && matrix(i, lastCol).getED() < currentBestED) {
            currentBestED = matrix(i, lastCol).getED();
            currentBestRow = i;
        }
    }
    // make a match in the text
    AppMatch match;
    match.range = Range(s, s + currentBestRow);
    match.editDist = currentBestED;
    return match;
}

bool BWT::goesOverWordDelimiting(length_t begin, length_t end,
                                 length_t& delimitPos) const {

    if (delimitingPos.size() == 0) {
        return false;
    }
    if (delimitingPos.size() == 1) {
        if (begin <= delimitingPos[0] && end > delimitingPos[0]) {
            delimitPos = delimitingPos[0];
            return true;
        }
    }

    // iterate over delimiters (should be small and sorted vector)

    for (const auto& dpos : delimitingPos) {
        // first beyond begin
        if (dpos >= begin) {
            // check if before end
            if (dpos < end) {
                delimitPos = dpos;
                return true;
            } else {
                // past begin and past end
                return false;
            }
        }
    }

    // no  worddelimiter in this range
    return false;
}

vector<AppMatch> BWT::convertToMatchesInText(const AppMatchSA& saMatch,
                                             const Substring& pat, int ED) {

    vector<AppMatch> textMatches;
    textMatches.reserve(saMatch.match.range.width());

    for (length_t i = saMatch.match.range.begin; i < saMatch.match.range.end;
         i++) {
        // find the startPosition in the text by looking at the SA
        length_t startPos = findSA(i);

        length_t endPos = startPos + saMatch.depth;

        AppMatch newMatch;
        newMatch.range = Range();

        length_t delimitP = 0;
        // check if goes over dollar
        if (endPos >= textLength) {
            length_t befDollar = textLength - startPos;
            // check after dollar
            if (saMatch.depth - befDollar >= pat.size() - ED) {
                newMatch = verifyInText(0, saMatch.depth - befDollar, ED, pat);
            }
            // check before dollar
            if (befDollar >= pat.size() - ED) {
                newMatch = verifyInText(startPos, textLength, ED, pat);
            }

        } else if (goesOverWordDelimiting(startPos, endPos, delimitP)) {
            // check if the piece before the wordDelimiting could form a
            // match
            if (delimitP - startPos >= pat.size() - ED) {
                newMatch = verifyInText(startPos, delimitP, ED, pat);

            } else if (endPos - (delimitP + 1) >= pat.size() - ED) {
                newMatch = verifyInText(delimitP + 1, endPos, ED, pat);
            }

        } else {
            // does not go over sentinel, create simple match
            newMatch.range = Range(startPos, endPos);
            newMatch.editDist = saMatch.match.editDist;
        }

        if (!newMatch.range.empty()) {
            textMatches.push_back(newMatch);
        }
    }
    return textMatches;
}