#include "BiBWT.h"

using namespace std;

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------

length_t BiBWT::findLF(length_t k, bool reversed) const {
    if (!reversed) {

        return counts[charToIndex[(unsigned char)bwt[k]]] +
               getNumberOfOcc(charToIndex[(unsigned char)bwt[k]], k);
    }
    length_t posInAlphabet = charToIndex[(unsigned char)bwt[k]];

    return counts[posInAlphabet] + getNumberOfOccRev(posInAlphabet, k);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR INITIALIZATION
// ----------------------------------------------------------------------------

SARangePair BiBWT::matchStringBidirectionally(const Substring& pattern,
                                              SARangePair rangesOfPrev) {

    for (length_t i = 0; i < pattern.size(); i++) {

        char c = pattern[i];
        int posInAlphabet = charToIndex[(unsigned char)c];
        if (posInAlphabet > -1) {

            if (!(this->*extraChar)(posInAlphabet, rangesOfPrev,
                                    rangesOfPrev)) {
                return SARangePair();
            }
            // each character that we look at is a new node
            nodeCounter++;
        } else {

            return SARangePair();
        }
    }

    return rangesOfPrev;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<AppMatch> BiBWT::approxMatches(const string& pattern, int maxED,
                                      bool forward) {

    resetCounters();
    vector<BiAppMatchSA> occurences;
    setDirection((forward) ? FORWARD : BACKWARD);

    // Intialize a band-diagonal matrix. The length of the diagonal should
    // be the (P.size + maxED + 1) because (P.size + maxED) is the largest
    // possible sequence we will align. The +1 is for the extra 0th row/col.

    EditMatrix matrix(pattern.size() + maxED + 1, maxED, 0, false, false);

    Substring pat = Substring(pattern);
    pat.setDirection(dir);

    improvedRecApproxMatches(getCompleteRange(), pat, matrix, maxED + 1,
                             occurences, 0);

    // in the occurence vector are ranges in the suffix array that have their
    // correspondent edit distance these indexes correspond to values that are
    // the positions in the original string
    return mapOccurencesInSAToOccurencesInText(occurences, maxED, pattern);
}

vector<AppMatch> BiBWT::mapOccurencesInSAToOccurencesInText(
    const vector<BiAppMatchSA>& occurences, int ED, const Substring& pattern) {

    vector<AppMatchSA> newOccurences = discardReverseRanges(occurences);
    return BWT::mapOccurencesInSAToOccurencesInText(newOccurences, ED, pattern);
}

void BiBWT::improvedRecApproxMatches(const SARangePair& ranges,
                                     const Substring& P, EditMatrix& M,
                                     const EditDistance& parentED,
                                     vector<BiAppMatchSA>& occ,
                                     length_t depth) {
    const length_t W = M.getWidth(); // W = maximum edit distance

    // store all character extensions ((characters before this character)
    vector<pair<SARangePair, char>> nextChar = getCharExtensions(ranges);
    nodeCounter += nextChar.size();

    // check the last entry in the matrix at column corresponding to first
    // character of P
    EditDistance ED = (depth + W >= P.size()) ? M(depth, P.size()) : W + 1;

    // Matrix range we're filling in for each character: M[i, firstCol:lastCol]
    // bwt = vertical (index i) ; pattern P = horizontal (index j for columns)
    // children are at row depth + 1
    length_t row = depth + 1;
    length_t rowMinW = (W >= row) ? 0 : row - W;
    length_t firstCol = max<int>(1, rowMinW);
    length_t lastCol = min<int>(P.size(), row + W);

    for (size_t ci = 0; ci < nextChar.size(); ci++) {
        char c = nextChar[ci].second;

        unsigned int minChdED = W + 1;
        unsigned int prevMin = W + 1;
        // check if the row is completely flagged, if one is not flagged then
        // row is not flagged
        bool flaggedRow = true;
        for (length_t j = firstCol; j <= lastCol; j++) {
            prevMin = minChdED;
            M.updateMatrix(P[j - 1] != c, row, j);
            matrixElementCounter++;

            minChdED = min(prevMin, M(row, j).getED());
            flaggedRow &= M(row, j).isFlagged();
        }

        if (minChdED > W || flaggedRow) {
            continue;
        }
        if (lastCol != P.size()) {
            // P is not fully aligned yet continue searching deeper
            improvedRecApproxMatches(nextChar[ci].first, P, M, ED, occ, row);
        } else {
            // P is fully alligned
            EditDistance chdED = M(row, P.size());

            if (lastCol != firstCol && chdED.getED() >= prevMin) {
                // if not the last row (no left neigbour) or all left
                // neighbours greater
                improvedRecApproxMatches(nextChar[ci].first, P, M, ED, occ,
                                         row);
            } else {
                // this is a non-redundant branch, check if it is good
                // enough
                if (M(row, P.size()).getED() <= parentED.getED()) {
                    // if this cell is better than the parent report
                    occ.push_back(makeBiAppMatchSA(nextChar[ci].first,
                                                   chdED.getED(), row));
                }
            }
        }
    }

    if (ED.getED() <= W && ED.getED() <= parentED.getED() &&
        (depth + W >= P.size()) &&
        M(depth, P.size() - 1).getED() >= ED.getED()) {
        // push if this is better than the parent and this is not redundant
        // with ending gap in pattern
        occ.push_back(makeBiAppMatchSA(ranges, ED.getED(), depth));
    }
}

vector<pair<SARangePair, char>>
BiBWT::getCharExtensions(const SARangePair& rangesOfParent) const {
    vector<pair<SARangePair, char>> nextChars;

    nextChars.reserve(indexToChar.size());

    // iterate over the entire alphabet
    for (length_t i = 0; i < indexToChar.size(); i++) {

        SARangePair pairForNewChar;

        // check if this character occurs in the specified range
        if ((this->*extraChar)(i, rangesOfParent, pairForNewChar)) {
            // push this range and character for the next iteration
            nextChars.emplace_back(make_pair(pairForNewChar, indexToChar[i]));
        }
    }

    return nextChars;
}

bool BiBWT::findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                            const SARangePair& rangesOfP,
                                            SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    Range trivialRange = rangesOfP.rangeSA;

    // find the new range by using the LF property
    length_t occBefore = getNumberOfOcc(positionInAlphabet, trivialRange.begin);
    length_t occAfter = getNumberOfOcc(positionInAlphabet, trivialRange.end);

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the sizes of the ranges
    // of (dP) using B

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.rangeSARev.begin;

    // get the start of the child within this range
    // find the number of occurences of chars smaller than c in the parent
    // range

    length_t x = getNumberOfPrefOcc(positionInAlphabet, trivialRange.end) -
                 getNumberOfPrefOcc(positionInAlphabet, trivialRange.begin);

    // make the new range with width equal to that of the trivial range
    Range range2 = Range(s + x, s + x + range1.width());

    rangesOfChild = SARangePair(range1, range2);
    return !rangesOfChild.empty();
}

bool BiBWT::findRangesWithExtraCharForward(length_t positionInAlphabet,
                                           const SARangePair& rangesOfP,
                                           SARangePair& childRanges) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    Range rangeForTrivial = rangesOfP.rangeSARev;

    // find the new range by using the LF property
    length_t occBefore =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.begin);
    length_t occAfter =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.end);

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.rangeSA.begin;

    // get the start of the child within this range
    // find the number of occurences of chars smaller than c in the parent
    // range
    Range prevRange = rangesOfP.rangeSARev;

    length_t x = getNumberOfPrefOccRev(positionInAlphabet, prevRange.end) -
                 getNumberOfPrefOccRev(positionInAlphabet, prevRange.begin);

    // make the new range
    Range range2 = Range(s + x, s + x + range1.width());

    childRanges = SARangePair(range2, range1);
    return !childRanges.empty();
}

void BiBWT::iteApproxMatch(const Substring& pattern, vector<BiAppMatchSA>& occ,
                           const BiAppMatchSA& match, length_t maxED,
                           bool leadingGapsAllowed, int maxEDNext) {

    length_t W = maxED - match.editDist;
    reserveSize(pattern.size());

    BandMatrix matrix(pattern.size() + W + 1, W, match.editDist);

    // get the initial nodes to check
    pushChildren(match.ranges);
    // the size of the pattern, in order not to recalculate this value every
    // iteration
    unsigned int patternSize = pattern.size();

    // the row of the previous parent that was checked, initialize on zero
    int lastRow = 0;

    while (!nodesToCheck.empty()) {
        // get the last node of the vector
        BidirectionalNode& currentNode = nodesToCheck.back();
        int row = currentNode.getRow();

        if (row >= lastRow) {
            // in case this is not the result of backtracking (= going a row
            // up) it is a child, calculate the value for this child and
            // update the matrix, also check the children

            // collumns of the matrix for this row
            length_t firstCol = max<int>(1, row - W);
            length_t lastCol = min<int>(patternSize, row + W);

            unsigned int minimalEDOfRow = maxED + 1;
            for (length_t j = firstCol; j <= lastCol; j++) {

                matrix.updateMatrix(
                    pattern[j - 1] != currentNode.getCharacter(), row, j);
                matrixElementCounter++;
                minimalEDOfRow = min(minimalEDOfRow, matrix(row, j));
            }

            // edit distance threshold exceeded: backtracking
            if (minimalEDOfRow > maxED) {
                nodesToCheck.pop_back();

            } else {
                // set the found EDEDForThis for this node if the pattern
                // has been completely alligned, else set max allowed ED + 1
                length_t EDForThisNode = (lastCol == patternSize)
                                             ? matrix(row, patternSize)
                                             : maxED + 1;
                currentNode.setFoundED(EDForThisNode);

                if (firstCol == patternSize) {
                    // will not go to children
                    lastRow = row + 1;
                    continue;
                }

                if (minimalEDOfRow == maxED && lastCol != patternSize) {
                    // all children should have exact matches

                    length_t start = max<int>(0, row - W);
                    for (length_t j = start; j <= lastCol; j++) {

                        if (matrix(row, j) == minimalEDOfRow) {

                            const Substring& remainingPiece =
                                pattern.getSubPiece(j);

                            SARangePair pair = matchStringBidirectionally(
                                remainingPiece, currentNode.getRanges());

                            if (!pair.empty()) {

                                const BiAppMatchSA newMatch = makeBiAppMatchSA(
                                    pair, minimalEDOfRow,
                                    row + remainingPiece.length() +
                                        match.depth);
                                occ.emplace_back(newMatch);
                            }
                        }
                    }
                    // as if this is backtrack
                    lastRow = row + 1;
                    continue;
                }

                // continue the search for children of this node
                pushChildren(currentNode);
            }

        } else {
            // this is the result of backtracking, so we can try and search
            // for cluster centrums and only report these centrums
            if (row + W >= patternSize) {
                // if fully alligned
                length_t cED = currentNode.getFoundED();
                int pIdx = currentNode.parentIdx;

                if (maxEDNext == 0) {

                    currentNode.report(occ, match.depth, maxED);

                } else if (pIdx == -1) {
                    if (currentNode.oneChildWorseOrEqual) {
                        // no parent and one child was worse
                        currentNode.report(occ, match.depth, maxED);
                    }
                } else {
                    // cluster time
                    BidirectionalNode& p = nodesToCheck[pIdx];

                    // parent has a child that got here
                    p.hasChildInThisStage |= true;
                    p.oneChildWorseOrEqual |= cED >= p.getFoundED();

                    if (p.getFoundED() >= cED) {
                        // current could be a centre
                        if (!currentNode.hasChildInThisStage ||
                            currentNode.oneChildWorseOrEqual) {
                            currentNode.report(occ, match.depth, maxED);
                        }
                    }
                }
            }

            // search for this node is completed so we can safely backtrack
            nodesToCheck.pop_back();
        }
        // update lastRow for next iteration
        lastRow = row;
    }
}

// assumes first piece is searched exactly
void BiBWT::depthFirstApproxMatch(const Search& s,
                                  const BiAppMatchSA& startMatch,
                                  vector<BiAppMatchSA>& occ, int pieceInSearch,
                                  length_t originalNodesSize) {
    const int pieceNumber = s.pieceOrder[pieceInSearch];
    const length_t& maxED = s.upperBounds[pieceInSearch];
    const length_t& W = maxED - startMatch.editDist;
    const Substring& piece = s.pieces[pieceNumber];
    const Direction& dir = s.directions[pieceInSearch - 1];
    const length_t& nextPiece = pieceInSearch + 1;

    setDirection(dir);

    if (W == 0) {
        const SARangePair newPair =
            matchStringBidirectionally(piece, startMatch.ranges);

        if (!newPair.empty()) {
            const BiAppMatchSA newMatch = makeBiAppMatchSA(
                newPair, startMatch.editDist, piece.size() + startMatch.depth);
            if (nextPiece == s.pieces.size()) {
                occ.emplace_back(newMatch);
            } else {
                depthFirstApproxMatch(s, newMatch, occ, nextPiece,
                                      originalNodesSize);
            }
        }
        return;
    }

    const length_t& lowerBound = s.lowerBounds[pieceInSearch];

    BandMatrix matrix(piece.size() + W + 1, W, startMatch.editDist);

    // get the initial nodes to check
    pushChildren(startMatch.ranges);

    // the size of the pattern, in order not to recalculate this value every
    // iteration
    const unsigned int& patternSize = piece.size();

    // the row of the previous parent that was checked, initialize on zero
    int lastRow = 0;
    int row;

    while (nodesToCheck.size() != originalNodesSize) {

        // get the last node of the vector
        BidirectionalNode& currentNode = nodesToCheck.back();
        row = currentNode.getRow();

        if (row >= lastRow) {
            // in case this is not the result of backtracking (= going a row
            // up) it is a child, calculate the value for this child and
            // update the matrix, also check the children

            // collumns of the matrix for this row
            const length_t firstCol = max<int>(1, row - W);
            const length_t lastCol = min<int>(patternSize, row + W);

            unsigned int minimalEDOfRow = maxED + 1;
            for (length_t j = firstCol; j <= lastCol; j++) {

                matrix.updateMatrix(piece[j - 1] != currentNode.getCharacter(),
                                    row, j);
                matrixElementCounter++;
                minimalEDOfRow = min(minimalEDOfRow, matrix(row, j));
            }

            // edit distance threshold exceeded: backtracking
            if (minimalEDOfRow > maxED) {
                nodesToCheck.pop_back();
                lastRow = row;
                continue;
            }

            // set the found EDEDForThis for this node if the pattern
            // has been completely alligned, else set max allowed ED + 1

            if (lastCol == patternSize) {
                // piece is fully alligned
                currentNode.setFoundED(matrix(row, patternSize));
                if (nextPiece == s.pieces.size()) {
                    // this is the final piece
                }
            }

            length_t EDForThisNode =
                (lastCol == patternSize) ? matrix(row, patternSize) : maxED + 1;
            currentNode.setFoundED(EDForThisNode);

            if (minimalEDOfRow == maxED && lastCol != patternSize) {
                // all children should be exact matches to the character at
                // position just after the minimalED find all j's for which
                // minimalEDofRow == matrix(row, j), each of these j's is a
                // startposition for exact search
                length_t start = max<int>(0, row - W);
                for (length_t j = start; j <= lastCol; j++) {

                    if (matrix(row, j) == minimalEDOfRow) {

                        const Substring& remainingPiece = piece.getSubPiece(j);
                        Substring(piece, piece.begin() + j, piece.end());
                        SARangePair pair = matchStringBidirectionally(
                            remainingPiece, currentNode.getRanges());

                        if (!pair.empty()) {

                            const BiAppMatchSA newMatch =
                                makeBiAppMatchSA(pair, minimalEDOfRow,
                                                 row + remainingPiece.length() +
                                                     startMatch.depth);
                            if (nextPiece == s.pieces.size()) {
                                occ.emplace_back(newMatch);
                            } else {

                                depthFirstApproxMatch(s, newMatch, occ,
                                                      nextPiece,
                                                      nodesToCheck.size());
                                // set direction correct again
                                setDirection(dir);
                            }
                        }
                    }
                }
                // as if this is backtrack
                lastRow = row + 1;
                continue;

            } else {

                // continue the search for children of this node
                pushChildren(currentNode);
            }

        } else {
            // this is the result of backtracking, so we can try and search
            // for cluster centrums and only report these centrums

            BiAppMatchSA newMatch;

            if (row + W >= patternSize) {
                // if fully alligned
                length_t cED = currentNode.getFoundED();
                int pIdx = currentNode.parentIdx;

                if (pIdx == -1) {
                    // presumably not going to happen as patterns are longer
                    // than 2*ED + 1
                    if (currentNode.oneChildWorseOrEqual) {
                        // no parent and one child was worse
                        currentNode.report(newMatch, startMatch.depth, maxED);
                    }
                } else {
                    // cluster time
                    BidirectionalNode& p = nodesToCheck[pIdx];

                    // parent has a child that got here
                    p.hasChildInThisStage |= true;
                    p.oneChildWorseOrEqual |= cED >= p.getFoundED();

                    if (p.getFoundED() >= cED) {
                        // current could be a centre
                        if (!currentNode.hasChildInThisStage ||
                            currentNode.oneChildWorseOrEqual) {
                            currentNode.report(newMatch, startMatch.depth,
                                               maxED);
                        }
                    }
                }

                if (!newMatch.ranges.empty()) {
                    // found a match
                    if (nextPiece == s.pieces.size()) {
                        // reached the end of the search
                        occ.push_back(newMatch);
                    } else if (newMatch.editDist >= lowerBound) {

                        // go deeper on next piece
                        depthFirstApproxMatch(s, newMatch, occ, nextPiece,
                                              nodesToCheck.size());
                        // set direction correct again upon return prev
                        setDirection(dir);
                    }
                }
            }
            // search for this node is completed so we can safely backtrack
            nodesToCheck.pop_back();
        }

        // update lastRow for next iteration
        lastRow = row;
    }
}

void BiBWT::pushChildren(const BidirectionalNode& parentNode) {
    const SARangePair& parentRanges = parentNode.getRanges();
    const length_t row = parentNode.getRow() + 1;
    // iterate over the entire alphabet
    const int parentIdx = nodesToCheck.size() - 1;
    for (length_t i = 0; i < indexToChar.size(); i++) {

        SARangePair pairForNewChar;

        // check if this character occurs in the specified range
        if ((this->*extraChar)(i, parentRanges, pairForNewChar)) {
            // push this range and character for the next iteration
            nodesToCheck.emplace_back(indexToChar[i], pairForNewChar, parentIdx,
                                      row);
            nodeCounter++;
        }
    }
}

void BiBWT::pushChildren(const SARangePair& parentRanges) {

    // iterate over the entire alphabet
    for (length_t i = 0; i < indexToChar.size(); i++) {

        SARangePair pairForNewChar;

        // check if this character occurs in the specified range
        if ((this->*extraChar)(i, parentRanges, pairForNewChar)) {
            // push this range and character for the next iteration
            nodesToCheck.emplace_back(indexToChar[i], pairForNewChar, -1, 1);
            nodeCounter++;
        }
    }
}