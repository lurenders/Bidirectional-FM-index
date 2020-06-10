#ifndef SEARCHSTRATEGY_H
#define SEARCHSTRATEGY_H

#include "BiBWT.h"

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// This is an abstract class. It has one member: a pointer to a bidirectional
// bwt. Every concrete derived class should overload the matchApprox
// memberfunction which matches a pattern approximately using the bidirectional
// bwt
class SearchStrategy {
  protected:
    int numberOfPieces; // the number of pieces the pattern will be split into
    std::vector<Search> searches; // the searches of this strategy

    BiBWT* biBWT; // pointer to the index of the text that is searched

    int lastMaxED =
        0; // the edit distance that was used for the last call to this strategy
           // if the max edit distance changes then the numberOfPieces and the
           // searches need to be recalculated

    std::string name;

    std::map<int, SARangePair> exactRangesPerIndex;

    SearchStrategy(BiBWT* index) {
        biBWT = index;
        lastMaxED = 0;
    }

    /**
     * Calculates the number of pieces for a certain max edit distance. This
     * calculation is strategy dependent
     * @param maxED the maximal allowed edit distance for the alligning
     */
    virtual void calculateNumberOfPieces(unsigned int maxED) = 0;

    /**
     * Creates all searches for this specific strategy
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual void createSearches(unsigned int maxED) = 0;

    /**
     * Splits the pattern into numberOfPieces pieces
     * @param pattern the pattern to be split
     * @param pieces the vector containing the substrings of this pattern, will
     * be cleared and filled during the execution of this method. If the
     * splitting fails for some reason, the vectro will be empty
     */
    void makePieces(const std::string& pattern, std::vector<Substring>& pieces,
                    int maxED) {
        pieces.clear();

        if (lastMaxED != maxED) {
            // calculate how many pieces there will be
            calculateNumberOfPieces(maxED);
            // create the searches
            searches.clear();
            createSearches(maxED);
            // update the last ED seen
            lastMaxED = maxED;
        }

        int numberOfChars = pattern.size();

        if (numberOfPieces >= numberOfChars || numberOfPieces == 1) {
            // no need of splitting up since all pieces would be one character
            // or less or there is only one piece
            return;
        }

        // the length of one piece
        int pieceLength = numberOfChars / numberOfPieces;

        for (int i = 0; i < numberOfPieces; i++) {
            pieces.emplace_back(pattern, i * pieceLength,
                                (i + 1) * pieceLength);
        }

        // change the end on the final piece to the patternsize
        pieces[pieces.size() - 1].setEnd(numberOfChars);
    }

  public:
    virtual void doSearch(const Search& s, std::vector<Substring>& pieces,
                          std::vector<BiAppMatchSA>& allMatches) {

        // the number of the previous piece, needed to know if the matching
        // for the next piece should be forward or backward
        int prevPieceNumber = s.pieceOrder[0];
        auto exactStart = exactRangesPerIndex.find(prevPieceNumber);

        bool exactFound = exactStart != exactRangesPerIndex.end();
        // the start of the search (if exact is found one can start at the
        // second piece)
        int start;
        // the vector with occurences for the previous piece in this search
        std::vector<BiAppMatchSA> prevReports;

        if (exactFound) {
            start = 1;
            // push the exact match to prevReports
            if (exactStart->second.rangeSA.empty()) {
                return;
            }

            prevReports.push_back(makeBiAppMatchSA(
                exactStart->second, 0, pieces[prevPieceNumber].size()));
        } else {
            start = 0;
            prevPieceNumber = numberOfPieces + 1;
            // initialize with the complete range of the suffix array, no depth
            // and no ED
            prevReports.push_back(
                makeBiAppMatchSA(biBWT->getCompleteRange(), 0, 0));
        }
        std::vector<BiAppMatchSA> newReports;

        // do the strategy in order of pieces according to this search
        for (int i = start; i < numberOfPieces; i++) {

            // reports for this piece
            std::vector<BiAppMatchSA>* vectorToAddTo;
            newReports.clear();

            // if this is the last piece than adding can be done directly to
            // allmatches vector
            vectorToAddTo =
                (i == numberOfPieces - 1) ? &allMatches : &newReports;

            // the number of this piece
            int pieceNumber = s.pieceOrder[i];

            // the current piece
            Substring& piece = pieces[pieceNumber];
            // the maximal allowed ED for this piece
            length_t maxEDPiece = s.upperBounds[i];
            // the maximal allowed ED for the next piece
            length_t maxEDNextPiece =
                (i == numberOfPieces - 1) ? maxEDPiece : s.upperBounds[i + 1];

            // set the directions correctly
            Direction d = (pieceNumber > prevPieceNumber) ? FORWARD : BACKWARD;
            biBWT->setDirection(d);
            piece.setDirection(d);

            // go over all the reports from the previous piece
            for (const auto& prevReport : prevReports) {
                if (maxEDPiece == prevReport.editDist) {
                    SARangePair exactRanges = biBWT->matchStringBidirectionally(
                        piece, prevReport.ranges);

                    if (!exactRanges.rangeSA.empty()) {
                        // push new not empty match with old edit distance and
                        // updated depth
                        vectorToAddTo->push_back(
                            makeBiAppMatchSA(exactRanges, prevReport.editDist,
                                             prevReport.depth + piece.size()));
                    }
                    continue;
                }
                biBWT->iteApproxMatch(piece, *vectorToAddTo, prevReport,
                                      maxEDPiece, true, maxEDNextPiece);
            }

            // remove the occurences with an ED lower than the lowerbound as
            // these should not be checked further in this search
            length_t lower = s.lowerBounds[i];
            if (i != numberOfPieces - 1 && lower != 0) {
                newReports.erase(
                    std::remove_if(newReports.begin(), newReports.end(),
                                   [lower, i](const BiAppMatchSA& occ) {
                                       return (occ.editDist < lower ||
                                               (0 == i && occ.depth == 0));
                                   }),
                    newReports.end());
            }

            // update the values to go to next piece
            prevPieceNumber = pieceNumber;
            prevReports.swap(newReports);
            if (prevReports.empty()) {
                return;
            }
        }
    }

    void doSearchDepthFirst(Search& s, std::vector<Substring>& pieces,
                            std::vector<BiAppMatchSA>& allMatches) {

        // first get the bidirectional match of first piece
        int first = s.pieceOrder[0];
        auto exactOfPiece = exactRangesPerIndex.find(first);

        const SARangePair& startRange = exactOfPiece->second;

        if (!startRange.empty()) {
            BiAppMatchSA startMatch =
                makeBiAppMatchSA(startRange, 0, pieces[first].size());
            setPieces(s, pieces);

            biBWT->depthFirstApproxMatch(s, startMatch, allMatches);
        }
    }
    /**
     * Matches the string approximately using a strategy. It does all the
     searches that are strategy-dependentely created.
     * @param pattern the string to match
     * @param ED the allowed edit distance from this pattern
     * @returns a vector containing all approximate matches in the text
     */
    virtual std::vector<AppMatch> matchApprox(const std::string& pattern,
                                              length_t maxED) {
        biBWT->resetCounters();

        if (maxED == 0) {
            biBWT->setDirection(BACKWARD);
            auto result = biBWT->exactMatches(pattern);
            std::vector<AppMatch> returnvalue;
            for (length_t startpos : result) {
                AppMatch m;
                m.editDist = 0;
                m.range = Range(startpos, startpos + pattern.size());
                returnvalue.push_back(m);
            }
            return returnvalue;
        }

        // create the pieces of the pattern
        std::vector<Substring> pieces;
        makePieces(pattern, pieces, maxED);

        if (pieces.empty() || numberOfPieces * maxED >= pattern.size()) {
            // splitting up was not viable just search the entire pattern
            std::cerr << "Warning: Normal bidirectional search was used as "
                         "entered pattern is too short"
                      << std::endl;

            return biBWT->approxMatches(pattern, maxED, true);
        }

        // the vector containing all matches in the sufffix array
        std::vector<BiAppMatchSA> allMatches;

        biBWT->setDirection(FORWARD);
        SARangePair initialRanges = biBWT->getCompleteRange();

        exactRangesPerIndex.clear();
        for (const Search& s : searches) {
            if (s.upperBounds[0] == 0) {
                int pieceNumber = s.pieceOrder[0];
                if (!exactRangesPerIndex.count(pieceNumber)) {
                    exactRangesPerIndex[pieceNumber] =
                        biBWT->matchStringBidirectionally(pieces[pieceNumber],
                                                          initialRanges);
                }
            }
        }

        // do all searches
        for (const Search& s : searches) {
            doSearch(s, pieces, allMatches);
        }

        // return all matches mapped to the text
        return biBWT->mapOccurencesInSAToOccurencesInText(allMatches, maxED,
                                                          pattern);
    }

    std::string getName() {
        return name;
    }

    std::vector<AppMatch> matchApproxDepthFirst(const std::string& pattern,
                                                length_t maxED) {
        biBWT->resetCounters();

        if (maxED == 0) {
            biBWT->setDirection(BACKWARD);
            auto result = biBWT->exactMatches(pattern);
            std::vector<AppMatch> returnvalue;
            for (length_t startpos : result) {
                AppMatch m;
                m.editDist = 0;
                m.range = Range(startpos, startpos + pattern.size());
                returnvalue.push_back(m);
            }
            return returnvalue;
        }
        // create the pieces of the pattern
        std::vector<Substring> pieces;
        makePieces(pattern, pieces, maxED);

        if (pieces.empty() || numberOfPieces * maxED >= pattern.size()) {
            // splitting up was not viable just search the entire pattern
            std::cerr << "Warning: Normal bidirectional search was used as "
                         "entered pattern is too short"
                      << std::endl;

            return biBWT->approxMatches(pattern, maxED, true);
        }

        // the vector containing all matches in the sufffix array
        std::vector<BiAppMatchSA> allMatches;

        biBWT->setDirection(FORWARD);
        SARangePair initialRanges = biBWT->getCompleteRange();
        exactRangesPerIndex.clear();
        for (const Search& s : searches) {
            if (s.upperBounds[0] == 0) {
                int pieceNumber = s.pieceOrder[0];
                if (!exactRangesPerIndex.count(pieceNumber)) {
                    exactRangesPerIndex[pieceNumber] =
                        biBWT->matchStringBidirectionally(pieces[pieceNumber],
                                                          initialRanges);
                }
            }
        }
        biBWT->reserveSize(pattern.length());
        /*for (auto it : exactRangesPerIndex) {
            std::cout << it.first << "\t" << it.second.rangeSA.begin << "\t"
                      << it.second.rangeSA.end << std::endl;
        }*/

        // do all searches
        for (Search& s : searches) {
            doSearchDepthFirst(s, pieces, allMatches);
        }

        // return all matches mapped to the text
        return biBWT->mapOccurencesInSAToOccurencesInText(allMatches, maxED,
                                                          pattern);
    }
};

class KucherovKplus1 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {makeSearch({0, 1}, {0, 0}, {0, 1}),
                                     makeSearch({1, 0}, {0, 0}, {0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({2, 1, 0}, {0, 0, 0}, {0, 1, 2}),
        makeSearch({1, 0, 2}, {0, 0, 1}, {0, 1, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}),
        makeSearch({1, 0, 2, 3}, {0, 0, 1, 1}, {0, 1, 3, 3}),
        makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 1, 3, 3}),
        makeSearch({3, 2, 1, 0}, {0, 0, 1, 1}, {0, 1, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 2, 2, 4, 4}),
        makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 4, 4}),
        makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}),
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}),
        makeSearch({3, 2, 4, 1, 0}, {0, 0, 0, 1, 1}, {0, 1, 2, 4, 4}),
        makeSearch({2, 1, 0, 3, 4}, {0, 0, 0, 1, 3}, {0, 1, 2, 4, 4}),
        makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 2, 4}, {0, 1, 2, 4, 4}),
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 3, 4}, {0, 0, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = maxED + 1;
    }
    void createSearches(unsigned int maxED) {
        if (maxED < 1 || maxED > 4) {
            throw std::invalid_argument("max ED should be between 1 and 4");
        }
        searches = schemePerED[maxED - 1];
    }

  public:
    KucherovKplus1(BiBWT* index) : SearchStrategy(index) {
        name = "Kucherov K + 1";
    };
};

class KucherovKplus2 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}),
        makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 1, 2}),
        makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({1, 2, 3, 0}, {0, 0, 0, 1}, {0, 0, 1, 2}),
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 2}, {0, 0, 2, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 3, 3}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 2, 2, 3}),
        makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 1, 3, 3}),
        makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 1, 2}, {0, 0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 4, 4}),
        makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 4, 4}),
        makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 1}, {0, 1, 2, 2, 4, 4}),
        makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 1, 2}, {0, 1, 1, 3, 4, 4}),
        makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 2, 3}, {0, 1, 1, 2, 4, 4}),
        makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 1, 3, 3}, {0, 0, 3, 3, 4, 4}),
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 3, 3, 3}, {0, 0, 3, 3, 4, 4}),
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 4, 4}, {0, 0, 2, 4, 4, 4}),
        makeSearch({2, 3, 1, 0, 4, 5}, {0, 0, 0, 1, 2, 4}, {0, 0, 2, 2, 4, 4}),
        makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 4, 4}, {0, 0, 1, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = maxED + 2;
    }
    void createSearches(unsigned int maxED) {
        searches = schemePerED[maxED - 1];
    }

  public:
    KucherovKplus2(BiBWT* index) : SearchStrategy(index) {
        name = "Kucherov K + 2";
    };
};

class OptimalKianfar : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {makeSearch({0, 1}, {0, 0}, {0, 1}),
                                     makeSearch({1, 0}, {0, 1}, {0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}),
        makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({1, 2, 0}, {0, 1, 1}, {0, 1, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 3}, {0, 2, 3, 3}),
        makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {1, 2, 3, 3}),
        makeSearch({2, 3, 1, 0}, {0, 0, 2, 2}, {0, 0, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 4}, {0, 3, 3, 4, 4}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {2, 2, 3, 3, 4}),
        makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 3, 3}, {0, 0, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = maxED + 1;
    }
    void createSearches(unsigned int maxED) {
        if (maxED < 1 || maxED > 4) {
            throw std::invalid_argument("max ED should be between 1 and 4");
        }
        searches = schemePerED[maxED - 1];
    }

  public:
    OptimalKianfar(BiBWT* index) : SearchStrategy(index) {
        name = "Optimal Kianfar";
    };
};

// ============================================================================
// CLASS O1StarSearchStrategy
// ============================================================================

// A concrete derived class of SearchStrategy. The strategy here is founded on
// this observation: if x errors are allowed and the pattern is divided up in (x
// + 2) pieces then every match with max x erros contains a seed consisting of n
// pieces, where the first and last piece of the seed contain no errors and all
// pieces inbetween these contain exacly one error. (2 <= n <= x + 2)

class O1StarSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}),
        makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 0, 2, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 4, 4, 4, 4}),
    };

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = maxED + 2;
    }
    void createSearches(unsigned int maxED) {
        searches = schemePerED[maxED - 1];
    }

  public:
    O1StarSearchStrategy(BiBWT* index) : SearchStrategy(index) {
        name = "01*0";
    };
};

// man best for ed 4 otherwise 01*0
class ManBestStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}),
        makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 0, 2, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 4}, {0, 3, 3, 3, 4, 4}),
        makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0}, {0, 2, 2, 3, 3, 4}),
        makeSearch({2, 1, 3, 4, 5, 0}, {0, 1, 1, 1, 1, 1}, {0, 2, 2, 3, 3, 4}),
        makeSearch({3, 2, 1, 4, 5, 0}, {0, 1, 2, 2, 2, 2}, {0, 1, 2, 3, 3, 4}),
        makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 3, 3}, {0, 0, 4, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = maxED + 2;
    }
    void createSearches(unsigned int maxED) {
        searches = schemePerED[maxED - 1];
    }

  public:
    ManBestStrategy(BiBWT* index) : SearchStrategy(index) {
        name = "manbest";
    };
};
// ============================================================================
// CLASS PIGEONHOLESEARCHSTRATEGY
// ============================================================================

// A concrete derived class of SearchStrategy. The strategy here is founded on
// this observation: if x errors are allowed and the pattern is divided up in (x
// + 1) sections then every approximate match has an exact match with at least
// one of the sections. The strategy iterates over the sections, it tries to
// exactly match the current section, then approximately match the pattern
// before this section and after the the pattern after this section with the
// remaining edit distance.
class PigeonHoleSearchStrategy : public SearchStrategy {

  private:
    const std::vector<Search> ED1 = {makeSearch({0, 1}, {0, 0}, {0, 1}),
                                     makeSearch({1, 0}, {0, 0}, {0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({1, 2, 0}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        makeSearch({1, 0, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = maxED + 1;
    }
    void createSearches(unsigned int maxED) {
        searches = schemePerED[maxED - 1];
    }

  public:
    PigeonHoleSearchStrategy(BiBWT* index) : SearchStrategy(index) {
        name = "Pigeon";
    };
};

class BackTrackStrategyNaive : public SearchStrategy {
  private:
    BWT* bwt;
    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = 1;
    }
    void createSearches(unsigned int maxED) {
    }

  public:
    virtual std::vector<AppMatch> matchApprox(const std::string& pattern,
                                              length_t maxED) {

        bwt->resetCounters();
        if (maxED == 0) {

            auto result = bwt->exactMatches(pattern);
            std::vector<AppMatch> returnvalue;
            for (length_t startpos : result) {
                AppMatch m;
                m.editDist = 0;
                m.range = Range(startpos, startpos + pattern.size());
                returnvalue.push_back(m);
            }
            return returnvalue;
        }

        return bwt->approxMatches(pattern, maxED);
    }

    BackTrackStrategyNaive(BiBWT* index) : SearchStrategy(index) {
        name = "Naive backtracking";
        bwt = (BWT*)index;
    };
};

class BackTrackStrategyImproved : public SearchStrategy {
  private:
    void calculateNumberOfPieces(unsigned int maxED) {
        numberOfPieces = 1;
    }
    void createSearches(unsigned int maxED) {
    }

  public:
    virtual std::vector<AppMatch> matchApprox(const std::string& pattern,
                                              length_t maxED) {

        biBWT->resetCounters();
        if (maxED == 0) {
            biBWT->setDirection(BACKWARD);
            auto result = biBWT->exactMatches(pattern);
            std::vector<AppMatch> returnvalue;
            for (length_t startpos : result) {
                AppMatch m;
                m.editDist = 0;
                m.range = Range(startpos, startpos + pattern.size());
                returnvalue.push_back(m);
            }
            return returnvalue;
        }

        return biBWT->approxMatches(pattern, maxED, false);
    }

    BackTrackStrategyImproved(BiBWT* index) : SearchStrategy(index) {
        name = "Improved backtracking";
    };
};

#endif