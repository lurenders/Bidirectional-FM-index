#ifndef BIBWT_H
#define BIBWT_H

#include "BWT.h"
#include <chrono>
class BiBWT;

typedef bool (BiBWT::*ExtraCharPtr)(length_t, const SARangePair&,
                                    SARangePair&) const;

class BiBWT : public BWT {
  private:
    class BidirectionalNode {
      public:
        char c;             // the character of this node
        SARangePair ranges; // the ranges of this node
        int row;            // the row of this node in this search

        length_t EDFound; // the found edit distance

        // BidirectionalNode* parent; // a pointer to the parent of this node
        int parentIdx;

        bool reported = false;
        bool hasChildInThisStage = false;
        bool oneChildWorseOrEqual = false;

        /**
         * Create a node of the search tree
         * @param character the character of this node
         * @param ranges the ranges over the suffix and reversed suffix array
         * that go to this node
         * @param parentidx the index of the parent in the stack (=
         * nodesToCheck)
         * @param row the row of this node in the allignment matrix = depth of
         * this node
         */
        BidirectionalNode(char character, const SARangePair& ranges,
                          const int& parentidx, const length_t& row) {
            this->c = character;
            this->ranges = ranges;
            this->parentIdx = parentidx;
            this->row = row;
            // intit all bools on false
            this->EDFound = 0;

            this->reported = false;
            this->hasChildInThisStage = false;
            this->oneChildWorseOrEqual = false;
        }

        /**
         * Sets the report flag to true
         */
        void report() {
            reported = true;
        }

        /**
         * Reports the match at this node (with added depth) to the occ vector
         * if it hasn't done this already
         * @param occ the vector of occurences, the match at this node might be
         * added if it hasn't been added before
         * @param startDepth the depth to add to this match
         *
         */
        void report(std::vector<BiAppMatchSA>& occ, const length_t startDepth,
                    length_t maxED) {
            if (!reported && EDFound <= maxED) {
                // make a match with the ranges, ED and depth = row of this node
                occ.push_back(
                    makeBiAppMatchSA(ranges, EDFound, row + startDepth));
                report();
            }
        }

        /**
         * Reports the match (with added depth) at this node, if this node
         * hasn't been reported yet
         * @param occ the match will be stored here
         * @param startDepth the depth to add to the match
         */
        void report(BiAppMatchSA& occ, const length_t startDepth,
                    length_t maxED) {
            if (!reported && EDFound <= maxED) {
                occ = makeBiAppMatchSA(ranges, EDFound, row + startDepth);
                report();
            }
        }

        /**
         * Update the value for the hasAWorseOrEqualChild member. If the
         * parameter is true than this bool will be set to true, else noting
         * will change
         * @param reportingValue the update value, if it is true the
         * hasAWorseOrEqualChild member will be true, if it is false the member
         * will stay as it was
         */
        void updateChildIsWorseOrEqual(bool reportingValue) {
            oneChildWorseOrEqual |= reportingValue;
        }

        void updateChildIsWorseOrEqual(const length_t childED) {
            updateChildIsWorseOrEqual(childED >= EDFound);
        }

        /**
         * Gets the ranges of this node
         */
        const SARangePair& getRanges() const {
            return ranges;
        }

        /**
         * Sets the found ED for this node
         * @param edit, the found edit distance
         */
        void setFoundED(length_t edit) {
            EDFound = edit;
        }

        const char getCharacter() const {
            return c;
        }

        unsigned int getRow() const {
            return row;
        }

        const length_t getFoundED() const {
            return EDFound;
        }
        /**
         * Sets the child was checked flag to true
         */
        void childWasChecked() {
            hasChildInThisStage = true;
        }
    };
    PrefixOccurrences4 prefOccRev; // the prefix occurences of the rev BWT
    Direction dir = BACKWARD;
    ExtraCharPtr extraChar; // pointer to extra char method (for direction)

    std::vector<BidirectionalNode> nodesToCheck; // stack of nodes to check

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Finds the LF mapping of the character at index k in the bwt string
     * @param k the index to find the LF mapping off
     * @param reversed boolean to indicate if the reversed bwt should be used
     * or the normal bwt
     * @returns the row that is the LF mapping of k. It is so that the entry in
     * the suffix array of this return value is one less than the entry in the
     * sufffix array at index k
     */
    length_t findLF(length_t k, bool reversed) const;

    /**
     * Function that returns the nummber of occurences before an index of the
     * symbol at symbolindex in the alphabet
     * @param symbolIndex the index of the the symbol in the alphabet to count
     * the occrences of at index index
     * @param index the index whose entry for symbol in the occurences table is
     * asked
     * @return the number of occurences of the symbol before index in the bwt
     */
    length_t getNumberOfOcc(length_t symbolIndex, length_t index) const {
        return prefixOccurrences.getNumberOfOccurences(symbolIndex, index);
    }
    /**
     * Same as BiBWT::getNumberOfOccurences, but now in the bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet to get the number of
     * occurences of
     * @param index the index in the occurences table, for which the number of
     * occurences of the symbol at symbolindex is asked.
     * @return the number of occurences at index index in the occurences table
     * of the bwt of the reversed text for the symbol at symbolIndex in the
     * alphabet
     */
    length_t getNumberOfOccRev(length_t symbolIndex, length_t index) const {
        return prefOccRev.getNumberOfOccurences(symbolIndex, index);
    }

    /**
     * Function that returns the nummber of occurences before the index of all
     * symbols smaller than the symbol at symbolindex in the  bwt
     * @param symbolIndex the index in the alphabet whose number of prefix
     * occurences is queried.
     * @param index the index whose entry for symbol in the prefixoccurences
     * table is asked
     * @return the number of occurences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfPrefOcc(length_t symbolIndex, length_t index) const {
        return prefixOccurrences.getNumberOfPrefixOccurences(symbolIndex,
                                                             index);
    }

    /**
     * Function that returns the nummber of occurences before the index of all
     * symbols smaller than the symbol at symbolindex in the  bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet whose number of prefix
     * occurences is queried.
     * @param index the index whose entry for symbol in the rprefixoccurences
     * table of the reverse text is asked
     * @return the number of occurences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfPrefOccRev(length_t symbolIndex, length_t index) const {
        return prefOccRev.getNumberOfPrefixOccurences(symbolIndex, index);
    }

    // ----------------------------------------------------------------------------
    // HELP ROUTINES FOR APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Private helper function that gets all possible next characters and their
     * ranges given the ranges of the current pattern
     * @param rangesOfParent the ranges of the parent in the bwt
     * @returns a vector with pairs of ranges and charachters, these are the
     * ranges of all possible previous/next characters
     */
    std::vector<std::pair<SARangePair, char>>
    getCharExtensions(const SARangePair& rangesOfParent) const;

    /**
     * Helper function for getCharExtensions. Finds the ranges of cP
     * using the principle explained in the paper of Lahm
     * @param positionInAlphabet the postition in alhabet of the character c
     * that is added in the front
     * @param rangesOfP the ranges of pattern P
     * @returns the ranges cP
     */
    bool findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                         const SARangePair& rangesOfP,
                                         SARangePair& childRanges) const;

    /**
     * Helper function for getCharExtensions. Finds the ranges of a character Pc
     * using the principle explained in the paper of Lahm
     * @param positionInAlphabet the postition in alhabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @returns the ranges of Pc
     */
    bool findRangesWithExtraCharForward(length_t positionInAlphabet,
                                        const SARangePair& rangesOfP,
                                        SARangePair& childRanges) const;

    // ----------------------------------------------------------------------------
    // HELPER ROUTINES FOR APPROXIMATE MATCHING (ITERATIVELY)
    // ----------------------------------------------------------------------------

    /**
     * Creates all children of a rangepair. And adds these nodes to the
     * nodesToCheck vector
     * @param nodesToCheck the vector of all nodes to check, the children will
     * be pushed to this vector
     * @param ranges the ranges to get the children of
     * @param parent pointer to the parent node of this childern. Default is
     * NULL
     */
    void pushChildren(const BidirectionalNode& parentNode);

    /**
     * Creates all children of a rangepair. And adds these nodes to the
     * nodesToCheck vector
     * @param nodesToCheck the vector of all nodes to check, the children will
     * be pushed to this vector
     * @param ranges the ranges to get the children of

     */
    void pushChildren(const SARangePair& ranges);

    /**
     * Reports the clusterCentrum of a specific branch. Each node will check up
     * to clusterWidth ancestors to find the minimum of the cluster. And this
     * ancestor will then report. A cluster always exists of a minimum and then
     * increasing values upwards and downwards. This minimum is the centrum of
     * the cluster.
     * @param occ the vector to report the matches to
     * @param maxED the maximal allowed ED
     * @param clusterWidth the max width of a cluster (both down and up)
     *
     */
    void reportClusterCentrum(std::vector<BiAppMatchSA>& occ,
                              const length_t maxED, const int clusterWidth,
                              const length_t startDepth) __attribute_noinline__;

    /**
     * Reports the clusterCentrum of a specific branch. Each node will check up
     * to clusterWidth ancestors to find the minimum of the cluster. And this
     * ancestor will then report. A cluster always exists of a minimum and then
     * increasing values upwards and downwards. This minimum is the centrum of
     * the cluster.
     * @param occ the vector to report the matches to
     * @param maxED the maximal allowed ED
     * @param clusterWidth the max width of a cluster (both down and up)
     *
     */
    void reportClusterCentrum(BiAppMatchSA& occ, const length_t maxED,
                              const int clusterWidth,
                              const length_t startDepth) __attribute_noinline__;

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor with default values for the sparseness factors
     * @param prefix prefix of the files that contain the info
     */
    BiBWT(const std::string& prefix) : BWT(prefix) {

        // read the reverse prefix occurrence table
        std::cout << "Reading " << prefix << ".rev.poc.dna" << std::endl;
        if (!prefOccRev.read(prefix + ".rev.poc.dna"))
            throw std::runtime_error("Cannot open file: " + prefix +
                                     ".rev.poc.dna");
        std::cout << "Done reading reverse prefix occurrence table"
                  << std::endl;

        extraChar = &BiBWT::findRangesWithExtraCharBackward;
    }

    /**
     * Constructor
     * @param prefix prefix of the files that contain the info
     * @param sa_spase sparseness factor of suffix array. It is assumed this is
     * a power of two
     */
    BiBWT(const std::string& prefix, int sa_sparse) : BWT(prefix, sa_sparse) {

        // read the reverse prefix occurrence table
        if (!prefOccRev.read(prefix + ".rev.poc.dna"))
            throw std::runtime_error("Cannot open file: " + prefix +
                                     ".rev.poc.dna");
        std::cout << "Done reading reverse prefix occurrence table"
                  << std::endl;
        // initializeBiBWT(prefix + ".rev.txt", prefix + ".rev.sa");
    }

    /**
     * Get the complete range of this index
     * @returns an SARangePair with both ranges the complete range of the index
     */
    SARangePair getCompleteRange() const {
        return SARangePair(Range(0, bwt.size()), Range(0, bwt.size()));
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR EXACT MATCHING
    // ----------------------------------------------------------------------------

    /**
     * This function matches a string exactly and returns the ranges in the sa
     * and saRev
     * @param pattern the string to match
     * @returns the pair of ranges of this pattern
     */
    SARangePair matchStringBidirectionally(const Substring& pattern) {
        return matchStringBidirectionally(pattern, getCompleteRange());
    }

    /**
     * This function matches a string exactly and returns the ranges in the sa
     * and saRev
     * @param pattern the string to match
     * @returns the pair of ranges
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           SARangePair startRange);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    void setDirection(Direction d) {
        dir = d;
        extraChar = (d == FORWARD) ? &BiBWT::findRangesWithExtraCharForward
                                   : &BiBWT::findRangesWithExtraCharBackward;
    }
    /**
     * recursively calculate which ranges in the suffix array correspond to
     * approximate matches to a pattern. Where the redundancy at the start and
     * end is removed
     * @param rangeW the initial range we are considering
     * @param P the pattern to match
     * @param M the banded matrix for keeping track of the edit distance
     * @param bestEDFromParent the bestED untill this position
     * @param occ a vector of found matches, these contain ranges in the suffix
     * array and their ED and depth
     * @param depth how many iterations (= length of matched substring untill
     * now)
     */
    void improvedRecApproxMatches(const SARangePair& ranges, const Substring& P,
                                  EditMatrix& M,
                                  const EditDistance& bestEDFromParent,
                                  std::vector<BiAppMatchSA>& occ,
                                  length_t depth);

    /**
     * Matches the pattern approximately. All matches are at most a certain edit
     * distance away from the pattern
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @returns a vector with matches which contain a range (the range of the
     * text that matched) and the edit distance this substring is away from the
     * pattern
     * @param forward a boolean to indicate if the search is forward or backward
     */
    std::vector<AppMatch> approxMatches(const std::string& pattern, int maxED,
                                        bool forward);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING (ITERATIVELY)
    // ----------------------------------------------------------------------------

    /**
     * Matches the parent approximately without recursive calls.
     * @param pattern the pattern to match
     * @param occ the vector wich stores all occurences found during the metoh
     * @param match the match to start searchin on, it contains ranges to search
     * in, a start edit distance and start depth
     * @param maxED the maximal allowed edit distance
     * @param leadingGapsAllowed bool indicating wheter the match can start with
     * leading gaps or not
     * @param maxEDNext the max allowed edit distance for the next piece that
     * might be matched after the occurences found for this pattern
     */
    void iteApproxMatch(const Substring& pattern,
                        std::vector<BiAppMatchSA>& occ,
                        const BiAppMatchSA& match, length_t maxED,
                        bool leadingGapsAllowed, int maxEDNext);

    void reserveSize(const length_t size) {
        nodesToCheck.reserve(size * indexToChar.size());
    }

    void depthFirstApproxMatch(const Search& search,
                               const BiAppMatchSA& startMatch,
                               std::vector<BiAppMatchSA>& occ,
                               int pieceInSearch = 1,
                               length_t originalNodesSize = 0);

    /**
     * Matches the parent approximately without recursive calls. It calls the
     * full method with maxEDNext = maxED
     */
    void iteApproxMatch(const Substring& pattern, const BiAppMatchSA& match,
                        std::vector<BiAppMatchSA>& occ, length_t maxED,
                        bool leadingGapsAllowed) {
        iteApproxMatch(pattern, occ, match, maxED, leadingGapsAllowed, maxED);
    }

    // ----------------------------------------------------------------------------
    // POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * This function maps matches in the sa to matches in the text. It takes the
     * ranges of the matches and together with the depth this is matched to a
     * range in the text (this new range has a width of depth). The edit
     * distance is also mapped to this new range in the text
     * @returns a vector of matches in the text containing a range and the edit
     * distance
     */
    std::vector<AppMatch> mapOccurencesInSAToOccurencesInText(
        const std::vector<BiAppMatchSA>& occurences, int ED,
        const Substring& pattern);
};

#endif