
#include "searchstrategy.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <string.h>

using namespace std;

// try to find the word that encompasses this range (for words seperated by %
// signs)
string findWordWithRange(Range& range, const string& text) {
    length_t start = 0;
    length_t end = text.size();
    for (length_t i = range.begin; i-- > 0;) {
        if (text[i] == '%') {
            start = i + 1;
            break;
        }
    }

    for (length_t i = range.end; i <= text.size(); i++) {
        if (text[i] == '%') {
            end = i;
            break;
        }
    }

    return text.substr(start, end - start);
}

void printMatches(vector<AppMatch> matches, string text, bool printLine,
                  string duration, BWT* mapper, string name) {

    cout << endl;

    tuple<length_t, length_t, length_t> counters = mapper->getCounters();

    cout << name << ":\tduration: " << duration
         << "µs\t nodes visited: " << get<0>(counters)
         << "\t matrix elements written: " << get<1>(counters)
         << "\t startpositions reported: " << get<2>(counters)
         << " #matches: " << matches.size() << endl;

    for (auto match : matches) {
        cout << "Found match at position " << match.range.begin << " with ED "
             << match.editDist << endl;

        cout << "\tCorresponding substring:\t"
             << text.substr(match.range.begin,
                            match.range.end - match.range.begin)
             << endl;
        if (printLine) {
            cout << "\tCorresponding word:\t"
                 << findWordWithRange(match.range, text) << endl;
        }
    }
}

void outputForMapper(BiBWT* mapper, bool print, bool depth) {
    cout << "Recompute the original text for displaying " << endl;
    // this can take a long time, but is not necessary for the program
    // this is only for displaying and debugging purposes. The BWT does not need
    // the text to do its magic
    string text = mapper->getText();
    if (print) {
        cout << "Original text:\t" << text << endl;
    }
    cout << "Length of text = " << text.size() << endl;
    string input;
    vector<AppMatch> matches;

    cout << "Type string to match or 0000 to stop" << endl;
    cin >> input;

    PigeonHoleSearchStrategy pigeon = PigeonHoleSearchStrategy(mapper);
    O1StarSearchStrategy O1star = O1StarSearchStrategy(mapper);
    KucherovKplus1 kucherov = KucherovKplus1(mapper);
    int ED = 2;
    while (input.compare("0000") != 0) {
        cout << "Approximate matches for " << input << " with max ED " << ED
             << endl;

        auto t1 = chrono::high_resolution_clock::now();
        matches = ((BWT*)mapper)->approxMatches(input, ED);
        auto t2 = chrono::high_resolution_clock::now();

        auto duration =
            chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
        string d = to_string(duration);
        printMatches(matches, text, print, d, mapper, "NAIVE BACKWARD");

        t1 = chrono::high_resolution_clock::now();
        matches = mapper->approxMatches(input, ED, false);
        t2 = chrono::high_resolution_clock::now();

        duration = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
        d = to_string(duration);
        printMatches(matches, text, print, d, mapper, "Bidirectional BACKWARD");

        t1 = chrono::high_resolution_clock::now();
        matches = mapper->approxMatches(input, ED, true);
        t2 = chrono::high_resolution_clock::now();

        duration = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
        d = to_string(duration);
        printMatches(matches, text, print, d, mapper, "Bidirectional FORWARD");

        t1 = chrono::high_resolution_clock::now();
        if (depth) {
            matches = pigeon.matchApproxDepthFirst(input, ED);
        } else {
            matches = pigeon.matchApprox(input, ED);
        }
        t2 = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
        d = to_string(duration);
        printMatches(matches, text, print, d, mapper, "PIGEON");

        /** t1 =  chrono::high_resolution_clock::now();
         matches = O1star.matchApprox(input, ED);
         t2 =  chrono::high_resolution_clock::now();
         duration =
              chrono::duration_cast< chrono::microseconds>(t2 - t1)
                 .count();
         d =  to_string(duration);
         printMatches(matches, text, print, d, mapper, "O1*0");

         t1 =  chrono::high_resolution_clock::now();
         matches = kucherov.matchApprox(input, ED);
         t2 =  chrono::high_resolution_clock::now();
         duration =
              chrono::duration_cast< chrono::microseconds>(t2 - t1)
                 .count();
         d =  to_string(duration);
         printMatches(matches, text, print, d, mapper, "Kucherov");*/

        cout << endl;

        cout << "Type string to match or 0000 to stop" << endl;
        cin >> input;
    }
}

void printHelp() {
    cout << "This program takes a few inputs. The first input is required: it "
            "is the "
            "prefix of the files that will be used, the extensions "
            "needed for this prefix are .txt .rev.txt .sa and .rev.sa"
         << endl;
    cout << endl;
    cout << "the following flags are optional" << endl;
    cout << "\t"
         << "--sparse"
         << "\t\t"
         << "to set a sparseness "
            "value for the suffix array add --sparse flag and the value"
         << endl;
    cout << "\t-depth\t\t"
         << "Add this flag to test the depth first search" << endl;

    cout << "\t--reads\t\t"
         << "Add this flag to specify the prefix of the file containing the "
            "reads (without the _counts_250.csv)"
         << endl;
}

map<int, vector<string>> getReads(const string& file) {
    map<int, vector<string>> reads;
    for (int i = 0; i < 5; i++) {
        reads[i] = {};
        reads[i].reserve(500000);
    }

    // putting the pos on 100
    reads[100] = {};
    reads[100].reserve(500000);

    ifstream ifile = getStream(file);
    string line;
    // get the first line we do not need this
    getline(ifile, line);

    while (getline(ifile, line)) {
        istringstream iss{line};

        vector<string> tokens;
        string token;

        while (getline(iss, token, ',')) {
            tokens.push_back(token);
        }
        for (auto token : tokens) {
            //  cout << token <<  endl;
        }
        // push the position on ED 100
        reads[100].push_back(tokens[1]);

        for (int ED = 0; ED < 5; ED++) {
            // the ED + 2nd column contains a read for ED
            reads[ED].push_back(tokens[ED + 2]);
        }

        if (reads[0].size() % 10000 == 0) {
            cout << reads[0].size() << " lines read" << endl;
        }
    }

    return reads;
}

double
avgVec(vector<length_t> const& v) // note: the average must not be an integer
{
    return v.empty() ? 0.0 : accumulate(v.begin(), v.end(), 0.0) / v.size();
    ;
}

void doBench(map<int, vector<string>>& reads, BiBWT* mapper,
             SearchStrategy* strategy, bool depth, bool positionsFromSameFile,
             length_t ED) {
    string text = mapper->getText();
    vector<length_t> durations;
    vector<length_t> nodes;
    vector<length_t> matrixElements;
    vector<length_t> uniqueMatches;
    vector<length_t> totalreportedmatches;

    vector<AppMatch> matches;
    int i = 1;
    length_t totalDuration = 0;
    length_t subDuration = 0;
    cout << "Benchmarking with " << strategy->getName() << " strategy for ED "
         << ED << endl;
    cout.precision(2);
    for (auto read : reads[ED]) {

        auto t1 = chrono::high_resolution_clock::now();
        if (depth) {

            matches = strategy->matchApproxDepthFirst(read, ED);

        } else {
            matches = strategy->matchApprox(read, ED);
        }
        auto t2 = chrono::high_resolution_clock::now();
        auto duration =
            chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
        totalDuration += duration;
        subDuration += duration;
        if (i % 10000 == 0) {
            double div = totalDuration / 1000;
            div /= 1000;

            cout << i << ": " << subDuration << " µs";
            cout << "\ttotal:  " << totalDuration << " µs (" << fixed << div
                 << " s)" << endl;

            subDuration = 0;
        }

        // stoi is the problem -> string to unsigned long should probably be
        length_t originalPos = (length_t)stoull(reads[100][i - 1]);

        bool originalFound = false;

        for (auto match : matches) {
            if (match.range.begin >= originalPos - (ED + 2) &&
                match.range.begin <= originalPos + (ED + 2)) {
                originalFound = true;
                break;
            }
        }

        i++;
        uniqueMatches.push_back(matches.size());
        durations.push_back(duration);
        auto counters = mapper->getCounters();

        nodes.push_back(get<0>(counters));
        matrixElements.push_back(get<1>(counters));
        totalreportedmatches.push_back(get<2>(counters));

        if (matches.size() == 0 || (positionsFromSameFile && !originalFound)) {
            printMatches(matches, text, false, " ", (BWT*)mapper, read);
            cout << i - 1 << " woops" << endl;
            cout << originalPos << endl;
            auto n = ((BWT*)mapper)->approxMatches(read, ED);
            cout << "size with original mapper " << n.size() << endl;
        }
    }
    cout << "Results for " << strategy->getName() << endl;

    double div = totalDuration / 1000;
    div /= 1000;

    cout << "Total duration: " << totalDuration << " µs (" << fixed << div
         << " s)" << endl;
    cout << "Average duration: " << avgVec(durations) << " µs" << endl;
    cout << "Average nodes: " << avgVec(nodes) << endl;
    cout << "Average matrix elements written: " << avgVec(matrixElements)
         << endl;
    cout << "Average number of unique matches: " << avgVec(uniqueMatches)
         << endl;
    cout << "Average total number of reported matches "
         << avgVec(totalreportedmatches) << endl;
}

void benchmark(map<int, vector<string>>& reads, BiBWT* mapper, bool depth,
               bool positionsFromSameFile, length_t ED) {

    vector<SearchStrategy*> strategies;
    PigeonHoleSearchStrategy s = PigeonHoleSearchStrategy(mapper);
    strategies.push_back(&s);
    O1StarSearchStrategy s2 = O1StarSearchStrategy(mapper);
    strategies.push_back(&s2);

    KucherovKplus1 s3 = KucherovKplus1(mapper);
    strategies.push_back(&s3);

    KucherovKplus2 s4 = KucherovKplus2(mapper);
    strategies.push_back(&s4);

    ManBestStrategy s8 = ManBestStrategy(mapper);
    if (ED == 4) {
        strategies.push_back(&s8);
    }

    OptimalKianfar s5 = OptimalKianfar(mapper);
    strategies.push_back(&s5);

    BackTrackStrategyImproved s7 = BackTrackStrategyImproved(mapper);
    strategies.push_back(&s7);

    BackTrackStrategyNaive s6 = BackTrackStrategyNaive(mapper);
    strategies.push_back(&s6);

    for (SearchStrategy* s : strategies) {
        doBench(reads, mapper, s, depth, positionsFromSameFile, ED);
    }
}

void showUsage() {
    cout << "Usage: ./BWT [options] basefilename patterns.csv\n\n";
    cout << " [options arg]\n";
    cout << "  -e  --max-ed\t\tmaximum edit distance [default = 0]\n";
    cout << "  -s  --sa-sparseness\tsuffix array sparseness factor [default = "
            "1]\n\n";
    cout << "  -d  --depth\tAdd flag to do the search depth first [default = "
            "false]\n";
    cout << "  -ri  --read-positions-incorrect\tAdd flag if the positions in "
            "the file containing the reads are NOT based on the basefiles "
            "[default = "
            "false]\n";

    cout << "Following input files are required:\n";
    cout << "\t<base filename>.txt: input text T\n";
    cout << "\t<base filename>.cct: charachter counts table\n";
    cout << "\t<base filename>.sa.[saSF]: suffix array sample every [saSF] "
            "elements\n";
    cout << "\t<base filename>.bwt: BWT of T\n";
    cout << "\t<base filename>.poc: Prefix occurrence table of T\n";
    cout << "\t<base filename>.rev.poc: Prefix occurrence table of the reverse "
            "of T\n";
}

int main(int argc, char* argv[]) {
    const int requiredArguments =
        2; // prefix of files and file containing reads

    if (argc < requiredArguments) {
        cerr << "Insufficient number of arguments" << endl;
        showUsage();
        return EXIT_FAILURE;
    }
    if (argc == 2 && strcmp("help", argv[1]) == 0) {
        showUsage();
        return EXIT_SUCCESS;
    }

    string saSparse = "1";
    string maxED = "0";

    bool depthFirst = false;
    bool readPosCorrect = true;

    // process optional arguments
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];

        if (arg == "-d" || arg == "--depth") {
            depthFirst = true;

        } else if (arg == "-s" || arg == "--sa-sparseness") {
            if (i + 1 < argc) {
                saSparse = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-e" || arg == "--max-ed") {
            if (i + 1 < argc) {
                maxED = argv[++i];
            }
        } else if (arg == "-ri" || arg == "--read-positions-incorrect") {
            readPosCorrect = false;
        }

        else {
            cerr << "Unknown argument: " << arg << " is not an option" << endl;
            return false;
        }
    }

    length_t ed = stoi(maxED);
    if (ed < 0 || ed > 4) {
        cerr << ed << " is not allowed as maxED should be in [0, 4]" << endl;

        return EXIT_FAILURE;
    }
    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
    }

    string depth = (depthFirst) ? "with" : "without";
    cout << "searching " << depth << " depth first" << endl;

    // check if the sparse factors are ints
    string prefix = argv[argc - 2];
    string readsFile = argv[argc - 1];

    cout << "Reading in reads from " << readsFile << endl;
    auto reads = getReads(readsFile);
    cout << "Start creation of BWT approximate matcher" << endl;

    BiBWT bwt = BiBWT(prefix, saSF);
    benchmark(reads, &bwt, depthFirst, readPosCorrect, ed);

    // outputForMapper(&bwt, true, depthFirst);
}
