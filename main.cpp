#include "utils/libdeflate_wrapper.h"

#include <iostream>
#include <fstream>
#include "utils/kseq.h"
#include "utils/helper.h"
#include "matching/CopMEMMatcher.h"
#include <omp.h>

KSEQ_INIT(gzFile, gzread)

using namespace std;

int k;

CopMEMMatcher* matcher;

string ref;
int rcStart;

void loadRef(kseq_t* seq) {
    while (kseq_read(seq) >= 0) {
        ref.append(seq->seq.s, seq->seq.l);
    }
    rcStart = ref.size();
    if (rcStart == 0) {
        fprintf(stderr, "ERROR: reference sequence is empty\n");
        exit(EXIT_FAILURE);
    }
    ref.resize(rcStart * 2);
    char *refPtr = (char *) ref.data();
    cout << "loaded reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    PgHelpers::reverseComplement(refPtr, rcStart, refPtr + rcStart);
    cout << "reversed reference - " << PgHelpers::time_millis() << " [ms]" << endl;
}

size_t unmatchedCharsAll = 0;
size_t totalMatchedAll = 0;
size_t totalDestOverlapAll = 0;
size_t totalDestLenAll = 0;

string getTotalMatchStat(uint32_t totalMatchLength, size_t destLen) {
    return PgHelpers::toString(totalMatchLength)
        + " (" + PgHelpers::toString((totalMatchLength * 100.0) / destLen, 1) + "%)";
}

size_t resCount;

void processMatches(vector<PgTools::TextMatch>& textMatches, char *destPtr, size_t destLen) {
    resCount += textMatches.size();
    uint32_t pos = 0;
    uint32_t nPos = 0;
    uint32_t totalDestOverlap = 0;
    uint32_t totalMatched = 0;
    for (TextMatch &match: textMatches) {
        if (match.posDestText < pos) {
            uint32_t overflow = pos - match.posDestText;
            if (overflow >= match.length) {
                totalDestOverlap += match.length;
                match.length = 0;
                continue;
            }
            totalDestOverlap += overflow;
            match.length -= overflow;
            match.posDestText += overflow;
            match.posSrcText += overflow;
        }
        if (match.length < k) {
            totalDestOverlap += match.length;
            continue;
        }
        totalMatched += match.length;
        uint64_t length = match.posDestText - pos;
//        memcpy(destPtr + nPos, destPtr + pos, length);
        nPos += length;
        //destPg[nPos++] = MATCH_MARK;
        //PgHelpers::writeValue<uint32_t>(pgMapOffDest, match.posSrcText);
        //PgHelpers::writeUIntByteFrugal(pgMapLenDest, match.length - minMatchLength);
        pos = match.endPosDestText();
    }
    uint64_t length = destLen - pos;
    //memcpy(destPtr + nPos, destPtr + pos, length);
    nPos += length;
    //destPg.resize(nPos);

    textMatches.clear();
//    resPgMapOff = pgMapOffDest.str();
//    pgMapOffDest.clear();
//    resPgMapLen = pgMapLenDest.str();
//    pgMapLenDest.clear();

//    *logout << "Preparing output time: " << time_millis(post_start) << " msec." << endl;
//    cout << "Final size of Pg: " << nPos << " (removed: " <<
//         getTotalMatchStat(totalMatched, destLen) << "; " << totalDestOverlap << " chars in overlapped dest symbol)" << endl;
    unmatchedCharsAll += nPos;
    totalMatchedAll += totalMatched;
    totalDestOverlapAll += totalDestOverlap;
    totalDestLenAll += destLen;
}


int main(int argc, char *argv[]) {

    gzFile fp;
    kseq_t *seq;
    int l;
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <k> <sequences.list.file>\n", argv[0]);
        return 1;
    }

    k = atoi(argv[1]);
    string seqListFileName = argv[2];

    ifstream listSrc(seqListFileName, ios_base::in | ios_base::binary);
    if (listSrc.fail()) {
        fprintf(stderr, "cannot open sequences list file %s\n", seqListFileName.c_str());
        exit(EXIT_FAILURE);
    }

    PgHelpers::numberOfThreads = 1;
    if (PgHelpers::numberOfThreads <= 0) {
        fprintf(stderr, "The number of threads must be positive.\n");
        exit(EXIT_FAILURE);
    }
    omp_set_num_threads(PgHelpers::numberOfThreads);

    PgHelpers::time_checkpoint();

    string seqName;
    std::getline(listSrc, seqName);
    if(!seqName.empty() && *seqName.rbegin() == '\r') seqName.resize(seqName.length() - 1);

    fp = gzopen(seqName.c_str(), "r");
    seq = kseq_init(fp);
    loadRef(seq);
    matcher = new CopMEMMatcher(ref.data(), ref.size(), k, false);
    //build_dict();
    cout << "build - " << PgHelpers::time_millis() << " [ms]" << endl;
    kseq_destroy(seq);
    gzclose(fp);

    vector<string> fnas;
    vector<PgTools::TextMatch> resMatches;
    while (std::getline(listSrc, seqName)) {
        if(!seqName.empty() && *seqName.rbegin() == '\r') seqName.resize(seqName.length() - 1);
        fp = gzopen(seqName.c_str(), "r");
        seq = kseq_init(fp);
        while ((l = kseq_read(seq)) >= 0) {
//            fnas.push_back(string(seq->seq.s, seq->seq.l));
            matcher->matchTexts(resMatches, seq->seq.s, seq->seq.l, false, false, k);
            processMatches(resMatches,  seq->seq.s, seq->seq.l);
//            printf(">%s", seq->name.s);
//            if (seq->comment.l) printf(" %s", seq->comment.s);
//            printf("\n%s\n", seq->seq.s);
//            if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        }
        kseq_destroy(seq);
        gzclose(fp);
    }
/*    cout << "targets loaded - " << PgHelpers::time_millis() << " [ms]" << endl;
    for(string fna : fnas) {
        matcher->matchTexts(resMatches, fna, false, false, k);
        processMatches(resMatches, (char*) fna.data(), fna.length());
    }*/
    cout << "reference length: " << rcStart << endl;
    cout << "targets total length: " << totalDestLenAll << endl;
    cout << "exact matches total: " << resCount << endl;
    cout << "Final unmatched chars: " << unmatchedCharsAll << " (removed: " <<
         getTotalMatchStat(totalMatchedAll, totalDestLenAll) << "; " << totalDestOverlapAll << " chars in overlapped target symbol)" << endl;

    cout << "matching finished - " << PgHelpers::time_millis() << " [ms]" << endl;
    return 0;
}
