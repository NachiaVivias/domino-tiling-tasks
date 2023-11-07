
#define NDEBUG
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <unordered_map>
using namespace std;

using CountingInt = uint32_t;
using MemoInt = uint8_t;
static constexpr const size_t GCOLS = 6;
static constexpr const size_t GROWS_U = 3;
static constexpr const size_t GROWS_L = 4;
static constexpr const size_t GROWS = GROWS_U + GROWS_L;
static constexpr const size_t EXPGCOLS = (size_t)1 << GCOLS;
static constexpr const size_t OUTPUT_STYLE = 2;

template<size_t ROWS>
using MemoStructure = array<MemoInt, ROWS>;

template<size_t COLS>
struct RowScanline {
    static constexpr const size_t EXPCOLS = (size_t)1 << COLS;
    array<CountingInt, EXPCOLS> state;
    RowScanline(){
        state.fill((CountingInt)0);
        state.back() = 1;
    }
    void scanARow() {
        for(size_t i=0; i+1<COLS; i++){
            const size_t mask = (size_t)3 << i;
            for(size_t f=mask; f<EXPCOLS; f=(f+1)|mask){
                state[f] += state[f-mask];
            }
        }
    }
    RowScanline<COLS-1> reduceAColumn(size_t pos) const {
        RowScanline<COLS-1> res;
        res.state.back() = 0;
        constexpr size_t x = res.EXPCOLS;
        size_t lobits = ((size_t)1 << pos) - 1;
        size_t hibits = x - 1 - lobits;
        for(size_t q=(size_t)1 << pos; q<EXPCOLS; q = (q+1) | ((size_t)1 << pos)){
            res.state[(q & lobits) | ((q >> 1) & hibits)] = state[q];
        }
        return res;
    }
};

template<size_t COLS>
bool operator<(const RowScanline<COLS>& l, const RowScanline<COLS>& r){
    return l.state < r.state;
}

template<size_t ROWS, size_t COLS, size_t REDUCED_INTO = COLS>
struct DpState {
    static constexpr const size_t SCAN_LINE_SIZE = (size_t)1 << REDUCED_INTO;
    using StateValueList = vector<pair<RowScanline<REDUCED_INTO>, MemoStructure<ROWS>>>;
    StateValueList A;
    static void unique(StateValueList& B){
        size_t p = 0;
        for(size_t i=1; i<B.size(); i++){
            if(B[p].first < B[i].first) if(++p != i) B[p] = B[i];
        }
        B.resize(p + 1);
    }
    DpState(){
        A.emplace_back();
        A.back().second.fill((MemoInt)0);
    }
    DpState(int){
    }
    void scanARow(){
        for(auto& a : A) a.first.scanARow();
    }
    void scanAMass(size_t row, size_t col){
        StateValueList B;
        MemoInt maskMemo = (MemoInt)1 << col;
        size_t mask = (size_t)1 << col;
        for(auto& a : A){
            {
                B.push_back(a);
                auto& b = B.back();
                for(size_t f=mask; f<SCAN_LINE_SIZE; f=(f+1)|mask){
                    b.first.state[f - mask] = 0;
                }
            }
            {
                B.push_back(a);
                auto& b = B.back();
                for(size_t f=mask; f<SCAN_LINE_SIZE; f=(f+1)|mask){
                    swap(b.first.state[f], b.first.state[f - mask]);
                }
                b.second[row] |= maskMemo;
            }
        }
        sort(B.begin(), B.end());
        unique(B);
        swap(A, B);
    }
    DpState<ROWS, COLS, REDUCED_INTO-1> reduceAColumn(size_t pos) const {
        using NextState = DpState<ROWS, COLS, REDUCED_INTO-1>;
        typename NextState::StateValueList B;
        for(auto& a : A) B.push_back({ a.first.reduceAColumn(pos), a.second });
        sort(B.begin(), B.end());
        NextState::unique(B);
        NextState res; res.A = move(B);
        return res;
    }
};

struct DpMergeResult {
    array<MemoInt, GROWS> grid;
};

bool operator<(const DpMergeResult& l, const DpMergeResult& r){
    return l.grid < r.grid;
}

template<size_t COLS, size_t ITERATION = 0, size_t REDUCED_INTO = COLS>
unordered_map<size_t, DpMergeResult> MeetInTheMiddle(
    const DpState<GROWS_U, COLS, REDUCED_INTO> &dpStateUpper,
    const DpState<GROWS_L, COLS, REDUCED_INTO> &dpStateLower
){
    static constexpr size_t SCAN_LINE_SIZE = dpStateUpper.SCAN_LINE_SIZE;
    cerr << "MeetInTheMiddle : " << endl;
    cerr << "   COLS = " << COLS << " , ITERATION = " << ITERATION << " , REDUCED_INTO = " << REDUCED_INTO << endl;
    cerr << "   upper size = " << dpStateUpper.A.size() << " , lower size = " << dpStateLower.A.size() << endl;
    if(dpStateUpper.A.empty()) return {};
    if(dpStateLower.A.empty()) return {};

    if constexpr(ITERATION == COLS){

        auto dpOutputLower = move(dpStateLower.A);
        sort(dpOutputLower.begin(), dpOutputLower.end(), [](auto l, auto r){ return l.second < r.second; });
        vector<size_t> lowerDominoCount;
        for(auto& a : dpOutputLower){
            size_t cnt = 0;
            for(size_t r=0; r<GROWS_L; r++) for(size_t c=0; c<GCOLS; c++){
                cnt += (a.second[r] >> c) & 1;
            }
            lowerDominoCount.push_back(cnt);
        }

        auto dpOutputUpper = move(dpStateUpper.A);
        sort(dpOutputUpper.begin(), dpOutputUpper.end(), [](auto l, auto r){ return l.second < r.second; });
        vector<size_t> upperDominoCount;
        for(auto& a : dpOutputUpper){
            size_t cnt = 0;
            for(size_t r=0; r<GROWS_U; r++) for(size_t c=0; c<GCOLS; c++){
                cnt += (a.second[r] >> c) & 1;
            }
            upperDominoCount.push_back(cnt);
        }

        unordered_map<size_t, DpMergeResult> ansEnum;
        for(size_t loi=0; loi<dpOutputLower.size(); loi++){
            for(size_t upi=0; upi<dpOutputUpper.size(); upi++){
                if((lowerDominoCount[loi] + upperDominoCount[upi]) % 2 != 0) continue;
                auto& lo = dpOutputLower[loi];
                auto& up = dpOutputUpper[upi];
                size_t f = 0;
                for(size_t q=0; q<SCAN_LINE_SIZE; q++){
                    f += (size_t)lo.first.state[q] * (size_t)up.first.state[q];
                }
                if(ansEnum.count(f)) continue;
                array<MemoInt, GROWS> resGrid;
                for(size_t r=0; r<GROWS_L; r++) resGrid[r] = lo.second[r];
                for(size_t r=0; r<GROWS_U; r++) resGrid[GROWS_L+r] = up.second[r];
                ansEnum[f].grid = resGrid;
            }
        }

        return ansEnum;
    }
    else{
        static constexpr size_t TARGET_BIT = ITERATION - (COLS - REDUCED_INTO);

        DpState<GROWS_U, COLS, REDUCED_INTO> dpStateUpper0(0);
        DpState<GROWS_U, COLS, REDUCED_INTO> dpStateUpper1(0);
        for(auto& a : dpStateUpper.A){
            ((((size_t)a.second[0] >> ITERATION) & 1) ? dpStateUpper1 : dpStateUpper0)
                .A.push_back(a);
        }
        
        DpState<GROWS_L, COLS, REDUCED_INTO> dpStateLower0(0);
        DpState<GROWS_L, COLS, REDUCED_INTO> dpStateLower1(0);
        for(auto& a : dpStateLower.A){
            ((((size_t)a.second[GROWS_L-1] >> ITERATION) & 1) ? dpStateLower1 : dpStateLower0)
                .A.push_back(a);
        }

        auto ans = MeetInTheMiddle<COLS, ITERATION+1, REDUCED_INTO-1>(
            dpStateUpper0.reduceAColumn(TARGET_BIT),
            dpStateLower.reduceAColumn(TARGET_BIT)
        );

        auto mergeIntoAns = [&ans](const unordered_map<size_t, DpMergeResult>& src){
            for(auto& f : src){
                if(ans.count(f.first)){
                    ans[f.first] = min(ans[f.first], f.second);
                } else {
                    ans[f.first] = f.second;
                }
            }
        };
        mergeIntoAns(
            MeetInTheMiddle<COLS, ITERATION+1, REDUCED_INTO-1>(
                dpStateUpper1.reduceAColumn(TARGET_BIT),
                dpStateLower0.reduceAColumn(TARGET_BIT)));
        mergeIntoAns(
            MeetInTheMiddle<COLS, ITERATION+1, REDUCED_INTO>(
                dpStateUpper1,
                dpStateLower1));

        return ans;
    }
}

int main(){

    DpState<GROWS_L, GCOLS> dpStateLower;
    for(size_t r=0; r<GROWS_L; r++){
        for(size_t c=0; c<GCOLS; c++){
            cerr << "r = " << r << " , c = " << c << endl;
            dpStateLower.scanAMass(r, c);
        }
        dpStateLower.scanARow();
    }
    cerr << "Lower # state = " << dpStateLower.A.size() << endl;

    DpState<GROWS_U, GCOLS> dpStateUpper;
    for(size_t r=GROWS_U-1; r<GROWS_U; r--){
        for(size_t c=0; c<GCOLS; c++){
            cerr << "r = " << r << " , c = " << c << endl;
            dpStateUpper.scanAMass(r, c);
        }
        dpStateUpper.scanARow();
    }
    cerr << "Upper # state = " << dpStateUpper.A.size() << endl;

    vector<pair<size_t, DpMergeResult>> ansEnum;
    {
        auto buf = MeetInTheMiddle(dpStateUpper, dpStateLower);
        array<MemoInt, GROWS> zeroMatch;
        zeroMatch.fill(0);
        zeroMatch[GROWS-1] = 1;
        buf[0].grid = zeroMatch;
        ansEnum.assign(buf.begin(), buf.end());
        sort(ansEnum.begin(), ansEnum.end(), [&](auto l, auto r){ return l.first < r.first; });
    }

    if(OUTPUT_STYLE == 0){
        bool comma = false;
        cout << "{\"counts\":{";
        for(auto ans : ansEnum){
            if(comma) cout << ",";
            cout << "\"" << ans.first << "\":[";
            for(size_t r=0; r<GROWS; r++){
                cout << (size_t)ans.second.grid[GROWS-1-r];
                if(r != GROWS-1) cout << ",";
            }
            cout << "]";
            comma = true;
        }
        cout << "}}" << endl;
    }
    
    if(OUTPUT_STYLE == 1){
        bool comma = false;
        cout << "{\"counts\":[";
        for(auto ans : ansEnum){
            if(comma) cout << ",";
            cout << ans.first;
            comma = true;
        }
        cout << "]}" << endl;
    }
    
    if(OUTPUT_STYLE == 2){
        cout << "BOARD : " << min(GCOLS, GROWS) << " x " << max(GCOLS, GROWS) << endl;
        cout << "    (ROWS_U, ROWS_L, COLS) : (" << GROWS_U << ", " << GROWS_L << ", " << GCOLS << ")" << endl;
        cout << "    COUNT : " << ansEnum.size() << endl;
        vector<size_t> excludedScores;
        for(size_t i=0; excludedScores.size()<20; i++){
            size_t targetIdx = i - excludedScores.size();
            if(targetIdx < ansEnum.size() && ansEnum[targetIdx].first == i) continue;
            excludedScores.push_back(i);
        }
        cout << "    MIN 20 EXCLUDED :";
        for(auto a : excludedScores){ cout << " " << a; }
        cout << endl;
        cout << "    MAX SCORE : " << ansEnum.back().first << endl;
    }
    return 0;
}
