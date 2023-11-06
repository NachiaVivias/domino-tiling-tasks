
#define NDEBUG
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <bitset>
using namespace std;

using CountingInt = uint32_t;
using MemoInt = uint8_t;
static constexpr const size_t COLS = 7;
static constexpr const size_t ROWS = 3;
static constexpr const size_t EXPCOLS = (size_t)1 << COLS;
using MemoStructure = array<MemoInt, ROWS>;

struct RowScanline {
    array<CountingInt, EXPCOLS> state;
    RowScanline(){
        state.fill((CountingInt)0);
        state[EXPCOLS - 1] = 1;
    }
    void scanARow() {
        for(size_t i=0; i+1<COLS; i++){
            const size_t mask = (size_t)3 << i;
            for(size_t f=mask; f<EXPCOLS; f=(f+1)|mask){
                state[f] += state[f-mask];
            }
        }
    }
};

bool operator<(const RowScanline& l, const RowScanline& r){
    return l.state < r.state;
}

struct DpState {
    vector<pair<RowScanline, MemoStructure>> A;
    DpState(){
        A.emplace_back();
        A.back().second.fill((MemoInt)0);
    }
    void scanARow(){
        for(auto& a : A) a.first.scanARow();
        // sort(A.begin(), A.end());
    }
    void scanAMass(size_t row, size_t col){
        vector<pair<RowScanline, MemoStructure>> B;
        MemoInt maskMemo = (MemoInt)1 << col;
        size_t mask = (size_t)1 << col;
        for(auto& a : A){
            {
                B.push_back(a);
                auto& b = B.back();
                for(size_t f=mask; f<EXPCOLS; f=(f+1)|mask){
                    b.first.state[f - mask] = 0;
                }
            }
            {
                B.push_back(a);
                auto& b = B.back();
                for(size_t f=mask; f<EXPCOLS; f=(f+1)|mask){
                    swap(b.first.state[f], b.first.state[f - mask]);
                }
                b.second[row] |= maskMemo;
            }
        }
        sort(B.begin(), B.end());
        size_t p = 0;
        for(size_t i=1; i<B.size(); i++){
            if(B[p].first < B[i].first) if(++p != i) B[p] = B[i];
        }
        B.resize(p + 1);
        swap(A, B);
    }
};

int main(){

    DpState dpStateLower;
    for(size_t r=0; r<ROWS; r++){
        for(size_t c=0; c<COLS; c++){
            cerr << "r = " << r << " , c = " << c << endl;
            dpStateLower.scanAMass(r, c);
        }
        dpStateLower.scanARow();
    }
    cerr << "Lower # state = " << dpStateLower.A.size() << endl;

    DpState dpStateUpper;
    for(size_t r=ROWS-1; r<ROWS; r--){
        for(size_t c=0; c<COLS; c++){
            cerr << "r = " << r << " , c = " << c << endl;
            dpStateUpper.scanAMass(r, c);
        }
        dpStateUpper.scanARow();
    }
    cerr << "Upper # state = " << dpStateUpper.A.size() << endl;

    auto dpOutputLower = move(dpStateLower.A);
    sort(dpOutputLower.begin(), dpOutputLower.end(), [](auto l, auto r){ return l.second < r.second; });
    vector<size_t> lowerDominoCount;
    for(auto& a : dpOutputLower){
        size_t cnt = 0;
        for(size_t r=0; r<ROWS; r++) for(size_t c=0; c<COLS; c++){
            cnt += (a.second[r] >> c) & 1;
        }
        lowerDominoCount.push_back(cnt);
    }

    auto dpOutputUpper = move(dpStateUpper.A);
    sort(dpOutputUpper.begin(), dpOutputUpper.end(), [](auto l, auto r){ return l.second < r.second; });
    vector<size_t> upperDominoCount;
    for(auto& a : dpOutputUpper){
        size_t cnt = 0;
        for(size_t r=0; r<ROWS; r++) for(size_t c=0; c<COLS; c++){
            cnt += (a.second[r] >> c) & 1;
        }
        upperDominoCount.push_back(cnt);
    }

    map<size_t, array<MemoInt, 6>> ansEnum;
    {
        array<MemoInt, 6> zeroMatch;
        zeroMatch.fill(0);
        zeroMatch[0] = 1;
        ansEnum[0] = zeroMatch;
    }

    for(size_t loi=0; loi<dpOutputLower.size(); loi++){
        for(size_t upi=0; upi<dpOutputUpper.size(); upi++){
            if((lowerDominoCount[loi] + upperDominoCount[upi]) % 2 != 0) continue;
            auto& lo = dpOutputLower[loi];
            auto& up = dpOutputUpper[upi];
            size_t border = (EXPCOLS - 1) - ((size_t)lo.second[2] & (size_t)up.second[0]);
            size_t f = 0;
            for(size_t q = border; q<EXPCOLS; q=(q+1)|border){
                f += (size_t)lo.first.state[q] * (size_t)up.first.state[q];
            }
            if(ansEnum.count(f)) continue;
            array<MemoInt, 6> resGrid;
            for(size_t r=0; r<3; r++) resGrid[5-r] = lo.second[r];
            for(size_t r=0; r<3; r++) resGrid[2-r] = up.second[r];
            ansEnum[f] = resGrid;
        }
    }

    bool comma = false;
    cout << "{\"counts\":{";
    for(auto ans : ansEnum){
        if(comma) cout << ",";
        cout << "\"" << ans.first << "\":[";
        for(size_t r=0; r<6; r++){
            cout << (size_t)ans.second[r];
            if(r != 5) cout << ",";
        }
        cout << "]";
        comma = true;
    }
    cout << "}}" << endl;
    return 0;
}
