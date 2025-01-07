#ifndef COMB_HPP
#define COMB_HPP

class comb {
    using ui = unsigned;
private://combinatorial number
    ui maxK = 30;
    ui maxSize = 30000;
    double ** CN = nullptr, *bf3 = nullptr;
    void computeC() {
        ui maxD = maxSize, maxD2 = maxK;
        CN = new double*[maxD];
        bf3 = new double[maxD2 * maxD];
        for(int i = 0; i < maxD; i++) {
            CN[i] = bf3 + i * maxD2;
        }
        CN[0][0] = 1;
        CN[1][0] = 1;
        CN[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            CN[i][0] = 1;
            if(i < maxD2) CN[i][i] = 1;
            for(int j = 1; j < i && j < maxD2; j++) {
                CN[i][j] = CN[i - 1][j - 1] + CN[i - 1][j];
            }
        }
    }
    

public:
    comb(ui k, ui sz) {
        setMaxKAndSize(k, sz);
        computeC();
    }
    comb() {
        computeC();
    }

    void setMaxKAndSize(ui k, ui sz) {
        maxK = k; maxSize = sz;
    }
    void reComputeC() {
        if(bf3 != nullptr) {
            delete [] bf3; delete [] CN; 
        }
        computeC();
    }

    ~comb() {
        if(bf3 != nullptr) {
            delete [] bf3; delete [] CN; 
        }
    }

    double operator () (ui n, ui m) const {
        return CN[n][m];
    }
};



#endif