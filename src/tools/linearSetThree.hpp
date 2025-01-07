#ifndef LINEARSETTHREE_HPP
#define LINEARSETTHREE_HPP

class linearSetThree {
private:
    uint32_t * vSet = nullptr;
    uint32_t * pSet = nullptr;
    uint32_t sz = 0;

public:
    uint32_t x, c, r;

    linearSetThree() {}
    linearSetThree(uint32_t size):sz(size), c(0), x(0), r(size) {
        vSet = new uint32_t[size];
        pSet = new uint32_t[size];

        for(uint32_t i = 0; i < size; i++) {
            vSet[i] = pSet[i] = i;
        }
    }
    void resize(uint32_t size) {
        sz = size;
        c = 0;
        x = 0;
        r = size;
        vSet = new uint32_t[size];
        pSet = new uint32_t[size];

        for(uint32_t i = 0; i < size; i++) {
            vSet[i] = pSet[i] = i;
        }
    }

    ~linearSetThree() {
        delete [] vSet;
        delete [] pSet;
    }

    void swapByPos(uint32_t i, uint32_t j) {
        std::swap(vSet[i], vSet[j]);
        std::swap(pSet[vSet[i]], pSet[vSet[j]]);
    }

    uint32_t operator [] (uint32_t i) {
        // if(i >= sz) {
        //     printf("%u %u\n", i, sz);fflush(stdout);
        //     exit(-1);
        // }
        return  vSet[i];
    }

    uint32_t pos(uint32_t v) {
        return pSet[v];
    }

    bool CIsEmpty() { return c == r; }
    bool XIsEmpty() { return x == c; }

    uint32_t cSize() {
        return r - c;
    }

    uint32_t * begin() { return vSet; }

    void print() {
        printf("X:");
        for(uint32_t i = x; i < c; i++) {
            printf("%u ", vSet[i]);
        }
        printf("C:");
        for(uint32_t i = c; i < r; i++) {
            printf("%u ", vSet[i]);
        }
        printf("R:");
        for(uint32_t i = r; i < sz; i++) {
            printf("%u ", vSet[i]);
        }
    }
};

#endif