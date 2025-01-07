#ifndef FASTIO_HPP
#define FASTIO_HPP

#include <cstdio>
#include <string>
#include "filesystem.hpp"
#include <cstdint>

constexpr int PAGESIZE = 1024*10000;


class fastIO {
private:
    char *buffer = nullptr, *S = nullptr, *T = nullptr;
    FILE * f;
	unsigned long bytes = 0;
    int sz;

public: 
    fastIO(FILE * f_) {
        f = f_;
        buffer = new char[PAGESIZE];
        sz = sizeof(char);
    }

    fastIO(const std::string &file, const std::string openStyle) {
        f = fopen(file.c_str(), openStyle.c_str());
        T = S = buffer = new char[PAGESIZE];
		if(openStyle[0] == 'r') {
			bytes = file_size(file);
		}
        sz = sizeof(char);
    }

    long load(unsigned * bf) {
        return fread(bf, 4, bytes / 4, f);
    }

    unsigned long leftBytes() {
        return bytes;
    }

    ~fastIO() {
        fclose(f);
        if(buffer != nullptr) {
            delete [] buffer;
            buffer = nullptr;
        }
    }

    char getChar()  
    {  
        if(S==T)  
        {   
            if(bytes == 0) return EOF;
			if(bytes < PAGESIZE) {
				T=(S=buffer) + fread(buffer, 1, bytes, f);
// printf("%c%c%c %d\n", S[0], S[1], S[2], T-S);
                bytes = 0;
			}
            else {
				T=(S=buffer) + fread(buffer, 1, PAGESIZE, f);
				bytes -= PAGESIZE;
			}
            if(S==T) return EOF;  
        }
// putchar(*S);
// putchar('\n');
        return *S++;  
    }

    bool empty() {
        if(S == T) {
            if(bytes == 0) return true;
        }
        return false;
    }
    
    void getInt(int & re)
    {
        char c;   
        for(c=getChar();c<'0'||c>'9';c=getChar());  
        while(c>='0'&&c<='9')  
            re=(re<<1)+(re<<3)+(c-'0'),c=getChar();
    }
    void getUInt(unsigned & re)
    {
        char c;
        re = 0;
        for(c=getChar();(c<'0'||c>'9') && c!=EOF;c=getChar());  
        while(c>='0'&&c<='9')  {
            re=(re<<1)+(re<<3)+(c-'0');
            if((c=getChar()) == EOF) break;
        }
    }
    void getUllInt(unsigned long long & re)
    {
        char c;
        re = 0;
        for(c=getChar();(c<'0'||c>'9') && c!=EOF;c=getChar());  
        while(c>='0'&&c<='9')  {
            re=(re<<1)+(re<<3)+(c-'0');
            if((c=getChar()) == EOF) break;
        }
    }

    void jump()
    {
        char c;
        unsigned re = 0;
        for(c=getChar();c<'0'||c>'9';c=getChar());  
        while(c>='0'&&c<='9')  
            re=(re<<1)+(re<<3)+(c-'0'),c=getChar();
        if(c == '.') {
            c = getChar();
            while(c>='0'&&c<='9')  
                re=(re<<1)+(re<<3)+(c-'0'),c=getChar();
        }
        
        // return re;
        // while((c = getChar()) != ' ' || c != '\n' || c != '\t')
        //     ;
    }

    void seek(uint32_t b) {
        assert(fseek(f, b, SEEK_SET) == 0);
    }
};

#endif