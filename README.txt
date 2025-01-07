# compile
    cmake -S ./
    make

We get bin/mce and bin/kcc now.

# data format and command

## each edge per line

    bin/mce noUVM -f_txt dataFile
    bin/kcc noUVM -f_txt dataFile -k k

## The first line is the number of nodes and edges, then each edge per line

    bin/mce -f_txt dataFile
    bin/kcc -f_txt dataFile -k k


