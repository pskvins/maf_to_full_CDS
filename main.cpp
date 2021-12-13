#include <iostream>
#include <vector>
#include <fstream>
#include <string>

std::string reverse_base(std::string seq) {
    std::string reversed = "";
    for (int i = 0; i < seq.length(); i++) {
        if (seq[i] == 'A') {
            reversed += "T";
        } else if (seq[i] == 'T') {
            reversed += "A";
        } else if (seq[i] == 'G') {
            reversed += "C";
        } else if (seq[i] == 'C') {
            reversed += "G";
        } else if (seq[i] == '-') {
            reversed += "-";
        } else if (seq[i] == 'N') {
            reversed += "N";
        } else {
            std::cout << "weird character in base : " << seq[i] << std::endl;
            exit(-1);
        }
    }
    return reversed;
}

char reverse_strand(char strand) {
    if (strand == '+') {
        return '-';
    } else if (strand == '-') {
        return '+';
    } else {
        std::cout << "strand is weird" << std::endl;
        exit(-1);
    }
}

int main() {
    std::string species[7] = {"sacCer3", "sacPar", "sacMik", "sacKud", "sacBay", "sacCas", "sacKlu"};
    std::string alignment[7];
    std::string gap = "";
    for (int i = 0; i < 25000000; i++) {
        gap += "-";
    }
    for (int i = 0; i < 7; i++) {
        alignment[i] = gap;
    }
    std::string line;
    int start;
    int end;
    char strand;
    size_t pos_sec = 0;

    std::ifstream maf;
    maf.open("/Users/sukhwanpark/Downloads/yeast_maf/chrM_tab_full.maf");
    getline(maf, line);
    std::string align_temp[7];
    size_t pos;
    std::ofstream maf_new;
    maf_new.open("/Users/sukhwanpark/Downloads/yeast_maf/chrM_CDS_full.maf");
    while (!line.empty()) {
        if (line[0] == 'a') {
            getline(maf, line);
            //first 's'
            pos = line.find('\t');
            pos_sec = line.find('\t', pos + 1);
            line.erase(0, pos_sec + 1);
            pos_sec = line.find('\t');
            start = std::stoi(line.substr(0, pos_sec));
            line.erase(0, pos_sec + 1);
            pos_sec = line.find('\t');
            line.erase(0, pos_sec + 1);
            strand = line[0];
            line.erase(0, 2);
            pos_sec = line.find('\t');
            line.erase(0, pos_sec + 1);
            for (int i = 0; i < line.length(); i++) {
                line[i] = std::toupper(line[i]);
            }
            align_temp[0] = line;
            getline(maf, line);

            int id = 1;
            while (line[0] == 's') {
                pos = line.find('\t');
                pos_sec = line.find('\t', pos + 1);
                line.erase(0, pos_sec + 1);
                pos_sec = line.find('\t');
                line.erase(0, pos_sec + 1);
                pos_sec = line.find('\t');
                line.erase(0, pos_sec + 1);
                line.erase(0, 2);
                pos_sec = line.find('\t');
                line.erase(0, pos_sec + 1);
                for (int i = 0; i < line.length(); i++) {
                    line[i] = std::toupper(line[i]);
                }
                align_temp[id] = line;
                getline(maf, line);
                id++;
            }
            while (align_temp[0].find('-') != std::string::npos) {
                size_t gap_pos = align_temp[0].find('-');
                for (int i = 0; i < 7; i++) {
                    align_temp[i].erase(gap_pos, 1);
                }
            }
            for (int i = 0; i < 7; i++) {
                for (int base = start, j = 0; base < start + align_temp[0].length(); base++, j++) {
                    alignment[i][base] = align_temp[i][j];
                }
            }
        } else {
            getline(maf, line);
        }
    }

    std::ifstream coor;
    coor.open("/Users/sukhwanpark/Downloads/chrM.txt");
    getline(coor, line);
    while (!line.empty()) {
        pos_sec = line.find('\t');
        start = std::stoi(line.substr(0, pos_sec)) - 1;
        line.erase(0, pos_sec + 1);
        pos_sec = line.find('\t');
        end = std::stoi(line.substr(0, pos_sec));
        line.erase(0, pos_sec + 1);
        strand = line[0];
        if (start >= end) {
            std::cout << "start & end coordinate is wierd" << std::endl;
            exit(-1);
        }
        if (strand == '+') {
            for (int i = 0; i < 7; i++) {
                align_temp[i] = alignment[i].substr(start, end - start);
            }
        } else {
            for (int i = 0; i < 7; i++) {
                align_temp[i] = reverse_base(alignment[i].substr(start, end - start));
            }
        }
        maf_new << "a\n";
        for (int i = 0; i < 7; i++) {
            maf_new << "s\t" << species[i] << ".chr\t" << std::to_string(start) << "\t" << std::to_string(end - start)
            << "\t" << strand << "\t1\t" << align_temp[i] << "\n";
        }
        getline(coor, line);
    }

    maf_new.close();
    maf.close();
    coor.close();
    return 0;
}
