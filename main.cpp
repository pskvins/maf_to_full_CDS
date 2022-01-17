#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <omp.h>
#define THREAD_NUM 2

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
            std::cout << "weird character in base : ";
            for (int j = 0; j < 95; j++) {
                std::cout << seq[j];
            }
            std::cout << std::endl;
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
    std::string species[100] = {"hg38", "panTro4", "gorGor3", "ponAbe2", "nomLeu3", "rheMac3", "macFas5", "papAnu2", "chlSab2", "calJac3", "saiBol1", "otoGar3", "tupChi1", "speTri2", "jacJac1", "micOch1", "criGri1", "mesAur1", "mm10", "rn6", "hetGla2", "cavPor3", "chiLan1", "octDeg1", "oryCun2", "ochPri3", "susScr3", "vicPac2", "camFer1", "turTru2", "orcOrc1", "panHod1", "bosTau8", "oviAri3", "capHir1", "equCab2", "cerSim1", "felCat8", "canFam3", "musFur1", "ailMel1", "odoRosDiv1", "lepWed1", "pteAle1", "pteVam1", "eptFus1", "myoDav1", "myoLuc2", "eriEur2", "sorAra2", "conCri1", "loxAfr3", "eleEdw1", "triMan1", "chrAsi1", "echTel2", "oryAfe1", "dasNov3", "monDom5", "sarHar1", "macEug2", "ornAna1", "colLiv1", "falChe1", "falPer1", "ficAlb2", "zonAlb1", "geoFor1", "taeGut2", "pseHum1", "melUnd1", "amaVit1", "araMac1", "anaPla1", "galGal4", "melGal1", "allMis1", "cheMyd1", "chrPic2", "pelSin1", "apaSpi1", "anoCar2", "xenTro7", "latCha1", "tetNig2", "fr3", "takFla1", "oreNil2", "neoBri1", "hapBur1", "mayZeb1", "punNye1", "oryLat2", "xipMac1", "gasAcu1", "gadMor1", "danRer10", "astMex1", "lepOcu1", "petMar2"};
    int species_size = sizeof(species) / sizeof(species[0]);
    std::string *alignment = new std::string[species_size];

    std::string directory = "/Users/sukhwanpark/Downloads/compare_100vertebrates/test/";
    std::string list = directory + "list.txt";
    std::vector<std::string> all_lines;
    std::ifstream lists;
    lists.open(list);
    std::string line_tmp;
    getline(lists, line_tmp);
    size_t pos;
    while(!line_tmp.empty()) {
        all_lines.push_back(line_tmp);
        getline(lists, line_tmp);
    }
#pragma omp parallel for
    for (int thread = 0; thread < all_lines.size(); thread++) {
        pos = all_lines[thread].find(' ');
        std::string chr = all_lines[thread].substr(0, pos);
        all_lines[thread].erase(0, pos);
        int length = std::stoi(all_lines[thread]);

        //now process other files
        std::string gap = "";
        for (int i = 0; i < length; i++) {
            gap += "-";
        }
        for (int i = 0; i < species_size; i++) {
            alignment[i] = gap;
        }
        int start;
        int end;
        char strand;
        size_t pos_sec = 0;

        std::string input_directory = directory + chr + "_space.maf";
        std::ifstream maf;
        maf.open(input_directory);
        std::string line;
        getline(maf, line);
        std::string *align_temp = new std::string[species_size];
        std::string name;

        std::string out_directory = directory + chr + "_CDS.maf";
        std::ofstream maf_new;
        maf_new.open(out_directory);
        while (!line.empty()) {
            if (line[0] == 'a') {
                for (int i = 0; i < species_size; i++) {
                    align_temp[i] = "";
                }
                getline(maf, line);
                //first 's'
                pos = line.find(' ');
                pos_sec = line.find(' ', pos + 1);
                line.erase(0, pos_sec + 1);
                pos_sec = line.find(' ');
                start = std::stoi(line.substr(0, pos_sec));
                line.erase(0, pos_sec + 1);
                pos_sec = line.find(' ');
                line.erase(0, pos_sec + 1);
                strand = line[0];
                line.erase(0, 2);
                pos_sec = line.find(' ');
                line.erase(0, pos_sec + 1);
                for (int i = 0; i < line.length(); i++) {
                    line[i] = std::toupper(line[i]);
                }
                align_temp[0] = line;
                getline(maf, line);

                int id = 1;
                while (line[0] == 's') {
                    pos = line.find(' ');
                    pos_sec = line.find('.', pos + 1);
                    name = line.substr(pos + 1, pos_sec - pos - 1);
                    line.erase(0, pos_sec + 1);
                    pos_sec = line.find_last_of(' ');
                    line.erase(0, pos_sec + 1);
                    for (int i = 0; i < line.length(); i++) {
                        line[i] = std::toupper(line[i]);
                    }
                    int index = 1;
                    while (index < species_size) {
                        if (species[index] != name) {
                            index++;
                        } else {
                            break;
                        }
                    }
                    align_temp[index] = line;
                    getline(maf, line);
                    id++;
                }

                int align_length = align_temp[0].length();
                std::string temp = "";
                for (int i = 0; i < align_length; i++) {
                    temp += '-';
                }
                for (int i = 0; i < species_size; i++) {
                    if (align_temp[i] == "") {
                        align_temp[i] = temp;
                    }
                }
                while (align_temp[0].find('-') != std::string::npos) {
                    size_t gap_pos = align_temp[0].find('-');
                    for (int i = 0; i < species_size; i++) {
                        align_temp[i].erase(gap_pos, 1);
                    }
                }
                for (int i = 0; i < species_size; i++) {
                    for (int base = start, j = 0; base < start + align_temp[0].length(); base++, j++) {
                        alignment[i][base] = align_temp[i][j];
                    }
                }
            } else {
                getline(maf, line);
            }
        }

        std::string info_directory = directory + chr + "_temp.txt";
        std::ifstream coor;
        coor.open(info_directory);
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
                std::cout << "start & end coordinate is wierd: " << start << ", " << end << std::endl;
                exit(-1);
            }
            if (strand == '+') {
                for (int i = 0; i < species_size; i++) {
                    align_temp[i] = alignment[i].substr(start, end - start);
                }
            } else {
                for (int i = 0; i < species_size; i++) {
                    align_temp[i] = reverse_base(alignment[i].substr(start, end - start));
                }
            }
            maf_new << "a\n";
            for (int i = 0; i < species_size; i++) {
                maf_new << "s\t" << species[i] << ".chr\t" << std::to_string(start) << "\t" << std::to_string(end - start)
                        << "\t" << strand << "\t1\t" << align_temp[i] << "\n";
            }
            getline(coor, line);
        }
        coor.close();
        maf_new.close();
        maf.close();
    }
    return 0;
}
