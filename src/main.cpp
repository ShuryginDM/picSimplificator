#include "simplificator.h"
#include <vector>
#include <fstream>

template<typename ValueType> ValueType read_value(std::string s) {
    std::stringstream ss(s);
    ValueType res;
    ss >> res;
    if (ss.fail() or not ss.eof())
        throw std::string("bad argument: ") + s;
    return res;
}

void getColors(std::vector<unsigned long long> &v, std::ifstream &inp) {
    int count, r, g, b;
    inp >> count;
    for(int i = 0; i < count; ++i) {
        inp >> r >> g >> b;
        v.push_back((r << 16) | (g << 8) | b);
    }
}

enum class Actions {
    HELP,
    IMAGE_TO_CROSSES,
    IMAGE_TO_BLOCKS,
    MEDIAN,
    AVERAGE,
    MAX,
    SLOW,
    FAST,
    COLOR_SIMPLIFICATE_PARAMETERS,
    MAX_COLOR_COUNT,
    DRAWING_LINES,
    BLOCK_SIZE,
    CROSS_SIZE,
    INPUT_FILE,
    OUTPUT_FILE,
    COLORS_LIST_FILE
};

Actions getAction(char *str) {
    if(std::string(str) == "--help") {
        return Actions::HELP;
    }
    if(std::string(str) == "--image-to-crosses") {
        return Actions::IMAGE_TO_CROSSES;
    }
    if(std::string(str) == "--image-to-blocks") {
        return Actions::IMAGE_TO_BLOCKS;
    }
    if(std::string(str) == "--median") {
        return Actions::MEDIAN;
    }
    if(std::string(str) == "--average") {
        return Actions::AVERAGE;
    }
    if(std::string(str) == "--max") {
        return Actions::MAX;
    }
    if(std::string(str) == "--slow") {
        return Actions::SLOW;
    }
    if(std::string(str) == "--fast") {
        return Actions::FAST;
    }
    if(std::string(str) == "--simple-colors") {
        return Actions::COLOR_SIMPLIFICATE_PARAMETERS;
    }
    if(std::string(str) == "--max-count") {
        return Actions::MAX_COLOR_COUNT;
    }
    if(std::string(str) == "--drawing-lines") {
        return Actions::DRAWING_LINES;
    }
    if(std::string(str) == "--block-size") {
        return Actions::BLOCK_SIZE;
    }
    if(std::string(str) == "--cross-size") {
        return Actions::CROSS_SIZE;
    }
    if(std::string(str) == "--input") {
        return Actions::INPUT_FILE;
    }
    if(std::string(str) == "--output") {
        return Actions::OUTPUT_FILE;
    }
    if(std::string(str) == "--color-list") {
        return Actions::COLORS_LIST_FILE;
    }
}

int main(int argc, char **argv) {
    try {
        bool mainActionChosen = false;
        bool simplificationTypeChosen = false;
        bool blockSimplificationTypeChosen = false;
        bool simplificationParametersChosen = false;
        bool colorsFileChosen = false;
        bool blockSizeChosen = false;
        bool crossSizeChosen = false;
        int mainAction, simplificationFunc, blockSimplificationFunc;
        int rParametr = 1, gParametr = 1, bParametr = 1;
        int maxCount = std::numeric_limits<int>::max();
        bool drawingLines = false;
        int sizeX = 1, sizeY = 1;
        int crossSize = 1, crossRadius = 1;
        std::string inputFileName, outputFileName, colorsFileName;
        check_argc(argc, 2);
        std::vector<unsigned long long> colors;
        for(int i = 1; i < argc; ++i) {
            Actions action = getAction(argv[i]);
            switch(action) {
                case Actions::HELP:
                    print_help(argv[0]);
                    return 0;
                case Actions::IMAGE_TO_BLOCKS:
                    if(mainActionChosen) {
                        throw "main action must be chosen just one time";
                    }
                    mainAction = 0;
                    mainActionChosen = true;
                    break;
                case Actions::IMAGE_TO_CROSSES:
                    if(mainActionChosen) {
                        throw "main action must be chosen just one time";
                    }
                    mainAction = 1;
                    mainActionChosen = true;
                    break;
                case Actions::AVERAGE:
                    blockSimplificationFunc = 0;
                    break;
                case Actions::MEDIAN:
                    blockSimplificationFunc = 1;
                    break;
                case Actions::MAX:
                    blockSimplificationFunc = 2;
                    break;
                case Actions::SLOW:
                    simplificationFunc = 0;
                    break;
                case Actions::FAST:
                    simplificationFunc = 1;
                    break;
                case Actions::COLOR_SIMPLIFICATE_PARAMETERS:
                    check_argc(argc, i + 4);
                    rParametr = read_value<int>(argv[i + 1]);
                    gParametr = read_value<int>(argv[i + 2]);
                    bParametr = read_value<int>(argv[i + 3]);
                    i += 3;
                    break;
                case Actions::MAX_COLOR_COUNT:
                    check_argc(argc, i + 2);
                    maxCount = read_value<int>(argv[i + 1]);
                    i++;
                    break;
                case Actions::DRAWING_LINES:
                    drawingLines = true;
                    break;
                case Actions::BLOCK_SIZE:
                    check_argc(argc, i + 3);
                    sizeX = read_value<int>(argv[i + 1]);
                    sizeY = read_value<int>(argv[i + 2]);
                    i += 2;
                    break;
                case Actions::CROSS_SIZE:
                    check_argc(argc, i + 3);
                    crossSize = read_value<int>(argv[i + 1]);
                    crossRadius = read_value<int>(argv[i + 2]);
                    i += 2;
                    break;
                case Actions::INPUT_FILE:
                    check_argc(argc, i + 2);
                    inputFileName = argv[i + 1];
                    i++;
                    break;
                case Actions::OUTPUT_FILE:
                    check_argc(argc, i + 2);
                    outputFileName = argv[i + 1];
                    i++;
                    break;
                case Actions::COLORS_LIST_FILE:
                    check_argc(argc, i + 2);
                    colorsFileName = argv[i + 1];
                    colorsFileChosen = true;
                    i++;
                    break;
                default:
                    throw "unknown action" + std::string(argv[i]);
            }
        }
        if(!mainActionChosen) {
            throw "main action must be chosen";
        }
        if(colorsFileChosen) {
            std::ifstream colorsFile;
            colorsFile.open(colorsFileName);
            getColors(colors, colorsFile);
            colorsFile.close();
        }
        Image src_image = load_image(inputFileName.c_str()), dst_image;
        switch(mainAction) {
            case 0:
                dst_image = blockSimplificator(src_image, simplificationFunc, blockSimplificationFunc, sizeX, sizeY, maxCount, colors);
                break;
            case 1:
                dst_image = crossSimplificator(src_image, simplificationFunc, blockSimplificationFunc, sizeX, sizeY, crossSize, crossRadius, drawingLines, maxCount, colors);
                break;
        }
        
        save_image(dst_image, outputFileName.c_str());

        // if (std::string(argv[1]) == "--help") {
        //     print_help(argv[0]);
        //     return 0;
        // }
        // std::string dst_image_name;
        // check_argc(argc, 3);
        // Image src_image = load_image(argv[2]), dst_image;
        // dst_image_name = argv[2];
        // std::string action(argv[1]);
        // if (action == "--image-to-crosses") {
        //     int colorSimplificatorType = 0;
        //     int sizeX = 15, sizeY = 15;
        //     int crossSize = 20, crossRadius = 1;
        //     check_argc(argc, 4, 14);
        //     int simplificatorFunc;
        //     if(std::string(argv[3]) == "--median") {
        //         simplificatorFunc = 1;
        //     } else if(std::string(argv[3]) == "--average") {
        //         simplificatorFunc = 0;
        //     } else if(std::string(argv[3]) == "--max") {
        //         simplificatorFunc = 2;
        //     } else {
        //         throw "3rd argument must be --median or --average";
        //     }
        //     switch(argc) {
        //         case 14:
        //             if(std::string(argv[13]) == "--fast") {
        //                 colorSimplificatorType = 1;
        //             } else {
        //                 throw "13th argument must be --fast";
        //             }
        //         case 13:
        //             MAX_COLOR_COUNT = read_value<int>(argv[12]);
        //             check_number("MAX_COLOR_COUNT", MAX_COLOR_COUNT, 2, 0x1000000);
        //         case 12:
        //             COLOR_SIMPLIFICATE_PARAMETR_B = read_value<int>(argv[11]);
        //             check_number("COLOR_SIMPLIFICATE_PARAMETER_B", COLOR_SIMPLIFICATE_PARAMETR_B, 1, 254);
        //         case 11:
        //             COLOR_SIMPLIFICATE_PARAMETR_G = read_value<int>(argv[10]);
        //             check_number("COLOR_SIMPLIFICATE_PARAMETER_G", COLOR_SIMPLIFICATE_PARAMETR_G, 1, 254);
        //         case 10:
        //             COLOR_SIMPLIFICATE_PARAMETR_R = read_value<int>(argv[9]);
        //             check_number("COLOR_SIMPLIFICATE_PARAMETER_R", COLOR_SIMPLIFICATE_PARAMETR_R, 1, 254);
        //         case 9:
        //             sizeY = read_value<int>(argv[8]);
        //         case 8:
        //             sizeX = read_value<int>(argv[7]);
        //         case 7:
        //             crossRadius = read_value<int>(argv[6]);
        //         case 6:
        //             crossSize = read_value<int>(argv[5]);
        //         case 5:
        //             dst_image_name = argv[4];
        //         case 4:
        //             dst_image = crossSimplificator(src_image, simplificatorFunc, sizeX, sizeY, crossSize, crossRadius, colorSimplificatorType);
        //     }
        // } else if (action == "--merge-blocks") {
        //     int colorSimplificatorType = 0;
        //     int sizeX = 15, sizeY = 15;
        //     int crossSize = 20, crossRadius = 1;
        //     check_argc(argc, 4, 12);
        //     int simplificatorFunc;
        //     if(std::string(argv[3]) == "--median") {
        //         simplificatorFunc = 1;
        //     } else if(std::string(argv[3]) == "--average") {
        //         simplificatorFunc = 0;
        //     } else if(std::string(argv[3]) == "--max") {
        //         simplificatorFunc = 2;
        //     } else {
        //         throw "3rd argument must be --median or --average";
        //     }
        //     switch(argc){
        //         case 12:
        //             if(std::string(argv[11]) == "--fast") {
        //                 colorSimplificatorType = 1;
        //             } else {
        //                 throw "13th argument must be --fast";
        //             }
        //         case 11:
        //             MAX_COLOR_COUNT = read_value<int>(argv[10]);
        //             check_number("MAX_COLOR_COUNT", MAX_COLOR_COUNT, 2, 0x1000000);
        //         case 10:
        //             COLOR_SIMPLIFICATE_PARAMETR_B = read_value<int>(argv[9]);
        //             check_number("COLOR_SIMPLIFICATE_PARAMETER_B", COLOR_SIMPLIFICATE_PARAMETR_B, 1, 254);
        //         case 9:
        //             COLOR_SIMPLIFICATE_PARAMETR_G = read_value<int>(argv[8]);
        //             check_number("COLOR_SIMPLIFICATE_PARAMETER_G", COLOR_SIMPLIFICATE_PARAMETR_G, 1, 254);
        //         case 8:
        //             COLOR_SIMPLIFICATE_PARAMETR_R = read_value<int>(argv[7]);
        //             check_number("COLOR_SIMPLIFICATE_PARAMETER_R", COLOR_SIMPLIFICATE_PARAMETR_R, 1, 254);
        //         case 7:
        //             sizeY = read_value<int>(argv[6]);
        //         case 6:
        //             sizeX = read_value<int>(argv[5]);
        //         case 5:
        //             dst_image_name = argv[4];
        //         case 4:
        //             dst_image = blockSimplificator(src_image, simplificatorFunc, sizeX, sizeY, colorSimplificatorType);
        //     }
        // } else if (action == "--change-color") {
        //     check_argc(argc, 9, 10);
        //     switch(argc){
        //         case 10:
        //             dst_image_name = argv[9];
        //         case 9:
        //             int r0 = read_value<int>(argv[3]);
        //             int g0 = read_value<int>(argv[4]);
        //             int b0 = read_value<int>(argv[5]);
        //             int r1 = read_value<int>(argv[6]);
        //             int g1 = read_value<int>(argv[7]);
        //             int b1 = read_value<int>(argv[8]);
        //             dst_image = changeColor(src_image, r0, g0, b0, r1, g1, b1);
        //     }
        // } else if (action == "--count-colors") {
        //     check_argc(argc, 4, 4);
        //     int crossSize = read_value<int>(argv[3]);
        //     std::set<unsigned long long> difColors;
        //     int r, g, b;
        //     for(int i = 0; i < src_image.n_rows / crossSize; ++i){
        //         for(int j = 0; j < src_image.n_cols / crossSize; ++j){
        //             r = std::get<0>(src_image(i * crossSize + crossSize / 2, j * crossSize + crossSize / 2));
        //             g = std::get<1>(src_image(i * crossSize + crossSize / 2, j * crossSize + crossSize / 2));
        //             b = std::get<2>(src_image(i * crossSize + crossSize / 2, j * crossSize + crossSize / 2));
        //             if(r + g + b < 255 * 3) { 
        //                 difColors.insert((r << 16) | (g << 8) | b);
        //             }
        //         }
        //     }
        //     std::cout << difColors.size() << std::endl;
        //     for(auto i : difColors) {
        //         std::cout << (i >> 16) << " " << ((i >> 8) & 0xFF) << " " << (i & 0xFF) << std::endl;
        //     }
        //     return 0;
        // } else {
        //     throw std::string("unknown action ") + action;
        // }
        // save_image(dst_image, dst_image_name.c_str());
    } catch (const std::string &s) {
        std::cerr << "Error: " << s << std::endl;
        std::cerr << "For help type: " << std::endl << argv[0] << " --help" << std::endl;
        return 1;
    }
}
