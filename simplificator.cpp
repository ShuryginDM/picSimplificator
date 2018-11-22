#include "include/io.h"
#include <map>
#include <string>
#include <cstdlib>
#include <sstream>
#include <limits>
#include <algorithm>
#include <set>
#include <functional>

std::map<unsigned long long, int> cols;
int crosses_count = 0;
int COLOR_SIMPLIFICATE_PARAMETR_R = 1;
int COLOR_SIMPLIFICATE_PARAMETR_G = 1;
int COLOR_SIMPLIFICATE_PARAMETR_B = 1;
int MAX_COLOR_COUNT = 0x1000000;



void print_help(const char *argv0) {
    const char *usage = R"(where PARAMS are from list:

    --merge-blocks <input_file_name> (--median || --average || --max) [ <output_file_name> ||
           
                <block_size_x> ||
               
                <block_size_y> ||

                <simplification_r> ||

                <simplification_g> ||

                <simplification_b> ||

                <max_color_count> ]
        makes an image as blocks <block_size_x> X <block_size_y> 
        by default, <block_size_x> == <block_size_y>, <block_size_x> == 15, <output_file_name> == <input_file_name>
        simplification_color means step between neighbour colors (and must be from 1 to 254)

    --image-to-crosses <input_file_name> (--median || --average || --max) [ <output_file_name> || 
            
                <cross_size> ||
                
                <cross_radius> ||
                
                <block_size_x> ||
                
                <block_size_y> ||

                <simplification_r> ||

                <simplification_g> ||

                <simplification_b> ||

                <max_color_count> ]
        makes an image as crosses <cross_size> X <cross_size> merging blocks size of <block_size_x> X <block_size_y>
        by default, <block_size_x> == <block_size_y>, <block_size_x> == 15, <output_file_name> == <input_file_name>,
        <cross_size> = 20, <cross_radius> = 1

    --change-color <input_file_name> <from_R> <from_G> <from_B> <to_R> <to_G> <to_B> [ <output_file_name> ]
        changes color (<from_R>, <from_G>, <from_B>) to (<to_R>, <to_G>, <to_B>) in picture <input_file_name>

    --count-colors <input_file_name> <crossSize>
        prints all cross colors in <input_file_name> (if <crossSize> is wrong answer will be also wrong)


    [<param>=default_val] means that parameter is optional.
    )";
    std::cout << "Usage: " << argv0 << " PARAMS" << std::endl;
    std::cout << usage;
}

template<typename ValueType> ValueType read_value(std::string s) {
    std::stringstream ss(s);
    ValueType res;
    ss >> res;
    if (ss.fail() or not ss.eof())
        throw std::string("bad argument: ") + s;
    return res;
}


void check_argc(int argc, int from, int to=std::numeric_limits<int>::max()) {
    if (argc < from)
        throw std::string("too few arguments for operation");

    if (argc > to)
        throw std::string("too many arguments for operation");
}

template<typename ValueT> void check_number(std::string val_name, ValueT val, ValueT from, ValueT to=std::numeric_limits<ValueT>::max()) {
    if (val < from)
        throw val_name + std::string(" is too small");
    if (val > to)
        throw val_name + std::string(" is too big");
}

Image changeColor(Image &src, int r0, int g0, int b0, int r1, int g1, int b1) {
    for(int i = 0; i < src.n_rows; ++i){
        for(int j = 0; j < src.n_cols; ++j){
            if(r0 == std::get<0> (src(i, j)) && g0 == std::get<1> (src(i, j)) && b0 == std::get<2> (src(i, j))) {
                src(i, j) = std::make_tuple(r1, g1, b1);
            }
        }    
    }
    return src;
}

std::tuple<int, int, int> linSimplification(int r, int g, int b) {
    r = r / COLOR_SIMPLIFICATE_PARAMETR_R * COLOR_SIMPLIFICATE_PARAMETR_R;
    g = g / COLOR_SIMPLIFICATE_PARAMETR_G * COLOR_SIMPLIFICATE_PARAMETR_G;
    b = b / COLOR_SIMPLIFICATE_PARAMETR_B * COLOR_SIMPLIFICATE_PARAMETR_B;
    r = r >= (255 / COLOR_SIMPLIFICATE_PARAMETR_R * COLOR_SIMPLIFICATE_PARAMETR_R) ? 255 : r;
    g = g >= (255 / COLOR_SIMPLIFICATE_PARAMETR_G * COLOR_SIMPLIFICATE_PARAMETR_G) ? 255 : g;
    b = b >= (255 / COLOR_SIMPLIFICATE_PARAMETR_B * COLOR_SIMPLIFICATE_PARAMETR_B) ? 255 : b;
    return std::make_tuple(r, g, b);
}

int colDistance(int r, int g, int b, int r1, int g1, int b1) {
    return std::abs(r - r1) + std::abs(g - g1) + std::abs(b - b1);
}

std::tuple<int, int, int> findNearestColorInBlock(Image &src, int x0, int y0, int x1, int y1, int r, int g, int b) {
    int rNear, gNear, bNear, rt, gt, bt, dist = colDistance(0, 0, 0, 255, 255, 255) + 1;
    for(int i = x0; i < x1; ++i) {
        for(int j = y0; j < y1; ++j) {
            std::tie(rt, gt, bt) = src(i, j);
            if(dist > colDistance(r, g, b, rt, gt, bt)) {
                dist = colDistance(r, g, b, rt, gt, bt);
                rNear = rt;
                gNear = gt;
                bNear = bt;
            }
        }
    }
    return std::make_tuple(rNear, gNear, bNear);
}

std::tuple<int, int, int> findColorAverage(Image &src, int x0, int y0, int x1, int y1) {
    long long sumR = 0, sumB = 0, sumG = 0;
    if(x0 >= x1 || y0 >= y1) {
        return std::make_tuple(0, 0, 0);
    }
    for(int i = x0; i < x1; ++i) {
        for(int j = y0; j < y1; ++j) {
            sumR += std::get<0>(src(i, j));
            sumG += std::get<1>(src(i, j));
            sumB += std::get<2>(src(i, j));
        }
    }
    int size = (x1 - x0) * (y1 - y0);
    sumR /= size;
    sumG /= size;
    sumB /= size;
    std::tie(sumR, sumG, sumB) = findNearestColorInBlock(src, x0, y0, x1, y1, sumR, sumG, sumB);
    return std::make_tuple(sumR, sumG, sumB);
}

std::tuple<int, int, int> findColorMedian(Image &src, int x0, int y0, int x1, int y1) {
    std::vector<int> rV(256), gV(256), bV(256);
    int r, g, b;
    for(int i = x0; i < x1; ++i) {
        for(int j = y0; j < y1; ++j) {
            std::tie(r, g, b) = src(i, j);
            rV[r]++;
            gV[g]++;
            bV[b]++;
        }
    }
    int rSize = 0, gSize = 0, bSize = 0;
    for(int i = 0; i < 256; ++i) {
        if(rV[i]) {
            rSize++;
        }
        if(gV[i]) {
            gSize++;
        }
        if(bV[i]) {
            bSize++;
        }
    }
    std::sort(rV.begin(), rV.end(), std::greater<int>());
    std::sort(gV.begin(), gV.end(), std::greater<int>());
    std::sort(bV.begin(), bV.end(), std::greater<int>());
    std::tie(r, g, b) = linSimplification(rV[rSize / 2], gV[gSize / 2], bV[bSize / 2]);
    std::tie(r, g, b) = findNearestColorInBlock(src, x0, y0, x1, y1, r, g, b);
    int rMed, gMed, bMed, rt, gt, bt, dist = colDistance(0, 0, 0, 255, 255, 255) + 1;
    // if(r > 240 && g > 240 && b > 240) {
    //     r = g = b = 255;
    // }
    return std::make_tuple(r, g, b);
}

std::tuple<int, int, int> findColorMax(Image &src, int x0, int y0, int x1, int y1) {
    std::vector<int> rCount(256), gCount(256), bCount(256);
    int r, g, b;
    int maxR, maxG, maxB;
    for(int i = x0; i < x1; ++i) {
        for(int j = y0; j < y1; ++j) {
            std::tie(r, g, b) = src(i, j);
            rCount[r]++;
            gCount[g]++;
            bCount[b]++;
        }
    }
    r = 0;
    g = 0;
    b = 0;
    for(int i = 0; i < 255; ++i) {
        if(rCount[i] > r) {
            r = rCount[i];
            maxR = i;
        }
        if(gCount[i] > g) {
            g = gCount[i];
            maxG = i;
        }
        if(bCount[i] > b) {
            b = bCount[i];
            maxB = i;
        }
    }
    std::tie(maxR, maxG, maxB) = findNearestColorInBlock(src, x0, y0, x1, y1, maxR, maxG, maxB);
    return std::make_tuple(maxR, maxG, maxB);
}

void simplificateBlock(Image &src, int simplificatorFunc, int x, int y, int sizeX, int sizeY) {
    int r, g, b;
    std::tie(r, g, b) = simplificatorFunc   ? findColorMedian(src, x, y, x + sizeX, y + sizeY)
                                            : findColorAverage(src, x, y, x + sizeX, y + sizeY);
    for(int i = 0; i < sizeX; ++i){
        for(int j = 0; j < sizeY; ++j){
            src(x + i, y + j) = std::make_tuple(r / COLOR_SIMPLIFICATE_PARAMETR_R * COLOR_SIMPLIFICATE_PARAMETR_R, g / COLOR_SIMPLIFICATE_PARAMETR_G * COLOR_SIMPLIFICATE_PARAMETR_G, b / COLOR_SIMPLIFICATE_PARAMETR_B * COLOR_SIMPLIFICATE_PARAMETR_B);
        }
    }
}

unsigned long long minColor(std::map<unsigned long long, int> &cols) {
    bool b = false;
    unsigned long long res;
    int minCount;
    for(auto &i : cols) {
        if(!b || i.second < minCount) {
            minCount = i.second;
            res = i.first;
            b = true;
        }
    }
    return res;
}

unsigned long long nearestColor(std::map<unsigned long long, int> &cols, unsigned long long color) {
    unsigned long long res;
    int dist = colDistance(0, 0, 0, 255, 255, 255) + 1, distT;
    for(auto &i : cols) {
        if(color == i.first) {
            continue;
        }
        distT = colDistance(i.first >> 16, (i.first >> 8) & 0xFF, i.first & 0xFF, color >> 16, (color >> 8) & 0xFF, color & 0xFF);
        if(distT < dist) {
            dist = distT;
            res = i.first;
        }
    }
    return res;
}



std::tuple<int, int, int> colorFromBlock(Image &src, int simplificatorFunc, int x, int y, int sizeX, int sizeY) {
    int r, g, b;
    switch(simplificatorFunc) {
        case 0:
            std::tie(r, g, b) = findColorAverage(src, x, y, x + sizeX, y + sizeY);
            break;
        case 1:
            std::tie(r, g, b) = findColorMedian(src, x, y, x + sizeX, y + sizeY);
            break;
        case 2:
            std::tie(r, g, b) = findColorMax(src, x, y, x + sizeX, y + sizeY);
            break;
    }
     if(r != g && r != b && g != b) {
        int min = std::min({r, g, b});
        int max = std::max({r, g, b});
        int mid = r + g + b - max - min; 
        if(r == min) {
            r = mid;
        }
        if(g == min) {
            g = mid;
        }
        if(b == min) {
            b = mid;
        }
    }
    if(r < 255 || g < 255 || b < 255) {
        unsigned long long color = (r << 16) | (g << 8) | b; 
        if(cols.find(color) == cols.end()) {
            cols[color] = 1;
        } else {
            cols[color]++;
        }
        crosses_count++;
    }
    return std::make_tuple(r, g, b);
}

void cross(Image &res, std::tuple<int, int, int> color, int x, int y, int size, int radius) {
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j) {
            if(!i || !j){
                res(x + i, y + j) = std::make_tuple(0, 0, 0);
                continue;
            }
            if((i - radius <= j && i + radius >= j) || (i - radius <= size - j && i + radius >= size - j)) {
                res(x + i, y + j) = color;
            } else {
                res(x + i, y + j) = std::make_tuple(255, 255, 255);
            }
        }
    }
}

Image blockSimplificator(Image src, int simplificatorFunc, int sizeX, int sizeY) {
    for(int i = 0; i < src.n_rows / sizeX; ++i) {
        for(int j = 0; j < src.n_cols / sizeY; ++j) {
            simplificateBlock(src, simplificatorFunc, i * sizeX, j * sizeY, sizeX, sizeY);
        }
    } 
    while(cols.size() > MAX_COLOR_COUNT) {
        unsigned long long prevCol = minColor(cols);
        unsigned long long nextCol = nearestColor(cols, prevCol);
        src = changeColor(src, prevCol >> 16, (prevCol >> 8) & 0xFF, prevCol & 0xFF, nextCol >> 16, (nextCol >> 8) & 0xFF, nextCol & 0xFF);
        cols[nextCol] += cols[prevCol];
        cols.erase(cols.find(prevCol));
    }
    return src;
}

Image crossSimplificator(Image& src, int simplificatorFunc, int sizeX, int sizeY, int crossSize, int crossRadius) {
    Image resT(src.n_rows / sizeX, src.n_cols / sizeY);
    Image res((src.n_rows / sizeX) * crossSize, (src.n_cols / sizeY) * crossSize);
    std::cout << "Image size (in crosses) : " << res.n_cols / crossSize << " x " << res.n_rows / crossSize << std::endl;
    for(int i = 0; i < src.n_rows / sizeX; ++i) {
        for(int j = 0; j < src.n_cols / sizeY; ++j) {
            //cross(res, colorFromBlock(src, simplificatorFunc, i * sizeX, j * sizeY, sizeX, sizeY), i * crossSize, j * crossSize, crossSize, crossRadius);
            resT(i, j) = colorFromBlock(src, simplificatorFunc, i * sizeX, j * sizeY, sizeX, sizeY);
        }
    }
    while(cols.size() > MAX_COLOR_COUNT) {
        unsigned long long prevCol = minColor(cols);
        unsigned long long nextCol = nearestColor(cols, prevCol);
        resT = changeColor(resT, prevCol >> 16, (prevCol >> 8) & 0xFF, prevCol & 0xFF, nextCol >> 16, (nextCol >> 8) & 0xFF, nextCol & 0xFF);
        cols[nextCol] += cols[prevCol];
        cols.erase(cols.find(prevCol));
    }
    std::cout << "Total : " << crosses_count << std::endl;
    for(auto &i : cols){
        std::cout << (i.first >> 16) << " " << ((i.first >> 8) & 0xFF) << " " << (i.first & 0xFF) << " X " << i.second << std::endl;
    }
    std::cout << cols.size() << std::endl;
    for(int i = 0; i < resT.n_rows; ++i) {
        for(int j = 0; j < resT.n_cols; ++j) {
            cross(res, resT(i, j), i * crossSize, j * crossSize, crossSize, crossRadius);
        }
    }
    return res;
}


int main(int argc, char **argv) {
    try {
        check_argc(argc, 2);
        if (std::string(argv[1]) == "--help") {
            print_help(argv[0]);
            return 0;
        }
        std::string dst_image_name;
        check_argc(argc, 3);
        Image src_image = load_image(argv[2]), dst_image;
        dst_image_name = argv[2];
        std::string action(argv[1]);

        if (action == "--image-to-crosses") {
            int sizeX = 15, sizeY = 15;
            int crossSize = 20, crossRadius = 1;
            check_argc(argc, 4, 14);
            int simplificatorFunc;
            if(std::string(argv[3]) == "--median") {
                simplificatorFunc = 1;
            } else if(std::string(argv[3]) == "--average") {
                simplificatorFunc = 0;
            } else if(std::string(argv[3]) == "--max") {
                simplificatorFunc = 2;
            } else {
                throw "3rd argument must be --median or --average";
            }
            switch(argc) {
                case 13:
                    MAX_COLOR_COUNT = read_value<int>(argv[12]);
                    check_number("MAX_COLOR_COUNT", MAX_COLOR_COUNT, 2, 0x1000000);
                case 12:
                    COLOR_SIMPLIFICATE_PARAMETR_B = read_value<int>(argv[11]);
                    check_number("COLOR_SIMPLIFICATE_PARAMETER_B", COLOR_SIMPLIFICATE_PARAMETR_B, 1, 254);
                case 11:
                    COLOR_SIMPLIFICATE_PARAMETR_G = read_value<int>(argv[10]);
                    check_number("COLOR_SIMPLIFICATE_PARAMETER_G", COLOR_SIMPLIFICATE_PARAMETR_G, 1, 254);
                case 10:
                    COLOR_SIMPLIFICATE_PARAMETR_R = read_value<int>(argv[9]);
                    check_number("COLOR_SIMPLIFICATE_PARAMETER_R", COLOR_SIMPLIFICATE_PARAMETR_R, 1, 254);
                case 9:
                    sizeY = read_value<int>(argv[8]);
                case 8:
                    sizeX = read_value<int>(argv[7]);
                case 7:
                    crossRadius = read_value<int>(argv[6]);
                case 6:
                    crossSize = read_value<int>(argv[5]);
                case 5:
                    dst_image_name = argv[4];
                case 4:
                    dst_image = crossSimplificator(src_image, simplificatorFunc, sizeX, sizeY, crossSize, crossRadius);
            }
        } else if (action == "--merge-blocks") {
            int sizeX = 15, sizeY = 15;
            int crossSize = 20, crossRadius = 1;
            check_argc(argc, 4, 10);
            int simplificatorFunc;
            if(std::string(argv[3]) == "--median") {
                simplificatorFunc = 1;
            } else if(std::string(argv[3]) == "--average") {
                simplificatorFunc = 0;
            } else if(std::string(argv[3]) == "--max") {
                simplificatorFunc = 2;
            } else {
                throw "3rd argument must be --median or --average";
            }
            switch(argc){
                case 10:
                    COLOR_SIMPLIFICATE_PARAMETR_B = read_value<int>(argv[9]);
                    check_number("COLOR_SIMPLIFICATE_PARAMETER_B", COLOR_SIMPLIFICATE_PARAMETR_B, 1, 254);
                case 9:
                    COLOR_SIMPLIFICATE_PARAMETR_G = read_value<int>(argv[8]);
                    check_number("COLOR_SIMPLIFICATE_PARAMETER_G", COLOR_SIMPLIFICATE_PARAMETR_G, 1, 254);
                case 8:
                    COLOR_SIMPLIFICATE_PARAMETR_R = read_value<int>(argv[7]);
                    check_number("COLOR_SIMPLIFICATE_PARAMETER_R", COLOR_SIMPLIFICATE_PARAMETR_R, 1, 254);
                case 7:
                    sizeY = read_value<int>(argv[6]);
                case 6:
                    sizeX = read_value<int>(argv[5]);
                case 5:
                    dst_image_name = argv[4];
                case 4:
                    dst_image = blockSimplificator(src_image, simplificatorFunc, sizeX, sizeY);
            }
        } else if (action == "--change-color") {
            check_argc(argc, 9, 10);
            switch(argc){
                case 10:
                    dst_image_name = argv[9];
                case 9:
                    int r0 = read_value<int>(argv[3]);
                    int g0 = read_value<int>(argv[4]);
                    int b0 = read_value<int>(argv[5]);
                    int r1 = read_value<int>(argv[6]);
                    int g1 = read_value<int>(argv[7]);
                    int b1 = read_value<int>(argv[8]);
                    dst_image = changeColor(src_image, r0, g0, b0, r1, g1, b1);
            }
        } else if (action == "--count-colors") {
            check_argc(argc, 4, 4);
            int crossSize = read_value<int>(argv[3]);
            std::set<unsigned long long> difColors;
            int r, g, b;
            for(int i = 0; i < src_image.n_rows / crossSize; ++i){
                for(int j = 0; j < src_image.n_cols / crossSize; ++j){
                    r = std::get<0>(src_image(i * crossSize + crossSize / 2, j * crossSize + crossSize / 2));
                    g = std::get<1>(src_image(i * crossSize + crossSize / 2, j * crossSize + crossSize / 2));
                    b = std::get<2>(src_image(i * crossSize + crossSize / 2, j * crossSize + crossSize / 2));
                    if(r + g + b < 255 * 3) { 
                        difColors.insert((r << 16) | (g << 8) | b);
                    }
                }
            }
            std::cout << difColors.size() << std::endl;
            for(auto i : difColors) {
                std::cout << (i >> 16) << " " << ((i >> 8) & 0xFF) << " " << (i & 0xFF) << std::endl;
            }
            return 0;
        } else {
            throw std::string("unknown action ") + action;
        }
        save_image(dst_image, dst_image_name.c_str());
    } catch (const std::string &s) {
        std::cerr << "Error: " << s << std::endl;
        std::cerr << "For help type: " << std::endl << argv[0] << " --help" << std::endl;
        return 1;
    }
}
