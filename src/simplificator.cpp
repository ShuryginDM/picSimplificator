#include <cstdlib>
#include <algorithm>
#include <functional>

#include "simplificator.h"

std::map<unsigned long long, int> cols;
int crosses_count = 0;

void changeColorsSimp(Image &src, std::map<unsigned long long, unsigned long long> &colChanges) {
    unsigned long long prevCol, nextCol, nextColT;
    for(int i = 0; i < src.n_rows; ++i){
        for(int j = 0; j < src.n_cols; ++j){
            prevCol = std::get<0>(src(i, j)) << 16 | std::get<1>(src(i, j)) << 8 | std::get<2>(src(i, j));
            if(colChanges.find(prevCol) != colChanges.end()) {
                nextColT = colChanges[prevCol];
                nextCol = nextColT;
                src(i, j) = std::make_tuple(nextCol >> 16, (nextCol >> 8) & 0xFF, nextCol & 0xFF);
            }
        }    
    }    
}


Image blockSimplificator(Image &src, int simplificatorFunc, int blockSimplificatorFunc, int sizeX, int sizeY, int maxCount, std::vector<unsigned long long> colors) {
    Image resT(src.n_rows / sizeX, src.n_cols / sizeY);
    for(int i = 0; i < src.n_rows / sizeX; ++i) {
        for(int j = 0; j < src.n_cols / sizeY; ++j) {
            resT(i, j) = colorFromBlock(src, simplificatorFunc, i * sizeX, j * sizeY, sizeX, sizeY);
        }
    }
    if(cols.size() > maxCount) {
        if(simplificatorFunc == 1){
            removeColorMinCount(resT, cols, maxCount);
        } else {
            removeColorNearest(resT, cols, maxCount);
        }
    }
    if(colors.size()) {
        chColorToNearestInList(resT, colors, cols);
    }
    Image res(resT.n_rows * sizeX, resT.n_cols * sizeY);
    for(int i = 0; i < resT.n_rows; ++i) {
        for(int j = 0; j < resT.n_cols; ++j) {
            colorBlock(res, resT, i, j, sizeX, sizeY);
        }
    }
    return res;
}

void changeColors(Image &src, std::map<unsigned long long, unsigned long long> &colChanges) {
    unsigned long long prevCol, nextCol, nextColT;
    for(int i = 0; i < src.n_rows; ++i){
        for(int j = 0; j < src.n_cols; ++j){
            prevCol = std::get<0>(src(i, j)) << 16 | std::get<1>(src(i, j)) << 8 | std::get<2>(src(i, j));
            if(colChanges.find(prevCol) != colChanges.end()) {
                nextColT = colChanges[prevCol];
                while(colChanges.find(nextColT) != colChanges.end()) {
                    nextColT = colChanges[nextColT];
                }
                nextCol = nextColT;
                src(i, j) = std::make_tuple(nextCol >> 16, (nextCol >> 8) & 0xFF, nextCol & 0xFF);
            }
        }    
    }    
}

void check_argc(int argc, int from, int to) {
    if (argc < from)
        throw std::string("too few arguments for operation");

    if (argc > to)
        throw std::string("too many arguments for operation");
}

template<typename ValueT> void check_number(std::string val_name, ValueT val, ValueT from, ValueT to) {
    if (val < from)
        throw val_name + std::string(" is too small");
    if (val > to)
        throw val_name + std::string(" is too big");
}

int colDistance(std::tuple<int, int, int> color1, std::tuple<int, int, int> color2) {
    int dr = std::abs(std::get<0>(color1) - std::get<0>(color2));
    int dg = std::abs(std::get<1>(color1) - std::get<1>(color2));
    int db = std::abs(std::get<2>(color1) - std::get<2>(color2));
    return std::pow(dr * dr * dr + dg * dg * dg + db * db * db, 1.0 / 3);
}

void colorBlock(Image &res, Image &src, int x0, int y0, int sizeX, int sizeY) {
    for(int i = 0; i < sizeX; ++i) {
        for(int j = 0; j < sizeY; ++j) {
            res(x0 * sizeX + i, y0 * sizeY + j) = src(x0, y0);
        }
    }
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
    unsigned long long color = (r << 16) | (g << 8) | b; 
    if(cols.find(color) == cols.end()) {
        cols[color] = 1;
    } else {
        cols[color]++;
    }
    crosses_count++;
    return std::make_tuple(r, g, b);
}

std::map<unsigned long long, unsigned long long>::iterator colorInList(std::map<unsigned long long, unsigned long long> &removal, unsigned long long color) {
    for(auto i = removal.begin(); i != removal.end(); ++i) {
        if(i->second == color) {
            return i;
        }
    }
    return removal.end();
}

void cross(Image &res, std::tuple<int, int, int> color, int x, int y, int size, int radius, bool isDrawingLines) {
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j) {
            if(isDrawingLines && (!i || !j)){
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

Image crossSimplificator(Image& src, int simplificatorFunc, int blockSimplificatorFunc, int sizeX, int sizeY, int crossSize, int crossRadius, bool isDrawindLines, int maxCount, std::vector<unsigned long long> colors) {
    Image resT(src.n_rows / sizeX, src.n_cols / sizeY);
    Image res((src.n_rows / sizeX) * crossSize, (src.n_cols / sizeY) * crossSize);
    std::cout << "Image size (in crosses) : " << res.n_cols / crossSize << " x " << res.n_rows / crossSize << std::endl;
    for(int i = 0; i < src.n_rows / sizeX; ++i) {
        for(int j = 0; j < src.n_cols / sizeY; ++j) {
            resT(i, j) = colorFromBlock(src, simplificatorFunc, i * sizeX, j * sizeY, sizeX, sizeY);
        }
    }
    if(cols.size() > maxCount) {
        if(simplificatorFunc == 1) {
            removeColorMinCount(resT, cols, maxCount);
        } else {
            removeColorNearest(resT, cols, maxCount);
        }
    }

    std::cout << "Total : " << crosses_count << std::endl;
    for(auto &i : cols){
        std::cout << (i.first >> 16) << " " << ((i.first >> 8) & 0xFF) << " " << (i.first & 0xFF) << " X " << i.second << std::endl;
    }
    std::cout << cols.size() << std::endl;
    
    if(colors.size()) {
        chColorToNearestInList(resT, colors, cols);
    }
    for(int i = 0; i < resT.n_rows; ++i) {
        for(int j = 0; j < resT.n_cols; ++j) {
            cross(res, resT(i, j), i * crossSize, j * crossSize, crossSize, crossRadius, isDrawindLines);
        }
    }

    return res;
}

std::tuple<int, int, int> findClosestColorInBlock(Image &src, int x0, int y0, int x1, int y1, std::tuple<int, int, int> color) {
    int rNear, gNear, bNear, rt, gt, bt;
    int dist = colDistance(std::make_tuple(0, 0, 0), std::make_tuple(255, 255, 255)) + 1;
    int r = std::get<0>(color);
    int g = std::get<1>(color);
    int b = std::get<2>(color);
    for(int i = x0; i < x1; ++i) {
        for(int j = y0; j < y1; ++j) {
            std::tie(rt, gt, bt) = src(i, j);
            if(dist > colDistance(std::make_tuple(r, g, b), std::make_tuple(rt, gt, bt))) {
                dist = colDistance(std::make_tuple(r, g, b), std::make_tuple(rt, gt, bt));
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
    std::tie(sumR, sumG, sumB) = findClosestColorInBlock(src, x0, y0, x1, y1, std::make_tuple(sumR, sumG, sumB));
    return std::make_tuple(sumR, sumG, sumB);
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
    std::tie(maxR, maxG, maxB) = findClosestColorInBlock(src, x0, y0, x1, y1, std::make_tuple(maxR, maxG, maxB));
    return std::make_tuple(maxR, maxG, maxB);
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
    std::tie(r, g, b) = findClosestColorInBlock(src, x0, y0, x1, y1, std::make_tuple(r, g, b));
    int rMed, gMed, bMed, rt, gt, bt, dist = colDistance(std::make_tuple(0, 0, 0), std::make_tuple(255, 255, 255)) + 1;

    return std::make_tuple(r, g, b);
}

std::tuple<int, int, int> linColorSimplification(std::tuple<int, int, int> &color, int rParam, int gParam, int bParam) {
    int r, g, b;
    r = std::get<0>(color) / rParam * rParam;
    g = std::get<1>(color) / gParam * gParam;
    b = std::get<2>(color) / bParam * bParam;
    r = r >= (255 / rParam * rParam) ? 255 : r;
    g = g >= (255 / gParam * gParam) ? 255 : g;
    b = b >= (255 / bParam * bParam) ? 255 : b;
    return std::make_tuple(r, g, b);
}

unsigned long long mostRareColor(std::map<unsigned long long, int> &cols) {
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
    int dist = colDistance(std::make_tuple(0, 0, 0), std::make_tuple(255, 255, 255)) + 1, distT;
    for(auto &i : cols) {
        if(color == i.first) {
            continue;
        }
        distT = colDistance(std::make_tuple(i.first >> 16, (i.first >> 8) & 0xFF, i.first & 0xFF), std::make_tuple(color >> 16, (color >> 8) & 0xFF, color & 0xFF));
        if(distT < dist) {
            dist = distT;
            res = i.first;
        }
    }
    return res;
}

void print_help(const char *argv0) {
    const char *usage = R"(where PARAMS are from list:

    --image-to-blocks --input <input_file_name> (--median || --average || --max) [ --output <output_file_name>]
           
                [--block-size <block_size_x> <block_size_y>]

                [--simple-colors <simplification_r> <simplification_g> <simplification_b>]

                [--max-count <max_color_count>]

                [--drawing-lines]

                [--fast || --slow]

                [--color-list <color_list_filename>]
                
        makes an image as blocks <block_size_x> X <block_size_y> 
        by default, <block_size_x> == <block_size_y>, <block_size_x> == 15, <output_file_name> == <input_file_name>
        simplification_color means step between neighbour colors (and must be from 1 to 254)

    --image-to-crosses <input_file_name> (--median || --average || --max) [ <output_file_name> || 
                
                [--cross_size <cross_size>]

                [--block-size <block_size_x> <block_size_y>]

                [--simple-colors <simplification_r> <simplification_g> <simplification_b>]

                [--max-count <max_color_count>]

                [--drawing-lines]

                [--fast || --slow]

                [--color-list <color_list_filename>]

        makes an image as crosses <cross_size> X <cross_size> merging blocks size of <block_size_x> X <block_size_y>
        by default, <block_size_x> == <block_size_y>, <block_size_x> == 15, <output_file_name> == <input_file_name>,
        <cross_size> = 20, <cross_radius> = 1

    --change-color <input_file_name> <from_R> <from_G> <from_B> <to_R> <to_G> <to_B> [ <output_file_name> ]
        changes color (<from_R>, <from_G>, <from_B>) to (<to_R>, <to_G>, <to_B>) in picture <input_file_name>

    --count-colors <input_file_name> <crossSize>
        prints all cross colors in <input_file_name> (if <crossSize> is wrong answer will be also wrong)

    


    [<param>=default_val] means that parameter is optional.
    )";
    std::cout << "Usage: " << argv0 << " <input_image_path> "
         << "PARAMS" << std::endl;
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

void removeColorMinCount(Image &src, std::map<unsigned long long, int> &cols, int maxCount) {
    unsigned long long prevCol, nextCol;
    std::map<unsigned long long, unsigned long long> removal;

    while(cols.size() > maxCount){
        prevCol = mostRareColor(cols);
        nextCol = nearestColor(cols, prevCol);
        auto index = colorInList(removal, nextCol);
        if(index != removal.end()) {
            index->second = nextCol;
        }
        removal[prevCol] = nextCol;
        cols[nextCol] += cols[prevCol];
        cols.erase(cols.find(prevCol));
    }
    changeColors(src, removal);
}

void chColorToNearestInList(Image &src, std::vector<unsigned long long> colList, std::map<unsigned long long, int> &cols) {
    unsigned long long prevCol, nextCol;
    std::map<unsigned long long, unsigned long long> removal;
    int dist;
    for(auto &i : cols) {
        dist = colDistance(std::make_tuple(0, 0, 0), std::make_tuple(255, 255, 255)) + 1;
        prevCol = i.first;
        for(auto j : colList) {
            if(dist > colDistance(std::make_tuple(prevCol >> 16, (prevCol >> 8) & 0xFF, prevCol & 0xFF), std::make_tuple(j >> 16, (j >> 8) & 0xFF, j & 0xFF))) {
                dist = colDistance(std::make_tuple(prevCol >> 16, (prevCol >> 8) & 0xFF, prevCol & 0xFF), std::make_tuple(j >> 16, (j >> 8) & 0xFF, j & 0xFF));
                nextCol = j;
            }
        }
        removal[prevCol] = nextCol;
    }
    changeColorsSimp(src, removal);
}

void removeColorNearest(Image &src, std::map<unsigned long long, int> &cols, int maxCount) {
    unsigned long long prevCol, nextCol;
    std::map<unsigned long long, unsigned long long> removal;
    int dist;
    while(cols.size() > maxCount){
        dist = colDistance(std::make_tuple(0, 0, 0), std::make_tuple(255, 255, 255)) + 1;
        for(auto &i : cols) {
            for(auto &j : cols) {
                if(i.first == j.first) {
                    continue;
                }
                unsigned long long c1 = i.first, c2 = j.first;
                if(dist > colDistance(std::make_tuple(c1 >> 16, (c1 >> 8) & 0xFF, c1 & 0xFF), std::make_tuple(c2 >> 16, (c2 >> 8) & 0xFF, c2 & 0xFF))) {
                    dist = colDistance(std::make_tuple(c1 >> 16, (c1 >> 8) & 0xFF, c1 & 0xFF), std::make_tuple(c2 >> 16, (c2 >> 8) & 0xFF, c2 & 0xFF));
                    prevCol = c1;
                    nextCol = c2;
                }
            }
        }
        if(cols[prevCol] > cols[nextCol]) {
            std::swap(prevCol, nextCol);
        }
        auto index = colorInList(removal, nextCol);
        if(index != removal.end()) {
            index->second = nextCol;
        }
        removal[prevCol] = nextCol;
        cols[nextCol] += cols[prevCol];
        cols.erase(cols.find(prevCol));
    }   
    changeColors(src, removal);
}

void simplificateBlock(Image &src, int simplificatorFunc, int x, int y, int sizeX, int sizeY) {
    int r, g, b;
    std::tuple<int, int, int> color;
    switch(simplificatorFunc) {
        case 0: 
            color = findColorAverage(src, x, y, x + sizeX, y + sizeY);
            break;
        case 1:
            color = findColorMedian(src, x, y, x + sizeX, y + sizeY);
            break;
        case 2:
            color = findColorMax(src, x, y, x + sizeX, y + sizeY);
    }
    for(int i = 0; i < sizeX; ++i){
        for(int j = 0; j < sizeY; ++j){
            src(x + i, y + j) = color;
        }
    }
}