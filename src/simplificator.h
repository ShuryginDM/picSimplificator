#pragma once

#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../include/io.h"

Image blockSimplificator(Image &src, int simplificatorFunc, int blockSimplificatorFunc, int sizeX, int sizeY, int maxCount, std::vector<unsigned long long> colors);

void changeColors(Image &src, std::map<unsigned long long, unsigned long long> &colChanges);

void chColorToNearestInList(Image &src, std::vector<unsigned long long> colList, std::map<unsigned long long, int> &cols);

void check_argc(int argc, int from, int to=std::numeric_limits<int>::max());

template<typename ValueT> void check_number(std::string val_name, ValueT val, ValueT from, ValueT to=std::numeric_limits<ValueT>::max());

int colDistance(std::tuple<int, int, int> &color1, std::tuple<int, int, int> &color2);

void colorBlock(Image &res, Image &src, int x0, int y0, int sizeX, int sizeY);

std::tuple<int, int, int> colorFromBlock(Image &src, int simplificatorFunc, int x, int y, int sizeX, int sizeY);

std::map<unsigned long long, unsigned long long>::iterator colorInList(std::map<unsigned long long, unsigned long long> &removal, unsigned long long color);

void cross(Image &res, std::tuple<int, int, int> &color, int x, int y, int size, int radius, bool isDrawingLines);

Image crossSimplificator(Image& src, int simplificatorFunc, int blockSimplificatorFunc, int sizeX, int sizeY, int crossSize, int crossRadius, bool isDrawindLines, int maxCount, std::vector<unsigned long long> colors);

std::tuple<int, int, int> findClosestColorInBlock(Image &src, int x0, int y0, int x1, int y1, std::tuple<int, int, int> color);

std::tuple<int, int, int> findColorAverage(Image &src, int x0, int y0, int x1, int y1);

std::tuple<int, int, int> findColorMax(Image &src, int x0, int y0, int x1, int y1);

std::tuple<int, int, int> findColorMedian(Image &src, int x0, int y0, int x1, int y1);

std::tuple<int, int, int> linColorSimplification(std::tuple<int, int, int> color, int rParam, int gParam, int bParam);

unsigned long long mostRareColor(std::map<unsigned long long, int> &cols);

void print_help(const char *argv0);

template<typename ValueType> ValueType read_value(std::string s);

void removeColorMinCount(Image &src, std::map<unsigned long long, int> &cols, int maxCount);

void removeColorNearest(Image &src, std::map<unsigned long long, int> &cols, int maxCount);

void simplificateBlock(Image &src, int simplificatorFunc, int x, int y, int sizeX, int sizeY);