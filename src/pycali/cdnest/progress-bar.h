/*
 * Progress Bar
 * A simple progress bar implementation in C
 * 
 * Adapted from https://github.com/ankddev/progress-bar-c
 */

#ifndef _PROGRESS_H
#define _PROGRESS_H

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
    char symbol;
    char startSymbol;
    char endSymbol;
    int length;
    int progress;
    int total;
    char* format;
    char* completedText;
    bool showPercent;
    bool showCount;
} ProgressBar;

ProgressBar* pb_alloc();

void pb_free(ProgressBar *pb);

void pb_init(ProgressBar *pb, char symbol, int length, int total);

ProgressBar pb_update(ProgressBar *pb, int progress);

ProgressBar showPercent(ProgressBar *pb, bool show);

char getSymbol(ProgressBar *pb);
ProgressBar setSymbol(ProgressBar *pb, char symbol);

bool getShowPercent(ProgressBar *pb);

ProgressBar showCount(ProgressBar *pb, bool show);

bool getShowCount(ProgressBar *pb);

ProgressBar setCompletedText(ProgressBar *pb, char* text);

char* getCompletedText(ProgressBar *pb);

ProgressBar setStartEndSymbols(ProgressBar *pb, char start, char end);

char getStartSymbol(ProgressBar *pb);

char getEndSymbol(ProgressBar *pb);

ProgressBar setCustomFormat(ProgressBar *pb, char* format);

char* getCustomFormat(ProgressBar *pb);

ProgressBar pb_tick(ProgressBar *pb);

void pb_print(ProgressBar *pb);
    
#ifdef __cplusplus
}
#endif

#endif
