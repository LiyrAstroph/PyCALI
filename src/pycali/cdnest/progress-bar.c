/*
 * Progress Bar
 * A simple progress bar implementation in C
 * 
 * Adapted from https://github.com/ankddev/progress-bar-c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "progress-bar.h"

ProgressBar* pb_alloc()
{
    ProgressBar* pb = (ProgressBar *)malloc(sizeof(ProgressBar));
    return pb;
}

void pb_free(ProgressBar *pb)
{
    free(pb);
}

void pb_init(ProgressBar *pb, char symbol, int length, int total) {
    pb->symbol = symbol;
    pb->length = length;
    pb->progress = 0;
    pb->showPercent = false;
    pb->total = total;
    pb->startSymbol = '[';
    pb->endSymbol = ']';
    pb->completedText = NULL;
    pb->format = "{bar} {percent} {count}";
}

ProgressBar pb_update(ProgressBar *pb, int progress) {
    pb->progress = progress;
    return *pb;
}

ProgressBar showPercent(ProgressBar *pb, bool show) {
    pb->showPercent = show;
    return *pb;
}

char getSymbol(ProgressBar *pb) {
    return pb->symbol;
}

ProgressBar setSymbol(ProgressBar *pb, char symbol) {
    pb->symbol = symbol;
    return *pb;
}

bool getShowPercent(ProgressBar *pb) {
    return pb->showPercent;
}

ProgressBar showCount(ProgressBar *pb, bool show) {
    pb->showCount = show;
    return *pb;
}

bool getShowCount(ProgressBar *pb) {
    return pb->showCount;
}

ProgressBar setCompletedText(ProgressBar *pb, char* text) {
    pb->completedText = text;
    return *pb;
}

char* getCompletedText(ProgressBar *pb) {
    return pb->completedText;
}

ProgressBar setStartEndSymbols(ProgressBar *pb, char start, char end) {
    pb->startSymbol = start;
    pb->endSymbol = end;
    return *pb;
}

char getStartSymbol(ProgressBar *pb) {
    return pb->startSymbol;
}

char getEndSymbol(ProgressBar *pb) {
    return pb->endSymbol;
}

ProgressBar setCustomFormat(ProgressBar *pb, char* format) {
    pb->format = format;
    return *pb;
}

char* getCustomFormat(ProgressBar *pb) {
    return pb->format;
}

ProgressBar pb_tick(ProgressBar *pb) {
    pb->progress++;
    return *pb;
}

void pb_print(ProgressBar *pb) {
    char bar[256] = "";
    char percent_str[32] = "";
    char count_str[32] = "";
    char result[512] = "";
    char *format = pb->format;
    
    // Generate bar component
    sprintf(bar, "%c", pb->startSymbol);
    int scaled_progress = (int)((float)pb->progress * pb->length / pb->total);
    
    for (int i = 0; i < pb->length; i++) {
        if (i < scaled_progress) {
            sprintf(bar + strlen(bar), "%c", pb->symbol);
        } else {
            strcat(bar, " ");
        }
    }
    sprintf(bar + strlen(bar), "%c", pb->endSymbol);
    
    // Generate percent component
    int percent = (pb->progress * 100) / pb->total;
    if (pb->showPercent) {
        sprintf(percent_str, "%d%%", percent);
    }
    
    // Generate count component
    if (pb->showCount) {
        sprintf(count_str, "%d/%d", pb->progress, pb->total);
    }
    
    // Process format string
    char *ptr = format;
    while (*ptr) {
        if (*ptr == '{') {
            if (strncmp(ptr, "{bar}", 5) == 0) {
                strcat(result, bar);
                ptr += 5;
            } else if (strncmp(ptr, "{percent}", 9) == 0) {
                if (pb->showPercent) {
                    strcat(result, percent_str);
                }
                ptr += 9;
            } else if (strncmp(ptr, "{count}", 7) == 0) {
                if (pb->showCount) {
                    strcat(result, count_str);
                }
                ptr += 7;
            } else {
                strncat(result, ptr, 1);
                ptr++;
            }
        } else {
            strncat(result, ptr, 1);
            ptr++;
        }
    }
    
    printf("%s", result);
    
    if (percent < 100) {
        printf("\b\r");
    } else {
        if (pb->completedText) {
            printf("\33[2K\r");
            printf("%s\n", pb->completedText);
        } else {
            printf("\n");
        }
    }
}