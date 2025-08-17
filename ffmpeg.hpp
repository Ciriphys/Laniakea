#pragma once

#include "solver.hpp"

#include <atomic>

extern std::atomic<bool> streaming_active;

void start_ffmpeg(settings config);
void stream_rgb_frame(const uint8_t* rgb_frame, int width, int height);
void apply_colormap_rgb(uint8_t* out_rgb, const double* in, int size, double vmin, double vmax);
void apply_colormap_rgb_log(uint8_t* out_rgb, const double* in, int size, double vmin, double vmax);
void apply_colormap_rgb_gamma(uint8_t* out_rgb, const double* in, int size, double vmin, double vmax, double gamma);
void draw_text(uint8_t* out_rgb, int width, int height, int x, int y, const char* text, uint8_t r, uint8_t g, uint8_t b);
void end_ffmpeg();
