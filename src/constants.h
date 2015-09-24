/*
 * Copyright
 * 2010 Martin Schreiber
 * 2013 Arash Bakhtiari
 * 2016 Christoph Riesinger, Ayman Saleem
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define TYPE double

#define FLAG_OBSTACLE           (1<<0)
#define FLAG_FLUID              (1<<1)
#define FLAG_VELOCITY_INJECTION (1<<2)
#define FLAG_GHOST_LAYER        (1<<3)
#define FLAG_GHOST_LAYER_BETA   (FLAG_GHOST_LAYER | (1<<4))

#define STORE_VELOCITY          1
#define STORE_DENSITY           1

#define BENCHMARK_OUTPUT_DIR    "output/benchmark"
#define PROFILE_OUTPUT_DIR      "output/profile"
#define VTK_OUTPUT_DIR          "output/vtk"
#define LOG_OUTPUT_DIR          "output/log"
#define LOG_OUTPUT_FILE_PREFIX  "log"

#endif
