/*
  GTfold: compute minimum free energy of RNA secondary structure
  Copyright (C) 2008  David A. Bader
  http://www.cc.gatech.edu/~bader

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _LOADER_H
#define _LOADER_H

#include <string>

#include "constants.h"
#include "data.h"
#include "helper-structs.h"

void readThermodynamicParameters(const char *userdatadir, ParameterVector params);

int initStackValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initMiscloopValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initDangleValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initLoopValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initTstkhValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initTstkiValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initTloopValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initInt21Values(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initInt22Values(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int initInt11Values(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int	initTstkmValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int	initTstkeValues(const std::string& fileName, const std::string& dirPath, ParameterVector params);
int	initTstk23Values(const std::string& fileName, const std::string& dirPath, ParameterVector params);

extern std::string EN_DATADIR;

#endif
