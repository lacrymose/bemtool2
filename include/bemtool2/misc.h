#ifndef MISC_H
#define MISC_H

//===================================================================
//
//  Copyright 2015 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================


#include <ctime>

class progress{
    
private:
    const char* title;
    const int length;
    int       it;
    clock_t   t0;
    
public:
    progress(const char* aff,const int& l): title(aff), length(l) {
        t0 = clock();
        it=0;

        std::cout << "\r";
        std::cout << title << ": \t";
        std::cout << 0 << "%";
        std::cout.flush(); }
    
    void operator++(int n){
        it++;

        if( it%int(length/100) == 0){
            std::cout << "\r";
            std::cout << title << ": \t";;
            std::cout << int(it*100./length) << "%";
            std::cout.flush();
        }
    }
    
    void end(){
        t0 = clock()-t0;
        time_t now; time(&now);
        std::cout << "\r";
        std::cout << title << ": \t";
        std::cout << ((float)t0)/CLOCKS_PER_SEC << " sec." << std::endl;
    }
    
};


#endif
