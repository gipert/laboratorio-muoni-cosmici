/* Progress Bar for c++ loops with internal running variable
 *
 * Author: Luigi Pertoldi
 * Created: 3 dic 2016
 *
 * Notes: The bar must be used when there's no other possible source of output
 *        inside the for loop
 *
 *
 */

#include "progressbar.h"
#include <iostream>

void ProgressBar::Update( int i ) {

	int perc = 0;

	// compute percentage
	perc = i*100. / (nCycles-1);

	// update percentage each unit
	if (perc == savedPerc + 1) {
		// erase the correct  number of characters
		if (perc  < 10)               std::cout << "\b\b"   << perc << "%";
		if (perc == 10)               std::cout << "\b\b"   << perc << "%";
		if (perc  > 10 && perc < 100) std::cout << "\b\b\b" << perc << "%";
		if (perc == 100)              std::cout << "\b\b\b" << perc << "%";
	}

	// update bar every ten units
	if (perc % 2 == 0) {
		// erase trailing percentage characters
		if (perc  < 10)               std::cout << "\b\b\b\b";
		if (perc == 10)               std::cout << "\b\b\b\b\b";
		if (perc  > 10 && perc < 100) std::cout << "\b\b\b\b\b";
		if (perc == 100)              std::cout << "\b\b\b\b\b\b";

		// erase "-"
		for (int j = 0; j < 50 - (perc-1) / 2; j++) std::cout << "\b";
		
        // add one additional "#"
		if (perc == 0) std::cout << "-" << std::flush;
		else           std::cout << "#" << std::flush;
		
        // refill with "-"
		for (int j = 0; j < 50-(perc-1)/2-1; j++) std::cout << "-";
		
        // readd trailing percentage characters
		std::cout << "] " << perc << "%";
	}
	savedPerc = perc;
	std::cout << std::flush;

    return;
}

ProgressBar::~ProgressBar() {}

void ProgressBar::Init() {
	std::cout << "[--------------------------------------------------] 0" << "%" << std::flush;
}
