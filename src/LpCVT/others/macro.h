//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-1.
//

#ifndef MACRO_H
#define MACRO_H

#define LPCVT_VERBOSE
#ifdef LPCVT_VERBOSE
#	define VERBOSE(x) std::cout << "[" << __FUNCTION__ << "]" << " " << x  << std::endl
#	define WARNING(x) std::cerr << "[File:" << __FILE__ << "]" << "[Line:" << __LINE__ << "]" << " " << x << std::endl
#else
#	define VERBOSE(x)
#	define WARNING(x)
#endif

#endif //MACRO_H
