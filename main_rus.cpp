#include <iostream>
#include <vector>
#include <complex>
#include <gmpxx.h>
#include <string>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <fstream>

#include "rus/rus.h"
namespace po = boost::program_options;
namespace btm = boost::timer;
using std::cout;
using std::endl;
using std::string;
static void print_about_message()
{
    cout << "Copyright (c) 2013 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com" << endl
         << endl;
    cout << "SQCT is free software: you can redistribute it and/or modify" << endl;
    cout << "it under the terms of the GNU Lesser General Public License as published by" << endl;
    cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
    cout << "(at your option) any later version." << endl;
    cout << "" << endl;
    cout << "SQCT is distributed in the hope that it will be useful," << endl;
    cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
    cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
    cout << "GNU Lesser General Public License for more details." << endl;
    cout << "" << endl;
    cout << "You should have received a copy of the GNU Lesser General Public License" << endl;
    cout << "along with SQCT.  If not, see <http://www.gnu.org/licenses/>." << endl;
    cout << "" << endl;
}

int main(int ac, char *av[])
{
    string help_topic;
    string output_file_name = "out.qasm";
    string theta_str;
    int debug_level;
    int effort;
    mpf_class epsilon;
    string crit;
    int iterations;
    try
    {

        po::options_description desc("Allowed options");
        desc.add_options()

            ("help,H", po::value<string>(&help_topic)->implicit_value(""),
             "Produce help message, see help <option name> for more details "
             "about specific option.")

                ("theta,T", po::value<string>(&(theta_str))->default_value("pi/16"),
                 "Theta value for the unitary approximation. Example: -T pi/16")

                    ("output,O", po::value<string>(&(output_file_name))->implicit_value("out.qasm"),
                     "Output file name. Example: -O out.qasm")

                        ("debuglevel,D", po::value<int>(&debug_level)->default_value(0), "Debug level")

                            ("effort,F", po::value<int>(&effort)->default_value(0), "Effort level(>=1)")

                                ("epsilon,E", po::value<mpf_class>(&epsilon)->default_value(1e-10), "Epsilon value for the unitary approximation. Example: -e 1e-10")

                                    ("criterion,C", po::value<string>(&crit)->default_value("t-count"), "Criterion for selecting the best circuit. Example: -c t-count")

                                        ("iterations,I", po::value<int>(&iterations)->implicit_value(10000), "Number of iterations in PSLQ stage. Default: 10000*effort. Example: -I 10000")

                                            ("about", "Information about the program.");

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("about"))
        {
            print_about_message();
            return 0;
        }

        // if (vm.count("help"))
        // {
        //     print_help(help_topic);
        //     return 0;
        // }

        {
            RUS::CRITERION criterion;
            if (crit == "g-count")
                criterion = RUS::CRITERION::G_COUNT;
            else
                criterion = RUS::CRITERION::T_COUNT;
            mpf_class theta;
            theta = RUS::parse_theta(theta_str);

            if (debug_level >= 2)
                cout << "Theta: " << theta << endl;

            if (vm.count("iterations") == 0)
                iterations = 10000 * effort;

            RUS rus(debug_level, epsilon, effort, iterations, criterion);

            circuit res;
            btm::cpu_timer timer;

            rus.run(theta, res);
            if (debug_level >= 1)
                cout << "Time: " << timer.format() << endl;
            std::ofstream ofs(output_file_name);
            QasmGenerator::to_qasm(res, ofs, "", "", false);
        }
    }
    catch (std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Exception of unknown type!" << std::endl;
        return 1;
    }
}