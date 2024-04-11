#include <iostream>
#include <vector>
#include <complex>
#include <gmpxx.h>
#include <string>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <fstream>

#if __cplusplus >= 201703L
#include <filesystem>
#else
#include <experimental/filesystem>
#endif

#include "rus/rus.h"
namespace po = boost::program_options;
namespace btm = boost::timer;
using std::cout;
using std::endl;
using std::string;
#if __cplusplus >= 201703L
namespace fs = std::filesystem;
#else
namespace fs = std::experimental::filesystem;
#endif
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
    string epsilon_str;
    string crit;
    int iterations;
    string precision_str;
    string qubit_name;
    string ancil_name;
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
                     "Output file(directory) name. Example: -O out.qasm")

                        ("debuglevel,D", po::value<int>(&debug_level)->default_value(0), "Debug level")

                            ("effort,F", po::value<int>(&effort)->default_value(0), "Effort level(>=1)")

                                ("epsilon,E", po::value<string>(&epsilon_str)->default_value("1e-10"), "Epsilon value for the unitary approximation. Example: -e 1e-10")

                                    ("criterion,C", po::value<string>(&crit)->default_value("t-count"), "Criterion for selecting the best circuit. Example: -c t-count")

                                        ("iterations,I", po::value<int>(&iterations)->implicit_value(10000), "Number of iterations in PSLQ stage. Default: 10000*effort. Example: -I 10000")

                                            ("precision,P", po::value<string>(&precision_str)->default_value("pi/128"), "Precision for the library generation. Example: -p 1e-2")

                                                ("qubit-name", po::value<string>(&qubit_name)->default_value(""), "Name of the qubit. Default: q[0]")

                                                    ("ancil-name", po::value<string>(&ancil_name)->default_value(""), "Name of the ancillary qubit. Default: a[0]")

                                                        ("database", "Produce database of gates. Based on precision and epsilon.")

                                                            ("about", "Information about the program.");

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("about"))
        {
            print_about_message();
            return 0;
        }

        if (vm.count("database"))
        {
            if (!vm.count("precision") || !vm.count("epsilon"))
            {
                std::cout << "Build database requires precision and epsilon" << std::endl;
                return 1;
            }
        }

        // if (vm.count("help"))
        // {
        //     print_help(help_topic);
        //     return 0;
        // }

        RUS::CRITERION criterion;
        if (crit == "g-count")
            criterion = RUS::CRITERION::G_COUNT;
        else
            criterion = RUS::CRITERION::T_COUNT;

        if (vm.count("iterations") == 0)
            iterations = 10000 * effort;

        mpf_class epsilon;
        epsilon = RUS::parse_theta(epsilon_str);
        if (!vm.count("database"))
        {
            mpf_class theta;
            theta = RUS::parse_theta(theta_str);

            if (debug_level >= 2)
                cout << "Theta: " << theta << endl;

            RUS rus(debug_level, epsilon, effort, iterations, criterion);

            circuit res;
            btm::cpu_timer timer;

            rus.run(theta, res);
            if (debug_level >= 1)
                cout << "Time: " << timer.format() << endl;
            std::ofstream ofs(output_file_name);
            QasmGenerator::to_rus_qasm(res, ofs, "", "", false);
        }
        else
        {
            string output_folder_name = output_file_name;
            fs::create_directories(output_folder_name);
            mpf_class precision;
            precision = RUS::parse_theta(precision_str);

            // ? change to +-pi/4 or pi/2~0 ?
            auto unit_count = mpf_class(ceil(M_PI / precision)).get_si(); // +-(pi-precision)
            {
                string t1, t2;
                t1 = precision_str.find("/") == string::npos ? precision_str : precision_str.replace(precision_str.find("/"), 1, "|");
                t2 = epsilon_str.find("/") == string::npos ? epsilon_str : epsilon_str.replace(epsilon_str.find("/"), 1, "|");
                output_file_name = t1 + "_" + t2;
            }

            RUS rus(debug_level, epsilon, effort, iterations, criterion);
            for (int i = -unit_count + 1; i < unit_count; i++)
            {
                mpf_class theta = i * precision;
                if (debug_level >= 2)
                    cout << "Theta: " << theta << endl;
                circuit res;
                cout << "output: " << (output_folder_name + "/out" + std::to_string(i) + "_" + output_file_name + ".qasm") << endl;
                std::ofstream ofs(output_folder_name + "/out" + std::to_string(i) + "_" + output_file_name + ".qasm");
                if (abs(theta - (M_PI / 4) * mpf_class(theta / (M_PI / 4)).get_si()) < epsilon)
                {
                    if (abs(theta) < epsilon)
                        res.push_back(gateLibrary::Id);
                    else if (abs(theta) - M_PI / 4 < epsilon)
                        res.push_back(theta > 0 ? gateLibrary::T : gateLibrary::Td);
                    else if (abs(theta) - M_PI / 2 < epsilon)
                        res.push_back(theta > 0 ? gateLibrary::P : gateLibrary::Pd);
                    else if (abs(theta) - 3 * M_PI / 4 < epsilon)
                    {
                        res.push_back(theta > 0 ? gateLibrary::T : gateLibrary::Td);
                        res.push_back(theta > 0 ? gateLibrary::P : gateLibrary::Pd);
                    }
                    else if (abs(theta) - M_PI < epsilon)
                        res.push_back(gateLibrary::Z);
                    else
                        throw std::runtime_error("Unknown Status");
                    QasmGenerator::to_qasm(res, ofs, qubit_name, false);
                    continue;
                }
                else
                    rus.run(theta, res);
                QasmGenerator::to_rus_qasm(res, ofs, qubit_name, ancil_name, false);
            }
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