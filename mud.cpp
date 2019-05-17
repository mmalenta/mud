#include <algorithm>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "boost/program_options.hpp"

#define TPACKET 128
#define FPACKET 64

namespace po = boost::program_options;

struct Config {

    std::string infile;
    std::string outfile;

    double netstart;
    double netlen;
    double netperc;

    double lvstart;
    double lvlen;
    double lvoff;
    double lvband;

    bool netrand;
    bool lvrand;
};

void read_config(int argc, char *argv[], Config &inconfig) {

    po::options_description information("Information");
    information.add_options()
                ("help,h", "prints this message")
                ("version", "prints the code version");

    po::options_description filenames("Files:");
    filenames.add_options()
                ("input, i", po::value<std::string>(&inconfig.infile), "input file name")
                ("output, o", po::value<std::string>(&inconfig.outfile), "output file name (<input>.muddy if not specified)");

    po::options_description network("Network");
    network.add_options()
                ("netstart", po::value<double>(&inconfig.netstart), "dropped packets start (seconds)")
                ("netlen", po::value<double>(&inconfig.netlen), "dropped packets length (seconds)")
                ("netperc", po::value<double>(&inconfig.netperc), "percentage of dropped packets")
                ("netrand", po::bool_switch(&inconfig.netrand) -> default_value(false), "random start, length and percantage - overrides any options provided explicitly");

    po::options_description levels("Levels");
    levels.add_options()
                ("lvstart", po::value<double>(&inconfig.lvstart), "levels mismatch start (seconds)")
                ("lvlen", po::value<double>(&inconfig.lvlen), "levels mismatch length (seconds)")
                ("lvoff", po::value<double>(&inconfig.lvoff), "levels offset")
                ("lvband", po::value<double>(&inconfig.lvband), "levels mismatch frequency bands")
                ("lvrand", po::bool_switch(&inconfig.lvrand) -> default_value(false), "random start, length, offset and bands - overrides any options provided explicitly");

    po::options_description combined("All options");
    combined.add(information).add(filenames).add(network).add(levels);

    po::variables_map variables;
    store(po::command_line_parser(argc, argv).options(combined).run(), variables);
    notify(variables);

    if (argc < 2) {
        std::ostringstream sstream("\n\nNot enough command line arguments!\n", std::ios_base::app);
        sstream << combined;
        throw std::runtime_error(sstream.str()); 
    }

    if (variables.count("help")) {
        std::cout << combined << std::endl;
    }

    if (variables.count("version")) {
        std::cout << "MUD, Messing Up the Data, version 0.1.whocares" << std::endl;
    }

    if (!variables.count("input")) {
        throw std::runtime_error("Did not specify the input type / filename");
    }

    if (!variables.count("output")) {
        inconfig.outfile = inconfig.infile + std::string(".muddy");
    }

}

int read_header(std::ifstream &infile, char *&header) {

    infile.seekg(0, infile.beg);


    const std::string header_end_str = "HEADER_END";
    int str_len = header_end_str.length();
    int rewind_len = -1 * (str_len - 1);
    int header_len = 0;
    char *read_chars = new char[str_len];

    while(true) {
        infile.read(read_chars, str_len);
        if (std::string(read_chars) == header_end_str) {
            header_len = infile.tellg();
            std::cout << "Header length: " << header_len << "B\n";
            break;
        }
        infile.seekg(rewind_len, infile.cur);
    }

    infile.seekg(0, infile.beg);
    header = new char[header_len];
    infile.read(reinterpret_cast<char*>(header), header_len);

    return header_len;

}

template <class ParamType>
ParamType read_value(char* header, std::string param) {

    ParamType param_value;

    int param_len = param.length();
    char *head_ptr = header;
    char *read_val = new char[sizeof(ParamType)];
    while(true) {
        if (std::string(head_ptr, head_ptr + param_len) == param) {
            std::copy(head_ptr + param_len, head_ptr + param_len + sizeof(ParamType), read_val);
            param_value = *(reinterpret_cast<ParamType*>(read_val));
            std::cout << "Found parameter " << param << ": " << param_value << "\n";
            break;
        }
        head_ptr++;
    }

    return param_value;

}

int main(int argc, char *argv[]) {

    try {
        
        // NOTE: Parse the command line options
        Config inconfig = {};
        read_config(argc, argv, inconfig);

        // NOTE: Open the input file
        std::string instring = inconfig.infile;
        std::string outstring = inconfig.outfile;
        std::cout << "Processing file " << instring << "...\n";
        std::ifstream infile(instring.c_str(), std::ios_base::binary);
        if (!infile) {
            throw std::invalid_argument("cannot open input file " + instring); 
        }

        // NOTE: Read the header and all the necessary values
        char *header;
        int headlen = read_header(infile, header);
        int nchans = read_value<int>(header, std::string("nchans"));
        int nbits = read_value<int>(header, std::string("nbits"));
        double tsamp = read_value<double>(header, std::string("tsamp"));

        infile.seekg(0, infile.end);
        size_t samples = (static_cast<size_t>(infile.tellg()) - static_cast<size_t>(headlen)) / static_cast<size_t>(nchans);
        if ((samples * nchans) != (static_cast<size_t>(infile.tellg()) - static_cast<size_t>(headlen))) {
            throw std::logic_error("did not read an integer number of samples. Is the input file corrupted?");
        }
        std::cout << "Read " << samples << " time samples\n";
        infile.seekg(headlen, infile.beg);

        // NOTE: Check if we have enough samples to corrupt what user has requested
        // NOTE: Make sure we start corrupting the data at the boundary of what is sent in a single packet
        size_t net_samples_start = static_cast<size_t>((inconfig.netstart / tsamp) / static_cast<float>(TPACKET)) * TPACKET;
        // NOTE: Make sure the length we corrupt is the multuiple of number of time samples we sent in a single packet
        size_t net_samples_length = static_cast<size_t>(((inconfig.netlen / tsamp) + static_cast<float>(TPACKET)) / static_cast<float>(TPACKET)) * TPACKET; 
        size_t net_samples_end = net_samples_start + net_samples_length;

        if (net_samples_end > samples) {
            throw std::runtime_error("Network packets dropping goes beyond the end of the file! Please review the start time and length");
        }

        size_t net_skip_bytes = net_samples_start * nchans * nbits / 8;
        size_t net_read_bytes = net_samples_length * nchans * nbits / 8;
        // NOTE: Read the chunk of the data to be network corrupted
        unsigned char *net_data = new unsigned char[net_read_bytes];
        infile.seekg(headlen + net_skip_bytes, infile.beg);
        infile.read(reinterpret_cast<char*>(net_data), net_read_bytes);

        // NOTE: Very rough estimate of how many packets were sent in order to produce what we have just read
        size_t time_packets = net_samples_length / TPACKET;
        size_t freq_packets = nchans / FPACKET;
        size_t total_packets = time_packets * freq_packets;
        size_t corrupted_packets = static_cast<size_t>(inconfig.netperc / 100.0 * total_packets);

        std::cout << "There are a total of " << total_packets << " packets in the requested chunk\n"
                    << "Will corrupt " << corrupted_packets << " out of those (" << inconfig.netperc << "%)\n";

        // NOTE: 'Choose' which packets to network corrupt
        std::random_device rd;
        std::mt19937 net_mt_gen(rd());

        std::vector<size_t> full_packet_id(total_packets);
        std::iota(full_packet_id.begin(), full_packet_id.end(), 0);
        std::shuffle(full_packet_id.begin(), full_packet_id.end(), net_mt_gen);
        std::vector<size_t> net_packet_id(full_packet_id.begin(), full_packet_id.begin() + corrupted_packets);
        std::sort(net_packet_id.begin(), net_packet_id.end());

        for (auto pid : net_packet_id) {
            std::cout << pid << " ";
        }

        std::cout << std::endl;

        size_t time_chunk = 0;
        size_t freq_chunk = 0;
        size_t time_skip = 0;
        size_t freq_skip = 0;
        size_t samp_idx = 0;

        // NOTE: Corrupting happens here
        for (size_t ipacket = 0; ipacket < corrupted_packets; ++ipacket) {
            time_chunk = ipacket / freq_packets; 
            freq_chunk = ipacket % freq_packets;

            freq_skip = freq_chunk * FPACKET;
            time_skip = time_chunk * TPACKET * nchans;

            for (int itime = 0; itime < TPACKET; ++itime) {

                samp_idx = time_skip + itime * nchans;

                for (int ichan = 0; itime < FPACKET; ++itime) {
                    net_data[samp_idx] = 0x00;
                    samp_idx++;
                }

            }

            std::cout << "Corrupted: " << (float)ipacket / (float)corrupted_packets * 100 << "% of data\r";
            std::cout.flush();

        }

        std::cout << std::endl;
        std::cout << "Saving the data into file " << outstring << "..." << std::endl;
        std::ofstream outfile(outstring.c_str(), std::ios_base::binary);
        if (!outfile) {
            throw std::runtime_error("cannot open output file " + outstring); 
        }

        unsigned char *tmp_buff = new unsigned char[TPACKET * FPACKET];
        infile.seekg(0, infile.beg);
        do {
            infile.read(reinterpret_cast<char*>(tmp_buff), TPACKET * FPACKET);
            outfile.write(reinterpret_cast<char*>(tmp_buff), infile.gcount());
        } while (infile.gcount() > 0);

        outfile.seekp(headlen + net_skip_bytes, infile.beg);
        outfile.write(reinterpret_cast<char*>(net_data), net_read_bytes);
        outfile.close();

        delete [] tmp_buff;


    } catch (const std::exception &exc) {
        std::cerr << "Something went wrong: " << exc.what() << std::endl;
    }

}