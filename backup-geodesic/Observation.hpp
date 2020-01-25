/**
 *  @2019 Kazuki IWANAGA, the University of Tsukuba.
 *  All rights reserved.
 */

#ifndef ___OBSERVATION___
#define ___OBSERVATION___

struct CSV_Observer {
    using state = BH_System::state;
    std::ofstream fout;
    CSV_Observer(const std::string& filename) : fout(filename, std::ios::app) {
        // fout << "t,r,theta,phi,p_t,p_r,p_theta,p_phi" << std::endl;
    }

    void data2csv(const std::array<double, 8>& data) {
        fout << data[0] << "," << data[1] << "," << data[2] << ","
             << data[3] << "," << data[4] << "," << data[5] << ","
             << data[6] << "," << data[7] << std::endl; 
    }

    void operator()(const state& data, double lambda){
        // if ((data[1] < 0.01) || (data[1] > 10000.0)) {
            fout << data[0] << "," << data[1] << "," << data[2] << ","
                << data[3] << "," << data[4] << "," << data[5] << ","
                << data[6] << "," << data[7] << std::endl;
        // }
    }
};

#endif // !___OBSERVATION___
