#include "Signal_Filter.h"
#include<fstream>
#include <vector>

vector<string> split(const string& str, const string& pattern)
{
    vector<string> ret;
    if(pattern.empty()) return ret;
    size_t start=0,index=str.find_first_of(pattern,0);
    while(index!=str.npos)
    {
        if(start!=index)
            ret.push_back(str.substr(start,index-start));
        start=index+1;
        index=str.find_first_of(pattern,start);
    }
    if(!str.substr(start).empty())
        ret.push_back(str.substr(start));
    return ret;
}

int main() {
    Signal_Filter sig_filter;
    //    std::string btype,double Wn,int order
    double fps = 30;
    float low = 3.0 / 60.0, high = 120.0 / 60.0;
    int order = 3;
    std::string type = "bandpass";

    std::string file_path = "/media/li/b806bc78-4cbd-4e31-ba5c-e0212d292e73/clion_project/ButterFilter/data.txt";
    fstream txt_file(file_path);
    string s;
    int data_len = 300;
    float x[data_len];
    while(getline(txt_file,s))
    {
        std::vector<string> temp= split(s,",");
        for (int i = 0; i < temp.size(); ++i) {
            x[i] = atof(temp[i].c_str());
        }
    }

    sig_filter.init(order,data_len,low,high,fps);
    float y[data_len];
    sig_filter.filtfilt(x,data_len,y);
    for (int i = 0; i < data_len; ++i) {
        printf("%f,",y[i]);
    }
    std::cout<<"\n";
    return 0;
}
