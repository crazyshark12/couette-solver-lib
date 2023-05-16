#include "datawriter.h"
DataWriter::DataWriter(std::string pathName_)
{

    directory = pathName_;
    fs::path dataDir = directory / "data";
    fs::remove_all(dataDir);
    fs::create_directories(dataDir);
}

fs::path DataWriter::createTimeDirectory(double time)
{
    fs::path localDir = directory / "data" / std::to_string(time);
    fs::create_directories(localDir);
    return localDir;
}

void DataWriter::writeData(vector<macroParam> data, double time)
{
    fs::path localDir = createTimeDirectory(time);

    ofstream pressure(localDir/"pressure.txt",std::ios::out);
    ofstream velocity(localDir/"velocity.txt",std::ios::out);
    ofstream temp(localDir/"temp.txt",std::ios::out);
    ofstream density(localDir/"density.txt",std::ios::out);
    for(size_t i = 0; i < data.size(); i++)
    {
        pressure<<dh*i<<" "<<data[i].pressure<<endl;
        velocity<<dh*i<<" "<<data[i].velocity<<endl;
        temp<<dh*i<<" "<<data[i].temp<<endl;
        density<<dh*i<<" "<<data[i].density<<endl;
    }
    pressure.close();
    velocity.close();
    temp.close();
}

void DataWriter::setDelta_h(double dh_)
{
    dh = dh_;
}
