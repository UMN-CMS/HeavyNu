#include <map>
#include <string>
#include <vector>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "systematics.h"
#include <algorithm>

static const int NMASSBIN = 20;

const char* SystematicsDB::GAMMASTATS = "GAMMASTATS";

static std::string to_upper(const std::string& item)
{
    std::string retval;
    for (std::string::const_iterator i = item.begin(); i != item.end(); i++)
        retval.push_back(toupper(*i));
    return retval;
}

SystematicsDB::SystematicsDB() : m_mode(0){ }

static std::string makeKey(const std::string& process, const std::string& systName)
{
    return to_upper(process + "-" + systName);
}

void SystematicsDB::load(const std::string& systdb)
{
    FILE* f = fopen(systdb.c_str(), "r");
    if (f == 0)
    {
        fprintf(stderr, "Unable to open systematics file '%s'\n", systdb.c_str());
        return;
    }

    char buffer[1024];
    float vals[NMASSBIN];
    char proc[100], syst[100];
    while (!feof(f))
    {
        buffer[0] = 0;
        fgets(buffer, 1000, f);
        if (strchr(buffer, '#') != 0) *(strchr(buffer, '#')) = 0;
        while (strchr(buffer, ',') != 0) *(strchr(buffer, ',')) = '\t';

        int matched = sscanf(buffer, "%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", proc, syst,
                             vals, vals + 1, vals + 2, vals + 3, vals + 4,
                             vals + 5, vals + 6, vals + 7, vals + 8, vals + 9, vals + 10,
                             vals + 11, vals + 12, vals + 13, vals + 14, vals + 15, vals + 16,
                             vals + 17, vals + 18, vals + 19);
        if (matched >= 3)
        {
            std::string key = makeKey(proc, syst);
            DBitem item;
            for (int i = 0; i < matched - 2; i++) item.values[i] = vals[i];
            for (int i = matched - 2; i < NMASSBIN; i++) item.values[i] = vals[matched - 3];
            m_db.insert(std::pair<std::string, DBitem>(key, item));
            m_processNames.insert(to_upper(proc));
        }
    }
    fclose(f);
}

void SystematicsDB::dump() const
{
    for (std::map<std::string, DBitem>::const_iterator i = m_db.begin(); i != m_db.end(); i++)
    {
        printf("'%s',%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", i->first.c_str(),
               i->second.values[0], i->second.values[1],
               i->second.values[2], i->second.values[3],
               i->second.values[4], i->second.values[5],
               i->second.values[6], i->second.values[7],
               i->second.values[8], i->second.values[9],
               i->second.values[10]);
    }
}

double SystematicsDB::getSystematic(const std::string& systName, const std::string& process, int imassbin) const
{
    std::string key = makeKey(process, systName);

    // No PDF error for signal?
    if ((m_mode & NO_PDF_FOR_SIGNAL) != 0 && key.find("PDF") != std::string::npos && key.find("SIGNAL") != std::string::npos) return 0;

    std::map<std::string, DBitem>::const_iterator i = m_db.find(key);
    if (i != m_db.end() && imassbin >= 0 && imassbin < NMASSBIN)
        return i->second.values[imassbin];
    else
    {
        // printf("No %s (%d)\n",key.c_str(),imassbin);
        return 0.0;
    }
}

std::vector<std::string> SystematicsDB::getSystematicsList() const
{
    return m_finalsystematics;
}

void SystematicsDB::setSimpleSystematic(const std::string& systName)
{
    m_finalsystematics.push_back(systName);
}

void SystematicsDB::defineSingleChannelSyst(const std::string& systName, const std::string& process, const std::vector<std::string>& contents)
{
    std::string key = makeKey(process, systName);
    DBitem dbi;
    // clear
    for (int i = 0; i < NMASSBIN; i++) dbi.values[i] = 0;
    // iterate over inputs and sum in quadrature
    for (int i = 0; i < NMASSBIN; i++)
    {
        for (std::vector<std::string>::const_iterator sourcesyst = contents.begin(); sourcesyst != contents.end(); sourcesyst++)
        {
            double asyst = getSystematic(*sourcesyst, process, i);
            //if(process == "ZJETS") printf("%s\t%s\t%i\t%f\n", sourcesyst->c_str(), process.c_str(), i, asyst);
            if (asyst < 1) continue; // no contribution
            asyst -= 1;
            dbi.values[i] += pow(asyst, 2);
        }
        dbi.values[i] = 1 + sqrt(dbi.values[i]);
        //if(process == "ZJETS") printf("%s\t%s\t%i\t%f\n", systName.c_str(), process.c_str(), i, dbi.values[i]);
    }

    m_db.insert(std::pair<std::string, DBitem>(key, dbi));
    if (std::find(m_finalsystematics.begin(), m_finalsystematics.end(), systName) == m_finalsystematics.end())
        m_finalsystematics.push_back(systName);
}

void SystematicsDB::defineCommonSyst(const std::string& systName, const std::vector<std::string>& contents)
{
    for (std::set<std::string>::const_iterator i = m_processNames.begin(); i != m_processNames.end(); i++)
        defineSingleChannelSyst(systName, *i, contents);
}

void SystematicsDB::defineSignalSyst(const std::string& systName, const std::vector<std::string>& contents)
{
    for (std::set<std::string>::const_iterator i = m_processNames.begin(); i != m_processNames.end(); i++)
    {
        if (i->find("SIGNAL") != std::string::npos)
            defineSingleChannelSyst(systName, *i, contents);
    }
}

void SystematicsDB::standardSystematics(const std::string& which)
{
    std::vector<std::string> systContents;

    if (which == "ELEC")
    {
        systContents.push_back("NORM");

        defineSingleChannelSyst("TTONLY", "TTJETS", systContents);
        defineSingleChannelSyst("ZJONLY", "ZJETS", systContents);
        defineSingleChannelSyst("OTHERONLY", "OTHER", systContents);

        systContents.clear();
        systContents.push_back("ZJSHAPE");
        defineSingleChannelSyst("ZJSHAPE", "ZJETS", systContents);

        systContents.clear();
        systContents.push_back("TTSHAPE");
        defineSingleChannelSyst("TTSHAPE", "TTJETS", systContents);

        systContents.clear();
        systContents.push_back("MCSTATS");
        defineSignalSyst("SIGONLY", systContents);

        systContents.clear();
        systContents.push_back("JES");
        systContents.push_back("JER");
        systContents.push_back("EES");
        systContents.push_back("EID");
        systContents.push_back("PU");
        systContents.push_back("TRIG");
        systContents.push_back("ER");
        //systContents.push_back("MODEL");
        defineCommonSyst("RECOID", systContents);

        systContents.clear();
        systContents.push_back("PDF");
        systContents.push_back("REN");
        systContents.push_back("FACT");
        systContents.push_back("ISRFSR");
        defineCommonSyst("PDFSCALE", systContents);

        setSimpleSystematic(GAMMASTATS);
        setSimpleSystematic("LUMI");
    }
    else if(which == "MUON")
    {
        systContents.push_back("NORM");

        defineSingleChannelSyst("TTONLY", "TTJETS", systContents);
        defineSingleChannelSyst("ZJONLY", "ZJETS", systContents);
        defineSingleChannelSyst("OTHERONLY", "OTHER", systContents);

        systContents.clear();
        systContents.push_back("ZJSHAPE");
        defineSingleChannelSyst("ZJSHAPE", "ZJETS", systContents);

        systContents.clear();
        systContents.push_back("TTSHAPE");
        defineSingleChannelSyst("TTSHAPE", "TTJETS", systContents);

        systContents.clear();
        systContents.push_back("MCSTATS");
        defineSignalSyst("SIGONLY", systContents);

        systContents.clear();
        systContents.push_back("JES");
        systContents.push_back("JER");
        systContents.push_back("MER");
        systContents.push_back("MID");
        systContents.push_back("PU");
        systContents.push_back("TRIG");
        defineCommonSyst("RECOID", systContents);

        systContents.clear();
        systContents.push_back("PDF");
        systContents.push_back("REN");
        systContents.push_back("FACT");
        systContents.push_back("ISRFSR");
        //systContents.push_back("MODEL");
        defineCommonSyst("PDFSCALE", systContents);

        setSimpleSystematic(GAMMASTATS);
        setSimpleSystematic("LUMI");
    }
    else // e-mu combination 
    {
        if(which == "EMUE")
        {
            // Electron specific systematics 
            systContents.push_back("NORM");

            defineSingleChannelSyst("TTONLYE", "TTJETS", systContents);
            defineSingleChannelSyst("ZJONLYE", "ZJETS", systContents);
            defineSingleChannelSyst("OTHERONLYE", "OTHER", systContents);

            systContents.clear();
            systContents.push_back("ZJSHAPE");
            defineSingleChannelSyst("ZJSHAPEE", "ZJETS", systContents);

            systContents.clear();
            systContents.push_back("MCSTATS");
            defineSignalSyst("SIGONLYE", systContents);

            systContents.clear();
            systContents.push_back("EES");
            systContents.push_back("EID");
            systContents.push_back("TRIG");
            systContents.push_back("ER");
            //systContents.push_back("MODEL");
            defineCommonSyst("RECOIDE", systContents);
        }

        if(which == "EMUM")
        {
            // Muon specific systematics 
            systContents.push_back("NORM");

            defineSingleChannelSyst("TTONLYM", "TTJETS", systContents);
            defineSingleChannelSyst("ZJONLYM", "ZJETS", systContents);
            defineSingleChannelSyst("OTHERONLYM", "OTHER", systContents);

            systContents.clear();
            systContents.push_back("ZJSHAPE");
            defineSingleChannelSyst("ZJSHAPEM", "ZJETS", systContents);

            systContents.clear();
            systContents.push_back("MCSTATS");
            defineSignalSyst("SIGONLYM", systContents);

            systContents.clear();
            systContents.push_back("MER");
            systContents.push_back("MID");
            systContents.push_back("TRIG");
            defineCommonSyst("RECOIDM", systContents);
        }

        // Common (ie corrolated) systematics)

        systContents.clear();
        systContents.push_back("TTSHAPE");
        defineSingleChannelSyst("TTSHAPE", "TTJETS", systContents);

        systContents.clear();
        systContents.push_back("JES");
        systContents.push_back("JER");
        systContents.push_back("PU");
        defineCommonSyst("RECOIDC", systContents);

        setSimpleSystematic(GAMMASTATS);
        setSimpleSystematic("LUMI");

        systContents.clear();
        systContents.push_back("PDF");
        systContents.push_back("REN");
        systContents.push_back("FACT");
        systContents.push_back("ISRFSR");
        defineCommonSyst("PDFSCALE", systContents);
    }
}

