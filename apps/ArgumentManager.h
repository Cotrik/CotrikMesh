#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

// This is a class that can parse the commnad line arguments we use in COSC 2430 homework.
class ArgumentManager {
private:
    std::map<std::string, std::string> m_argumentMap;
public:
    ArgumentManager() { }
    ArgumentManager(int argc, char *argv[], char delimiter=';');
    ArgumentManager(std::string rawArguments, char delimiter=';');
    void parse(int argc, char *argv[], char delimiter=';');
    void parse(std::string rawArguments, char delimiter=';');
    std::string get(std::string argumentName);
    void get(const std::string key, bool& value);
    void get(const std::string key, int& value);
    void get(const std::string key, float& value);
    void get(const std::string key, double& value);
    void get(const std::string key, std::vector<size_t>& value);
    void get(const std::string key, std::vector<int>& value);
    void get(const std::string key, std::vector<float>& value);
    void get(const std::string key, std::vector<double>& value);
    std::string toString();
    friend std::ostream& operator << (std::ostream &out, ArgumentManager &am);
};

void ArgumentManager::parse(std::string rawArguments, char delimiter) {
    std::stringstream currentArgumentName;
    std::stringstream currentArgumentValue;
    bool argumentNameFinished = false;
    
    for (unsigned int i = 0; i <= rawArguments.length(); i++) {
        if (i == rawArguments.length() || rawArguments[i] == delimiter) {
            if (currentArgumentName.str() != "") m_argumentMap[currentArgumentName.str()] = currentArgumentValue.str();
            // reset
            currentArgumentName.str("");
            currentArgumentValue.str("");
            argumentNameFinished = false;
        } else if (rawArguments[i] == '=') argumentNameFinished = true;
        else {
            if (argumentNameFinished) currentArgumentValue << rawArguments[i];
            else {
                // ignore any spaces in argument names. 
                if (rawArguments[i] == ' ') continue;
                currentArgumentName << rawArguments[i];
            }
        }
    }
}

void ArgumentManager::parse(int argc, char *argv[], char delimiter) {
    if (argc > 1) {
        for (int i=1; i<argc; i++) {
            parse(argv[i], delimiter);
        }
    }
}

ArgumentManager::ArgumentManager(int argc, char *argv[], char delimiter) {
    parse(argc, argv, delimiter);
}

ArgumentManager::ArgumentManager(std::string rawArguments, char delimiter) {
    parse(rawArguments, delimiter);
}

std::string ArgumentManager::get(std::string argumentName) {
    std::map<std::string, std::string>::iterator iter = m_argumentMap.find(argumentName);
    if (iter == m_argumentMap.end()) return "";
    else return iter->second;
}

void ArgumentManager::get(const std::string key, bool& value) {
    auto str = get(key);
    if (!str.empty()) value = str == "false" ? false : true;
}

void ArgumentManager::get(const std::string key, int& value) {
    auto str = get(key);
    if (!str.empty()) value = std::stoi(str);
}

void ArgumentManager::get(const std::string key, float& value) {
    auto str = get(key);
    if (!str.empty()) value = std::stof(str);
}

void ArgumentManager::get(const std::string key, double& value) {
    auto str = get(key);
    if (!str.empty()) value = std::stod(str);
}

void ArgumentManager::get(const std::string key, std::vector<size_t>& value) {
    auto str = get(key);
    std::stringstream ss(str);
    size_t id;
    while (ss >> id) value.push_back(id);
}

void ArgumentManager::get(const std::string key, std::vector<int>& value) {
    auto str = get(key);
    std::stringstream ss(str);
    int id;
    while (ss >> id) value.push_back(id);
}

void ArgumentManager::get(const std::string key, std::vector<float>& value) {
    auto str = get(key);
    std::stringstream ss(str);
    float id;
    while (ss >> id) value.push_back(id);
}

void ArgumentManager::get(const std::string key, std::vector<double>& value) {
    auto str = get(key);
    std::stringstream ss(str);
    double id;
    while (ss >> id) value.push_back(id);
}

std::string ArgumentManager::toString() {
    std::stringstream ss;
    for (std::map<std::string, std::string>::iterator iter = m_argumentMap.begin(); iter != m_argumentMap.end(); iter++) {
        ss << "Argument name: " << iter->first << std::endl;
        ss << "Argument value: " << iter->second << std::endl;
    }
    return ss.str();
}

std::ostream& operator << (std::ostream &out, ArgumentManager &am) {
    out << am.toString();
    return out;
}
