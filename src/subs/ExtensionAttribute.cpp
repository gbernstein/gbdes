#include "ExtensionAttribute.h"

using namespace img;

ExtensionAttributeBase::ExtensionAttributeBase(string columnName_, Status s)
        : columnName(columnName_), isInput(false), isOutput(false) {
    switch (s) {
        case ReadOnly:
            isInput = true;
            break;
        case WriteOnly:
            isOutput = true;
            break;
        case ReadWrite:
            isInput = true;
            isOutput = true;
            break;
        default:
            break;
    }
}

template <class T>
ExtensionAttribute<T>::ExtensionAttribute(string columnName_, ExtensionAttributeBase::Status s,
                                          const T &defaultValue_)
        : ExtensionAttributeBase(columnName_, s), defaultValue(defaultValue_) {
    value = defaultValue;
}

// Specialize the reading to particular types so we can check for mismatched
// column types
template <>
bool ExtensionAttribute<string>::readInputTable(const FTable &inTable, long row) {
    if (!isInput) return false;
    lookupInHeader = false;
    try {
        inTable.readCell(value, columnName, row);
        if (value.empty()) {
            value = defaultValue;
            return true;
        } else if (value[0] == '@') {
            lookupInHeader = true;
            headerKeyword = value.substr(1);
            value = defaultValue;
        }
        return true;
    } catch (FTableNonExistentColumn &f) {
        value = defaultValue;
        return false;
    }
}

template <>
bool ExtensionAttribute<int>::readInputTable(const FTable &inTable, long row) {
    if (!isInput) return false;
    string vstring;
    lookupInHeader = false;
    try {
        inTable.readCell(value, columnName, row);
        return true;
    } catch (FTableNonExistentColumn &f) {
        value = defaultValue;
        return false;
    } catch (FTableError &f) {
        // Try long column instead of int
        try {
            long j;
            inTable.readCell(j, columnName, row);
            value = j;
            return true;
        } catch (FTableError &f2) {
            // Try reading a string (and just throw if it fails this time
            inTable.readCell(vstring, columnName, row);
        }
    }
    // If we've made it here, we have a string.
    if (vstring.empty()) {
        value = defaultValue;
    } else if (vstring[0] == '@') {
        lookupInHeader = true;
        headerKeyword = vstring.substr(1);
        value = defaultValue;
    } else {
        istringstream iss(vstring);
        iss >> value;
        if (iss.fail()) throw std::runtime_error("Could not convert attribute <" + vstring + ">to int");
    }
    return true;
}

template <>
bool ExtensionAttribute<double>::readInputTable(const FTable &inTable, long row) {
    if (!isInput) return false;
    string vstring;
    lookupInHeader = false;
    try {
        inTable.readCell(value, columnName, row);
        return true;
    } catch (FTableNonExistentColumn &f) {
        value = defaultValue;
        return false;
    } catch (FTableError &f) {
        // Try float column instead of double
        try {
            float x;
            inTable.readCell(x, columnName, row);
            value = x;
            return true;
        } catch (FTableError &f2) {
            // Try reading a string (and just throw if it fails this time
            inTable.readCell(vstring, columnName, row);
        }
    }
    // If we've made it here, we have a string.
    if (vstring.empty()) {
        value = defaultValue;
    } else if (vstring[0] == '@') {
        lookupInHeader = true;
        headerKeyword = vstring.substr(1);
        value = defaultValue;
    } else {
        istringstream iss(vstring);
        iss >> value;
        if (iss.fail()) throw std::runtime_error("Could not convert attribute <" + vstring + ">to double");
    }
    return true;
}

template <class T>
void ExtensionAttribute<T>::makeOutputColumn(FTable &outTable) const {
    if (!isOutput) return;
    vector<T> v;
    /**/ cerr << "Making output column for <" << columnName << ">" << endl;
    outTable.addColumn(v, columnName);
}
template <class T>
void ExtensionAttribute<T>::writeOutputTable(FTable &outTable, long row) const {
    if (!isOutput) return;
    outTable.writeCell(value, columnName, row);
}

template <class T>
bool ExtensionAttribute<T>::checkHeader(const Header &h) {
    if (!isInput) return false;
    if (!lookupInHeader) return true;
    // Return true if a keyword of the correct type exists.
    if (h.getValue(headerKeyword, value)) {
        return true;
    } else {
        value = defaultValue;
        return false;
    }
}

// Specialize the reading from headers, again to look for alternate types of data
template <>
bool ExtensionAttribute<string>::checkHeader(const Header &h) {
    if (!isInput) return false;
    if (!lookupInHeader) return true;
    if (h.getValue(headerKeyword, value)) {
        return true;
    }
    // Try once more, reading an integer and converting to a string
    int i;
    if (!h.getValue(headerKeyword, i)) {
        value = defaultValue;
        return false;
    }
    ostringstream oss;
    oss << i;
    value = oss.str();
    return true;
}

// Specialize the reading from headers, again to look for alternate types of data
template <>
bool ExtensionAttribute<double>::checkHeader(const Header &h) {
    if (!isInput) return false;
    if (!lookupInHeader) return true;
    if (h.getValue(headerKeyword, value)) return true;
    // Try once more, reading a float and converting
    float f;
    if (!h.getValue(headerKeyword, f)) {
        value = defaultValue;
        return false;
    }
    value = f;
    return true;
}

template class ExtensionAttribute<string>;
template class ExtensionAttribute<int>;
template class ExtensionAttribute<double>;
