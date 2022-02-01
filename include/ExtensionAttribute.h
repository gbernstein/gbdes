// Classes used by WCSFoF.h to propagate attributes of input catalog files to
// the table describing every catalog extension.
#ifndef EXTENSIONATTRIBUTE_H
#define EXTENSIONATTRIBUTE_H

#include "FTable.h"

class ExtensionAttributeBase {
public:
    enum Status { ReadOnly, WriteOnly, ReadWrite };
    ExtensionAttributeBase(string columnName_, Status s);
    virtual ~ExtensionAttributeBase() {}
    // Get input value for an attribute, return value is whether its column exists
    virtual bool readInputTable(const img::FTable &inTable, long row) = 0;
    // Add an appropriate column to output table
    virtual void makeOutputColumn(img::FTable &outTable) const = 0;
    // Search header, if needed, for value of attribute
    // Returns false if header keyword is missing (in which case default value is assigned)
    virtual bool checkHeader(const img::Header &h) = 0;
    // Save value into output table
    virtual void writeOutputTable(img::FTable &outTable, long row) const = 0;
    string getName() const { return columnName; }

    bool doRead() const { return isInput; }
    bool doWrite() const { return isOutput; }

protected:
    string columnName;
    bool isInput;
    bool isOutput;
};

template <class T>
class ExtensionAttribute : public ExtensionAttributeBase {
public:
    ExtensionAttribute(string columnName_, ExtensionAttributeBase::Status s, const T &defaultValue_ = T());
    virtual ~ExtensionAttribute() {}
    bool readInputTable(const img::FTable &inTable, long row);
    void makeOutputColumn(img::FTable &outTable) const;
    bool checkHeader(const img::Header &h);
    T getValue() const { return value; }
    void setValue(const T &v) { value = v; }
    void writeOutputTable(img::FTable &outTable, long row) const;

private:
    T defaultValue;
    T value;
    bool lookupInHeader;
    string headerKeyword;
};

#endif  // EXTENSIONATTRIBUTE_H
