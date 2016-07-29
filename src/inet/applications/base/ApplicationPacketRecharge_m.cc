//
// Generated file, do not edit! Created by nedtool 5.0 from inet/applications/base/ApplicationPacketRecharge.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "ApplicationPacketRecharge_m.h"

namespace omnetpp {

// Template pack/unpack rules. They are declared *after* a1l type-specific pack functions for multiple reasons.
// They are in the omnetpp namespace, to allow them to be found by argument-dependent lookup via the cCommBuffer argument

// Packing/unpacking an std::vector
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::vector<T,A>& v)
{
    int n = v.size();
    doParsimPacking(buffer, n);
    for (int i = 0; i < n; i++)
        doParsimPacking(buffer, v[i]);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::vector<T,A>& v)
{
    int n;
    doParsimUnpacking(buffer, n);
    v.resize(n);
    for (int i = 0; i < n; i++)
        doParsimUnpacking(buffer, v[i]);
}

// Packing/unpacking an std::list
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::list<T,A>& l)
{
    doParsimPacking(buffer, (int)l.size());
    for (typename std::list<T,A>::const_iterator it = l.begin(); it != l.end(); ++it)
        doParsimPacking(buffer, (T&)*it);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::list<T,A>& l)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i=0; i<n; i++) {
        l.push_back(T());
        doParsimUnpacking(buffer, l.back());
    }
}

// Packing/unpacking an std::set
template<typename T, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::set<T,Tr,A>& s)
{
    doParsimPacking(buffer, (int)s.size());
    for (typename std::set<T,Tr,A>::const_iterator it = s.begin(); it != s.end(); ++it)
        doParsimPacking(buffer, *it);
}

template<typename T, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::set<T,Tr,A>& s)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i=0; i<n; i++) {
        T x;
        doParsimUnpacking(buffer, x);
        s.insert(x);
    }
}

// Packing/unpacking an std::map
template<typename K, typename V, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::map<K,V,Tr,A>& m)
{
    doParsimPacking(buffer, (int)m.size());
    for (typename std::map<K,V,Tr,A>::const_iterator it = m.begin(); it != m.end(); ++it) {
        doParsimPacking(buffer, it->first);
        doParsimPacking(buffer, it->second);
    }
}

template<typename K, typename V, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::map<K,V,Tr,A>& m)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i=0; i<n; i++) {
        K k; V v;
        doParsimUnpacking(buffer, k);
        doParsimUnpacking(buffer, v);
        m[k] = v;
    }
}

// Default pack/unpack function for arrays
template<typename T>
void doParsimArrayPacking(omnetpp::cCommBuffer *b, const T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimPacking(b, t[i]);
}

template<typename T>
void doParsimArrayUnpacking(omnetpp::cCommBuffer *b, T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimUnpacking(b, t[i]);
}

// Default rule to prevent compiler from choosing base class' doParsimPacking() function
template<typename T>
void doParsimPacking(omnetpp::cCommBuffer *, const T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: no doParsimPacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

template<typename T>
void doParsimUnpacking(omnetpp::cCommBuffer *, T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: no doParsimUnpacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

}  // namespace omnetpp

namespace inet {

// forward
template<typename T, typename A>
std::ostream& operator<<(std::ostream& out, const std::vector<T,A>& vec);

// Template rule which fires if a struct or class doesn't have operator<<
template<typename T>
inline std::ostream& operator<<(std::ostream& out,const T&) {return out;}

// operator<< for std::vector<T>
template<typename T, typename A>
inline std::ostream& operator<<(std::ostream& out, const std::vector<T,A>& vec)
{
    out.put('{');
    for(typename std::vector<T,A>::const_iterator it = vec.begin(); it != vec.end(); ++it)
    {
        if (it != vec.begin()) {
            out.put(','); out.put(' ');
        }
        out << *it;
    }
    out.put('}');
    
    char buf[32];
    sprintf(buf, " (size=%u)", (unsigned int)vec.size());
    out.write(buf, strlen(buf));
    return out;
}

Register_Class(ApplicationPacketRecharge);

ApplicationPacketRecharge::ApplicationPacketRecharge(const char *name, int kind) : ::omnetpp::cPacket(name,kind)
{
    this->sequenceNumber = 0;
    this->appAddr = 0;
    this->batteryLevelAbs = 0;
    this->batteryLevelPerc = 0;
}

ApplicationPacketRecharge::ApplicationPacketRecharge(const ApplicationPacketRecharge& other) : ::omnetpp::cPacket(other)
{
    copy(other);
}

ApplicationPacketRecharge::~ApplicationPacketRecharge()
{
}

ApplicationPacketRecharge& ApplicationPacketRecharge::operator=(const ApplicationPacketRecharge& other)
{
    if (this==&other) return *this;
    ::omnetpp::cPacket::operator=(other);
    copy(other);
    return *this;
}

void ApplicationPacketRecharge::copy(const ApplicationPacketRecharge& other)
{
    this->sequenceNumber = other.sequenceNumber;
    this->addr = other.addr;
    this->appAddr = other.appAddr;
    this->pos = other.pos;
    this->batteryLevelAbs = other.batteryLevelAbs;
    this->batteryLevelPerc = other.batteryLevelPerc;
}

void ApplicationPacketRecharge::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cPacket::parsimPack(b);
    doParsimPacking(b,this->sequenceNumber);
    doParsimPacking(b,this->addr);
    doParsimPacking(b,this->appAddr);
    doParsimPacking(b,this->pos);
    doParsimPacking(b,this->batteryLevelAbs);
    doParsimPacking(b,this->batteryLevelPerc);
}

void ApplicationPacketRecharge::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->sequenceNumber);
    doParsimUnpacking(b,this->addr);
    doParsimUnpacking(b,this->appAddr);
    doParsimUnpacking(b,this->pos);
    doParsimUnpacking(b,this->batteryLevelAbs);
    doParsimUnpacking(b,this->batteryLevelPerc);
}

long ApplicationPacketRecharge::getSequenceNumber() const
{
    return this->sequenceNumber;
}

void ApplicationPacketRecharge::setSequenceNumber(long sequenceNumber)
{
    this->sequenceNumber = sequenceNumber;
}

L3Address& ApplicationPacketRecharge::getAddr()
{
    return this->addr;
}

void ApplicationPacketRecharge::setAddr(const L3Address& addr)
{
    this->addr = addr;
}

int ApplicationPacketRecharge::getAppAddr() const
{
    return this->appAddr;
}

void ApplicationPacketRecharge::setAppAddr(int appAddr)
{
    this->appAddr = appAddr;
}

Coord& ApplicationPacketRecharge::getPos()
{
    return this->pos;
}

void ApplicationPacketRecharge::setPos(const Coord& pos)
{
    this->pos = pos;
}

double ApplicationPacketRecharge::getBatteryLevelAbs() const
{
    return this->batteryLevelAbs;
}

void ApplicationPacketRecharge::setBatteryLevelAbs(double batteryLevelAbs)
{
    this->batteryLevelAbs = batteryLevelAbs;
}

double ApplicationPacketRecharge::getBatteryLevelPerc() const
{
    return this->batteryLevelPerc;
}

void ApplicationPacketRecharge::setBatteryLevelPerc(double batteryLevelPerc)
{
    this->batteryLevelPerc = batteryLevelPerc;
}

class ApplicationPacketRechargeDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertynames;
  public:
    ApplicationPacketRechargeDescriptor();
    virtual ~ApplicationPacketRechargeDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyname) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyname) const override;
    virtual int getFieldArraySize(void *object, int field) const override;

    virtual std::string getFieldValueAsString(void *object, int field, int i) const override;
    virtual bool setFieldValueAsString(void *object, int field, int i, const char *value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual void *getFieldStructValuePointer(void *object, int field, int i) const override;
};

Register_ClassDescriptor(ApplicationPacketRechargeDescriptor);

ApplicationPacketRechargeDescriptor::ApplicationPacketRechargeDescriptor() : omnetpp::cClassDescriptor("inet::ApplicationPacketRecharge", "omnetpp::cPacket")
{
    propertynames = nullptr;
}

ApplicationPacketRechargeDescriptor::~ApplicationPacketRechargeDescriptor()
{
    delete[] propertynames;
}

bool ApplicationPacketRechargeDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<ApplicationPacketRecharge *>(obj)!=nullptr;
}

const char **ApplicationPacketRechargeDescriptor::getPropertyNames() const
{
    if (!propertynames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
        const char **basenames = basedesc ? basedesc->getPropertyNames() : nullptr;
        propertynames = mergeLists(basenames, names);
    }
    return propertynames;
}

const char *ApplicationPacketRechargeDescriptor::getProperty(const char *propertyname) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : nullptr;
}

int ApplicationPacketRechargeDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 6+basedesc->getFieldCount() : 6;
}

unsigned int ApplicationPacketRechargeDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldTypeFlags(field);
        field -= basedesc->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,
        FD_ISCOMPOUND,
        FD_ISEDITABLE,
        FD_ISCOMPOUND,
        FD_ISEDITABLE,
        FD_ISEDITABLE,
    };
    return (field>=0 && field<6) ? fieldTypeFlags[field] : 0;
}

const char *ApplicationPacketRechargeDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldName(field);
        field -= basedesc->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sequenceNumber",
        "addr",
        "appAddr",
        "pos",
        "batteryLevelAbs",
        "batteryLevelPerc",
    };
    return (field>=0 && field<6) ? fieldNames[field] : nullptr;
}

int ApplicationPacketRechargeDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount() : 0;
    if (fieldName[0]=='s' && strcmp(fieldName, "sequenceNumber")==0) return base+0;
    if (fieldName[0]=='a' && strcmp(fieldName, "addr")==0) return base+1;
    if (fieldName[0]=='a' && strcmp(fieldName, "appAddr")==0) return base+2;
    if (fieldName[0]=='p' && strcmp(fieldName, "pos")==0) return base+3;
    if (fieldName[0]=='b' && strcmp(fieldName, "batteryLevelAbs")==0) return base+4;
    if (fieldName[0]=='b' && strcmp(fieldName, "batteryLevelPerc")==0) return base+5;
    return basedesc ? basedesc->findField(fieldName) : -1;
}

const char *ApplicationPacketRechargeDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldTypeString(field);
        field -= basedesc->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "long",
        "L3Address",
        "int",
        "Coord",
        "double",
        "double",
    };
    return (field>=0 && field<6) ? fieldTypeStrings[field] : nullptr;
}

const char **ApplicationPacketRechargeDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldPropertyNames(field);
        field -= basedesc->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *ApplicationPacketRechargeDescriptor::getFieldProperty(int field, const char *propertyname) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldProperty(field, propertyname);
        field -= basedesc->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int ApplicationPacketRechargeDescriptor::getFieldArraySize(void *object, int field) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldArraySize(object, field);
        field -= basedesc->getFieldCount();
    }
    ApplicationPacketRecharge *pp = (ApplicationPacketRecharge *)object; (void)pp;
    switch (field) {
        default: return 0;
    }
}

std::string ApplicationPacketRechargeDescriptor::getFieldValueAsString(void *object, int field, int i) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldValueAsString(object,field,i);
        field -= basedesc->getFieldCount();
    }
    ApplicationPacketRecharge *pp = (ApplicationPacketRecharge *)object; (void)pp;
    switch (field) {
        case 0: return long2string(pp->getSequenceNumber());
        case 1: {std::stringstream out; out << pp->getAddr(); return out.str();}
        case 2: return long2string(pp->getAppAddr());
        case 3: {std::stringstream out; out << pp->getPos(); return out.str();}
        case 4: return double2string(pp->getBatteryLevelAbs());
        case 5: return double2string(pp->getBatteryLevelPerc());
        default: return "";
    }
}

bool ApplicationPacketRechargeDescriptor::setFieldValueAsString(void *object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->setFieldValueAsString(object,field,i,value);
        field -= basedesc->getFieldCount();
    }
    ApplicationPacketRecharge *pp = (ApplicationPacketRecharge *)object; (void)pp;
    switch (field) {
        case 0: pp->setSequenceNumber(string2long(value)); return true;
        case 2: pp->setAppAddr(string2long(value)); return true;
        case 4: pp->setBatteryLevelAbs(string2double(value)); return true;
        case 5: pp->setBatteryLevelPerc(string2double(value)); return true;
        default: return false;
    }
}

const char *ApplicationPacketRechargeDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldStructName(field);
        field -= basedesc->getFieldCount();
    }
    switch (field) {
        case 1: return omnetpp::opp_typename(typeid(L3Address));
        case 3: return omnetpp::opp_typename(typeid(Coord));
        default: return nullptr;
    };
}

void *ApplicationPacketRechargeDescriptor::getFieldStructValuePointer(void *object, int field, int i) const
{
    omnetpp::cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount())
            return basedesc->getFieldStructValuePointer(object, field, i);
        field -= basedesc->getFieldCount();
    }
    ApplicationPacketRecharge *pp = (ApplicationPacketRecharge *)object; (void)pp;
    switch (field) {
        case 1: return (void *)(&pp->getAddr()); break;
        case 3: return (void *)(&pp->getPos()); break;
        default: return nullptr;
    }
}

} // namespace inet

