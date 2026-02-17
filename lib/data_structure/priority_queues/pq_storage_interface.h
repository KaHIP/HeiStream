#ifndef PQ_STORAGE_INTERFACE_H
#define PQ_STORAGE_INTERFACE_H

#include <functional>
#include <optional>

template <typename Key, typename Value>
class PQStorageBase {
   public:
    virtual void insert(const Key& key, const Value& value) = 0;
    virtual void emplace(const Key& key, Value&& value) = 0;
    virtual std::optional<std::reference_wrapper<Value>> get(const Key& key) = 0;
    virtual std::optional<std::reference_wrapper<const Value>> get(const Key& key) const = 0;
    virtual Value* find_ptr(const Key& key) = 0;
    virtual const Value* find_ptr(const Key& key) const = 0;
    virtual void erase(const Key& key) = 0;
    virtual bool contains(const Key& key) const = 0;
    virtual Value& operator[](const Key& key) = 0;
    virtual void reserve(size_t size) = 0;
    virtual void set_max_load_factor(float factor) = 0;
    virtual ~PQStorageBase() = default;
};

#endif
