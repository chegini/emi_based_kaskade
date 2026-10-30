#ifndef EMI_MATERIAL_TEST_MATRIX_HH
#define EMI_MATERIAL_TEST_MATRIX_HH

#include <array>
#include <map>
#include <stdexcept>
#include <vector>

// Small sparse matrix with the NumaBCRS scalar-block/proxy interface. Its
// iterators deliberately offer only prefix increment, as NumaBCRS does.
// Production codec/audit tests therefore need no Dune or Kaskade libraries.
class MaterialTestMatrix
{
public:
  using Block = std::array<std::array<double,1>,1>;
  using Entries = std::map<size_t,Block>;
  struct Column
  {
    Entries::const_iterator position;
    size_t index() const { return position->first; }
    Block const& operator*() const { return position->second; }
    void operator++() { ++position; }
    bool operator!=(Column const& other) const { return position != other.position; }
    bool operator==(Column const& other) const { return position == other.position; }
  };
  struct ConstRow
  {
    Entries const& entries;
    Column begin() const { return {entries.begin()}; }
    Column end() const { return {entries.end()}; }
  };
  struct Row
  {
    Entries& entries;
    Block& operator[](size_t column)
    {
      auto const found = entries.find(column);
      if (found == entries.end())
        throw std::runtime_error("Material test matrix write outside pattern");
      return found->second;
    }
  };
  struct Rows
  {
    MaterialTestMatrix const& matrix;
    size_t row;
    size_t index() const { return row; }
    Column begin() const { return {matrix.rows.at(row).begin()}; }
    Column end() const { return {matrix.rows.at(row).end()}; }
    void operator++() { ++row; }
    bool operator!=(Rows const& other) const { return row != other.row; }
  };

  MaterialTestMatrix() = default;
  explicit MaterialTestMatrix(size_t n) : rows(n)
  {
    for (auto& row : rows)
      for (size_t col = 0; col < n; ++col)
        row[col] = {{{0.0}}};
  }
  MaterialTestMatrix& operator=(double value)
  {
    for (auto& row : rows)
      for (auto& [col,block] : row)
      {
        (void)col;
        block[0][0] = value;
      }
    return *this;
  }
  size_t N() const { return rows.size(); }
  size_t M() const { return rows.size(); }
  Rows begin() const { return {*this,0}; }
  Rows end() const { return {*this,rows.size()}; }
  ConstRow operator[](size_t row) const { return {rows.at(row)}; }
  Row operator[](size_t row) { return {rows.at(row)}; }
  void erase(size_t row, size_t col) { rows.at(row).erase(col); }

private:
  std::vector<Entries> rows;
};

#endif
