/*
  Copyright (c) 2015, 2016 Genome Research Ltd.

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <unistd.h>

#include "gcsa.h"

using namespace gcsa;

//------------------------------------------------------------------------------

void identifyGCSA(const std::string& input_name);
void compressGCSA(const std::string& input_name, const std::string& output_name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: gcsa_format [options] input [output]" << std::endl;
    std::cerr << "  -c    Compress the GCSA in an experimental file format" << std::endl;
    std::cerr << "  -i    Identify the file format version (default)" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  int c = 0;
  bool compress = false, identify = false;
  bool output_needed = false;
  std::string input_name, output_name;
  while((c = getopt(argc, argv, "ci")) != -1)
  {
    switch(c)
    {
    case 'c':
      compress = true; output_needed = true; break;
    case 'i':
      identify = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(!compress) { identify = true; }
  if(optind < argc) { input_name = argv[optind]; optind++; }
  else
  {
    std::cerr << "gcsa_format: Input file not specified" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(output_needed)
  {
    if(optind < argc) { output_name = argv[optind]; optind++; }
    else
    {
      std::cerr << "gcsa_format: Output file not specified" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if(input_name == output_name)
  {
    std::cerr << "gcsa_format: Input and output files must be different" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(identify) { identifyGCSA(input_name); }
  if(compress) { compressGCSA(input_name, output_name); }

  return 0;
}

//------------------------------------------------------------------------------

/*
  This is GCSA file format header version 0. It was used in pre-releases 0.1 to 0.4.
*/

struct GCSAHeader_0
{
  uint64_t path_nodes;
  uint64_t order;

  const static uint32_t VERSION = 0;

  GCSAHeader_0();

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check() const;
};

GCSAHeader_0::GCSAHeader_0() :
  path_nodes(0), order(0)
{
}

size_type
GCSAHeader_0::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(this->path_nodes, out, child, "path_nodes");
  written_bytes += sdsl::write_member(this->order, out, child, "order");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSAHeader_0::load(std::istream& in)
{
  sdsl::read_member(this->path_nodes, in);
  sdsl::read_member(this->order, in);
}

bool
GCSAHeader_0::check() const
{
  return true;
}

std::ostream& operator<<(std::ostream& stream, const GCSAHeader_0& header)
{
  return stream << "GCSA header version " << GCSAHeader_0::VERSION << ": "
                << header.path_nodes << " path nodes, order " << header.order;
}

//------------------------------------------------------------------------------

/*
  This is a generic encoder for the GCSA body. With the default template parameters, it is
  identical to file format version 0 to 1.

  serialize() does not create a new structure tree node. Calling it is equivalent to
  serializing each of the fields separately.
*/

template<class BWTType = GCSA::bwt_type, class BitVector = GCSA::bit_vector>
struct GCSABody
{
  explicit GCSABody(const GCSA& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;

  BWTType                           bwt;
  Alphabet                          alpha;

  // The last BWT position in each path is marked with an 1-bit.
  BitVector                         path_nodes;
  typename BitVector::rank_1_type   path_rank;
  typename BitVector::select_1_type path_select;

  // The last outgoing edge from each path is marked with an 1-bit.
  BitVector                         edges;
  typename BitVector::rank_1_type   edge_rank;
  typename BitVector::select_1_type edge_select;

  // Paths containing samples are marked with an 1-bit.
  BitVector                         sampled_paths;
  typename BitVector::rank_1_type   sampled_path_rank;

  // The last sample belonging to the same path is marked with an 1-bit.
  sdsl::int_vector<0>               stored_samples;
  BitVector                         samples;
  typename BitVector::select_1_type sample_select;

  static sdsl::bit_vector extract(const BitVector& source);
};

template<class BWTType, class BitVector>
GCSABody<BWTType, BitVector>::GCSABody(const GCSA& source)
{
  // Extract plain BWT.
  std::string filename = TempFile::getName(".gcsa_format");
  sdsl::int_vector_buffer<8> bwt_buffer(filename, std::ios::out);
  for(size_type i = 0; i < source.bwt.size(); i++) { bwt_buffer[i] = source.bwt[i]; }
  bwt_buffer.close();

  bwt_buffer = sdsl::int_vector_buffer<8>(filename);
  this->bwt = BWTType(bwt_buffer, source.bwt.size());
  bwt_buffer.close(); remove(filename.c_str());
  this->alpha = source.alpha;

  this->path_nodes = extract(source.path_nodes);
  sdsl::util::init_support(this->path_rank, &(this->path_nodes));
  sdsl::util::init_support(this->path_select, &(this->path_nodes));

  this->edges = extract(source.edges);
  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->edge_select, &(this->edges));

  this->sampled_paths = extract(source.sampled_paths);
  sdsl::util::init_support(this->sampled_path_rank, &(this->sampled_paths));

  this->stored_samples = source.stored_samples;
  this->samples = extract(source.samples);
  sdsl::util::init_support(this->sample_select, &(this->samples));
}

template<class BWTType, class BitVector>
sdsl::bit_vector
GCSABody<BWTType, BitVector>::extract(const BitVector& source)
{
  sdsl::bit_vector buffer(source.size());
  for(size_type i = 0; i < source.size(); i++) { buffer[i] = source[i]; }
  return buffer;
}

template<class BWTType, class BitVector>
size_type
GCSABody<BWTType, BitVector>::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string) const
{
  sdsl::structure_tree_node* child = v; // We do not consider the body a real object.
  size_type written_bytes = 0;

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->alpha.serialize(out, child, "alpha");

  written_bytes += this->path_nodes.serialize(out, child, "path_nodes");
  written_bytes += this->path_rank.serialize(out, child, "path_rank");
  written_bytes += this->path_select.serialize(out, child, "path_select");

  written_bytes += this->edges.serialize(out, child, "edges");
  written_bytes += this->edge_rank.serialize(out, child, "edge_rank");
  written_bytes += this->edge_select.serialize(out, child, "edge_select");

  written_bytes += this->sampled_paths.serialize(out, child, "sampled_paths");
  written_bytes += this->sampled_path_rank.serialize(out, child, "sampled_path_rank");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  return written_bytes;
}

//------------------------------------------------------------------------------

struct SadaSparseSparse
{
  typedef sdsl::sd_vector<> sd_vector;

  explicit SadaSparseSparse(const SadaCount& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;

  // Values equal to 1 are marked with an 1-bit.
  sd_vector                ones;
  sd_vector::rank_1_type   one_rank;

  // Values greater than 1 are marked with an 1-bit.
  sd_vector                filter;
  sd_vector::rank_1_type   filter_rank;

  // Non-zero values encoded in unary: k becomes 0^{k-1} 1.
  sd_vector                values;
  sd_vector::select_1_type value_select;
};

SadaSparseSparse::SadaSparseSparse(const SadaCount& source)
{
  size_type values = 0, ones = 0, large = 0, large_total = 0;
  for(size_type i = 0, curr_val = 0; i < source.data.size(); i++)
  {
    if(source.data[i] == 1)
    {
      values++;
      if(curr_val == 1) { ones++; }
      else if(curr_val > 1) { large++; large_total += curr_val; }
      curr_val = 0;
    }
    else { curr_val++; }
  }

  sdsl::sd_vector_builder one_builder(values, ones);
  sdsl::sd_vector_builder filter_builder(values - ones, large);
  sdsl::sd_vector_builder value_builder(large_total, large);
  for(size_type i = 0, curr_pos = 0, curr_val = 0, curr_ones = 0, tail = 0; i < source.data.size(); i++)
  {
    if(source.data[i] == 1)
    {
      if(curr_val == 1) { one_builder.set(curr_pos); curr_ones++; }
      else if(curr_val > 1)
      {
        filter_builder.set(curr_pos - curr_ones);
        tail += curr_val; value_builder.set(tail - 1);
      }
      curr_pos++; curr_val = 0;
    }
    else { curr_val++; }
  }
  this->ones = sd_vector(one_builder);
  this->filter = sd_vector(filter_builder);
  this->values = sd_vector(value_builder);

  sdsl::util::init_support(this->one_rank, &(this->ones));
  sdsl::util::init_support(this->filter_rank, &(this->filter));
  sdsl::util::init_support(this->value_select, &(this->values));
}

size_type
SadaSparseSparse::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->ones.serialize(out, child, "ones");
  written_bytes += this->one_rank.serialize(out, child, "one_rank");

  written_bytes += this->filter.serialize(out, child, "filter");
  written_bytes += this->filter_rank.serialize(out, child, "filter_rank");

  written_bytes += this->values.serialize(out, child, "values");
  written_bytes += this->value_select.serialize(out, child, "value_select");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

//------------------------------------------------------------------------------

/*
  An experiment for compressing the GCSA. Bitvectors are entropy-compressed, while
  samples are stored as their ranks in the set of sampled nodes.
*/

struct ExperimentalGCSA
{
  typedef GCSA::size_type           size_type;
  typedef sdsl::rrr_vector<>        bit_vector;
  typedef sdsl::sd_vector<>         sd_vector;
  typedef sdsl::wt_huff<bit_vector> bwt_type;

  explicit ExperimentalGCSA(const GCSA& source);
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;

  GCSAHeader                        header;
  GCSABody<bwt_type, bit_vector>    body;

  // Structures used for compressing the samples.
  sd_vector                      sample_values;
  sd_vector::select_1_type       sample_value_select;

  // Structures used for counting queries.
  SadaSparse                        extra_pointers;
  SadaSparseSparse                  redundant_pointers;
};

ExperimentalGCSA::ExperimentalGCSA(const GCSA& source) :
  body(source),
  extra_pointers(source.extra_pointers), redundant_pointers(source.redundant_pointers)
{
  this->header.path_nodes = source.size();
  this->header.edges = source.edgeCount();
  this->header.order = source.order();
  this->header.flags = GCSAHeader::COMPRESSED;

  // Find the distinct sampled nodes.
  std::vector<node_type> values(source.stored_samples.begin(), source.stored_samples.end());
  removeDuplicates(values, true);
  size_type distinct = values.size();
  std::cout << distinct << " distinct values in " << source.stored_samples.size() << " samples" << std::endl;

  // Build a mapping between sampled nodes and their ranks.
  this->sample_values = sd_vector(values.begin(), values.end());
  sd_vector::rank_1_type value_rank(&(this->sample_values));
  sdsl::util::init_support(this->sample_value_select, &(this->sample_values));

  // Convert sampled nodes to their ranks.
  values.resize(source.stored_samples.size());
  #pragma omp parallel for schedule(static)
  for(size_type i = 0; i < values.size(); i++)
  {
    values[i] = value_rank(source.stored_samples[i]);
  }

  // Replace the samples copied from source with compressed samples.
  this->body.stored_samples = sdsl::int_vector<0>(source.stored_samples.size(), 0, bit_length(distinct - 1));
  for(size_type i = 0; i < this->body.stored_samples.size(); i++) { this->body.stored_samples[i] = values[i]; }
}

ExperimentalGCSA::size_type
ExperimentalGCSA::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->body.serialize(out, child, "body");

  written_bytes += this->sample_values.serialize(out, child, "sample_values");
  written_bytes += this->sample_value_select.serialize(out, child, "sample_value_select");

  written_bytes += this->extra_pointers.serialize(out, child, "extra_pointers");
  written_bytes += this->redundant_pointers.serialize(out, child, "redundant_pointers");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

//------------------------------------------------------------------------------

std::string
stripSuffix(const std::string& value, const std::string& suffix)
{
  if(value.length() < suffix.length()) { return value; }

  bool found = true;
  for(size_type i = 0, j = value.length() - suffix.length(); i < suffix.length(); i++, j++)
  {
    if(value[j] != suffix[i]) { found = false; break; }
  }

  if(found) { return value.substr(0, value.length() - suffix.length()); }
  return value;
}

template<class IndexType>
void
componentSizes(const IndexType& index, const std::string& file_name)
{
  std::string base_name = stripSuffix(file_name, GCSA::EXTENSION);

  std::string html_name = base_name + ".html";
  std::ofstream html_output(html_name.c_str());
  if(!html_output)
  {
    std::cerr << "componentSizes(): Cannot open HTML file " << html_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  sdsl::write_structure<sdsl::HTML_FORMAT>(index, html_output);
  html_output.close();
  std::cout << "Size breakdown written to " << html_name << std::endl;

  std::string json_name = base_name + ".json";
  std::ofstream json_output(json_name.c_str());
  if(!json_output)
  {
    std::cerr << "componentSizes(): Cannot open JSON file " << json_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  sdsl::write_structure<sdsl::JSON_FORMAT>(index, json_output);
  json_output.close();
  std::cout << "Size breakdown written to " << json_name << std::endl;
}

void
compressGCSA(const std::string& input_name, const std::string& output_name)
{
  std::cout << "Experimental GCSA compressor" << std::endl;
  std::cout << std::endl;

  std::cout << "Input: " << input_name << std::endl;
  std::cout << "Output: " << output_name << std::endl;
  std::cout << std::endl;

  GCSA index;
  sdsl::load_from_file(index, input_name);
  std::cout << "GCSA size: " << inMegabytes(sdsl::size_in_bytes(index)) << " MB" << std::endl;
  componentSizes(index, input_name);
  std::cout << std::endl;

  ExperimentalGCSA alt(index);
  std::cout << "Alternative encoding: " << inMegabytes(sdsl::size_in_bytes(alt)) << " MB" << std::endl;
  sdsl::store_to_file(alt, output_name);
  componentSizes(alt, output_name);
  std::cout << std::endl;

  std::cout << std::endl;
}

//------------------------------------------------------------------------------

bool
tryVersion0(std::ifstream& input)
{
  std::streampos pos = input.tellg();
  GCSAHeader_0 header;
  header.load(input);
  input.seekg(pos);

  if(header.check())
  {
    std::cout << header << std::endl;
    return true;
  }
  return false;
}

bool
tryCurrentVersion(std::ifstream& input)
{
  std::streampos pos = input.tellg();
  GCSAHeader header;
  header.load(input);
  input.seekg(pos);

  for(uint32_t version = GCSAHeader::VERSION; version >= GCSAHeader::MIN_VERSION; version--)
  {
    if(header.check(version))
    {
      std::cout << header << std::endl;
      return true;
    }
  }
  if(header.checkNew())
  {
    std::cout << header << std::endl;
    std::cout << "Note: The file is newer than this version of GCSA" << std::endl;
    return true;
  }

  return false;
}

void
identifyGCSA(const std::string& input_name)
{
  std::cout << "GCSA file format inspector" << std::endl;
  std::cout << std::endl;

  std::cout << "Input: " << input_name << std::endl;
  std::cout << std::endl;

  std::ifstream input(input_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "identifyGCSA(): Cannot open input file " << input_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(!tryCurrentVersion(input))
  {
    std::cout << "File format cannot be identified, trying version 0" << std::endl;
    tryVersion0(input);
  }
  input.close();

  std::cout << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
