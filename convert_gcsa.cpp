/*
  Copyright (c) 2015 Genome Research Ltd.

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

/*
  This will eventually become an utility for converting GCSA to the current format.
  Right now it is used for experimenting with different encodings.
*/

//------------------------------------------------------------------------------

void identifyGCSA(const std::string& input_name);
void compressGCSA(const std::string& input_name, const std::string& output_name);
void convertGCSA(const std::string& input_name, const std::string& output_name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: convert_gcsa [options] input [output]" << std::endl;
    std::cerr << "  -c    Convert the GCSA to the current file format (default)" << std::endl;
    std::cerr << "  -C    Compress the GCSA in an experimental file format" << std::endl;
    std::cerr << "  -i    Identify the file format version without converting" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  int c = 0;
  bool convert = false, compress = false, identify = false;
  bool output_needed = false;
  std::string input_name, output_name;
  while((c = getopt(argc, argv, "cCi")) != -1)
  {
    switch(c)
    {
    case 'c':
      convert = true; compress = false; output_needed = true; break;
    case 'C':
      compress = true; convert = false; output_needed = true; break;
    case 'i':
      identify = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(!(compress | identify)) { convert = output_needed = true; }
  if(optind < argc) { input_name = argv[optind]; optind++; }
  else
  {
    std::cerr << "convert_gcsa: Input file not specified" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(output_needed)
  {
    if(optind < argc) { output_name = argv[optind]; optind++; }
    else
    {
      std::cerr << "convert_gcsa: Output file not specified" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if(identify) { identifyGCSA(input_name); }
  if(convert)  { convertGCSA(input_name, output_name); }
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
  identical to file format versions 0 and 1.

  serialize() does not create a new structure tree node. Calling it is equivalent to
  serializing each of the fields separately.
*/

template<class BWTType = GCSA::bwt_type, class BitVector = GCSA::bit_vector>
struct GCSABody
{
  GCSABody();
  explicit GCSABody(const GCSA& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

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
GCSABody<BWTType, BitVector>::GCSABody()
{
}

template<class BWTType, class BitVector>
GCSABody<BWTType, BitVector>::GCSABody(const GCSA& source)
{
  // Extract plain BWT.
  std::string filename = TempFile::getName(".convert_gcsa");
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
//  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
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

//  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

template<class BWTType, class BitVector>
void
GCSABody<BWTType, BitVector>::load(std::istream& in)
{
  this->bwt.load(in);
  this->alpha.load(in);

  this->path_nodes.load(in);
  this->path_rank.load(in, &(this->path_nodes));
  this->path_select.load(in, &(this->path_nodes));

  this->edges.load(in);
  this->edge_rank.load(in, &(this->edges));
  this->edge_select.load(in, &(this->edges));

  this->sampled_paths.load(in);
  this->sampled_path_rank.load(in, &(this->sampled_paths));

  this->stored_samples.load(in);
  this->samples.load(in);
  this->sample_select.load(in, &(this->samples));
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

  GCSAHeader                     header;
  GCSABody<bwt_type, bit_vector> body;

  sd_vector                      sample_values;
  sd_vector::select_1_type       sample_value_select;
};

ExperimentalGCSA::ExperimentalGCSA(const GCSA& source) :
  body(source)
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

std::ostream& operator<<(std::ostream& stream, const GCSAHeader_0& header);
template<class Header>
bool tryVersion(std::ifstream& input)
{
  std::streampos pos = input.tellg();
  Header header;
  header.load(input);
  input.seekg(pos);

  if(header.check())
  {
    std::cout << header << std::endl;
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

  if(!tryVersion<GCSAHeader>(input))
  {
    std::cout << "File format cannot be identified, trying version 0" << std::endl;
    tryVersion<GCSAHeader_0>(input);
  }
  input.close();

  std::cout << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

/*
  A simple GCSA encoding consisting of a header and a body. File format versions 0 and 1
  are of that form.
*/

template<class Header, class Body>
struct GenericGCSA
{
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  Header header;
  Body   body;
};

template<class Header, class Body>
size_type
GenericGCSA<Header, Body>::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->body.serialize(out, child, "body");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

template<class Header, class Body>
void
GenericGCSA<Header, Body>::load(std::istream& in)
{
  this->header.load(in);
  this->body.load(in);
}

//------------------------------------------------------------------------------

void
convertGCSA(const std::string& input_name, const std::string& output_name)
{
  std::cout << "GCSA converter" << std::endl;
  std::cout << std::endl;

  std::cout << "Input: " << input_name << std::endl;
  std::cout << "Output: " << output_name << std::endl;
  std::cout << std::endl;

  std::ifstream input(input_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "convertGCSA(): Cannot open input file " << input_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::ofstream output(output_name.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "convertGCSA(): Cannot open output file " << output_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(tryVersion<GCSAHeader>(input))
  {
    GenericGCSA<GCSAHeader, GCSABody<>> source;
    source.load(input);
    source.serialize(output);
  }
  else
  {
    std::cout << "File format cannot be identified, trying version 0" << std::endl;
    tryVersion<GCSAHeader_0>(input);

    // Load GCSA version 0.
    GenericGCSA<GCSAHeader_0, GCSABody<sdsl::wt_huff<>, sdsl::bit_vector>> source;
    source.load(input);

    // Build the current version of GCSA. At the moment the body is identical.
    GCSAHeader header;
    header.path_nodes = source.header.path_nodes;
    header.edges = source.body.bwt.size();
    header.order = source.header.order;
    header.serialize(output);
    source.body.serialize(output);
  }
  input.close(); output.close();

  std::cout << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
