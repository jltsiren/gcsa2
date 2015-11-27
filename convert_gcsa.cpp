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

#include "gcsa.h"

using namespace gcsa;

/*
  This will eventually become an utility for converting GCSA to the current format.
  Right now it is used for experimenting with different encodings.
*/

//------------------------------------------------------------------------------

struct ExperimentalGCSA
{
  typedef GCSA::size_type    size_type;
//  typedef GCSA::bit_vector   bit_vector;
  typedef sdsl::rrr_vector<> bit_vector;
  typedef sdsl::sd_vector<>  sd_vector;
  typedef GCSA::bwt_type     bwt_type;

  explicit ExperimentalGCSA(const GCSA& source);
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;

  size_type                 path_node_count;
  size_type                 max_query_length;

  bwt_type                  bwt;
  Alphabet                  alpha;

  // The last BWT position in each path is marked with an 1-bit.
  bit_vector                path_nodes;
  bit_vector::rank_1_type   path_rank;
  bit_vector::select_1_type path_select;

  // The last outgoing edge from each path is marked with an 1-bit.
  bit_vector                edges;
  bit_vector::rank_1_type   edge_rank;
  bit_vector::select_1_type edge_select;

  // Paths containing samples are marked with an 1-bit.
  bit_vector                sampled_paths;
  bit_vector::rank_1_type   sampled_path_rank;

  sd_vector                 sample_values;
  sd_vector::select_1_type  sample_value_select;

  // The last sample belonging to the same path is marked with an 1-bit.
  sdsl::int_vector<0>       stored_samples;
  bit_vector                samples;
  bit_vector::select_1_type sample_select;
};

//------------------------------------------------------------------------------

template<class IndexType>
void componentSizes(const IndexType& index, const std::string& base_name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: convert_gcsa base_name" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "GCSA format converter" << std::endl;
  std::cout << std::endl;

  std::string base_name = argv[1];
  std::string input_name = base_name + GCSA::EXTENSION;
  std::cout << "Input: " << input_name << std::endl;
  std::cout << std::endl;

  GCSA index;
  sdsl::load_from_file(index, input_name);
  std::cout << "GCSA size: " << inMegabytes(sdsl::size_in_bytes(index)) << " MB" << std::endl;
  componentSizes(index, base_name);
  std::cout << std::endl;

  std::string alt_name = base_name + "_alt";
  ExperimentalGCSA alt(index);
  std::cout << "Alternative encoding: " << inMegabytes(sdsl::size_in_bytes(alt)) << " MB" << std::endl;
  componentSizes(alt, alt_name);
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

ExperimentalGCSA::ExperimentalGCSA(const GCSA& source)
{
  this->path_node_count = source.path_node_count;
  this->max_query_length = source.max_query_length;

  this->bwt = source.bwt;
  this->alpha = source.alpha;

  this->path_nodes = source.path_nodes;
  sdsl::util::init_support(this->path_rank, &(this->path_nodes));
  sdsl::util::init_support(this->path_select, &(this->path_nodes));

  this->edges = source.edges;
  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->edge_select, &(this->edges));

  this->sampled_paths = source.sampled_paths;
  sdsl::util::init_support(this->sampled_path_rank, &(this->sampled_paths));

  // Find the distinct nodes.
  std::vector<node_type> values(source.stored_samples.begin(), source.stored_samples.end());
  removeDuplicates(values, true);
  size_type distinct = values.size();
  std::cout << distinct << " distinct values in " << source.stored_samples.size() << " samples" << std::endl;

  this->sample_values = sd_vector(values.begin(), values.end());
  sdsl::util::init_support(this->sample_value_select, &(this->sample_values));

  // Convert sampled nodes to their ranks.
  sd_vector::rank_1_type value_rank(&(this->sample_values));
  values.resize(source.stored_samples.size());
  #pragma omp parallel for schedule(static)
  for(size_type i = 0; i < values.size(); i++)
  {
    values[i] = value_rank(source.stored_samples[i]);
  }

  this->stored_samples = sdsl::int_vector<0>(source.stored_samples.size(), 0, bit_length(distinct - 1));
  for(size_type i = 0; i < this->stored_samples.size(); i++) { this->stored_samples[i] = values[i]; }
  this->samples = source.samples;
  sdsl::util::init_support(this->sample_select, &(this->samples));
}

ExperimentalGCSA::size_type
ExperimentalGCSA::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->path_node_count, out, child, "path_node_count");
  written_bytes += sdsl::write_member(this->max_query_length, out, child, "max_query_length");

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

  written_bytes += this->sample_values.serialize(out, child, "sample_values");
  written_bytes += this->sample_value_select.serialize(out, child, "sample_value_select");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

//------------------------------------------------------------------------------

template<class IndexType>
void
componentSizes(const IndexType& index, const std::string& base_name)
{
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

//------------------------------------------------------------------------------
