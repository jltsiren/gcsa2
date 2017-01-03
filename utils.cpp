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

#include <cstdio>
#include <cstdlib>

#include <sys/resource.h>
#include <unistd.h>

#include <gcsa/internal.h>

namespace gcsa
{

//------------------------------------------------------------------------------

size_type Verbosity::level = Verbosity::DEFAULT;

void
Verbosity::set(size_type new_level)
{
  level = Range::bound(new_level, SILENT, FULL);
}

std::string
Verbosity::levelName()
{
  switch(level)
  {
    case SILENT:
      return "silent"; break;
    case BASIC:
      return "basic"; break;
    case EXTENDED:
      return "extended"; break;
    case FULL:
      return "full"; break;
  }
  return "unknown";
}

//------------------------------------------------------------------------------

void
printHeader(const std::string& header, size_type indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }
  std::cout << header << ":" << padding;
}

void
printTime(const std::string& header, size_type queries, double seconds, size_type indent)
{
  printHeader(header, indent);
  std::cout << queries << " queries in " << seconds << " seconds ("
            << inMicroseconds(seconds / queries) << " µs/query)" << std::endl;
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return omp_get_wtime();
}

size_type
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#ifdef RUSAGE_IN_BYTES
  return usage.ru_maxrss;
#else
  return KILOBYTE * usage.ru_maxrss;
#endif
}

size_type
readVolume()
{
  return DiskIO::read_volume;
}

size_type
writeVolume()
{
  return DiskIO::write_volume;
}

//------------------------------------------------------------------------------

std::atomic<size_type> TempFile::counter(0);

const std::string TempFile::DEFAULT_TEMP_DIR = ".";
std::string TempFile::temp_dir = TempFile::DEFAULT_TEMP_DIR;

void
TempFile::setDirectory(const std::string& directory)
{
  if(directory.empty()) { temp_dir = DEFAULT_TEMP_DIR; }
  else if(directory[directory.length() - 1] != '/') { temp_dir = directory; }
  else { temp_dir = directory.substr(0, directory.length() - 1); }
}

std::string
TempFile::getName(const std::string& name_part)
{
  char hostname[32];
  gethostname(hostname, 32); hostname[31] = 0;

  return temp_dir + '/' + name_part + '_'
    + std::string(hostname) + '_'
    + sdsl::util::to_string(sdsl::util::pid()) + '_'
    + sdsl::util::to_string(counter++);
}

void
TempFile::remove(std::string& filename)
{
  if(!(filename.empty()))
  {
    std::remove(filename.c_str());
    filename.clear();
  }
}

size_type
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  size_type chars = 0;
  while(input)
  {
    std::string buf;
    std::getline(input, buf);
    if(skip_empty_rows && buf.empty()) { continue; }
    rows.push_back(buf);
    chars += buf.length();
  }

  input.close();
  return chars;
}

size_type
fileSize(std::ifstream& file)
{
  std::streamoff curr = file.tellg();

  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  size -= file.tellg();

  file.seekg(curr, std::ios::beg);
  return size;
}

size_type
fileSize(std::ofstream& file)
{
  std::streamoff curr = file.tellp();

  file.seekp(0, std::ios::end);
  std::streamoff size = file.tellp();
  file.seekp(0, std::ios::beg);
  size -= file.tellp();

  file.seekp(curr, std::ios::beg);
  return size;
}

//------------------------------------------------------------------------------

} // namespace gcsa
