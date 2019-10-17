#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
#include <cctype>
using namespace std;

class GenomeImpl {
public:
    GenomeImpl(const string& nm, const string& sequence);

    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;

private:
	string m_name;
	string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
	: m_name(nm), m_sequence(sequence) {}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) {
	/*char c;
	string name, sequence;
	while (genomeSource.get(c)) {
		if (c == '\n') {
			continue;
		}

		// beginning of a name
		if (c == '>') {
			// if there is another name but no sequence belonging to the previous name
			if (name != "" && sequence == "") {
				return false;
			}

			// if sequence isn't empty, then we know this is a new genome in the text file,
			// so create a new genome and push onto vector
			if (sequence != "") {
				Genome creation(name, sequence);
				genomes.push_back(creation);
				sequence = "";
			}

			// read name line
			name = "";
			while (genomeSource.get(c) && c != '\n') {
				name += c;
			}

			// if nothing follows '>' (so name is nothing)
			if (name == "") {
				return false;
			}
		}

		// everything else belongs in the sequence
		else {
			// can't read in a sequence without a name
			if (name == "") {
				return false;
			}

			// char must be valid
			if (c != 'a' && c != 'A' &&
				c != 'c' && c != 'C' &&
				c != 't' && c != 'T' &&
				c != 'g' && c != 'G' &&
				c != 'n' && c != 'N') {
				return false;
			}

			// every base letter must be upper case, even if lower case in the file
			sequence += toupper(c);
		}
	}

	// end of loop, create last genome

	// if there is just a name without a following sequence
	if (sequence == "" || name == "") {
		return false;
	}

	Genome creation(name, sequence);
	genomes.push_back(creation);
	return true;*/

	string name, sequence, temp;
	while (getline(genomeSource, temp)) {
		// beginning of a name
		if (temp[0] == '>') {
			// if nothing follows '>' (name is nothing)
			if (temp == ">") {
				return false;
			}

			// if there is another name but no sequence belonging to the previous name
			if (name != "" && sequence == "") {
				return false;
			}

			// if sequence isn't empty, then we know this is a new genome in the text file,
			// so create a new genome and push onto vector
			if (sequence != "") {
				Genome creation(name, sequence);
				genomes.push_back(creation);
				sequence = "";
			}

			// set name as everything after '>'
			name = temp.substr(1, temp.size() - 1);
		}

		// everything else belongs in the sequence
		else {
			// can't read in a sequence without a name
			if (name == "") {
				return false;
			}

			for (int i = 0; i < temp.size(); i++) {
				temp[i] = toupper(temp[i]);
				// char must be valid
				if (temp[i] != 'A' &&
					temp[i] != 'C' &&
					temp[i] != 'T' &&
					temp[i] != 'G' &&
					temp[i] != 'N') {
					return false;
				}
			}

			// append to sequence
			sequence += temp;
		}
	}

	// end of loop, create last genome

	// if there is just a name without a following sequence
	if (sequence == "" || name == "") {
		return false;
	}

	Genome creation(name, sequence);
	genomes.push_back(creation);
	return true;
}

int GenomeImpl::length() const {
	return m_sequence.size();
}

string GenomeImpl::name() const {
	return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const {
	// can't read nonexistent positions
	if (position < 0 || position >= m_sequence.size() ||
		length <= 0 || position + length > m_sequence.size()) {
		return false;
	}

	fragment = m_sequence.substr(position, length);
	return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
