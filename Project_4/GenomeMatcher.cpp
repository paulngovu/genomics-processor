#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
using namespace std;

class GenomeMatcherImpl {
public:
    GenomeMatcherImpl(int minSearchLength);

    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;

private:
	int m_minSearchLength;
	std::vector<Genome> m_genomes;		// stores genomes
	Trie<DNAMatch> m_genomeTrie;
	unordered_map<string, int> gMap;	// maps genome name and index

	// member functions
	bool findDNAHelper(const string& cmp1, const string& cmp2, bool exactMatchOnly) const;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength) {
	m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome) {
	// add genome to collection of genomes
	m_genomes.push_back(genome);

	// index DNA sequences of minSearchLength into trie
	for (int i = 0; i + m_minSearchLength <= genome.length(); i++) {
		string temp;
		genome.extract(i, m_minSearchLength, temp);

		DNAMatch dna;
		dna.genomeName = genome.name();
		dna.position = i;
		dna.length = m_minSearchLength;
		m_genomeTrie.insert(temp, dna);
	}

	// maps genome name and index
	gMap.insert(pair<string, int>(genome.name(), m_genomes.size() - 1));
}

int GenomeMatcherImpl::minimumSearchLength() const {
	return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const {
	if (fragment.size() < minimumLength || minimumLength < m_minSearchLength) {
		return false;
	}

	matches.clear();

	std::vector<DNAMatch> candidates = m_genomeTrie.find(fragment.substr(0, m_minSearchLength), exactMatchOnly);

	// remove candidates that are invalid (if index is too far back)
	std::vector<DNAMatch>::iterator ci;
	for (ci = candidates.begin(); ci != candidates.end();) {
		if ((*ci).position > m_genomes[gMap.at((*ci).genomeName)].length() - fragment.size()) {
			ci = candidates.erase(ci);
		}
		else {
			ci++;
		}
	}
	
	// verify that we can match minimumLength or more characters
	if (!candidates.empty()) {
		string tempFrag;
		int tempMax = minimumLength;

		for (ci = candidates.begin(); ci != candidates.end(); ci++) {
			for (int i = 0; i + minimumLength <= fragment.size(); i++) {
				// use map to access proper genome in vector
				m_genomes[gMap.at((*ci).genomeName)].extract((*ci).position, i + minimumLength, tempFrag);

				// compare extraction with fragment
				if (findDNAHelper(tempFrag, fragment.substr(0, i + minimumLength), exactMatchOnly)) {
					(*ci).length = tempFrag.size();

					if (tempMax < tempFrag.size()) {
						tempMax = tempFrag.size();
					}
				}
			}

			if ((*ci).length == tempMax) {
				matches.push_back(*ci);
			}
		}
	}

	// there must be at least one match to return true
	return (matches.size() > 0);
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const {
	if (fragmentMatchLength < m_minSearchLength) {
		return false;
	}

	double S = query.length() / fragmentMatchLength;	// used to calculate percentage later
	std::vector<DNAMatch> tempMatches;
	std::vector<DNAMatch> matches;

	int i;
	for (i = 0; i * fragmentMatchLength < query.length(); i++) {
		// extract fragment from query genome
		string frag;
		query.extract(i * fragmentMatchLength, fragmentMatchLength, frag);

		// search across all genomes in library
		findGenomesWithThisDNA(frag, fragmentMatchLength, exactMatchOnly, tempMatches);
		matches.insert(matches.end(), tempMatches.begin(), tempMatches.end());
	}

	// search for number of sequences that match a certain genome
	int* nMatchingSequences = new int[m_genomes.size()];
	for (i = 0; i < m_genomes.size(); i++) {
		nMatchingSequences[i] = 0;			
	}										
	// array corresponds to index of genome vector. Every time we reach a match, we increment
	// at the proper index.
	for (i = 0; i < matches.size(); i++) {
		nMatchingSequences[gMap.at(matches[i].genomeName)]++;		// use the map to find the proper index to increment
	}

	// calculate percentage
	for (i = 0; i < m_genomes.size(); i++) {
		double percentage = (nMatchingSequences[i] * 1.0) / S;

		if (percentage * 100.0 >= matchPercentThreshold) {	// matching genome
			GenomeMatch creation;
			creation.genomeName = m_genomes[i].name();
			creation.percentMatch = percentage * 100.0;
			results.push_back(creation);
		}
	}

	delete[] nMatchingSequences;
	// there must be at least one match to return true
	return (results.size() > 0);
}

bool GenomeMatcherImpl::findDNAHelper(const string & cmp1, const string & cmp2, bool exactMatchOnly) const {
	// cmp1 and cmp2 should be the same size

	if (exactMatchOnly) {
		for (int i = 0; i < cmp1.size(); i++) {
			if (cmp1[i] != cmp2[i]) {
				return false;
			}
		}
	}

	// SNiP
	else {
		bool alreadyMismatch = false;
		for (int i = 0; i < cmp1.size(); i++) {
			if (cmp1[i] != cmp2[i]) {
				// can only have one mismatch
				if (!alreadyMismatch) {
					alreadyMismatch = true;
					continue;
				}

				return false;
			}
		}
	}

	return true;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
