#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

#include <iostream>

template<typename ValueType>
class Trie {
public:
	Trie();
	~Trie();
	
	void reset();
	void insert(const std::string& key, const ValueType& value);
	std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

	// C++11 syntax for preventing copying and assignment
	Trie(const Trie&) = delete;
	Trie& operator=(const Trie&) = delete;
	
private:
	struct Node {
		char m_label;						// each node has a label
		std::vector<ValueType> m_values;	// each node can have multiple values
		std::vector<Node*> m_children;		// each node has children nodes
	};
	
	Node* m_root;
	
	// member functions
	void cleanUp(Node* head);
	Node* insertHelper(const std::string& key, const int& index, Node* head);
	Node* findExactHelper(const std::string& key, const int& index, Node* head) const;
	void findSNiPHelper(const std::string& key, const int& index, Node* head, std::vector<ValueType>& v, bool alreadyMismatch) const;
	void printHelper(Node* root);			// testing
};

template<typename ValueType>
Trie<ValueType>::Trie() {	
	m_root = new Node;		// create a root trie node with no children
	// label doesn't matter because we're not going to analyze the root node
}

template<typename ValueType>
Trie<ValueType>::~Trie() {
	cleanUp(m_root);
}

template<typename ValueType>
void Trie<ValueType>::reset() {
	cleanUp(m_root);		// free all allocated memory

	m_root = new Node;		// create a root trie node with no children
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value) {
	Node* n = insertHelper(key, 0, m_root);		// searches for node that satisfies key
	n->m_values.push_back(value);
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const {
	// if we want the vector with exact key
	if (exactMatchOnly) {
		Node* n = findExactHelper(key, 0, m_root);	// searches for node of entire key

		// if no values associated with specific key, return empty vector
		if (n == nullptr) {
			return std::vector<ValueType>();
		}

		// values found, so return vector of values
		return n->m_values;
	}

	// if we want a SNiP
	else {
		std::vector<ValueType> vv;						// empty vector
		findSNiPHelper(key, 0, m_root, vv, false);		// adds values to vector
		return vv;
	}
}

template<typename ValueType>
void Trie<ValueType>::cleanUp(Node* head) {
	// base case: head == nullptr
	if (head != nullptr) {
		// iterator
		for (auto i = head->m_children.begin(); i != head->m_children.end(); i++) {
			cleanUp(*i);
		}

		delete head;
	}
}

template<typename ValueType>
typename Trie<ValueType>::Node* Trie<ValueType>::insertHelper(const std::string& key, const int& index, Node* head) {
	// base case / if recursion reaches end of string
	if (head == nullptr || index == key.size()) {
		return head;
	}

	for (auto i = head->m_children.begin(); i != head->m_children.end(); i++) {
		// node label matches char in key string
		if (*i != nullptr && (*i)->m_label == key[index]) {
			return insertHelper(key, index + 1, *i);
		}
	}

	// if we exited the loop, no matching label was found, so we create a new child
	Node* creation = new Node;
	creation->m_label = key[index];
	head->m_children.push_back(creation);
	return insertHelper(key, index + 1, creation);
}

template<typename ValueType>
typename Trie<ValueType>::Node* Trie<ValueType>::findExactHelper(const std::string& key, const int& index, Node* head) const {
	// base case / if recursion reaches end of string
	if (head == nullptr || index == key.size()) {
		return head;
	}

	for (auto i = head->m_children.begin(); i != head->m_children.end(); i++) {
		// // node label matches char in key string
		if (*i != nullptr && (*i)->m_label == key[index]) {
			return findExactHelper(key, index + 1, *i);
		}
	}

	// if we exited the loop, no matching label was found, so return nullptr
	return nullptr;
}

template<typename ValueType>
void Trie<ValueType>::findSNiPHelper(const std::string& key, const int& index, Node* head, std::vector<ValueType>& v, bool alreadyMismatch) const {
	// base case / if recursion reaches end of string
	if (head == nullptr || index == key.size()) {
		return;
	}

	// search children
	for (auto i = head->m_children.begin(); i != head->m_children.end(); i++) {
		// first character must match exactly
		if (index == 0) {
			if ((*i)->m_label == key[0]) {
				findSNiPHelper(key, index + 1, *i, v, false);
			}
		}

		// all succeeding characters can only have one mismatch
		else {
			if ((*i)->m_label == key[index]) {
				// if function has reached full length of key, copy every value into our vector v
				if (index == key.size() - 1) {
					for (int j = 0; j < (*i)->m_values.size(); j++) {
						v.push_back((*i)->m_values[j]);
					}
				}

				findSNiPHelper(key, index + 1, *i, v, alreadyMismatch);
			}

			// if the function hasn't met a mismatch yet, search through all children,
			// but if there is already a mismatch, don't do anything
			else if (!alreadyMismatch) {
				// if function has reached full length of key, copy every value into our vector v
				if (index == key.size() - 1) {
					for (int j = 0; j < (*i)->m_values.size(); j++) {
						v.push_back((*i)->m_values[j]);
					}
				}

				// now we pass in true for boolean value
				findSNiPHelper(key, index + 1, *i, v, true);
			}
		}
	}
}

// TESTING ONLY
template<typename ValueType>
void Trie<ValueType>::printHelper(Node* root) {
	if (root != nullptr) {
		for (auto i = root->m_children.begin(); i != root->m_children.end(); i++) {
			std::cout << (*i)->m_label;
			for (int j = 0; j < root->m_values.size(); j++) {
				std::cout << (*i)->m_values[j];
			}

			printHelper(*i);
		}
	}
}

#endif // TRIE_INCLUDED