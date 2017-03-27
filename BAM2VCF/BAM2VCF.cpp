//============================================================================
// Name        : BAM2VCF.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <map>
#include <assert.h>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <set>
#include <tuple>
#include <utility>
#include <algorithm>

using namespace std;

std::vector<std::string> split(std::string input, std::string delimiter);
std::string join(std::vector<std::string> parts, std::string delim);
void eraseNL(std::string& s);
int StrtoI(std::string s);
std::string ItoStr(int i);
unsigned int StrtoUI(std::string s);
std::string removeGaps(std::string in);

class startingHaplotype;
void produceVCF(const std::string referenceSequenceID, const std::string& referenceSequence, const std::map<unsigned int, std::vector<startingHaplotype*>>& alignments_starting_at, std::string outputFn);
void printHaplotypesAroundPosition(const std::string& referenceSequence, const std::map<unsigned int, std::vector<startingHaplotype*>>& alignments_starting_at, int posI);

class startingHaplotype
{
public:
	std::string ref;
	std::string query;
	std::string query_name;
	unsigned int aligment_start_pos;
	unsigned int alignment_last_pos;
};

int main(int argc, char *argv[]) {
	std::vector<std::string> ARG (argv + 1, argv + argc + !argc);
	std::map<std::string, std::string> arguments;

	arguments["input"] = "C:\\Users\\diltheyat\\Desktop\\Temp\\chr21";
	arguments["referenceSequenceID"] = "chr21";

	for(unsigned int i = 0; i < ARG.size(); i++)
	{
		if((ARG.at(i).length() > 2) && (ARG.at(i).substr(0, 2) == "--"))
		{
			std::string argname = ARG.at(i).substr(2);
			std::string argvalue = ARG.at(i+1);
			arguments[argname] = argvalue;
		}
	}

	assert(arguments.count("input"));
	assert(arguments.count("referenceSequenceID"));

	std::ifstream inputStream;
	inputStream.open(arguments.at("input").c_str());
	if(! inputStream.is_open())
	{
		throw std::runtime_error("Could not open file "+arguments.at("input"));
	}
	assert(inputStream.good());
	std::string referenceSequence;
	std::getline(inputStream, referenceSequence);
	eraseNL(referenceSequence);

	int max_gap_length = 100;
	int n_alignments_loaded = 0;
	int n_alignments_skipped = 0;
	std::map<unsigned int, std::vector<startingHaplotype*>> alignments_starting_at;
	std::string line;
	while(inputStream.good())
	{
		std::getline(inputStream, line);
		eraseNL(line);
		if(line.length())
		{
			std::vector<std::string> line_fields = split(line, "\t");
			assert(line_fields.size() == 5);
			startingHaplotype* h = new startingHaplotype();
			h->ref = line_fields.at(0);
			h->query = line_fields.at(1);
			h->query_name = line_fields.at(2);
			h->aligment_start_pos = StrtoUI(line_fields.at(3));
			h->alignment_last_pos = StrtoUI(line_fields.at(4))+1;

			int running_gap_length = 0;
			int max_running_gap_length = 0;
			for(unsigned int i = 0; i < h->ref.length(); i++)
			{
				unsigned char c_ref = h->ref.at(i);
				if((c_ref == '-') or (c_ref == '*'))
				{
					running_gap_length++;
				}
				else
				{
					if(running_gap_length)
					{
						if(running_gap_length > max_running_gap_length)
							max_running_gap_length = running_gap_length;
					}
					running_gap_length = 0;
				}
			}
			assert(running_gap_length == 0);

			// std::cout << "Alignment " << h->query_name << " max running gap: " << max_running_gap_length << "\n" << std::flush;
			if(max_running_gap_length <= max_gap_length)
			{
				alignments_starting_at[h->aligment_start_pos].push_back(h);
				n_alignments_loaded++;
			}
			else
			{
				delete(h);
				n_alignments_skipped++;
			}
		}
	}
	std::cout << "For max. gap length " << max_gap_length << "\n";
	std::cout << "\t" << "n_alignments_loaded" << ": " << n_alignments_loaded << "\n";
	std::cout << "\t" << "n_alignments_skipped" << ": " << n_alignments_skipped << "\n";
	std::cout << std::flush;

	produceVCF(arguments.at("referenceSequenceID"), referenceSequence, alignments_starting_at, arguments.at("input")+".VCF");

	return 0;
}

void produceVCF(const std::string referenceSequenceID, const std::string& referenceSequence, const std::map<unsigned int, std::vector<startingHaplotype*>>& alignments_starting_at, std::string outputFn)
{
	std::ofstream outputStream;
	outputStream.open(outputFn.c_str());
	if(! outputStream.is_open())
	{
		throw std::runtime_error("Cannot open " + outputFn + " for writing!");
	}

	int n_alignments = 0;
	for(auto startPos : alignments_starting_at)
	{
			n_alignments += startPos.second.size();
	}

	std::vector<int> gap_structure;
	std::vector<int> coverage_structure;
	gap_structure.resize(referenceSequence.length(), -1);
	coverage_structure.resize(referenceSequence.length(), 0);
	int examine_gaps_n_alignment = 0;
	for(auto startPos : alignments_starting_at)
	{
		for(startingHaplotype* alignment : startPos.second)
		{
			assert(startPos.first == alignment->aligment_start_pos);

			int start_pos = (int)startPos.first - 1;
			int ref_pos = start_pos;
			int running_gaps = 0;

			for(unsigned int i = 0; i < alignment->ref.length(); i++)
			{
				unsigned char c_ref = alignment->ref.at(i);
				if((c_ref == '-') or (c_ref == '*'))
				{
					running_gaps++;
				}
				else
				{
					if(ref_pos != start_pos)
					{
						if(gap_structure.at(ref_pos) == -1)
						{
							gap_structure.at(ref_pos) = running_gaps;
						}
						else
						{
							if(gap_structure.at(ref_pos) != running_gaps)
							{
								std::cerr << "Gap structure mismatch at position " << ref_pos << " - this is alignment " << examine_gaps_n_alignment << " / " << alignment->query_name << ", have existing value " << gap_structure.at(ref_pos) << ", want to set " << running_gaps << "\n" << std::flush;
								std::cerr << "Alignment start " << alignment->aligment_start_pos << "\n";
								std::cerr << "Alignment stop " << alignment->alignment_last_pos << "\n";
								std::cerr << std::flush;
								throw std::runtime_error("Gap structure mismatch");
							}

						}
					}
					ref_pos++;
					running_gaps = 0;
					coverage_structure.at(ref_pos)++;
				}
			}

			examine_gaps_n_alignment++;
		}
	}

	std::cout << "Loaded " << examine_gaps_n_alignment << " alignments.\n";

	std::cout << "Coverage structure:\n";
	int coverage_window_length = 10000;
	for(unsigned int pI = 0; pI < coverage_structure.size(); pI += coverage_window_length)
	{
		int coverage_in_window = 0;
		unsigned int last_window_pos = (pI+coverage_window_length) - 1;
		if(last_window_pos > (coverage_structure.size() - 1))
			last_window_pos = coverage_structure.size() - 1;

		for(unsigned int j = pI; j <= last_window_pos; j++)
		{
			assert(j < coverage_structure.size());
			coverage_in_window += coverage_structure.at(j);
		}
		double avg_coverage = (double) coverage_in_window / (double)(last_window_pos - pI + 1);

		if((pI >= 15000000) && (pI <= 17000000))
			std::cout << "\t" << "Window starting at pI = " << pI << " => avg. coverage " << avg_coverage << "\n";
	}
	std::cout << std::flush;

	std::set<const startingHaplotype*> known_haplotype_pointers;
	for(auto startingPos : alignments_starting_at)
	{
		for(startingHaplotype* sH :  startingPos.second)
		{
			known_haplotype_pointers.insert(sH);
		}
	}
	using openHaplotype = std::tuple<std::string, const startingHaplotype*, int>;
	std::vector<openHaplotype> open_haplotypes;
	openHaplotype initH = std::make_tuple("", (const startingHaplotype*)0, -1);
	open_haplotypes.push_back(initH);

	int start_open_haplotypes = 0;
	int opened_alignments = 0;

	for(int posI = 0; posI < (int)referenceSequence.length(); posI++)
	{
		if(((posI % 10000) == 0) or (0 && open_haplotypes.size() > 100))
		{
			std::cout << posI << ", open haplotypes: " << open_haplotypes.size() << "\n";
		}

		for(openHaplotype& haplotype : open_haplotypes)
		{
			//consume all gaps "before" the current reference position
			if(std::get<1>(haplotype) != 0)
			{
				if(std::get<2>(haplotype) == ((int)std::get<1>(haplotype)->ref.length() - 1))
				{
					int n_gaps = gap_structure.at(posI-1);
					if(n_gaps == -1)
						n_gaps = 0;

					std::string gaps;
					gaps.resize(n_gaps, '-');
					assert((int)gaps.length() == n_gaps);
					std::get<0>(haplotype) += gaps;
				}
				else
				{
					if(((std::get<1>(haplotype)->ref.at(std::get<2>(haplotype))) == '-') || (std::get<1>(haplotype)->ref.at(std::get<2>(haplotype)) == '*'))
					{
						std::cerr << "Position " << std::get<2>(haplotype) << "is gap in one of our haplotypes!";
					}

					int nextPos = std::get<2>(haplotype)+1;
					std::string additionalExtension;
					while((nextPos < (int)std::get<1>(haplotype)->ref.length()) && ((std::get<1>(haplotype)->ref.at(nextPos) == '-') || (std::get<1>(haplotype)->ref.at(nextPos) == '*')))
					{
						additionalExtension += std::get<1>(haplotype)->query.substr(nextPos, 1);
						nextPos++;
					}
					int consumedUntil = nextPos - 1;

					std::get<0>(haplotype).append(additionalExtension);
					std::get<2>(haplotype) = consumedUntil;
				}
			}
			else
			{
				if(posI > 0)
				{

					int n_gaps = gap_structure.at(posI-1);
					if(n_gaps == -1)
						n_gaps = 0;

					std::string gaps;
					gaps.resize(n_gaps, '-');
					assert((int)gaps.length() == n_gaps);
					std::get<0>(haplotype) += gaps;
				}
			}
		}

		// check that all extensions - i.e. up to the current reference position - have the same length
		int assembled_h_length = -1;
		for(openHaplotype haplotype : open_haplotypes)
		{
			if(assembled_h_length == -1)
			{
				assembled_h_length = std::get<0>(haplotype).length();
			}

			if(assembled_h_length != (int)std::get<0>(haplotype).length())
			{
				std::cerr << "Initial II length mismatch " << posI << " " << assembled_h_length << "\n"; // [@gap_structure[(posI-3) .. (posI+1)]]
				for(openHaplotype oH2 : open_haplotypes)
				{
					std::cerr << "\t" << std::get<0>(oH2).length() << "\tconsumed until: " << std::get<2>(oH2) << ", of length " << ((std::get<1>(oH2) == 0) ? "REF" : ("nonRef " + std::get<1>(oH2)->query_name + " / length " + ItoStr(std::get<1>(oH2)->ref.length()))) << "\n";
				}
				printHaplotypesAroundPosition(referenceSequence, alignments_starting_at, posI);
				assert(2 == 4);
			}
		}

		// copy in new haplotypes
		unsigned char refC = referenceSequence.at(posI);
		std::vector<const startingHaplotype*> new_haplotypes;
		if(alignments_starting_at.count(posI))
		{
			for(startingHaplotype* sH :  alignments_starting_at.at(posI))
			{
				assert(known_haplotype_pointers.count(sH));
				new_haplotypes.push_back(sH);
			}
		}
		unsigned int open_haplotypes_size = open_haplotypes.size();
		for(const startingHaplotype* new_haplotype : new_haplotypes)
		{
			if(open_haplotypes_size > 0)
			{
				opened_alignments++;

				for(int existingHaploI = 0; existingHaploI < (int)open_haplotypes_size; existingHaploI++)
				{
					openHaplotype new_haplotype_copy_this = std::make_tuple(std::get<0>(open_haplotypes.at(existingHaploI)), new_haplotype, -1);
					open_haplotypes.push_back(new_haplotype_copy_this);
				}

				int open_span = posI - start_open_haplotypes;
				int start_reference_extraction = start_open_haplotypes;
				int stop_reference_extraction = posI - 1;
				assert(stop_reference_extraction >= start_reference_extraction); // die Dumper("Weird", start_reference_extraction, stop_reference_extraction) unless(stop_reference_extraction >= start_reference_extraction);
				std::string referenceExtraction;
				referenceExtraction.reserve(stop_reference_extraction - start_reference_extraction + 1);
				for(int refI = start_reference_extraction; refI <= stop_reference_extraction; refI++)
				{
					referenceExtraction.push_back(referenceSequence.at(refI));

					int n_gaps = gap_structure.at(refI);
					if(n_gaps == -1)
						n_gaps = 0;
					std::string gaps;
					gaps.resize(n_gaps, '-');
					assert((int)gaps.length() == n_gaps);
					referenceExtraction.append(gaps);
				}

				//new_haplotype_referenceSequence = [substr(referenceSequence, start_open_haplotypes, open_span), new_haplotype, -1];
				openHaplotype new_haplotype_referenceSequence = std::make_tuple(referenceExtraction, new_haplotype, -1);

				//missing = assembled_h_length - open_span;
				//die Dumper(posI, start_open_haplotypes, missing, open_span, assembled_h_length) unless(missing >= 0);
				//missingStr = '*' x missing;
				//die unless(length(missingStr) == missing);
				//new_haplotype_referenceSequence->[0] .= missingStr;

				open_haplotypes.push_back(new_haplotype_referenceSequence);

				std::cout << "Position " << posI << ", enter new haplotype " << new_haplotype->query_name << " --> " << open_haplotypes.size() << " haplotypes.\n" << std::flush;

			}
		}


		if(posI == 7652900)
		{
			std::cerr << "Pre-exit haplotype lengths " << posI << "\n"; // [@gap_structure[(posI-3) .. (posI+1)]]
			for(openHaplotype oH2 : open_haplotypes)
			{
				std::cerr << "\t" << std::get<0>(oH2).length() << "\tconsumed until: " << std::get<2>(oH2) << ", of length " << ((std::get<1>(oH2) == 0) ? "REF" : ("nonRef " + std::get<1>(oH2)->query_name + " / length " + ItoStr(std::get<1>(oH2)->ref.length()))) << "\n";
			}
		}

		// we exit these haplotypes be making them ref / another haplotype
		open_haplotypes_size = open_haplotypes.size();
		std::set<unsigned int> exitedHaplotype;
		for(unsigned int outer_haplotype_I = 0; outer_haplotype_I < open_haplotypes_size; outer_haplotype_I++)
		{
			openHaplotype& haplotype = open_haplotypes.at(outer_haplotype_I);

			if(std::get<1>(haplotype) != 0) // i.e. non-ref
			{
				if(std::get<2>(haplotype) == ((int)std::get<1>(haplotype)->ref.length() - 1))
				{
					std::cerr << "Position " << posI << ", exit haplotype " << std::get<1>(haplotype)->query_name << " length " << std::get<0>(haplotype).length() << "\n" << std::flush;
					// print "exit one\n";
					std::get<1>(haplotype) = 0;
					std::get<2>(haplotype) = -1;
					exitedHaplotype.insert(outer_haplotype_I);

					size_t expected_haplotype_length = std::get<0>(haplotype).length();
					std::cerr << "\texpected_haplotype_length: " << expected_haplotype_length << "\n";

					for(int existingHaploI = 0; existingHaploI < (int)open_haplotypes_size; existingHaploI++)
					{
						if(existingHaploI == existingHaploI)
						{
							continue;
						}

						if(exitedHaplotype.count(existingHaploI))
							continue;

						assert(std::get<0>(haplotype).length() == expected_haplotype_length);
						openHaplotype new_haplotype_copy_this = std::make_tuple(std::string(std::get<0>(haplotype)), std::get<1>(open_haplotypes.at(existingHaploI)), std::get<2>(open_haplotypes.at(existingHaploI)));
						assert(std::get<0>(haplotype).length() == expected_haplotype_length);

						if((std::get<1>(new_haplotype_copy_this) == 0) || (std::get<2>(new_haplotype_copy_this) != ((int)std::get<1>(new_haplotype_copy_this)->ref.length() - 1)))
						{
							assert((std::get<1>(haplotype) == 0) || (std::get<2>(haplotype) != ((int)std::get<1>(haplotype)->ref.length() - 1)));
							assert(std::get<0>(haplotype).length() == std::get<0>(new_haplotype_copy_this).length());
							if(posI == 7652900)
							{
								std::cerr << "Position " << posI << " add of length " << std::get<0>(new_haplotype_copy_this).length() << "\n"; // [@gap_structure[(posI-3) .. (posI+1)]]
							}
							assert(std::get<0>(haplotype).length() == expected_haplotype_length);
							assert(std::get<0>(new_haplotype_copy_this).length() == expected_haplotype_length);
							open_haplotypes.push_back(new_haplotype_copy_this);
						}
					}

					assert((std::get<1>(haplotype) == 0) || known_haplotype_pointers.count(std::get<1>(haplotype)));

					//assert(std::get<1>(haplotype) != 0);
					//std::cout << "Position " << posI << ", exit haplotype " << std::get<1>(haplotype)->query_name << " --> " << open_haplotypes.size() << " haplotypes.\n" << std::flush;
				}
			}
		}


		if(posI == 7652900)
		{
			std::cerr << "Post-exit haplotype lengths " << posI << "\n"; // [@gap_structure[(posI-3) .. (posI+1)]]
			for(openHaplotype oH2 : open_haplotypes)
			{
				std::cerr << "\t" << std::get<0>(oH2).length() << "\tconsumed until: " << std::get<2>(oH2) << ", of length " << ((std::get<1>(oH2) == 0) ? "REF" : ("nonRef " + std::get<1>(oH2)->query_name + " / length " + ItoStr(std::get<1>(oH2)->ref.length()))) << "\n";
			}
		}


		// print "\tLength ", assembled_h_length, "\n";

		/*
		if(1 == 0)
		{
			print "Haplotype info:\n";
			for(existingHaploI = 0; existingHaploI <= #open_haplotypes; existingHaploI++)
			{
				print "\t", existingHaploI, "\n";
				print "\t\t", open_haplotypes[existingHaploI][0], "\n";
				print "\t\t", open_haplotypes[existingHaploI][2], "\n";
				if(open_haplotypes[existingHaploI][1])
				{
					ref_str = open_haplotypes[existingHaploI][1][0];
					haplo_str = open_haplotypes[existingHaploI][1][1];
					print "\t\t", open_haplotypes[existingHaploI][1][2], "\n";
					printFrom = open_haplotypes[existingHaploI][2];
					printFrom = 0 if(printFrom < 0);
					print "\t\t", substr(ref_str, printFrom, 10), "\n";
					print "\t\t", substr(haplo_str, printFrom, 10), "\n";
				}
				else
				{
					print "\t\tREF\n";
				}
			}
			print "\n";
		}
		*/

		// we consume characters up to (but not including) the next reference non-gap
		std::set<std::string> extensions_nonRef;
		int extensions_nonRef_length = -1;
		for(openHaplotype haplotype : open_haplotypes)
		{
			std::string extension;
			int consumed_ref_start = -1;
			int consumed_ref = 0;
			std::string consumed_ref_sequence;

			if(std::get<1>(haplotype) == 0)
			{

			}
			else
			{
				int addIndex = 0;
				consumed_ref_start = std::get<2>(haplotype)+addIndex+1;
				do {
					int nextPosToConsume = std::get<2>(haplotype)+addIndex+1;
					if(!(nextPosToConsume < (int)std::get<1>(haplotype)->ref.length()))
					{
						std::cerr << "nextPosToConsume" << ": " << nextPosToConsume << "\n";
						std::cerr << "std::get<1>(haplotype)->ref.length()" << ": " << std::get<1>(haplotype)->ref.length() << "\n";
						std::cerr << "extension" << ": " << extension << "\n";
						std::cerr << std::flush;
					}
					assert(nextPosToConsume < (int)std::get<1>(haplotype)->ref.length());
					unsigned char refC = std::get<1>(haplotype)->ref.at(nextPosToConsume);
					consumed_ref_sequence.push_back(refC);
					if((refC != '-') && (refC != '*'))
					{
						consumed_ref++;
					}
					unsigned char hapC = std::get<1>(haplotype)->query.at(nextPosToConsume);
					extension.push_back(hapC);
					addIndex++;
				} while(consumed_ref < 1);
			}

			if(extension.length())
			{
				// push(@{extensions_nonRef{extension}}, [consumed_ref_start, consumed_ref, consumed_ref_sequence]);
				extensions_nonRef.insert(extension);
				if(extensions_nonRef_length == -1)
				{
					extensions_nonRef_length = extension.length();
				}
				assert((int)extension.length() == extensions_nonRef_length);
				/*
				unless(defined extensions_nonRef_length)
				{
					extensions_nonRef_length = length(extension);
				}
				die Dumper("Length mismatch", extension, \%extensions_nonRef, posI, "Length mismatch") unless(length(extension) == extensions_nonRef_length);
				*/
			}
		}

		// now carry out the actual extension
		std::set<std::string> extensions;
		for(openHaplotype& haplotype : open_haplotypes)
		{
			std::string extension;
			if(std::get<1>(haplotype) == 0)
			{
				std::string refExt = {(char)refC};
				if(extensions_nonRef_length != -1)
				{
					int missing = extensions_nonRef_length - refExt.length();
					assert(missing >= 0);
					std::string missingStr;
					missingStr.resize(missing, '*');
					assert((int)missingStr.length() == missing);
					refExt.append(missingStr);
				}
				extension.append(refExt);
			}
			else
			{
				int consumed_ref = 0;
				do {
					int nextPosToConsume = std::get<2>(haplotype)+1;
					assert(nextPosToConsume < (int)std::get<1>(haplotype)->ref.length());
					unsigned char refC = std::get<1>(haplotype)->ref.at(nextPosToConsume);
					if((refC != '-') && (refC != '*'))
					{
						consumed_ref++;
					}
					unsigned char hapC = std::get<1>(haplotype)->query.at(nextPosToConsume);
					extension.push_back(hapC);
					std::get<2>(haplotype)++;
				} while(consumed_ref < 1);
			}
			assert(extension.length());
			std::get<0>(haplotype).append(extension);
			extensions.insert(extension);
		}
		assert(extensions.size());

		// print "Extensions:\n", join("\n", map {"\t'"._."'"} keys %extensions), "\n\n";

		//#this_all_equal = ( (scalar(keys %extensions) == 0) or ((scalar(keys %extensions) == 1) and (exists extensions{refC})) );
		std::string refC_string = {(char)refC};
		bool this_all_equal = ((extensions.size() == 1) && (extensions.count(refC_string)));
		if(posI == 0)
		{
			assert(this_all_equal);
		}
		/*
		if(open_haplotypes.size() > 100)
		{
			std::cout << "Open haplotypes position " << posI << "\n";
			for(auto e : extensions)
			{
				std::cout << "\t" << e << "\n" << std::flush;
			}
 		}*/

		if(this_all_equal && (posI > 0))
		{
			// close
			int ref_span = posI - start_open_haplotypes;
			assert(ref_span > 0);
			std::string reference_sequence = referenceSequence.substr(start_open_haplotypes, ref_span);
			std::set<std::string> alternativeSequences;
			std::set<std::string> uniqueRemainers;

			std::vector<openHaplotype> new_open_haplotypes;
			int open_haplotypes_before = open_haplotypes.size();
			for(openHaplotype& haplotype : open_haplotypes)
			{
				assert((int)std::get<0>(haplotype).length() >= (ref_span + 1)); //die Dumper("Length mismatch II", ref_span+1, length(std::get<0>(haplotype)), "Length mismatch II") unless(length(std::get<0>(haplotype)) >= (ref_span + 1));
				std::string haplotype_coveredSequence = std::get<0>(haplotype).substr(0, std::get<0>(haplotype).length()-1);
				haplotype_coveredSequence = removeGaps(haplotype_coveredSequence);
				if(haplotype_coveredSequence != reference_sequence)
				{
					alternativeSequences.insert(haplotype_coveredSequence);
				}
				std::get<0>(haplotype) = std::get<0>(haplotype).substr(std::get<0>(haplotype).length()-1);
				assert(std::get<0>(haplotype).length() == 1);

				std::stringstream uniqueRemainerKey;
				uniqueRemainerKey << std::get<0>(haplotype) << "//" << (void*)std::get<1>(haplotype) << "//" << std::get<2>(haplotype);

				std::string k = uniqueRemainerKey.str();
				if(open_haplotypes.size() > 100)
				{
					std::cout << k << "\n";
				}

				if(uniqueRemainers.count(k) == 0)
				{
					new_open_haplotypes.push_back(haplotype);
					uniqueRemainers.insert(k);
				}
			}

			open_haplotypes = new_open_haplotypes;
			int open_haplotypes_after = open_haplotypes.size();

			if(alternativeSequences.size())
			{
				bool all_alternativeAlleles_length_2 = true;
				for(auto a : alternativeSequences)
				{
					if(a.length() != 2)
						all_alternativeAlleles_length_2 = false;
				}

				// print "Starting at position start_open_haplotypes, have REF reference_sequence and alternative sequences " . join(' / ', @alternativeAlleles) . "\n";

				if((reference_sequence.length() == 2) and all_alternativeAlleles_length_2)
				{
					for(auto a : alternativeSequences)
					{
						assert(a.substr(0, 1) == reference_sequence.substr(0, 1)); // die Dumper("Some problem with supposed SNP", posI, reference_sequence, \@alternativeAlleles, "Some problem with supposed SNP") unless(substr(alt, 0, 1) eq firstRefChar);
					}

					std::vector<std::string> alternativeAlleles_reduced;
					for(auto a : alternativeSequences)
					{
						alternativeAlleles_reduced.push_back(a.substr(1,1));
					}

					outputStream <<
							referenceSequenceID << "\t" <<
							start_open_haplotypes+2 << "\t" <<
							"." << "\t" <<
							reference_sequence.substr(1,1) << "\t" <<
							join(alternativeAlleles_reduced, ",") << "\t" <<
							'.' << "\t" <<
							"PASS" << "\t" <<
							'.'
					<< "\n";
				}
				else
				{
					outputStream <<
							referenceSequenceID << "\t" <<
							start_open_haplotypes+1 << "\t" <<
							"." << "\t" <<
							reference_sequence << "\t" <<
							join(std::vector<std::string>(alternativeSequences.begin(), alternativeSequences.end()), ",") << "\t" <<
							'.' << "\t" <<
							"PASS" << "\t" <<
							'.'
					<< "\n";
				}
			}
			start_open_haplotypes = posI;

			if((n_alignments == opened_alignments) and (open_haplotypes_after == 1))
			{
				//
			}
			// std::cout << "Went from " << open_haplotypes_before << " to " << open_haplotypes_after << "\n";
		}

		if(posI == 7652900)
		{
			std::cerr << "Haplotype lengths " << posI << "\n"; // [@gap_structure[(posI-3) .. (posI+1)]]
			for(openHaplotype oH2 : open_haplotypes)
			{
				std::cerr << "\t" << std::get<0>(oH2).length() << "\tconsumed until: " << std::get<2>(oH2) << ", of length " << ((std::get<1>(oH2) == 0) ? "REF" : ("nonRef " + std::get<1>(oH2)->query_name + " / length " + ItoStr(std::get<1>(oH2)->ref.length()))) << "\n";
			}
		}

		// last_all_equal = this_all_equal;
	}

	std::cout << "Done.\n" << std::flush;
}



void printHaplotypesAroundPosition(const std::string& referenceSequence, const std::map<unsigned int, std::vector<startingHaplotype*>>& alignments_starting_at, int posI)
{
	std::cout << "Positions plot around " << posI << "\n" << std::flush;

	std::vector<int> positions;
	for(int i = posI - 2; i <= posI + 2; i++)
	{
		if(i >= 0)
			positions.push_back(i);
	}

	for(auto startPos : alignments_starting_at)
	{
		for(startingHaplotype* alignment : startPos.second)
		{
			assert(alignment->aligment_start_pos == startPos.first);
			int stopPos = alignment->alignment_last_pos;
			bool interesting = false;
			for(auto interestingPos : positions)
			{
				if((interestingPos >= (int)startPos.first) and (interestingPos <= stopPos))
				{
					interesting = true;
				}
			}

			if(interesting)
			{
				std::map<int, std::string> gt_per_position;

				int ref_pos = startPos.first - 1;
				std::string running_allele;

				for(int i = 0; i < (int)alignment->ref.length(); i++)
				{
					unsigned char c_ref = alignment->ref.at(i);
					unsigned char c_query = alignment->query.at(i);

					if((c_ref == '-') or (c_ref == '*'))
					{
						running_allele.push_back(c_query);
					}
					else
					{
						if(running_allele.length())
						{
							gt_per_position[ref_pos] = running_allele;
						}

						running_allele.clear();
						running_allele.push_back(c_query);
						ref_pos++;
					}
				}
				if(running_allele.length())
				{
					gt_per_position[ref_pos] = running_allele;
				}

				std::cout << "Positions " << alignment->query_name << "\n";
				for(auto interestingPos : positions)
				{
					if(gt_per_position.count(interestingPos))
					{
						std::cout << "\t" << interestingPos << "\t" << gt_per_position.at(interestingPos) << "\n";
					}
				}
			}
		}
	}

	std::cout << " -- end positions plot.\n" << std::flush;
}

vector<string> split(string input, string delimiter)
{
	vector<string> output;
	if(input.length() == 0)
	{
		return output;
	}

	if(delimiter == "")
	{
		output.reserve(input.size());
		for(unsigned int i = 0; i < input.length(); i++)
		{
			output.push_back(input.substr(i, 1));
		}
	}
	else
	{
		if(input.find(delimiter) == string::npos)
		{
			output.push_back(input);
		}
		else
		{
			int s = 0;
			int p = input.find(delimiter);

			do {
				output.push_back(input.substr(s, p - s));
				s = p + delimiter.size();
				p = input.find(delimiter, s);
			} while (p != (int)string::npos);
			output.push_back(input.substr(s));
		}
	}

	return output;
}

void eraseNL(string& s)
{
	if (!s.empty() && s[s.length()-1] == '\r') {
	    s.erase(s.length()-1);
	}
	if (!s.empty() && s[s.length()-1] == '\n') {
	    s.erase(s.length()-1);
	}
}

int StrtoI(string s)
{
	  stringstream ss(s);
	  int i;
	  ss >> i;
	  return i;
}

unsigned int StrtoUI(string s)
{
	  stringstream ss(s);
	  unsigned int i;
	  ss >> i;
	  return i;
}

string ItoStr(int i)
{
	std::stringstream sstm;
	sstm << i;
	return sstm.str();
}

string join(vector<string> parts, string delim)
{
	if(parts.size() == 0)
		return "";

	string ret = parts.at(0);

	for(unsigned int i = 1; i < parts.size(); i++)
	{
		ret.append(delim);
		ret.append(parts.at(i));
	}

	return ret;
}




std::string removeGaps(std::string in)
{
	std::string out;
	out.reserve(in.size());
	for(size_t i = 0; i < in.size(); i++)
	{
		if((in.at(i) != '_') && (in.at(i) != '-') && (in.at(i) != '*'))
		{
			out.push_back(in.at(i));
		}
	}
	return out;
}



