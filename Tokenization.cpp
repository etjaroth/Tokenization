#include <codecvt>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <fcntl.h>
#include <io.h>

#include <chrono>
#include <vector>

////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

template<typename T, unsigned int BUFFER_SIZE>
class CircularBuffer {
	T buffer[BUFFER_SIZE] = { 0 };
	unsigned int start = 0;
	unsigned int length = 0;

public:
	CircularBuffer() {
		for (int i = 0; i < BUFFER_SIZE; ++i) {
			buffer[i] = 0;
		}
	}

	void Take(T wc) {
		buffer[(start + length) % BUFFER_SIZE] = wc;

		if (length == BUFFER_SIZE) {
			start = (start + 1) % BUFFER_SIZE;
		}
		else {
			length = std::min(length + 1, BUFFER_SIZE);
		}
	}

	T Get(unsigned int i) const {
		return buffer[(start + i) % BUFFER_SIZE];
	}

	unsigned int GetLength() const {
		return length;
	}

	std::wstring GetDataAsWString() const {
		std::wstringstream wss{};
		for (unsigned int i = 0; i < length; ++i) {
			wss << std::wstring{Get(i)};
		}
		return wss.str();
	}

	void RemoveFirst(unsigned int n) {
		if (n >= length) {
			start = 0;
			length = 0;
			return;
		}

		length -= n;
		start = (start + n) % BUFFER_SIZE;
	}
};


class PermutationCounter
{
public:
	const static unsigned int WindowSize = 2;

private:
	using score_data = double;
	//const static unsigned int window_size = 15; // 16 - 1 for null char
	//const static unsigned int window_size = 7; // 8 - 1 for null char

	CircularBuffer<wchar_t, WindowSize> buffer{};

	unsigned int totalFrequencyCount[WindowSize + 1];
	std::unordered_map<std::wstring, unsigned int> frequencyCount[WindowSize + 1];
	std::unordered_map<std::wstring, score_data> chi2score[WindowSize + 1];
	std::unordered_map<std::wstring, score_data> totalScore;

public:
	PermutationCounter() {
		for (int i = 0; i < WindowSize + 1; ++i) {
			totalFrequencyCount[i] = 0;
			frequencyCount[i] = {};
			chi2score[i] = {};
		}
	}

	void Take(wchar_t wc) {
		buffer.Take(wc);
		if (buffer.GetLength() < WindowSize) {
			return;
		}

		for (unsigned int pLen = 1; pLen <= buffer.GetLength(); ++pLen) {
			std::wstring wstr = GetPermutationFromBuffer(pLen);
			if (wstr.length() != pLen) {
				break; // Flush is done
			}

			std::unordered_map<std::wstring, unsigned int>::iterator itr = frequencyCount[wstr.length()].find(wstr);

			if (itr == frequencyCount[wstr.length()].end()) {
				frequencyCount[wstr.length()].insert({ wstr, 1 });
			}
			else {
				itr->second += 1;
			}

			totalFrequencyCount[pLen] += 1;
		}
	}

	void Flush() {
		for (int i = 0; i < WindowSize; ++i) {
			Take(L'\0');
		}
	}

	void CalculateWeights() {
		const static score_data zero = 0;
		const static score_data one = 1;

		for (unsigned int i = 0; i <= WindowSize; ++i) {
			for (auto itr = frequencyCount[i].begin(); itr != frequencyCount[i].end(); ++itr) {
				const std::wstring& token = itr->first;
				const unsigned int y = frequencyCount[i][token];
				const unsigned int n = totalFrequencyCount[i];

				//const score_data nullTheta = std::log(static_cast<score_data>(y) / static_cast<score_data>(n));
				//score_data altTheta = 0;
				//for (wchar_t wc : token) { // ln(ab) = ln(a) + ln(b)
				//	altTheta = std::log(static_cast<score_data>(frequencyCount[1][std::wstring{wc}]) / static_cast<score_data>(totalFrequencyCount[1]));
				//}
				//altTheta = -altTheta;
				
				//const score_data score = abs(nullTheta - altTheta);
				//totalScore.insert({ token, score });
				totalScore.insert({ token, static_cast<score_data>(y) / static_cast<score_data>(n) });
			}
		}

		//// Score Monograms
		//for (auto itr = frequencyCount[1].begin(); itr != frequencyCount[1].end(); ++itr) {
		//	const std::wstring& token = itr->first;
		//	const unsigned int y = frequencyCount[1][token];
		//	const unsigned int n = totalFrequencyCount[1];

		//	const score_data nullTheta = static_cast<score_data>(y) / static_cast<score_data>(n);
		//	const score_data altTheta = one / static_cast<score_data>(frequencyCount[1].size());

		//	score_data score = (altTheta == nullTheta);
		//	if (!(
		//		(altTheta == zero) ||
		//		(nullTheta == zero) ||
		//		(altTheta == one) ||
		//		(nullTheta == one)
		//		)) {
		//		score = ((score_data)2) * (BinomRLF(y, n, nullTheta) - BinomRLF(y, n, altTheta));
		//	}

		//	chi2score[1].insert({ token, score });
		//}

		//// Score N-Grams
		//for (unsigned int i = 2; i <= window_size; ++i) {
		//	for (auto itr = frequencyCount[i].begin(); itr != frequencyCount[i].end(); ++itr) {
		//		const std::wstring& token = itr->first;
		//		const std::wstring lesserToken = token.substr(0, token.length() - 1);
		//		const std::wstring lastWChar{ token.back() };

		//		const unsigned int y = itr->second;
		//		const unsigned int n = frequencyCount[lesserToken.length()][lesserToken];

		//		const score_data nullTheta = static_cast<score_data>(y) / static_cast<score_data>(n);
		//		const score_data altTheta = static_cast<score_data>(frequencyCount[1][std::wstring{ token[token.length() - 1] }])
		//			/ static_cast<score_data>(totalFrequencyCount[1]);

		//		score_data score = 0;
		//		if (!(
		//			(altTheta == zero) ||
		//			(nullTheta == zero) ||
		//			(altTheta == one) ||
		//			(nullTheta == one)
		//			)) {
		//			score = ((score_data)2) * (BinomRLF(y, n, nullTheta) - BinomRLF(y, n, altTheta));
		//		}

		//		chi2score[i].insert({ token, score });
		//	}
		//}

		//// Free memory
		//for (int i = 0; i <= window_size; ++i) {
		//	frequencyCount[i].clear();
		//}
	}

	void SumWeights() {
		// Copy Monogram scores
		for (auto itr = chi2score[1].begin(); itr != chi2score[1].end(); ++itr) {
			totalScore.insert(*itr);
		}

		// N-Grams
		for (int i = 2; i <= WindowSize; ++i) {
			for (auto itr = chi2score[i].begin(); itr != chi2score[i].end(); ++itr) {
				const std::wstring token = itr->first;

				score_data score = itr->second;
				score += totalScore[token.substr(0, token.length() - 1)];

				for (wchar_t wc : token) {
					//score -= totalScore[std::wstring{ wc }];
					score -= chi2score[1][std::wstring{ wc }];
				}

				totalScore.insert({ token, score });
			}
		}

		// Free memory
		for (int i = 0; i <= WindowSize; ++i) {
			chi2score[i].clear();
		}
	}

	void OutputTokensToFile(std::string filename) {
		std::vector<std::pair<std::wstring, score_data>> tokens;
		for (auto itr = totalScore.begin(); itr != totalScore.end(); ++itr) {
			tokens.push_back(*itr);
		}

		std::sort(tokens.begin(), tokens.end(),
			[](std::pair<std::wstring, score_data>& a, std::pair<std::wstring, score_data>& b) {
				return a.second < b.second;
			});

		std::wofstream wof(filename);
		wof.imbue(std::locale(std::locale::empty(), new std::codecvt_utf8<wchar_t>));

		unsigned int i = 0;
		for (auto itr = tokens.rbegin(); (itr != tokens.rend()); ++itr, ++i) {
			wof << std::setw(5) << i << " | " <<
				"\"" << itr->first << "\"" <<
				std::setw(WindowSize + 3 - itr->first.length()) << std::wstring{ L" | " } <<
				std::setw(2) << itr->first.length() << std::wstring(L" | ") << itr->second << std::endl;
		}
		wof << std::endl;
		wof.close();
	}

	std::wstring BestNextToken(std::wstring window) const {
		score_data bestScore = std::numeric_limits<score_data>::lowest();
		std::wstring bestToken = L"";

		for (int i = 1; i <= window.length(); ++i) { // token[0:i]
			const std::wstring token = window.substr(0, i); // abcde
			auto itr = totalScore.find(token);
			if (itr == totalScore.end()) {
				std::wcout << "Couldn't find " << token << " (token)" << std::endl;
				continue;
			}
			score_data score = itr->second;

			//for (int j = 0; j < window.length(); ++j) {
				//score += (j < i ? 1.0 : -1.0) * totalScore.find(std::wstring{ window[j] })->second;
			//}

			/*const std::wstring wsubstr = window.substr(i);
			if (wsubstr.length() > 0) {
				itr = totalScore.find(wsubstr);
				if (itr == totalScore.end()) {
					std::wcout << "Couldn't find \"" << wsubstr << "\" (wsubstr)" << std::endl;
					continue;
				}
				score -= itr->second;
			}*/

			//for (int j = 1; j < i; ++j) { // a bcde, ab cde, abc de
			//	const std::wstring wsubstr = window.substr(j, i - j);
			//	score -= totalScore.find(wsubstr)->second;
			//}

			if (score > bestScore) {
				bestScore = score;
				bestToken = token;
			}
		}

		return bestToken;
	}

private:
	std::wstring GetPermutationFromBuffer(unsigned int pLen) {
		wchar_t* data = new wchar_t[pLen + 1];

		data[pLen] = 0;
		for (unsigned int i = 0; i < pLen; ++i) {
			data[i] = buffer.Get(i);
		}

		std::wstring wstr(data);
		delete[] data;
		return wstr;
	}

	score_data BinomRLF(unsigned int y, unsigned int n, score_data theta) {
		const static score_data one = 1;
		return static_cast<score_data>(y) * std::log(theta) + static_cast<score_data>(n - y) * std::log(one - theta);
	}

} pCounter{};

////////////////////////////////////////////////////////////////////////

int main()
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
	{
		std::wifstream wif("Data\\Worm.txt");
		//std::wifstream wif("Data\\Perm.txt");
		wif.imbue(std::locale(std::locale::empty(), new std::codecvt_utf8<wchar_t>));

#pragma warning(push)
#pragma warning(disable:6031)
		_setmode(_fileno(stdout), _O_U16TEXT); // Allows the console to print unicode
#pragma warning(pop)

		wchar_t wc;
		while (wif >> std::noskipws >> wc) {
			//std::wcout << wc;
			pCounter.Take(wc);
		}
		wif.close();

		pCounter.Flush();

		//std::wcout << std::endl;
	}
	std::chrono::time_point<std::chrono::high_resolution_clock> t2 = std::chrono::high_resolution_clock::now();
	std::wcout << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << "s" << std::endl;

	pCounter.CalculateWeights();

	std::chrono::time_point<std::chrono::high_resolution_clock> t3 = std::chrono::high_resolution_clock::now();
	std::wcout << std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count() << "s" << std::endl;

	//pCounter.SumWeights();
	pCounter.OutputTokensToFile("Data\\tokens.txt");

	std::chrono::time_point<std::chrono::high_resolution_clock> t4 = std::chrono::high_resolution_clock::now();
	std::wcout << std::chrono::duration_cast<std::chrono::seconds>(t4 - t3).count() << "s" << std::endl;

	// Tokenize
	{
		unsigned int i = 0;
		std::wifstream wif("Data\\Worm.txt");
		wif.imbue(std::locale(std::locale::empty(), new std::codecvt_utf8<wchar_t>));

		std::wofstream wof("Data/Tokenization.txt");
		wof.imbue(std::locale(std::locale::empty(), new std::codecvt_utf8<wchar_t>));

		CircularBuffer<wchar_t, pCounter.WindowSize> buffer{};

		while (wif && (++i < 100000)) {
			while (buffer.GetLength() < pCounter.WindowSize) {
				wchar_t wc;
				wif >> std::noskipws >> wc;
				buffer.Take(wc);
			}

			std::wstring wstr = buffer.GetDataAsWString();
			std::wstring token = pCounter.BestNextToken(wstr);
			std::wcout << wstr.length() << " : [" << token << "]" << std::endl;
			wof << "[[" << token << "]]";
			buffer.RemoveFirst(token.length());
		}

		wof << std::endl;
		wof.close();
		wif.close();
	}
}
