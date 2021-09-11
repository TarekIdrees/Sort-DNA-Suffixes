# Sort-DNA-Suffixes
* C++ program to sort DNA suffixes from a given fasta file  using two method (Prefix Doubling , Quick Sort) with complexity (n log(n) , n^2 log(n) ) .
# Functions
* save_Genome
```
Function to save the whole genom as a type of char* .
```
* save_Sequence_TienK_Char
```
Function to save first 10,000 charachter of sequence as a type of char*.
```
* Build_Suffixes
```
Function to build array of suffixes that contain index of each suffix.
```
* Print_Suffixes
```
Function to print each suffix of string using the index of the suffix only.
```
* comp_nSeqLog
```
Function to compare the two char until find the smallest or when the length of one of them ends.
```
* Sort_nSeqLog
```
Function to sort suffixes of the sequence with complexity n^2 log(n).
```
* comp_prefix_Seq , comp_prefix_Genom
```
Function to compare the letters according to their ASSCI code.
```
* Prefix_Doubling_Seq , Prefix_Doubling_Genom
```
Function to sort suffixes of the sequence with complexity n log(n).
``` 
* Print_TienKChar_Suffixes_nSeqLogn
```
Function to print the sorted suffixes of complexity n^2 log(n) and time duration in file format.
```
* Print_TienKChar_Suffixes_nLogSeqn
```
Function to print the sorted suffixes of complexity n log(n) and time duration in file format.
```
* Print_Genom_File
```
Function to print the sorted suffixes of complexity n log(n) and time duration of whole genome in file format.
```
* Delete_Sequence
```
Function to delete the sequence and the suffixes.
