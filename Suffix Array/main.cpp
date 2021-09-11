#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <chrono>
using namespace std::chrono;
using namespace std;

////////////////////////////////////////////////////////////////////////
#define DIR_DNA "C://Users//sfahg//OneDrive//Desktop//Assignment 1//"
////////////////////////////////////////////////////////////////////////

// Time variable  //
long  int Time;


/// Global Variables for first two Files  ///
int SizeOfSequence=0;
char* Sequence=new char[10000];
int* Suffixes_Seq=new int[10000];
int* Order_Seq=new int[10000];
int* NewOrder_Seq=new int[10000];
int CurLength_Seq=0;

/// Global Variable for whole Genome ///
char* Genome= new char[6000000];
int* Suffixes_G=new int[6000000];
int* Order_Genom=new int[6000000];
int* NewOrder_Genom=new int[6000000];
int CurLength_Genom=0;

/// Global Variable for test ///
char* Sequence_test=new char[11];
int* Suffixes_test=new int[11];


char* save_Genome()
{
    char* buf=new char[6000000];
    FILE* file_Genom=fopen(DIR_DNA "genome.fasta","r");
    int num_reads=0;
    int start=0;
    bool flag=true;
    while(!feof(file_Genom)) // While the file does not empty
    {
        fscanf(file_Genom, "%[^\n\r] ", buf); // Read the first line in the buffer
        if(num_reads==0) // ignore first line ID of Genome
        {
            num_reads++;
            buf[0]=0;
            continue;
        }
        if(buf[0]=='>') // the end of first Genome where the start of second Genome start with '>'
        {
            break;
        }
        for(int i=0; i<strlen(buf); i++)
        {
            Genome[start]=buf[i]; // Append into Genome Array
            start++;
        }
    }
    Genome[start]='$'; // Append $ at the end of Genome
    delete [] buf;
    return Genome;
}
char* save_Sequence_TienK_Char()
{
    char* buf=new char [10000];
    FILE* file=fopen(DIR_DNA "genome.fasta", "r");
    int num_reads=0;
    int start=0;
    int Max_Length=0;
    int Last_iteration_Length=0;
    while(!feof(file)) // Read the whole Sequence from the file
    {
        if(num_reads==0) // Ignore this line (ID)
        {
            fscanf(file, "%[^\n\r] ", buf);
            num_reads++;
            buf[0]=0;
            continue;
        }
        else if(Max_Length<9999) // Read first 9999 characters from the file and 1 for $
        {
            fscanf(file, "%[^\n\r] ", buf);
            Max_Length+=strlen(buf); //Calculate Size of all buffers

            /* if last length of Sequence was 9990 and last buffer size = 20 so final size will be 10010
               to avoid this condition we well check if Max Length has exceeds max size 10000 and subtract it from last size 10000
               in this case last iteration will loop only for 10 instead of 20 and Max Length will be 10000 */

            if(Max_Length>=10000)
            {
                Last_iteration_Length=9999-start; // to keep space for $ at the last of sequence (9999)
                Max_Length=Last_iteration_Length+start;

                for(int i=0; i<Last_iteration_Length; i++)
                {
                    Sequence[start]=buf[i];
                    start++;
                }
                Sequence[9999]='$'; // put $ at the end of Sequence in case of sequence length 9999
                break;
            }
            for(int i=0; i<strlen(buf); i++)// Append each character of the buffer in the array of the sequence
            {
                Sequence[start]=buf[i];
                start++;
            }
            start=strlen(Sequence);
        }
        else // Exit after 10000 characters
            break;
    }
    if(Max_Length!=9999) // if Length of Sequence less than  9999 we have to put $ at the end of this Sequence
        Sequence[start]='$';

    fclose(file);
    fflush(NULL);
    delete [] buf; // delete the buffer array
    return Sequence;
}
int *Build_Suffixes(char* Sequence,int length) // Build array of suffixes that contain index of each suffix
{
    int * Suffixes=new int [length];
    for(int i=0; i<length; i++)
    {
        Suffixes[i]=i;
    }
    return Suffixes;
}
void Print_Suffixes(int * Suffixes,char* Sequence,int length) // Print each suffix of string using the index of the suffix only
{
    for(int i=0; i<length; i++)
    {
        cout<<"Suffix  ID = "<<Suffixes[i]<<endl;
        cout<<"Sequence = ";
        for(int j=Suffixes[i]; j<length; j++)
        {
            cout<<Sequence[j];
        }
        cout<<endl;
    }
}
bool comp_nSeqLog(int a,int b)
{
    int Length_comp=strlen(Sequence);
    while(Sequence[a]==Sequence[b]) // Compare the two char until find the smallest or when the length of one of them ends
    {
        if(b>=Length_comp && a>=Length_comp)
            break;
        ++a;
        ++b;

    }
    return Sequence[a]<Sequence[b];
}
void Sort_nSeqLog(int* Suffixes,char * Sequence)
{
    auto start = high_resolution_clock::now(); // Start of Sorting
    int Length=strlen(Sequence);
    sort(Suffixes,Suffixes+Length,comp_nSeqLog);
    auto stop = high_resolution_clock::now(); // end of Sorting
    auto duration = duration_cast<microseconds>(stop - start); // Duration time needed t evaluate sort function
    Time=duration.count();
}
bool comp_prefix_Seq(int a,int b) //Compare the letters according to their ASSCI code
{
    return Order_Seq[a] < Order_Seq[b] || Order_Seq[a]==Order_Seq[b] && Order_Seq[a+CurLength_Seq] <Order_Seq[b+CurLength_Seq];
}
void Prefix_Doubling_Seq(char* Sequence,int *Suffixes)
{
    auto start = high_resolution_clock::now(); //  Start of Sorting
    int Length=strlen(Sequence);
    for(int i=0; i<Length; i++)
    {
        Order_Seq[i]=Sequence[i]; // Save the ASSCI of the letter
    }
    memset(NewOrder_Seq,0,10000); // Remove old and make new version of NewOrder
    for(int len=1; ; len*=2) // Determine how many letters will Compare (2^len)
    {
        CurLength_Seq=len;
        sort(Suffixes,Suffixes+Length,comp_prefix_Seq);
        for(int i=1; i<Length; i++)
        {
            NewOrder_Seq[i]=NewOrder_Seq[i-1]+comp_prefix_Seq(Suffixes[i-1],Suffixes[i]); //New Order after compare last Order
        }
        for(int i=0; i<Length; i++)
        {
            Order_Seq[Suffixes[i]]=NewOrder_Seq[i]; //Final Order
        }
        if(NewOrder_Seq[Length-1]==Length-1)
            break;
    }
    auto stop = high_resolution_clock::now(); // Stop of Sorting
    auto duration = duration_cast<microseconds>(stop - start); // Duration time taken to evaluate sort function
    Time=(duration.count());
    memset(NewOrder_Seq,0,10000);
    memset(Order_Seq,0,10000);
}
bool comp_prefix_Genom(int a,int b) //Compare the letters according to their ASSCI code
{
    return Order_Genom[a] < Order_Genom[b] || Order_Genom[a]==Order_Genom[b] && Order_Genom[a+CurLength_Genom] <Order_Genom[b+CurLength_Genom];
}
void Prefix_Doubling_Genom(char* Genome,int *Suffixes_G)
{
    auto start = high_resolution_clock::now(); //  Start of Sorting
    int Length=strlen(Genome);
    for(int i=0; i<Length; i++)
    {
        Order_Genom[i]=Genome[i]; // Save the ASSCI of the letter
    }
    memset(NewOrder_Genom,0,6000000); // Remove old and make new version of NewOrder
    for(int len=1; ; len*=2) // Determine how many letters will Compare (2^len)
    {
        CurLength_Genom=len;
        sort(Suffixes_G,Suffixes_G+Length,comp_prefix_Genom);
        for(int i=1; i<Length; i++)
        {
            NewOrder_Genom[i]=NewOrder_Genom[i-1]+comp_prefix_Genom(Suffixes_G[i-1],Suffixes_G[i]); //New Order after compare last Order
        }
        for(int i=0; i<Length; i++)
        {
            Order_Genom[Suffixes_G[i]]=NewOrder_Genom[i]; //Final Order
        }
        if(NewOrder_Genom[Length-1]==Length-1)
            break;
    }
    auto stop = high_resolution_clock::now(); // Stop of Sorting
    auto duration = duration_cast<seconds>(stop - start); // Duration time taken to evaluate sort function
    Time=(duration.count());
}
void Print_TienKChar_Suffixes_nSeqLogn(int *Suffixes)
{
    int Length=strlen(Sequence);
    FILE *in_file  = fopen("SortedSuffixes_1.txt", "w"); // Open new file
    fprintf(in_file,"Sorted Suffixes with Complexity = O(n^2log(n))\n"); // print first line in the file
    fprintf(in_file,"Time taken by Sort Function in Mircroseconds = %d\n",Time); // print second line in the file
    for(int i=0; i<Length; i++)
    {
        fprintf(in_file,"Suffix ID: %d\n",Suffixes[i]); // print the suffixes ID's
    }
    fclose(in_file);
}
void Print_TienKChar_Suffixes_nLogSeqn(int* Suffixes)
{
    int Length=strlen(Sequence);
    FILE *in_file  = fopen("SortedSuffixes_2.txt", "w"); // Open new file
    fprintf(in_file,"Sorted Suffixes with Complexity = O(nlog^2(n))\n"); // print first line in the file
    fprintf(in_file,"Time taken by Sort Function in mircoseconds = %d\n",Time); // print second line in the file
    for(int i=0; i<Length; i++)
    {
        fprintf(in_file,"Suffix ID: %d\n",Suffixes[i]);
    }
    fclose(in_file);

}
void Print_Genom_File(int* Suffixes_G)
{
    int Length=strlen(Genome);
    FILE *in_file  = fopen("Sorted_Suffix_Genom.txt", "w"); // Open new file
    fprintf(in_file,"Sorted Suffixes with Complexity = O(nlog^2(n))\n"); // print first line in the file
    fprintf(in_file,"Time taken by Sort Function in seconds = %d\n",Time); // print second line in the file
    for(int i=0; i<Length; i++)
    {
        fprintf(in_file,"Suffix ID: %d\n",Suffixes_G[i]);
    }
    fclose(in_file);
}
void Delete_Sequence(char* Sequence,int* Sufiixes)
{
    delete [] Sufiixes;
    delete [] Sequence;
}
char* TestCase(int Line)
{
    char *buf_Test=new char[90];   // temporary array to save lines
    FILE* file=fopen(DIR_DNA "genome.fasta", "r");
    int Start=0;  // variable to determine number of  line that will be ignored
    int End=0;  // variable to determine the last index of Sequence
    int flag=0;   // boolean variable to stop reading from file
    while(!feof(file)&& flag!=1)
    {

        while(Start<Line) // Determine the number of lines that will be ignored using start
        {
            fscanf(file, "%[^\n\r] ", buf_Test);
            buf_Test[0]=0;
            Start++;
            continue;
        }
        fscanf(file, "%[^\n\r] ", buf_Test);
        for(int i=0; i<strlen(buf_Test); i++)
        {
            if(i==10)
            {
                flag=1;
                break;
            }
            Sequence_test[End]=buf_Test[i];
            End++;
        }
    }
    Sequence_test[End]='$';       // Add $ to the end of Sequence
    fclose(file);
    fflush(NULL);
    delete[]buf_Test;
    return Sequence_test;
}
bool comp_nSeqLog_test(int a,int b)
{
    int Length_comp=strlen(Sequence_test);
    while(Sequence_test[a]==Sequence_test[b]) // Compare the two char until find the smallest or when the length of one of them ends
    {
        if(b>=Length_comp && a>=Length_comp)
            break;
        ++a;
        ++b;

    }
    return Sequence_test[a]<Sequence_test[b];
}
void Sort_nSeqLog_test(int* Suffixes,char * Sequence)
{
    auto start = high_resolution_clock::now(); // Start of Sorting
    int Length=strlen(Sequence);
    sort(Suffixes,Suffixes+Length,comp_nSeqLog_test);
    auto stop = high_resolution_clock::now(); // end of Sorting
    auto duration = duration_cast<microseconds>(stop - start); // Duration time needed t evaluate sort function
    Time=duration.count();
}
int main()
{
    ///    Write First 10,000 char  Sorted with complexity n^2Log(n) on File Name :SortedSuffixes_1.txt ///

    // Fill the array of Sequence
    Sequence=save_Sequence_TienK_Char();

    // Length of The Sequence
    int Len_Sequence=strlen(Sequence);

    // Initiate Suffixes
    Suffixes_Seq=Build_Suffixes(Sequence,Len_Sequence);

    //Sort Suffixes (n^2log(n))
    Sort_nSeqLog(Suffixes_Seq,Sequence);

    //Write on file
    Print_TienKChar_Suffixes_nSeqLogn(Suffixes_Seq);


    ///    Write First 10,000 char  Sorted with complexity nLog^2(n) on File Name :SortedSuffixes_2.txt ///


    //Sort Suffixes (nlog^2(n))
    Prefix_Doubling_Seq(Sequence,Suffixes_Seq);

    //Write on file
    Print_TienKChar_Suffixes_nLogSeqn(Suffixes_Seq);

    // Delete the Sequence & Suffixes
    Delete_Sequence(Sequence,Suffixes_Seq);


    ///  Write Whole Genome Suffixes sorted with Complexity nLog^2(n) on File Name : Sorted_Suffix_Genome.txt  ///

    // Fill the array of Sequence
    Genome=save_Genome();

    // Length of The Sequence
    int Len_Genome=strlen(Genome);

    // Initiate Suffixes
    Suffixes_G=Build_Suffixes(Genome,Len_Genome);

    //Sort Suffixes (nlog^2(n))
    Prefix_Doubling_Genom(Genome,Suffixes_G);

    //Write on file
    Print_Genom_File(Suffixes_G);

    //Delete Suffixes and Genome
    Delete_Sequence(Genome,Suffixes_G);


    /// 30 Test cases each test case contain at least 10 char ///


    for ( int i=1; i<=30; i++)
    {
        if(i%2==0) // Line sort using Prefix sort
        {

            Sequence_test=TestCase(i);
            cout<<"Test Case : "<<i<<endl<<endl;;
            cout<<"Sequence : "<<Sequence_test<<endl<<endl;
            Suffixes_test=Build_Suffixes(Sequence_test,strlen(Sequence_test));
            cout<<"Suffixes before sort :"<<endl<<endl;;
            Print_Suffixes(Suffixes_test,Sequence_test,strlen(Sequence_test));
            cout<<endl;
            Prefix_Doubling_Seq(Sequence_test,Suffixes_test);
            cout<<"Suffixes after sort using prefix doubling :"<<endl<<endl;;
            Print_Suffixes(Suffixes_test,Sequence_test,strlen(Sequence_test));
            cout<<endl<<"*************************"<<endl;
        }
        else // Line sort using naive sort
        {
            Sequence_test=TestCase(i);
            cout<<"Test Case : "<<i<<endl<<endl;;
            cout<<"Sequence : "<<Sequence_test<<endl<<endl;;
            Suffixes_test=Build_Suffixes(Sequence_test,strlen(Sequence_test));
            cout<<"Suffixes before sort :"<<endl<<endl;
            Print_Suffixes(Suffixes_test,Sequence_test,strlen(Sequence_test));
            cout<<endl;
            Sort_nSeqLog_test(Suffixes_test,Sequence_test);
            cout<<endl<<"Suffixes after sort using naive sort :"<<endl<<endl;;
            Print_Suffixes(Suffixes_test,Sequence_test,strlen(Sequence_test));
            cout<<endl<<"*************************"<<endl;

        }
    }
    Delete_Sequence(Sequence_test,Suffixes_test); // Delete Sequence and Suffixes

    return 0;
}
