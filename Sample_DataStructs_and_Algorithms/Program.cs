using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HelloWorld
{

    class Program
    {
        static void Main(string[] args)
        {
            
            Console.WriteLine("Starting Binary tree creation. Binary trees are data structures in which each node has at" +
                " most a reference to two other nodes, which are referred to as the left child and the right child. as a" +
                " means of accessing nodes based on some value or label associated with each node. Binary trees labelled" +
                " this way are used to implement binary search trees and binary heaps, and are used for efficient" +
                " searching and sorting.");

            Console.WriteLine();
            Console.WriteLine("Enter the number of elements that you would like the Binary tree to have");
            int treeSize = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Enter a series of {0} integers, in any order pressing enter after each entry:",treeSize);
            int initValue = Convert.ToInt32(Console.ReadLine());
            BinaryTree<int> binTree = new BinaryTree<int>(initValue);
            for (int i = 0; i < treeSize - 1; i++)
            {
                binTree.Add(Convert.ToInt32(Console.ReadLine()));
            }
                        
            Console.WriteLine("Enter a node value to check the values of its predecessor and successor:");
            BinaryTree<int> Node = new BinaryTree<int>(Convert.ToInt32(Console.ReadLine()));
            BinaryTree<int> pred = binTree.Predecessor(Node);
            Console.WriteLine("Value of predecessor to {0} is: {1}", Node.Value, pred.Value);
            BinaryTree<int> succ= binTree.Successor(Node);
            Console.WriteLine("Value of successor to {0} is: {1}", Node.Value, succ.Value );


            Console.WriteLine();
            BinaryTree<int> MaxNode = binTree.Max();
            BinaryTree<int> MinNode = binTree.Min();
            int max = MaxNode.Value;
            int min = MinNode.Value;
            Console.WriteLine("Min value is {0}, Max value is {1}", min, max);

            Console.WriteLine("Please enter two numbers to check if they are in the Binary Tree");
            int checkTree2 = Convert.ToInt32(Console.ReadLine());
            bool IsIn2 = binTree.Search(checkTree2);
            int checkTree1 = Convert.ToInt32(Console.ReadLine());
            bool IsIn1 = binTree.Search(checkTree1);
            Console.WriteLine("{0} is in Binary Tree: {1}",checkTree1,IsIn1);
            Console.WriteLine("{0} is in Binary Tree: {1}", checkTree2, IsIn2);

            Console.WriteLine();
            Console.WriteLine("Press Enter to move on to Linked Lists");
            var moveOn2 = Console.ReadLine();

            Console.WriteLine();
            Console.WriteLine("Starting Linked Lists");
            Console.WriteLine("a Linked List is a linear collection of data elements whose order is not given by" +
                " their physical placement in memory. Instead, each element points to the next. It is a data structure" +
                " consisting of a collection of nodes which together represent a sequence. They are dynamic data structures" +
                " with efficient memory utilization. They can be used to implement stacks, queues, and other abstract data types. ");
            Console.WriteLine();
            Console.WriteLine("Enter the number of elements that you would like the Linked List to have");
            int listSize = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Enter a series of {0} numbers in any order, pressing enter after each entry:", listSize);

            int initListValue = Convert.ToInt32(Console.ReadLine());
            LinkedList<int> aList = LinkedList<int>.CreateWithInitialValue(initListValue);
            for (int i = 0; i < listSize - 1; i++)
            {
                aList.Add(Convert.ToInt32(Console.ReadLine()));
            }

            Console.WriteLine();

            Console.WriteLine("Checking if certain values are in int Linkedlist");
            Console.WriteLine("Please enter a number to check if in list");
            int check = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Linkedlist contains {0} is {1}: ", check ,aList.Contains(check));

            Console.WriteLine("Please enter a number to delete from list");
            int delFromList = Convert.ToInt32(Console.ReadLine());

            Console.WriteLine("Deleting {0}", delFromList);
            aList.TryDelete(delFromList);
            Console.WriteLine("Linkedlist contains {0} is {1}: ", delFromList, aList.Contains(check));
            Console.WriteLine();

            Console.WriteLine();
            Console.WriteLine("Creating enumeration method for Linked list and iterating over each number");
            foreach(var listItem in aList)
            {
                Console.Write("{0}, ", listItem);
            }

            Console.WriteLine();
            Console.WriteLine("Press Enter key to move on to Stacks");
            var moveOn3 = Console.ReadLine();
            

            Console.WriteLine();
            Console.WriteLine("Starting stack creation");
            Console.WriteLine();
            Console.WriteLine("Enter the number of elements that you would like the Stack to have");
            int stackSize = Convert.ToInt32(Console.ReadLine());

            Console.WriteLine("Enter the initial number at the base of the Stack");
            int initStackValue = Convert.ToInt32(Console.ReadLine()); Stack<int> intStack = new Stack<int>(stackSize);
            Console.WriteLine("Starting with an empty stack");
                        for (int i = 0; i < stackSize; i++)
                        {
                            if (i==0)
                            {
                                Console.WriteLine("Pointer Before Pushing: {0}", intStack.GetPointer());
                                Console.WriteLine("Starting Stack with {0}", i + initStackValue);
                                intStack.Push(i + initStackValue);
                            }
                            else
                            {
                                Console.WriteLine("Pointer Before Pushing: {0}", intStack.GetPointer());
                                Console.WriteLine("Adding {0} to top of stack", i + initStackValue);
                                intStack.Push(i + initStackValue);
                            }
                                  
                        }

                        Console.WriteLine();

                        int popped;
                        for (int i = 0; i < stackSize; i++)
                        {
                            popped = intStack.Pop();
                            Console.WriteLine("Popped is: {0}", popped);
                            Console.WriteLine("Pointer After Popping: {0}", intStack.GetPointer());
                        }
            Console.WriteLine();
            Console.WriteLine("Press Enter key to move on to Sets");
            var moveOn4 = Console.ReadLine();
            Console.WriteLine();
    
            Console.WriteLine("Starting set manipulation");
            Console.WriteLine();

            string[] tools = { "knife", "pliers", "saw", "hammer" };
            string[] weapons = { "hammer", "sword", "axe", "knife" };
            Console.WriteLine("Our tools in the toolbox are:");
            tools.ToList().ForEach(x => Console.Write("{0}, ",x));
            Console.WriteLine();
            Console.WriteLine("Our weapons in the weapons chest are:");
            weapons.ToList().ForEach(x => Console.Write("{0}, ", x));
            Console.WriteLine();

            Console.WriteLine();
            Console.WriteLine("Union of sets:");
            string[] union = Union(tools, weapons);
            union.ToList().ForEach(x => Console.Write("{0}, ",x));
            Console.WriteLine();

            Console.WriteLine();
            string[] commonElements = getCommonElements(tools, weapons);
            Console.WriteLine("The common elements in these sets are: ");
            List<string> commonElementsList = commonElements.ToList();
            commonElementsList.ForEach(x => Console.Write(x+","));
            Console.WriteLine();

            Console.WriteLine();
            int numberCommonElements = getNumberCommonElements(tools, weapons);
            Console.WriteLine("Number of common elements: {0}", numberCommonElements);

            Console.WriteLine();
            Console.WriteLine("Press Enter key to move on to 'Make DNA Compliment'");
            var moveOn5 = Console.ReadLine();
            Console.WriteLine();

            string DNA = "AGTCCGATATCGCTGGATCAT";
            Console.WriteLine("Our DNA string of nucleoases is: {0}",DNA);
            Console.WriteLine("It requires a compliemnt string to form a full DNA helix");
            string DNA_compliment = Make_Compliment(DNA);
            Console.WriteLine("Our compliment function gives: {0}", DNA_compliment);
            Console.WriteLine("Which is correct!");

            /*
            Console.WriteLine();
            Console.WriteLine("Press Enter key to move on to '1 odd number out'");
            var moveOn6 = Console.ReadLine();
            Console.WriteLine();
            string numString = "2 15 4 6 8";
            int oddNumberOut = Test(numString);
            Console.WriteLine("Position of odd number out is: {0}", oddNumberOut);

            int fizzbuzz = Solution(10);
            Console.WriteLine("FizzBuzz is: {0}", fizzbuzz);

            string XOString = "ooxxo";
            bool XOequal = XO(XOString);

            Console.WriteLine(XOequal);

            Console.WriteLine("Starting WUBWUB Encoding. Write out some song lyrics ");
            string lyrics = Convert.ToString(Console.ReadLine());
            string wubString = WUBEncoder(lyrics);
            Console.WriteLine("{0}", wubString);

            int number = 919;
            int squaredDigits=SquareDigits(number);
            Console.WriteLine(squaredDigits);
            */


        }
        public static string WUBEncoder(string input)
        {
            Random rnd = new Random();
            string[] sepString = input.Split(' ');
            int strLen = sepString.Length;
            string[] wubString = new string[2 * strLen + 1];
            for (int i = 0; i < strLen; i++)
            {
                wubString[2 * i + 1] = sepString[i];
                int rand = rnd.Next(1, 3);
                string wubs = String.Concat(Enumerable.Repeat("WUB", rand));
                wubString[2 * i] = wubs;
            }
            wubString[2 * strLen] = "WUB";
            return string.Join("", wubString);
        }
        static bool Contains(string[] array, string key)
        {
            foreach (string item in array)
            {
                if (item == key)
                {
                    return true;
                }
            }
            return false;
        }
        static int getNumberCommonElements(string[] A, string[] B)
        {
            int commonCounter = 0;
            foreach (string item in A)
            {
                if (Contains(B, item))
                {
                    commonCounter++;
                }
            }
            return commonCounter;
        }
        static string[] Union(string[] A, string[] B)
        {
            int UnionSize = A.Length + B.Length - getNumberCommonElements(A, B);
            string[] Union = new string[UnionSize];
            int index = 0;
            foreach (string item in A)
            {
                Union[index++] = item;
            }
            foreach (string item in B)
            {
                if (!Contains(A, item))
                {
                    Union[index++] = item;
                }
            }
            return Union;
        }
        static string[] getCommonElements(string[] A, string[] B)
        {
            string[] commonElements = new string[getNumberCommonElements(A, B)];
            int gotCounter = 0;
            for (int i = 0; i < A.Length; i++)
            {
                if (Contains(B, A[i]))
                {
                    commonElements[gotCounter] = A[i];
                    gotCounter++;
                }
            }
            return commonElements;
        }
        public static string Make_Compliment(string DNA)
        {
            char[] DNAarray = DNA.ToCharArray();
            int arrLen = DNAarray.Length;
            Dictionary<char,char> DNA_Compliment = new Dictionary<char, char>();
            DNA_Compliment.Add('A','C');
            DNA_Compliment.Add('C','A');
            DNA_Compliment.Add('T', 'G');
            DNA_Compliment.Add('G', 'T');

            char[] complimentArray = new char[arrLen];
            for(int i = 0; i<arrLen; i++)
            {
                complimentArray[i] = DNA_Compliment[DNAarray[i]];
            }
            string DNACompString = new string(complimentArray);
            return DNACompString;
        }

        public static int SquareDigits(int n)
        {

            int numLen = (int)(Math.Floor(Math.Log10(n))) + 1;
            int[] numList = new int[numLen];
            int changeNum = n;
            for(int i =0;i<numLen;i++)
            {
                double num = Math.Floor(changeNum / Math.Pow(10, numLen - i - 1));
                numList[i] = (int)Math.Pow(num,2);
                changeNum = changeNum-(int)(num * Math.Pow(10, numLen - i - 1));
            }
            string numListString= string.Join("", numList);
            return Int32.Parse(numListString);
        }
        public static bool XO(string input)
        {
            int xNum = 0;
            int oNum = 0;
            char[] charInput = input.ToCharArray();

            foreach(char c in input)
            {
                if (Char.ToLower(c) == 'x')
                {
                    xNum++;
                }
                else if(Char.ToLower(c) == 'o')
                {
                    oNum++;
                }
                else
                {
                    ;
                }
            }
            return xNum == oNum;
        }
        public static int Solution(int value)
        {
            int sum = 0;
            if (value <= 0)
            {
                return 0;
            }
            else
            {
                for (int i = 1; i < value; i++)
                {
                    if (i%3==0 | i%5==0)
                    {
                        sum += i;
                    }
                }
            }
            return sum;
        }

        public static int Test(string numbers)
        {
            string[] numString = numbers.Split(new char[] {' '},
            StringSplitOptions.RemoveEmptyEntries);
            int[] numArray = Array.ConvertAll(numString, c =>int.Parse(c));
            int oddCount = 0;
            int evenCount = 0;
            foreach(int num in numArray)
            {

                if (num%2==0)
                {
                    evenCount++;
                }
                else
                {
                    oddCount++;
                }
            }
            bool outlierOdd = evenCount > oddCount;
            int PosOddNumOut = 0;
            for(int num = 0; num<numArray.Length;num++)
            {
                if (numArray[num]%2==1 == outlierOdd)
                {
                    PosOddNumOut = num+1;
                }
            }
            return PosOddNumOut;
        }
    }
}


