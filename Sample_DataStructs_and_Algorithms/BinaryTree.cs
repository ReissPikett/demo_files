using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HelloWorld
{
    class BinaryTree<T> : IEnumerable<T> where T : IComparable
    {
        // properties
        public T Value;
        public BinaryTree<T> Right { get; set; }
        public BinaryTree<T> Left { get; set; }

        // constructors
        public BinaryTree(T initValue)
        {
            Value = initValue;
        }

        // methods

        public BinaryTree<T> Predecessor(BinaryTree<T> tree)
        {
            if (tree.Left != null)
            {
                return tree.Left.Max();
            }
            else
            {
                
                BinaryTree<T> hub = tree;
                BinaryTree<T> current = this;

                if (tree.Value.Equals(current.Min().Value))
                {
                    Console.WriteLine("{0} does not have a predecessor",tree.Value);
                    return tree;
                }

                List <KeyValuePair<BinaryTree<T>, bool>> nodesRight = new List<KeyValuePair<BinaryTree<T>, bool>>();

                while (!current.Value.Equals(tree.Value))
                {
                    if (current.Value.CompareTo(tree.Value) > 0)
                    {
                        nodesRight.Add(new KeyValuePair<BinaryTree<T>, bool>(current, false));
                        current = current.Left;
                    }
                    else if (current.Value.CompareTo(tree.Value) < 0)
                    {
                        nodesRight.Add(new KeyValuePair<BinaryTree<T>, bool>(current, true));
                        current = current.Right;
                    }
                }
                if (current.Left!=null)
                {
                    return current.Left;
                }

                var nodesArray = nodesRight.ToArray();

                for (int i = nodesArray.Length - 1; i >= 0; i--)
                {
                    if (nodesArray[i].Value == true)
                    {
                        hub = nodesArray[i].Key;
                        return hub;
                    }
                }
                return hub;
            }
        }
        public BinaryTree<T> Successor(BinaryTree<T> tree)
        {
            if(tree.Right!=null)
            {
                return tree.Right.Min();
            }
            else
            {
                BinaryTree<T> hub = tree;
                BinaryTree<T> current = this;
                if (tree.Value.Equals(current.Max().Value))
                {
                    Console.WriteLine("{0} does not have a Successor", tree.Value);
                    return tree;
                }

                List<KeyValuePair<BinaryTree<T>, bool>> nodesLeft = new List<KeyValuePair<BinaryTree<T>, bool>>();

                while (!current.Value.Equals(tree.Value))
                {
                    if (current.Value.CompareTo(tree.Value) > 0)
                    {
                        nodesLeft.Add(new KeyValuePair<BinaryTree<T>, bool>(current, true));
                        current = current.Left;

                    }
                    else if (current.Value.CompareTo(tree.Value) < 0)
                    {
                        nodesLeft.Add(new KeyValuePair<BinaryTree<T>, bool>(current, false));
                        current = current.Right;
                    }
                }
                if (current.Right != null)
                {
                    return current.Right;
                }
                var nodesArray = nodesLeft.ToArray();
                for (int i = nodesArray.Length - 1; i >= 0; i--)
                {
                    if (nodesArray[i].Value == true)
                    {
                        hub = nodesArray[i].Key;
                        return hub;
                    }
                }
                return hub;
            }
        }

        public BinaryTree<T> MaxBranch()
        {
            return FindMaxBranch(this);
        }
        public BinaryTree<T> MinBranch()
        {
            return FindMinBranch(this);
        }
        public BinaryTree<T> FindMaxBranch(BinaryTree<T> tree)
        {
            BinaryTree<T> current = tree;
            while(current.Right!=null)
            {
                current = current.Right;
            }
            return current;
        }

        public BinaryTree<T> FindMinBranch(BinaryTree<T> tree)
        {
            BinaryTree<T> current = tree;
            while (current.Left != null)
            {
                current = current.Left;
            }
            return current;
        }
            public BinaryTree<T> Min()
        {
            return FindMinimum(this);
        }

        public BinaryTree<T> Max()
        {
            return FindMaxBranch(this);
        }

        public BinaryTree<T> FindMinimum(BinaryTree<T> tree)
        {
            if(tree.Left==null)
            {
                return tree;
            }
            else
            {
                return FindMinBranch(tree.Left);
            }
        }

        public BinaryTree<T> FindMaximum(BinaryTree<T> tree)
        {
            if (tree.Right == null)
            {
                return tree;
            }
            else
            {
                return FindMinimum(tree.Right);
            }
        }
        public bool Search(T value)
        {
            return SearchInNodes(this,value);
        }

        private bool SearchInNodes(BinaryTree<T> tree, T value)
        {
            if (tree == null)
            {
                return false;
            }
            else if (tree.Value.Equals(value))
            {
                return true;
            }
            else if(tree.Right == null && tree.Left==null)
            {
                return false;
            }
            else if(tree.Value.CompareTo(value) <0)
            {
                return SearchInNodes(tree.Right, value);
            }
            else if (tree.Value.CompareTo(value) > 0)
            {
                return SearchInNodes(tree.Left, value);
            }
            else
            {
                return false;
            }
        }

        public void Add(T newValue)
        {
            AddtoNode(this, newValue);
        }

        private void AddtoNode(BinaryTree<T> tree, T newValue)
        {
            if (tree.Value.CompareTo(newValue) < 0 )
            {
                if(tree.Right==null)
                {
                    tree.Right = new BinaryTree<T>(newValue);
                }
                else
                {
                    AddtoNode(tree.Right,newValue);
                }
            }
            if(tree.Value.CompareTo(newValue)>0)
            {
                if (tree.Left==null)
                {
                    tree.Left = new BinaryTree<T>(newValue);
                }
                else
                {
                    AddtoNode(tree.Left, newValue);
                }
            }
        }

        // Enumerator 
        public IEnumerator<T> GetEnumerator()
        {
            throw new NotImplementedException();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            throw new NotImplementedException();
        }
    }
}
