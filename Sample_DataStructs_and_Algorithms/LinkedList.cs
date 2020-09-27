using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HelloWorld
{
    class LinkedList<T> : IEnumerable<T>
    {
        // properties
        public ListItem<T> First { get; set; }

        // constructors
        public LinkedList(T initialValue)
        {
            First = new ListItem<T> { Value = initialValue };
        }
        public static LinkedList<T> CreateWithInitialValue(T initialValue)
        {
            return new LinkedList<T>(initialValue);
        }

        //methods
        public void Add(T addingValue)
        {
            ListItem<T> current = First;

            while (current.Next != null)
            {
                current = current.Next;
            }
            current.Next = new ListItem<T> { Value = addingValue };
        }
        public bool Contains(T containValue)
        {
            ListItem<T> current = First;
            do
            {
                if (current.Value.Equals(containValue))
                {
                    return true;
                }
                current = current.Next;
            } while (current != null);
            return false;
        }
        public void TryDelete(T deleteValue)
        {
            if (Contains(deleteValue))
            {
                if (First.Value.Equals(deleteValue))
                {
                    First = First.Next;
                }
                else
                {
                    ListItem<T> current = First;
                    while (current.Next != null)
                    {
                        if (current.Next.Value.Equals(deleteValue))
                        {
                            current.Next = current.Next.Next;
                            break;
                        }
                        current = current.Next;
                    }
                }
            }
        }

        public IEnumerator<T> GetEnumerator()
        {
            ListItem<T> current = First;
            while (current!=null)
            {
                yield return current.Value;
                current = current.Next;
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            throw new NotImplementedException();
        }
    }
}
