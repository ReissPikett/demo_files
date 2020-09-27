using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HelloWorld
{
    class Stack<T>
    {
        // private
        private T[] items;
        private int pointer=0;

        // public
        public Stack(int size)
        {
            items = new T[size];
        }
        public void Push(T value)
        {
            if (pointer < items.Length)
            {
                items[pointer++] = value;
            }
            else
            {
                throw new Exception("This Stack is full");
            }
        }
        public T Pop()
        {
            pointer -= 1;
            if (pointer >= 0)
            {
                return items[pointer];
            }
            else
            {
                throw new Exception("Stack is Empty");
            }
        }
        public int GetPointer()
        {
            return pointer;
        }

    }
}
