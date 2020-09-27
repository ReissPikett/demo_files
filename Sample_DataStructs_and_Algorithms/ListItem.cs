using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HelloWorld
{
    class ListItem<T>
    {
        public T Value { get; set; }
        public ListItem<T> Next { get; set; }

    }
}
