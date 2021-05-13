"""Ordered set module.

Based on oset by Raymond Hettinger, with the following license:
Copyright (c) 2009, Raymond Hettinger, and others All rights reserved.
Package structured based on the one developed to odict Copyright (c) 2010,
    BlueDynamics Alliance, Austria
 - Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 - Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 - Neither the name of the BlueDynamics Alliance nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY BlueDynamics Alliance AS IS AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL BlueDynamics Alliance BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from collections.abc import MutableSet
import itertools as it
from typing import (
    Any,
    Dict,
    Iterable,
    Iterator,
    List,
    MutableSet,
    Optional,
    Sequence,
    Set,
    TypeVar,
    Union,
    overload,
)
KEY, PREV, NEXT = range(3)
T = TypeVar("T")

class OrderedSet(MutableSet):
    """An ordered set object with additional convenience methods.

    Based on oset by Raymond Hettinger (see module for license)

    Methods
    -------
    add
    discard
    index
    list
    pop
    union
    """

    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]  # sentinel node for doubly linked list
        self.map = {}  # key --> [key, prev, next]
        self._list = None
        if iterable is not None:
            self |= iterable

    def __len__(self):
        """Return the length."""
        return len(self.map)

    def __contains__(self, key):
        """Determine whether the set contains a key."""
        return key in self.map

    def __iter__(self):
        """Return an iterable of the set."""
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        """Iterate through the set in reverse."""
        end = self.end
        curr = end[PREV]
        while curr is not end:
            yield curr[KEY]
            curr = curr[PREV]

    def __getitem__(self, value):
        """Get an item."""
        return list(self)[value]

    def __repr__(self):
        """Get the set representation."""
        if not self:
            return f"{self.__class__.__name__}"
        return f"{self.__class__.__name__}({list(self)})"

    def __eq__(self, other):
        """Determine if set is equal to other."""
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

    def __del__(self):
        """Delete the set."""
        self.clear()

    def union(self, *sets: Union[Sequence[T], Set[T]]):
        """Combine all unique items.

        Each items order is defined by its first appearance.

        Example:
        -------
            >>> oset = OrderedSet.union(OrderedSet([3, 1, 4, 1, 5]), [1, 3], [2, 0])
            >>> print(oset)
            OrderedSet([3, 1, 4, 5, 2, 0])
            >>> oset.union([8, 9])
            OrderedSet([3, 1, 4, 5, 2, 0, 8, 9])
            >>> oset | {10}
            OrderedSet([3, 1, 4, 5, 2, 0, 10])
        """
        cls = self.__class__ if isinstance(self, OrderedSet) else OrderedSet
        containers = map(list, it.chain([self], sets))
        items = it.chain.from_iterable(containers)
        return cls(items)

    def add(self, key):
        """Add the key to the set."""
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        """Delete the key from the set."""
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def pop(self, last=True):
        """Remove a value from the set and return the value."""
        if not self:
            raise KeyError("set is empty")
        key = next(reversed(self)) if last else next(iter(self))
        self.discard(key)
        return key

    def index(self, value):
        """Get the index of a value."""
        return list(self).index(value)
