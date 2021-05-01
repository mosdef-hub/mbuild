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

KEY, PREV, NEXT = range(3)


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

    # @property
    # def list(self):
    #    """Access a stored instance of the set as a list.

    #    List cannot be set. Attempting to set it will clear the list.
    #    """
    #    if self._list is None:
    #        self._list = list(self)
    #    return self._list

    # @list.setter
    # def list(self, value):
    #    self._list = None

    def __getitem__(self, value):
        """Get an item."""
        return list(self)[value]

    def add(self, key):
        """Add a value."""
        if key not in self.map:
            end = self.end
            curr = end[PREV]
            curr[NEXT] = end[PREV] = self.map[key] = [key, curr, end]

    def discard(self, key):
        """Remove a value."""
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[NEXT] = next
            next[PREV] = prev

    def __iter__(self):
        """Iterate through the set."""
        end = self.end
        curr = end[NEXT]
        while curr is not end:
            yield curr[KEY]
            curr = curr[NEXT]

    def __reversed__(self):
        """Iterate through the set in reverse."""
        end = self.end
        curr = end[PREV]
        while curr is not end:
            yield curr[KEY]
            curr = curr[PREV]

    def pop(self, last=True):
        """Remove a value from the set and return the value."""
        if not self:
            raise KeyError("set is empty")
        key = next(reversed(self)) if last else next(iter(self))
        self.discard(key)
        return key

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

    def index(self, value):
        """Get the index of a value."""
        return list(self).index(value)
