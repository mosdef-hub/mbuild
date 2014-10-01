# import os
# import glob
import Tkinter
import ttk
from compound import *

class TreeView(object):

    def __init__(self, compound, **kwargs):
        self.compound = compound
        self.nodemap = {"root": compound}
        self.kwargs = kwargs

    def populate_tree(self, tree, node):
        if tree.set(node, "type") != 'compound':
            return

        # path = tree.set(node, "fullpath")
        # tree.delete(*tree.get_children(node))

        parent = tree.parent(node)

        compound = self.nodemap[node]


        for child in tree.get_children(node):
            # print tree.item(child)
            if tree.item(child)['text'] == 'dummy':
                tree.delete(child)

        for v in compound.parts:
            ptype = None
            if isinstance(v,Compound):
                ptype = "compound"
            if isinstance(v,Atom):
                ptype = "atom"
            if isinstance(v,Bond):
                ptype = "bond"

            if not v in self.nodemap.values():
                # find a reference to the object
                k = 'na'
                for rk,rv in compound.labels.iteritems():
                    if v is rv:
                        k = rk
                        break

                id = tree.insert(node, "end", text=k, values=[k, ptype, v.__str__(), str(v.__class__)])
                self.nodemap[id] = v

                if ptype == 'compound':
                    tree.insert(id, "end", text="dummy")


    def populate_roots(self, tree):
        k = ""
        ptype = "compound"
        v = self.compound
        node = tree.insert("", "end", text="root", values=[k, ptype, v.kind, str(v.__class__)])

        # node = tree.insert('', 'end', text=str(self.compound.__class__), values=[str(self.compound.__class__), "compound"])
        self.nodemap[node] = self.compound
        self.populate_tree(tree, node)

    @staticmethod
    def update_tree(event):
        tree = event.widget
        tree.compoundTreeView.populate_tree(tree, tree.focus())

    @staticmethod
    def change_dir(event):
        tree = event.widget
        node = tree.focus()
        # if tree.parent(node):
        compound = tree.compoundTreeView.nodemap[node]

    @staticmethod
    def autoscroll(sbar, first, last):
        """Hide and show scrollbar as needed."""
        first, last = float(first), float(last)
        if first <= 0 and last >= 1:
            sbar.grid_remove()
        else:
            sbar.grid()
        sbar.set(first, last)

    def show(self):
        root = Tkinter.Tk()

        vsb = ttk.Scrollbar(orient="vertical")
        hsb = ttk.Scrollbar(orient="horizontal")

        tree = ttk.Treeview(columns=("label", "type", "kind", "class"),
            displaycolumns=("label", "kind", "class"), yscrollcommand=lambda f, l: self.autoscroll(vsb, f, l),
            xscrollcommand=lambda f, l: self.autoscroll(hsb, f, l))

        tree.compoundTreeView = self


        vsb['command'] = tree.yview
        hsb['command'] = tree.xview

        tree.heading("#0", text="Component hierarchy", anchor='w')
        tree.heading("label", text="Label", anchor='w')
        tree.heading("type", text="Type", anchor='w')
        tree.heading("kind", text="Kind", anchor='w')
        tree.heading("class", text="Class", anchor='w')
        tree.column("class", stretch=0, width=500)

        tree.bind('<<TreeviewOpen>>', self.update_tree)
        tree.bind('<Double-Button-1>', self.change_dir)


        self.populate_roots(tree)

        # Arrange the tree and its scrollbars in the toplevel
        tree.grid(column=0, row=0, sticky='nswe')
        vsb.grid(column=1, row=0, sticky='ns')
        hsb.grid(column=0, row=1, sticky='ew')
        root.grid_columnconfigure(0, weight=1)
        root.grid_rowconfigure(0, weight=1)

        root.mainloop()

if __name__ == "__main__":
    from ethane import Ethane
    ethane = Ethane()
    TreeView(ethane).show()