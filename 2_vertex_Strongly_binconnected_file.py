import io
import sys
sys.setrecursionlimit(10**6)

class Graph:
    def __init__(self, num_vertices):
        self.num_vertices = num_vertices
        self.adj_list = [[] for _ in range(num_vertices)]

    def add_edge(self, u, v):
        self.adj_list[u].append(v)

    def cheriyan_mehlhorn_gabow(self):
        index = 0
        stack = []
        low_links = [float('inf')] * self.num_vertices
        indices = [float('inf')] * self.num_vertices
        sccs = []

        def strong_connect(v):
            nonlocal index
            indices[v] = index
            low_links[v] = index
            index += 1
            stack.append(v)

            for w in self.adj_list[v]:
                if indices[w] == float('inf'):
                    strong_connect(w)
                    low_links[v] = min(low_links[v], low_links[w])
                elif w in stack:
                    low_links[v] = min(low_links[v], indices[w])

            if low_links[v] == indices[v]:
                scc = []
                while True:
                    w = stack.pop()
                    scc.append(w)
                    if w == v:
                        break
                sccs.append(scc)

        for v in range(self.num_vertices):
            if indices[v] == float('inf'):
                strong_connect(v)

        return sccs

    def jens_schmidt(self):
        articulation_points = set()
        low = [float('inf')] * self.num_vertices
        disc = [float('inf')] * self.num_vertices
        time = 0
        parent = [-1] * self.num_vertices

        def dfs(v):
            nonlocal time
            low[v] = time #lowest discovery time
            disc[v] = time  # discovery time
            time += 1
            children = 0  # (neighbors) 

            for w in self.adj_list[v]:
                if disc[w] == float('inf'):
                    parent[w] = v
                    children += 1
                    dfs(w)
                    low[v] = min(low[v], low[w])
                    if (parent[v] == -1 and children > 1) or (parent[v]!= -1 and low[w] >= disc[v]):
                        articulation_points.add(v)
                elif w!= parent[v]:
                    low[v] = min(low[v], disc[w])

        for v in range(self.num_vertices):
            if disc[v] == float('inf'):
                dfs(v)

        return len(articulation_points) == 0

    def is_2_vertex_strongly_biconnected(self):
        low = [float('inf')] * self.num_vertices
        disc = [float('inf')] * self.num_vertices
        time = 0
        parent = [-1] * self.num_vertices
        articulation_points = set()

        def dfs(v):
            nonlocal time
            low[v] = time  #lowest discovery time
            disc[v] = time #discovery time
            time += 1
            children = 0

            for w in self.adj_list[v]:
                if disc[w] == float('inf'):
                    parent[w] = v
                    children += 1
                    dfs(w)
                    low[v] = min(low[v], low[w])
                    if (parent[v] == -1 and children > 1) or (parent[v]!= -1 and low[w] >= disc[v]):
                        articulation_points.add(v)
                elif w!= parent[v]:
                    low[v] = min(low[v], disc[w])

        for v in range(self.num_vertices):
            if disc[v] == float('inf'):
                dfs(v)

        for i in range(self.num_vertices - 1):
            for j in range(i + 1, self.num_vertices):
                if i!= j and (i in articulation_points or j in articulation_points):
                    g_copy = self.copy_graph()
                    g_copy.remove_vertex(i)
                    g_copy.remove_vertex(j)
                    if g_copy.num_vertices > 0 and not g_copy.is_connected():
                        return False
        return True
    def construct_undirected_graph(self, sccs):
        undirected_graph = Graph(len(sccs))
        for i in range(len(sccs)):
            for j in range(i + 1, len(sccs)):
                for v in sccs[i]:
                    for w in sccs[j]:
                        if w in self.adj_list[v] or v in self.adj_list[w]:
                            undirected_graph.add_edge(i, j)
                            break
        return undirected_graph

    def copy_graph(self):
        g_copy = Graph(self.num_vertices)
        for i in range(self.num_vertices):
            for j in self.adj_list[i]:
                g_copy.add_edge(i, j)
        return g_copy

    def remove_vertex(self, v):
        if v < len(self.adj_list):
            del self.adj_list[v]
        for i in range(len(self.adj_list)):
            self.adj_list[i] = [j - 1 if j > v else j for j in self.adj_list[i] if j!= v]
        self.num_vertices -= 1
        if self.num_vertices == 0:
            self.adj_list = []

    def is_connected(self):
        if not self.adj_list:
            return True
        visited = [False] * len(self.adj_list)
        self.dfs(0, visited)
        return all(visited)

    def dfs(self, v, visited):
        if v < len(visited):
            visited[v] = True
            for w in self.adj_list[v]:
                if w < len(visited) and not visited[w]:
                    self.dfs(w, visited)



with open('gemsec_dezzer(HR).txt', 'r') as f: #we must to replace the file name with "dataset_file"
    num_vertices = int(f.readline().strip())
    g = Graph(num_vertices)

    num_edges = int(f.readline().strip())

    for _ in range(num_edges):
        line = f.readline().strip()
        u, v = line.split()
        g.add_edge(int(u), int(v))

    if g.is_2_vertex_strongly_biconnected():
        print("The graph is 2-vertex strongly biconnected.")
    else:
        print("The graph is not 2-vertex strongly biconnected.")