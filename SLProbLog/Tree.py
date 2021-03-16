import mpmath
import numpy

mpmath.dps = 200


class Node():
    def __init__(self):
        self._children = set([])
        self._parents = set([])
        self._name = None
        self._weight = None
        self._negated = False
        self._qchildren = set([])
        self._alias = None

    def add_child(self, n):
        self._children.add(n)
        self._qchildren.add(n)

    def children(self):
        return self._children

    def qchildren(self):
        return self._qchildren

    def add_parent(self, p):
        self._parents.add(p)

    def parents(self):
        return self._parents

    def name(self, n = None):
        if n is not None:
            self._name = n
        return self._name

    def type(self):
        return None

    def weight(self):
        return self._weight

    def negated(self, val=None):
        if val is not None:
            self._negated = val
        return self._negated

    def alias(self, al=None):
        if al is not None:
            self._alias = al
        return self._alias

    def swap_qchildren(self, orig, dest):
        if orig in self._qchildren:
            self._qchildren.remove(orig)
            self._qchildren.add(dest)

class AndNode(Node):
    def type(self):
        return "AND"

    def __str__(self):
        return self.type()

class OrNode(Node):
    def type(self):
        return "OR"

    def __str__(self):
        return self.type()

class Leaf(Node):
    def __init__(self, zero=0, one=1):
        super(Leaf, self).__init__()
        self._lambda = 1
        self._theta = zero
        self._zero = zero
        self._one = one

    def theta(self, t=None):
        if t is not None:
            self._theta = t
        return self._theta

    def llambda(self, l=None):
        if l is not None:
            self._lambda = l
        return self._lambda

    def type(self):
        return "LEAF"

    def weight(self):
        if self.llambda() == 1:
            return self.theta()
        else:
            return self._zero

    def __str__(self):
        sign = ""
        if self.negated():
            sign = "-"
        return "l: "  + str(self.llambda()) + "; t: " + sign +  str(self.theta())

    def set_neutral(self):
        self.negated(False)
        self.llambda(1)
        self.theta(1)

    def neutral(self):
        if not self.negated() and self.llambda() == 1 and self.theta() == 1:
            return True

class Query(Leaf):
    def type(self):
        return "QUERY"

class NegQuery(Leaf):
    def type(self):
        return "NEGQUERY"

class Evidence(Leaf):
    def __init__(self, zero, one):
        super().__init__(zero, one)
        self.theta(one)

    def type(self):
        return "EVIDENCE"


class NegEvidence(Leaf):
    def __init__(self, zero, one):
        super().__init__(zero, one)
        self.llambda(0)

    def type(self):
        return "NEGEVIDENCE"



class AliasNode():
    def __init__(self, ref):
        self._reference = ref

    def type(self):
        return "ALIAS"

    def reference(self):
        return self._reference


class Tree():
    def __init__(self, betas):
        self._nodes = {}
        self._leaves = []
        self._root = None
        self._evidences = []
        self._queries = []
        self._negqueries = []
        self._betas = betas
        self._cov = {}
        self._m = {}
        self._zero = 0
        self._one = 1
        self._lenorig = None
        self._debug = False

    def add_node(self, k, n):
        if k not in self._nodes:
            self._nodes[k] = n

            if n.type() == "LEAF" or n.type() == "NEGEVIDENCE":
                self._leaves.append(k)
            elif n.type() == "QUERY":
                self._leaves.append(k)
                self._queries.append(k)
            elif n.type() == "EVIDENCE":
                self._leaves.append(k)
                self._evidences.append(k)
            elif n.type() == "NEGQUERY":
                self._leaves.append(k)
                self._negqueries.append(k)

            if n.type() != "ALIAS" and not n.parents():
                self._root = k

    def node(self, k):
        return self._nodes[k]

    def cov(self, i, j, val = None):
        id = None
        if i <= j:
            id = (i,j)
        else:
            id = (j, i)

        if val is not None:
            self._cov[id] = val

        return self._cov[id]

    def m(self, i, val = None):
        if val is not None:
            self._m[i] = val
        return self._m[i]

    def create_shadow_circuit(self, k):
        idx = max(k for k in self._nodes) + 1

        if k > self._lenorig:
            negk = k - self._lenorig
        else:
            negk = k + self._lenorig

        links = []
        for p in self.node(negk).parents():
            links.append((negk, p))

        while links:
            (c, p) = links.pop()
            if self.node(c).alias() is None:
                self.add_node(idx, AliasNode(c))
                self.node(c).alias(idx)
                idx += 1

            self.node(p).swap_qchildren(c, self.node(c).alias())

            if self.node(p).parents():
                for pp in self.node(p).parents():
                    links.append((p, pp))
            elif self.node(p).alias() is None:  # this is the root
                self.add_node(idx, AliasNode(p))
                self.node(p).alias(idx)
                idx += 1

        if self._debug:
            for k in range(1, idx + 1):
                if k in self._nodes and self._nodes[k].type() != "ALIAS":
                    v = self._nodes[k]

                    add = ""
                    if v.type() != "AND" and v.type() != "OR":
                        add += "; " + str(v.llambda()) + "; " + str(self._betas[v.theta()])

                    print(str(k) + "; " + str(v.parents()) + "; " + str(v.children()) + "; " + str(v.qchildren()) + add)

    def compute_query_monte(self, querynode, samples=100, progression = False):
        self.create_shadow_circuit(querynode)
        results = []
        betasamples = {}
        cumulative = {}

        for k,v in self._betas.items():
            if k == 0:
                betasamples[k] = numpy.zeros(samples)
            elif k == 1:
                betasamples[k] = numpy.ones(samples)
            else:
                betasamples[k] = v.samples(samples)

        assoc_key_beta = {}
        for k in self._leaves:
            assoc_key_beta[k] = self.node(k).theta()

        # ## from https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists on 9 Jul 2019
        # # for s in (dict(zip(betas.keys(), values)) for values in itertools.product(*betas.values())):
        # # from https://stackoverflow.com/questions/5558418/list-of-dicts-to-from-dict-of-lists on 9 Jul 2019
        # for s in [dict(zip(betas, t)) for t in zip(*betas.values())]:
        #     for k, v in s.items():
        #         self.node(k).theta(v)
        for s in range(samples):

            for k,v in assoc_key_beta.items():
                self.node(k).theta(betasamples[v][s])

            rsample = self._compute_query_prob(querynode)

            results.append(rsample)

            if progression:
                cumulative[s+1] = [numpy.mean(results), numpy.var(results)]

        if progression:
            return cumulative
        else:
            return [numpy.mean(results), numpy.var(results)]

    def compute_query_prob(self, k):
        self.create_shadow_circuit(k)

        if k > self._lenorig:
            negk = k - self._lenorig
        else:
            negk = k + self._lenorig

        return self._compute_query_prob(k)


    def _compute_query_prob(self, k):

        if k > self._lenorig:
            negk = k - self._lenorig
        else:
            negk = k + self._lenorig


        means = numpy.zeros(len(self._leaves))

        visited = list(self._leaves)

        for i in range(len(visited)):
            n = visited[i]

            if self.node(n).negated():
                means[i] = 1 - self.node(n).weight()
            else:
                means[i] = self.node(n).weight()

        means = numpy.append(means, [0])
        visited.append(self.node(negk).alias())
        tovisit = set([k for k in self._nodes if self.node(k).type() != "ALIAS"]) - set(visited)

        while tovisit:
            nodei = None
            for nodei in tovisit:
                if self.node(nodei).children().issubset(set(visited)):
                    break
            ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).children()]

            if self.node(nodei).type() == "OR":
                visited.append(nodei)
                o = numpy.zeros(len(means))
                for c in ichildren:
                    o[c] = 1
                means = numpy.append(means, o @ means)

                if self.node(nodei).children().difference(self.node(nodei).qchildren()):
                    nnodei = self.node(nodei).alias()
                    visited.append(nnodei)
                    ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).qchildren()]
                    o = numpy.zeros(len(means))
                    for c in ichildren:
                        o[c] = 1
                    means = numpy.append(means, o @ means)

            elif self.node(nodei).type() == "AND":
                visited.append(nodei)
                ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).children()]
                child1 = ichildren[0]
                means = numpy.append(means, means[child1])
                for child2 in ichildren[1:]:
                    means[-1] *= means[child2]

                if self.node(nodei).children().difference(self.node(nodei).qchildren()):
                    nnodei = self.node(nodei).alias()
                    visited.append(nnodei)
                    ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).qchildren()]
                    child1 = ichildren[0]
                    means = numpy.append(means, means[child1])
                    for child2 in ichildren[1:]:
                        means[-1] *= means[child2]

            tovisit.remove(nodei)


        iroot = [i for i in range(len(visited)) if visited[i] == self._root][0]
        ialias = [i for i in range(len(visited)) if visited[i] == self.node(self._root).alias()][0]
        mralias = means[ialias]
        mroot = means[iroot]
        return mralias / mroot



    def compute_query(self, k):
        self.create_shadow_circuit(k)

        if k > self._lenorig:
            negk = k - self._lenorig
        else:
            negk = k + self._lenorig


        cov = numpy.zeros(shape=(len(self._leaves), len(self._leaves))) #(len(self._leaves),len(self._leaves)))
        means = numpy.zeros(len(self._leaves))

        visited = list(self._leaves)

        for i in range(len(visited)):
            n = visited[i]

            if self.node(n).negated():
                means[i]= self._betas[self.node(n).weight()].negate().mean()
            else:
                means[i] = self._betas[self.node(n).weight()].mean()
            pn = 1
            if self.node(n).negated():
                pn = -1
            for j in range(len(visited)):
                m = visited[j]
                pm = 1
                if self.node(m).negated():
                    pm = -1

                if self.node(n).weight() == self.node(m).weight():
                    if self.node(n).neutral():
                        cov[i, j] = 0.0
                        cov[j, i] = 0.0
                    else:
                        cov[i,j] = self._betas[self.node(n).weight()].variance() * pn * pm
                        cov[j, i] = cov[i, j]

        # add the negated query
        # self.cov(self.node(negk).alias(), self.node(-1*k).alias(), mpmath.mpf(1e-9))
        # self.m(self.node(negk).alias(), mpmath.mpf(1e-9))
        # visited.add(self.node(negk).alias())
        means = numpy.append(means, [0])
        visited.append(self.node(negk).alias())
        cov = numpy.concatenate((cov, numpy.zeros((len(visited)-1,1))), 1)
        cov = numpy.concatenate((cov, numpy.zeros((1, len(visited)))))


        tovisit = set([k for k in self._nodes if self.node(k).type() != "ALIAS"]) - set(visited)

        while tovisit:
            nodei = None
            for nodei in tovisit:
                if self.node(nodei).children().issubset(set(visited)):
                    break
            ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).children()]

            if self.node(nodei).type() == "OR":
                visited.append(nodei)
                o = numpy.zeros(len(means))
                for c in ichildren:
                    o[c] = 1
                means = numpy.append(means, o @ means)
                c = cov @ o
                v = o @ cov @ o
                cov = numpy.concatenate((cov, numpy.zeros((len(visited) - 1, 1))), 1)
                cov = numpy.concatenate((cov, numpy.zeros((1, len(visited)))))
                cov[-1,:-1] = c
                cov[:-1,-1] = c
                cov[-1, -1] = v

                if self.node(nodei).children().difference(self.node(nodei).qchildren()):
                    nnodei = self.node(nodei).alias()
                    visited.append(nnodei)
                    ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).qchildren()]
                    o = numpy.zeros(len(means))
                    for c in ichildren:
                        o[c] = 1
                    means = numpy.append(means, o @ means)
                    c = cov @ o
                    v = o @ cov @ o
                    cov = numpy.concatenate((cov, numpy.zeros((len(visited) - 1, 1))), 1)
                    cov = numpy.concatenate((cov, numpy.zeros((1, len(visited)))))
                    cov[- 1, :-1] = c
                    cov[:-1, - 1] = c
                    cov[- 1, - 1] = v

            elif self.node(nodei).type() == "AND":
                visited.append(nodei)
                ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).children()]
                child1 = ichildren[0]
                means = numpy.append(means, means[child1])
                c = cov[:, child1]
                c2 = c[child1]
                cov = numpy.concatenate((cov, numpy.zeros((len(visited) - 1, 1))), 1)
                cov = numpy.concatenate((cov, numpy.zeros((1, len(visited)))))
                cov[- 1, :-1] = c
                cov[:-1, - 1] = c
                cov[- 1, - 1] = c2
                for child2 in ichildren[1:]:
                    c = means[child2] * cov[:-1, -1] + means[-1] * cov[:-1, child2]
                    c2 = means[child2]**2 * cov[-1, -1] + means[-1]**2 * cov[child2, child2] + 2 * means[child2] * means[-1] *cov[child2, -1]
                    cov[-1, :-1] = c
                    cov[:-1, -1] = c
                    cov[-1, -1] = c2
                    means[-1] *= means[child2]

                if self.node(nodei).children().difference(self.node(nodei).qchildren()):
                    nnodei = self.node(nodei).alias()
                    visited.append(nnodei)
                    ichildren = [c for c in range(len(visited)) if visited[c] in self.node(nodei).qchildren()]
                    child1 = ichildren[0]
                    means = numpy.append(means, means[child1])
                    c = cov[:, child1]
                    c2 = c[child1]
                    cov = numpy.concatenate((cov, numpy.zeros((len(visited) - 1, 1))), 1)
                    cov = numpy.concatenate((cov, numpy.zeros((1, len(visited)))))
                    cov[- 1, :-1] = c
                    cov[:-1, - 1] = c
                    cov[- 1, - 1] = c2
                    for child2 in ichildren[1:]:
                        c = means[child2] * cov[:-1, -1] + means[-1] * cov[:-1, child2]
                        c2 = means[child2] ** 2 * cov[-1, -1] + means[-1] ** 2 * cov[child2, child2] + 2 * means[
                            child2] * means[-1] * cov[child2, -1]
                        cov[-1, :-1] = c
                        cov[:-1, -1] = c
                        cov[-1, -1] = c2
                        means[-1] *= means[child2]

            # removelist = set()
            # for i in range(len(visited)):
            #     if self.node(visited[i]).type() == "ALIAS":
            #         orig = self.node(visited[i]).reference()
            #         if self.node(orig).parents() and self.node(orig).parents().issubset(set(visited)):
            #             removelist.add(i)
            #     elif self.node(visited[i]).parents() and self.node(visited[i]).parents().issubset(set(visited)):
            #         removelist.add(i)
            # if removelist:
            #     selector = [i for i in range(len(visited)) if i not in removelist]
            #     cov = cov[selector, :][:, selector]
            #     means = means[selector]
            #     visited = [visited[i] for i in range(len(visited)) if i not in removelist]


            tovisit.remove(nodei)


        iroot = [i for i in range(len(visited)) if visited[i] == self._root][0]
        ialias = [i for i in range(len(visited)) if visited[i] == self.node(self._root).alias()][0]
        mralias = means[ialias]
        mroot = means[iroot]
        o = numpy.array([-mralias/mroot**2, 1/mroot])
        var = o @ cov[[iroot,ialias],:][:, [iroot,ialias]] @ o
        ev = mralias / mroot
        return [ev, var]







    def from_formula(self, formula, zero = 0, one = 1):
        negative = set([])
        associatedprobs = {}
        self._zero = zero
        self._one = one

        parents = {}

        idx = len(self._betas)

        self._lenorig = max([k for k,_,_ in formula])

        for index, _, _ in formula:
            node = formula.get_node(index)
            nodetype = type(node).__name__

            if nodetype == 'conj' or nodetype == 'disj':
                newnode = None
                if nodetype == 'conj':
                    newnode = AndNode()
                elif nodetype == 'disj':
                    newnode = OrNode()

                for c in node.children:
                    cindex = c
                    if c < 0: #and c not in negative:
                        cindex = abs(c) + self._lenorig
                        negative.add(cindex)

                    if c != 0:
                        newnode.add_child(cindex)
                        if cindex not in parents:
                            parents[cindex] = []
                        parents[cindex].append(index)
                self.add_node(index, newnode)

            elif nodetype == 'atom':

                newnode = None
                newnodeneg = None
                if index in dict(formula.evidence()).values(): #map(abs, dict(formula.evidence()).values()):
                    newnode = Evidence(zero, one)
                    newnodeneg = NegEvidence(zero, one)
                elif -index in dict(formula.evidence()).values():
                    newnode = NegEvidence(zero, one)
                    newnodeneg = Evidence(zero, one)
                elif index in dict(formula.queries()).values(): #map(abs, dict(formula.queries()).values()):
                    newnode = Query(zero)
                    newnodeneg = NegQuery(zero)
                elif -index in dict(formula.queries()).values(): #map(abs, dict(formula.queries()).values()):
                    newnode = NegQuery(zero)
                    newnodeneg = Query(zero)
                else:
                    newnode = Leaf(zero)
                    newnodeneg = Leaf(zero)

                if index not in map(abs, dict(formula.evidence()).values()):
                    if node.probability == formula.WEIGHT_NEUTRAL:
                        newnode.theta(one)
                        newnodeneg.theta(one)
                    elif node.group is None:
                        newnode.theta(int((node.probability)))
                        newnodeneg.theta(int(node.probability))
                        if index < 0:
                            newnode.negated(True)
                        else:
                            newnodeneg.negated(True)
                # else:
                #     raise Exception("Uninmplemented (yet)")
                    # clusters[node.group].append('%s [ shape="ellipse", label="%s", '
                    #                             'style="filled", fillcolor="white" ];\n'
                    #                             % (index, node.probability))
                if index > 0:
                    self.add_node(index, newnode)
                    self.add_node(index + self._lenorig, newnodeneg)
                else:
                    self.add_node(abs(index) + self._lenorig, newnode)
                    self.add_node(abs(index), newnodeneg)

            else:
                raise TypeError("Unexpected node type: '%s'" % nodetype)


        for k,v in parents.items():
            for p in v:
                self.node(k).add_parent(p)


    def _str_id(self, n):
        id = "n"
        if n > 0:
            id += str(n)
        else:
            id += "_" + str(abs(n))
        return id

    def __str__(self):
        ret = "digraph{\n"
        for n in self._nodes:
            if self.node(n).type() == "ALIAS":
                continue
            ret += self._str_id(n) + ' [label="' + str(n) + ": " + str(self.node(n)) + '"];\n'

            for c in self.node(n).children():
                ret += self._str_id(n) + " -> " + self._str_id(c) + ";\n"
        ret += "}"
        return ret