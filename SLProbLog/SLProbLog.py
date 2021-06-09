"""
Copyright (c) 2018 Federico Cerutti <CeruttiF@cardiff.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""
from problog.errors import InconsistentEvidenceError
from problog.evaluator import Semiring
from problog.engine import DefaultEngine
from problog.program import PrologString
from problog import get_evaluatable
from problog.formula import atom
from problog.formula import conj
from problog.logic import Constant
from problog import ddnnf_formula
import pandas as pd
import os
from collections import namedtuple, defaultdict, OrderedDict
from scipy.stats import beta
import numpy as np
import re
import itertools
from problog.logic import Term, Or, Clause, And, is_ground
from SLProbLog.Tree import Tree
import copy





import mpmath

EPSILON = 10e-100
mpmath.dps = 200

def from_sl_opinion(wb, W = 2):
    prior = mpmath.mpf(W)

    [belief, disbelief, uncertainty, base] = wb[0:4]

    if mpmath.almosteq(mpmath.mpf("0"), uncertainty, EPSILON):
        uncertainty = mpmath.mpf(EPSILON)

    if mpmath.almosteq(mpmath.mpf("0"), disbelief, EPSILON):
        disbelief = mpmath.mpf(EPSILON)

    if mpmath.almosteq(mpmath.mpf("0"), belief, EPSILON):
        belief = mpmath.mpf(EPSILON)

    mean = belief + uncertainty * base
    sx = prior / uncertainty
    variance = mean * (1- mean) / (sx + 1)

    return BetaDistribution(mean, variance)


def moment_matching(b):
    m = b.mean()
    v = b.variance()
    if v == 0:
        var = mpmath.mpf(1e-10)
    else:
        var = mpmath.mpf(v)

    mean = min(mpmath.mpf(1), max(mpmath.mpf(0), mpmath.mpf(m)))
    var = min(var, mean ** 2 * (1.0 - mean) / (1.0 + mean), (1.0 - mean) ** 2 * mean / (2 - mean))
    #var = min(var, mean ** 2 * (1.0 - mean) / (1.0 + mean), (1.0 - mean) ** 2 * mean / (2 - mean))

    #sx = ((mean * (1 - mean)) / var - 1)

    #return BetaDistribution(mean, (mean * (1-mean) / (sx + 1) ))
    return BetaDistribution(mean, var)



def from_alpha_beta(a, b):
    return BetaDistribution(a / (a + b), a * b / ((a + b)**2 * (a + b + 1)))


class BetaDistribution():

    def __init__(self, m, v):
        self._epsilon = EPSILON
        self._ZERO = mpmath.mpf("0")
        self._ONE = mpmath.mpf("1")
        self._mu = mpmath.mpf(m)
        self._var = mpmath.mpf(v)

    def is_complete_belief(self):
        if mpmath.almosteq(self.mean(), self._ONE, self._epsilon):
            return True
        return False

    def mean(self):
        return self._mu

    def variance(self):
        return self._var

    def strength(self):
        var = self.variance()
        if mpmath.almosteq(var, 0, self._epsilon):
            var = mpmath.mpf(self._epsilon)

        return (self.mean() * (1 - self.mean())) / var - 1

    def alpha(self):
        if self.mean() == 1.0:
            return mpmath.inf
        return max(mpmath.mpf(self._epsilon), self.mean() * self.strength())

    def beta(self):
        if self.mean() == 0.0:
            return mpmath.inf
        return max(mpmath.mpf(self._epsilon), (1 - self.mean()) * self.strength())

    def sum(self, Y):
        mean = self.mean() + Y.mean()
        var = self.variance() + Y.variance()

        var = min(var, mean ** 2 * (1.0 - mean) / (1.0 + mean), (1.0 - mean) ** 2 * mean / (2 - mean))

        return BetaDistribution(mean,var)

    def product(self, Y):
        mean = self.mean() * Y.mean()
        var = self.variance() * Y.variance() + \
              self.variance() * (Y.mean())**2 + Y.variance() * (self.mean()) ** 2

        var = min(var, mean ** 2 * (1.0 - mean) / (1.0 + mean), (1.0 - mean) ** 2 * mean / (2 - mean))

        return BetaDistribution(mean, var)

    def negate(self):
        if not 0 <= self.mean() <= 1:
            raise Exception("Error with negation: [%f, %f]", (self.mean(), self.variance()))
        return BetaDistribution(1.0 - self.mean(), self.variance())

    def conditioning(self, Y):
        mean = min(1.0-1e-6, self.mean() / Y.mean())

        muneg = Y.mean() - self.mean() #+ Y.mean() * self.mean()
        varsum = Y.variance() + self.variance()



        if self._mu <= 0:
            self._mu = mpmath.mpf(1e-10)

        if muneg <= 0:
            muneg = mpmath.mpf(1e-10)

        var = mean**2 * (1.0-mean)**2 * ((self.variance() / (self._mu ** 2)) + (varsum / (muneg ** 2)) - 2 * (self.variance() / (self._mu * muneg)))

        if var < 0:
            var = min(mean ** 2 * (1.0 - mean) / (1.0 + mean), (1.0 - mean) ** 2 * mean / (2 - mean))

        var = min(var, mean ** 2 * (1.0 - mean) / (1.0 + mean), (1.0 - mean) ** 2 * mean / (2 - mean))

        return BetaDistribution(mean, var)

    def __repr__(self):
        return "b(%s,%s)" % (mpmath.nstr(self.mean(), mpmath.mp.dps), mpmath.nstr(self.variance(), mpmath.mp.dps))

    def mean_str(self):
        return mpmath.nstr(self.mean(), mpmath.mp.dps)

    def variance_str(self):
        return mpmath.nstr(self.variance(), mpmath.mp.dps)

    def to_sl_opinion(self, a = 1/2, W=2):
        if self.alpha() == mpmath.inf:
            return (1, 0, 0, mpmath.mpf(a))

        if self.beta() == mpmath.inf:
            return (0, 1, 0, mpmath.mpf(a))

        rx = max(mpmath.mpf(0), self.alpha() - a * W)
        sx = max(mpmath.mpf(0), self.beta() - (1-a) * W)
        return ((rx / (rx + sx + W)), (sx / (rx + sx + W)), (W / (rx + sx + W)), mpmath.mpf(a))

    def betavariate(self):
        return "beta(%s,%s)" % (mpmath.nstr(self.alpha(), mpmath.mp.dps), mpmath.nstr(self.beta(), mpmath.mp.dps))

    def samples(self, numsamples = 100):
        samps = beta.rvs(float(self.alpha()), float(self.beta()), size=numsamples)
        # if abs(self.mean() - np.mean(samps)) > 0.01 or abs(self.variance() - np.var(samps)) > 0.01:
        #     print("Beta(%s, %s); samples mean: %s; samples var: %s" % (self.mean(), self.variance(), np.mean(samps), np.var(samps)))
        return samps

class SLSemiring(Semiring):

    def parse(self, w):
        start = w.find('(') + 1
        end = w.find(')')
        ret = [mpmath.mpf(x) for x in w[start:end].replace(" ","").split(',')]

        return ret

    def one(self):
        return "w(1.0, 0.0, 0.0, 0.99999999)"

    def zero(self):
        return "w(0.0, 1.0, 0.0, 0.00000001)"

    def plus(self, x, y):
        [b1,d1,u1,a1] = self.parse(x)[0:4]
        [b2,d2,u2,a2] = self.parse(y)[0:4]
        u = (a1 * u1 + a2 * u2) / (a1 + a2)
        d = max(0.0, (a1 * (d1 - b2) + a2 * (d2 - b1)) / (a1 + a2))
        b = min(b1 + b2, 1.0)
        a = min(a1 + a2, 1.0)
        return "w(%s,%s,%s,%s)" % (str(b), str(d), str(u), str(a))

    def times(self, x, y):
        [b1, d1, u1, a1] = self.parse(x)[0:4]
        [b2, d2, u2, a2] = self.parse(y)[0:4]
        a = a1 * a2
        b = b1 * b2 + ((1 - a1) * a2 * b1 * u2 + a1 * (1 - a2) * u1 * b2) / (1 - a1 * a2)
        u = u1 * u2 + ((1 - a2) * b1 * u2 + (1 - a1) * u1 * b2) / (1 - a1 * a2)
        d = min(1, d1 + d2 - d1 * d2)
        return "w(%s,%s,%s,%s)" % (str(b), str(d), str(u), str(a))

    def negate(self, a):
        [b1, d1, u1, a1] = self.parse(a)[0:4]
        return "w(%s,%s,%s,%s)" % (str(d1), str(b1), str(u1), str(1 - a1))

    def value(self, a):
        return str(a)

    def normalize(self, x, z):

        if z == self.one():
            return x

        [b1, d1, u1, a1] = self.parse(x)[0:4]
        [b2, d2, u2, a2] = self.parse(z)[0:4]
        e1 = b1 + u1*a1
        e2 = b2+ u2 * a2

        if not ((a1<=a2) and (d1>=d2) and (b1*(1-a1)*a2*(1-d2) >= a1*(1-a2)*(1-d1)*b2) and (u1*(1-a1)*(1-d2)>=u2*(1-a2)*(1-d1)) and a2!=0 ):
            return "w(%s,%s,%s,%s)" % (str(0.0), str(0.0), str(1.0), str(0.5))
        else:
            a = a1/a2
            b = 0.0
            d = 0.0
            u = 0.0
            if e1 == 0:
                d = 1.0
            elif a==1:
                b = 1.0
            else:
                e = e1 / e2
                d = min(max(0, (d1 - d2) / (1 - d2)), 1)
                u = min(max(0, (1 - d - e) / (1 - a)), 1)
                b = min(max(0, (1 - d - u)), 1)
            return "w(%s,%s,%s,%s)" % (str(b), str(d), str(u), str(a))

    def is_dsp(self):
        return True

class BetaSemiring(Semiring):

    def parse(self, w):
        start = str(w).find('(') + 1
        end = str(w).find(')')
        parsed = [mpmath.mpf(x) for x in str(w)[start:end].replace(" ","").split(',')]
        return BetaDistribution(parsed[0], parsed[1])

    def one(self):
        return "b(1.0,0.000000001)"

    def zero(self):
        return "b(0.0,0.000000001)"

    def plus(self, a, b):
        wa = self.parse(a)
        wb = self.parse(b)
        return self._to_str(wa.sum(wb))

    def times(self, a, b):
        wa = self.parse(a)
        wb = self.parse(b)
        return self._to_str(wa.product(wb))

    def negate(self, a):
        wa = self.parse(a)
        return self._to_str(wa.negate())

    def value(self, a):
        return str(a)

    def _to_str(self, r):
        wr = self.parse(r)
        return wr.__repr__()

    def normalize(self, a, z):
        wa = self.parse(a)
        wz = self.parse(z)

        if wz.is_complete_belief():
            return a
        return self._to_str(wa.conditioning(wz))

    def is_dsp(self):
        return True



class SingletonBetas(object):
    class __SingletonBetas:
        def __init__(self):
            self.betas = None

    instance = None
    def __new__(cls):
        if not SingletonBetas.instance:
            SingletonBetas.instance = SingletonBetas.__SingletonBetas()
        return SingletonBetas.instance

    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)


class SLProbLog:

    def __init__(self, program, sloutput = False):
        self._slproblog_program = program
        self._slout = sloutput
        self._lance_assoc_nodes = {}
        self._lance_new_nodeid = 0
        self._lance_assoc_f = {}
        self._lance_new_f_id = 0
        self._lance_assoc_p = {}
        self._lance_new_p_id = 2 #p1 is reserved to identity #p2 for zero #p0 for True meaning that must not be negated
        self.debug = False
        SingletonBetas.instance = {}
        SingletonBetas.instance[0] = BetaDistribution(0.0,mpmath.mpf(1e-200))
        SingletonBetas.instance[1] = BetaDistribution(1.0,mpmath.mpf(1e-200))


    def _convert_input(self, to_sl = False, to_beta = False, to_singleton = False):
        if to_sl and to_beta:
            raise Exception("Cannot convert both to SL and to Beta")

        lines = self._slproblog_program.splitlines()
        newlines = []
        idx = len(SingletonBetas.instance)
        for l in lines:
            label_prolog = l.split("::")
            if len(label_prolog) == 1:
                newlines.append(l)
            elif len(label_prolog) == 2:
                if to_sl:
                    wb = None
                    if label_prolog[0][0] == "b":
                        w = label_prolog[0]
                        start = str(w).find('(') + 1
                        end = str(w).find(')')
                        parsed = [mpmath.mpf(x) for x in str(w)[start:end].replace(" ", "").split(',')]
                        b = BetaDistribution(parsed[0], parsed[1])
                        wb = b.to_sl_opinion()

                    elif label_prolog[0][0] == "w":
                        wb = SLSemiring().parse(label_prolog[0])

                    else:
                        raise Exception("Problem with this line: %s" % (l))

                    newlines.append("w(%s,%s,%s,%s,%s)::%s" % (mpmath.nstr(wb[0]),
                                                                mpmath.nstr(wb[1]),
                                                                mpmath.nstr(wb[2]),
                                                                mpmath.nstr(wb[3]),
                                                                mpmath.nstr(idx),
                                                                label_prolog[1]))
                    idx += 1


                elif to_beta or to_singleton:
                    bw = None
                    if label_prolog[0][0] == "w":
                        w = SLSemiring().parse(label_prolog[0])
                        bw = from_sl_opinion(w)


                    elif label_prolog[0][0] == "b":
                        w = label_prolog[0]
                        start = str(w).find('(') + 1
                        end = str(w).find(')')
                        parsed = [mpmath.mpf(x) for x in str(w)[start:end].replace(" ", "").split(',')]
                        bw = BetaDistribution(parsed[0], parsed[1])
                        #bw = BetaSemiring().parse(label_prolog[0])

                    else:
                        raise Exception("Problem with this line: %s" % (l))

                    if to_singleton:
                        newlines.append("%d::%s" % (idx, label_prolog[1]))
                        SingletonBetas.instance[idx] = bw
                        idx += 1
                    else:
                        newlines.append("%s::%s" % (repr(bw), label_prolog[1]))

            else:
                raise Exception("Problem with this line: %s" % (l))

        return "\n".join(newlines)

    def _convert_output(self, res, to_sl = False, to_beta = False):
        if to_sl and to_beta:
            raise Exception("Cannot convert both to SL and to Beta")

        ret = {}
        for k, v in res.items():
            if to_sl and isinstance(v, BetaDistribution):
                ret[k] = v.to_sl_opinion()
            elif to_beta and isinstance(v, list):
                ret[k] = from_sl_opinion(v)

        return ret


    def run_SL(self):
        res = self._run_sl_operators_on_semiring(SLSemiring(), self._convert_input(to_sl = True))
        if self._slout:
            return res
        return self._convert_output(res, to_beta = True)

    def run_beta(self):
        res = self._run_sl_operators_on_semiring(BetaSemiring(), self._convert_input(to_beta = True))
        if self._slout:
            return self._convert_output(res, to_sl=True)
        return res

    def iter_sample(self, minsample = 20, maxsample = 200, step = 20):
        engine = DefaultEngine(label_all=True)
        db = engine.prepare(PrologString(self._convert_input(to_beta=True)))
        lf = engine.ground_all(db)
        kc_class = get_evaluatable()  # name='ddnnf')
        kc = kc_class.create_from(lf)
        if self.debug:
            print(kc)

        results = {}
        listindexes = []
        allbetas = {}
        for index, node, nodetype in kc:
            if nodetype == 'atom' and node.probability != kc.WEIGHT_NEUTRAL:
                listindexes.append(index)
                allbetas[index] = BetaSemiring().parse(str(node.probability)).samples(maxsample)

        idx = 0
        res = {}
        while minsample + idx * step <= maxsample:

            itersample = minsample + idx * step
            idx += 1

            betas = {}
            for index in listindexes:
                betas[index] = allbetas[index][:itersample]

            ## from https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists on 9 Jul 2019
            # for s in (dict(zip(betas.keys(), values)) for values in itertools.product(*betas.values())):
            # from https://stackoverflow.com/questions/5558418/list-of-dicts-to-from-dict-of-lists on 9 Jul 2019
            for s in [dict(zip(betas, t)) for t in zip(*betas.values())]:
                for k, v in s.items():
                    p = kc.get_node(k)
                    newconstant = Constant(v)
                    newconstant._cache_is_ground = True
                    kc._update(k, p._replace(probability=newconstant))
                    kc._weights[k] = newconstant

                if self.debug:
                    print(kc)
                    print(kc.evaluate())

                rsample = kc.evaluate()

                if not results:
                    for k, v in rsample.items():
                        results[str(k)] = [v]
                else:
                    for k, v in rsample.items():
                        results[str(k)].append(v)

            resbeta = {}
            for k, v in results.items():
                resbeta[str(k)] = moment_matching(BetaDistribution(mpmath.mpf(np.mean(v)), mpmath.mpf(np.var(v))))

            resbeta = self._order_dicts(resbeta)

            if self._slout:
                resbeta = self._convert_output(resbeta, to_sl=True)

            res[itersample] = resbeta
        return res

    def run_sample(self, samples=100, fitbeta = True, progression = False):
        program = self._convert_input(to_sl=False, to_beta=False, to_singleton=True)
        # program = self._convert_input(to_sl=True, to_beta=False)
        engine = DefaultEngine()
        if program == None:
            program = self._slproblog_program

        db = engine.prepare(PrologString(program))
        knowledge = get_evaluatable(name='ddnnf')

        kc = knowledge.create_from(db, engine=engine, database=db)

        results = {}
        resprogression = {}
        listindexes = []
        betas = {}

        betasamples = {}
        cumulative = {}

        for k, v in SingletonBetas.instance.items():
            if k == 0:
                betasamples[k] = np.zeros(samples)
            elif k == 1:
                betasamples[k] = np.ones(samples)
            else:
                betasamples[k] = v.samples(samples)

        assoc_key_beta = {}

        for index, node, nodetype in kc:
            if nodetype == 'atom' and node.probability != kc.WEIGHT_NEUTRAL:
                # listindexes.append(index)
                # betas[index] = BetaSemiring().parse(str(node.probability)).samples(samples)
                assoc_key_beta[index] = node.probability

        ## from https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists on 9 Jul 2019
        #for s in (dict(zip(betas.keys(), values)) for values in itertools.product(*betas.values())):
        # from https://stackoverflow.com/questions/5558418/list-of-dicts-to-from-dict-of-lists on 9 Jul 2019
        # for s in [dict(zip(betas,t)) for t in zip(*betas.values())]:
        #     for k, v in s.items():
        for s in range(samples):
            for k, v in assoc_key_beta.items():
                p = kc.get_node(k)
                newconstant = Constant(betasamples[v][s])
                newconstant._cache_is_ground = True
                kc._update(k, p._replace(probability=newconstant))
                kc._weights[k] = newconstant


            if self.debug:
                print(kc)
                print(kc.evaluate())

            rsample = kc.evaluate()


            if not results:
                for k,v in rsample.items():
                    results[str(k)] = [v]
            else:
                for k,v in rsample.items():
                    results[str(k)].append(v)

            if progression:
                resprogression[s+1] = {}
                toadd = {}
                for k,v in rsample.items():
                    toadd[str(k)] = moment_matching(BetaDistribution(mpmath.mpf(np.mean(results[str(k)])), mpmath.mpf(np.var(results[str(k)]))))

                if self._slout:
                    toadd= self._convert_output(toadd, to_sl=True)

                resprogression[s+1] = toadd

        if not fitbeta:
            return results

        resbeta = {}
        for k,v in results.items():
            resbeta[str(k)] = moment_matching(BetaDistribution(mpmath.mpf(np.mean(v)), mpmath.mpf(np.var(v))))

        resbeta = self._order_dicts(resbeta)

        if progression:
            return resprogression


        if self._slout:
            return self._convert_output(resbeta, to_sl=True)
        return resbeta



    def _lance_get_node_id(self, index, pos=True, sizeformula=None):
        if index not in self._lance_assoc_nodes and pos and index > 0:
            #            self._lance_new_nodeid += 1
            self._lance_assoc_nodes[index] = index
        if index not in self._lance_assoc_nodes and (not pos or index < 0):
            self._lance_assoc_nodes[index] = abs(index) + max([k for k, _, _ in self._formula])

        return self._lance_assoc_nodes[index]


    # def _lance_get_node_id(self, index): ### ORIGINAL
    #     if index not in self._lance_assoc_nodes:
    #         self._lance_new_nodeid += 1
    #         self._lance_assoc_nodes[index] = self._lance_new_nodeid
    #
    #     return self._lance_assoc_nodes[index]

    def _lance_get_f_id(self, index):
        if isinstance(index, int) and index < 0:
            index = -1 * index

        if index not in self._lance_assoc_f:
            self._lance_new_f_id += 1
            self._lance_assoc_f[index] = self._lance_new_f_id

        return self._lance_assoc_f[index]


    def _lance_get_p_id(self, index):
        if isinstance(index, int) and index < 0:
            index = -1 * index

        if index not in self._lance_assoc_p:
            self._lance_new_p_id += 1
            self._lance_assoc_p[index] = self._lance_new_p_id

        return self._lance_assoc_p[index]




    def run_beta_cov(self):
        program = self._convert_input(to_sl=False, to_beta=False, to_singleton=True)
        #program = self._convert_input(to_sl=True, to_beta=False)
        engine = DefaultEngine()
        if program == None:
            program = self._slproblog_program

        db = engine.prepare(PrologString(program))
        knowledge = get_evaluatable(name='ddnnf')

        formula = knowledge.create_from(db, engine=engine, database=db)

        mytree = Tree(SingletonBetas.instance)
        try:
            mytree.from_formula(formula)
        except ValueError as e:
            pass

        resunordered = {}
        for k, v in formula.queries():
            if v is None:
                resunordered[str(k)] = from_sl_opinion([0, 0, 1, 0.5])
            else:
                mytreecopy = copy.deepcopy(mytree)
                resmytree = mytreecopy.compute_query(v)
                meanvar = BetaDistribution(resmytree[0], resmytree[1])
                resunordered[str(k)] = moment_matching(meanvar)

        res = OrderedDict(sorted(resunordered.items()))

        if self._slout:
            return self._convert_output(res, to_sl=True)
        return res

    def run_sample_cov(self, samples=100, fitbeta=True, progression=False):
        program = self._convert_input(to_sl=False, to_beta=False, to_singleton=True)
        # program = self._convert_input(to_sl=True, to_beta=False)

        resprogression = {}

        engine = DefaultEngine()
        if program == None:
            program = self._slproblog_program

        db = engine.prepare(PrologString(program))
        knowledge = get_evaluatable(name='ddnnf')

        formula = knowledge.create_from(db, engine=engine, database=db)

        mytree = Tree(SingletonBetas.instance)
        mytree.from_formula(formula)

        resunordered = {}
        for k, v in formula.queries():
            if v is None:
                resunordered[str(k)] = BetaSemiring().parse(BetaSemiring().zero())
            else:
                mytreecopy = copy.deepcopy(mytree)

                if progression:
                    progr =  mytreecopy.compute_query_monte(v, samples, progression)
                    for pk,pv in progr.items():
                        resprogression[pk] = {str(k): moment_matching(BetaDistribution(pv[0], pv[1]))}
                else:
                    resmytree = mytreecopy.compute_query_monte(v, samples)
                    meanvar = BetaDistribution(resmytree[0], resmytree[1])
                    resunordered[str(k)] = moment_matching(meanvar)

        if progression:
            for pk,pv in resprogression.items():
                toorder = pv
                if self._slout:
                     toorder = self._convert_output(pv, to_sl=True)
                resprogression[pk] = OrderedDict(sorted(toorder.items()))

            return resprogression

        else:
            res = OrderedDict(sorted(resunordered.items()))

            if self._slout:
                return self._convert_output(res, to_sl=True)
            return res

    # Obsolete. Use run_beta_cov(). Kept for back-compatibility
    def run_KL(self, MatlabPath=None):
        return self.run_beta_cov()

    # Old code for Beta Cov
    def run_KL_obsolete(self, MatlabPath=None):
        program = self._convert_input(to_sl=True,to_beta=False)
        engine = DefaultEngine()
        if program == None:
            program = self._slproblog_program

        db = engine.prepare(PrologString(program))
        knowledge = get_evaluatable(name='ddnnf')

        formula = knowledge.create_from(db, engine=engine, database=db)

        #
        if self.debug:
            print(formula.to_dot())

        if MatlabPath is None:
            return self._to_KL(formula)
        else:
            return self._to_KL(formula, MATLABPATH=MatlabPath)

    def _run_sl_operators_on_semiring(self, givensemiring, program = None):
        engine = DefaultEngine()
        if program == None:
            program = self._slproblog_program

        db = engine.prepare(PrologString(program))
        semiring = givensemiring
        knowledge = get_evaluatable(name='ddnnf', semiring=semiring)

        formula = knowledge.create_from(db, engine=engine, database=db)

        #print(self.to_dot(formula))

        # myformula = MyDDNNF()
        # myformula.__dict__.update(formula.__dict__)

        res = formula.evaluate(semiring=semiring)

        ret = {}
        for k, v in res.items():
            if isinstance(semiring, BetaSemiring):
                ret[k] = moment_matching(semiring.parse(v))
            else:
                ret[k] = semiring.parse(v)

        return self._order_dicts(ret)


    def _order_dicts(self, dicinput):
        res = {}
        for k, v in dicinput.items():
            res[str(k)] = v

        ret = {}
        for k in sorted(res):
            ret[k] = res[k]

        return ret