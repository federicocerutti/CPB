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

import sys
from SLProbLog.SLProbLog import SLProbLog
from itertools import product
import numpy.random
import pickle
import datetime
import math
from string import Formatter,Template
from problog import get_evaluatable
from problog.program import PrologString
import mpmath
from oct2py import octave
import matplotlib.pyplot as plt
import matplotlib
plt.switch_backend('agg')
import os
from SLProbLog.SLProbLog import from_sl_opinion
import numpy as np
import matplotlib.patches as mpatches
import csv
import time
import scipy.special as sc

MATLABPATH = '/libmatlab'

class Node:
    """
    Node in a Bayesian network
    """

    def __init__(self, name):
        self._name = name
        self._parents = []
        self._children = []

    def get_name(self):
        return self._name

    def get_parents(self):
        return self._parents

    def add_parent(self, node):
        self._parents.append(node)

    def add_children(self, node):
        self._children.append(node)

    def get_children(self):
        return self._children

    def get_problog_name(self):
        return "n" + self._name

    def __repr__(self):
        return "n" + self._name

class Graph:
    """
    Data structure to represent a Bayesian network
    """

    def __init__(self):
        self._nodes = {}

    def storeFilename(self, fn):
        self.filename = fn

    def getNodes(self):
        return self._nodes

    def getFilename(self):
        return self.filename

    def add_edge(self, edge):
        for n in edge:
            if n not in self._nodes:
                self._nodes[n] = Node(n)

        self._nodes[edge[1]].add_parent(self._nodes[edge[0]])
        self._nodes[edge[0]].add_children(self._nodes[edge[1]])

    def __repr__(self):
        return str(self.__dict__)

    def get_problog_string(self):
        ret = ""
        end_loc = []
        query_loc = []
        p_index = 0

        self.semanticsprobs = {}
        self.semanticsevidences = {}
        for n in self._nodes:

            if len(self._nodes[n].get_parents()) + len(self._nodes[n].get_children()) == 1:
                end_loc.append(n)
            else:
                query_loc.append(n)

            if not self._nodes[n].get_parents():
                ret += "${p" + str(p_index) + "}::" + str(self._nodes[n].get_problog_name()) + ".\n"
                self.semanticsprobs[self._nodes[n]] = "p"+str(p_index)
                p_index += 1
            else:
                self.semanticsprobs[self._nodes[n]] = []
                for parents in product((False, True), repeat=len(self._nodes[n].get_parents())):
                    ret += "${p" + str(p_index) + "}::" + str(self._nodes[n].get_problog_name()) + " :- "
                    self.semanticsprobs[self._nodes[n]].append("p"+str(p_index))
                    p_index += 1
                    for i in range(len(parents)):
                        if not parents[i]:
                            ret += "\+"
                        ret += str(self._nodes[n].get_parents()[i].get_problog_name())

                        if i < len(parents) - 1:
                            ret += ", "
                    ret += ".\n"

        evidencenum = 0
        for n in end_loc:
            ret += "evidence("+str(self._nodes[n].get_problog_name()) + ", ${e" + str(evidencenum) + "}).\n"
            self.semanticsevidences[self._nodes[n]] = "e"+str(evidencenum)
            evidencenum += 1

        self.querynames = []
        for n in query_loc:
            ret += "query(" + str(self._nodes[n].get_problog_name()) + ").\n"
            self.querynames.append(str(self._nodes[n].get_problog_name()))
        return ret

class ProbProblog:
    """
    Given a template problog string, e.g.
    ${p1}::stress(X) :- person(X).
    evidence(..., ${e1}).
    randomly allocates values for those parameters.

    In addition, it servers as wrapper for ProbLog (function run)
    """
    def __init__(self, problogstring, network = None):
        self.network = network
        self.problogstring = problogstring
        self.keys = [ele[1] for ele in Formatter().parse(self.problogstring) if ele[1]]

        self.probabilities = {}
        self.evidences = {}
        for k in self.keys:
            if "e" in k:
                self.evidences[k] = ("true" if numpy.random.uniform(0, 1) < 0.5 else "false")
            else:
                self.probabilities[k] = numpy.random.uniform(0, 1)

    def getProbabilities(self):
        return self.probabilities

    def getEvidences(self):
        return self.evidences

    def getStructure(self):
        return self.problogstring

    def getProblogProgram(self):
        substitiutions = self.probabilities.copy()
        substitiutions.update(self.evidences)
        s = Template(self.problogstring).safe_substitute(substitiutions)
        return s

    def run(self):
        """
        Problog wrapper
        """
        res = {}
        for k, v in (get_evaluatable().create_from(PrologString(self.getProblogProgram())).evaluate()).items():
            res[str(k)] = v

        ret = {}
        for k in sorted(res):
            ret[k] = res[k]

        return ret

    def dataToSendToMatlab(self):
        tosend = []


        for k, v in sorted(self.network.getNodes().items(), key=lambda x: float(str(x[0]))):

            nodedesc = {}

            if len(v.get_children()) == 0:
                nodedesc["children"] = []
            elif len(v.get_children()) == 1:
                nodedesc["children"] = int(v.get_children()[0]._name)
            else:
                kids = []
                for k in v.get_children():
                    kids.append(int(k._name))
                nodedesc["children"] = kids

            if len(v.get_parents()) == 0:
                nodedesc["parents"] = []
            elif len(v.get_parents()) == 1:
                nodedesc["parents"] = int(v.get_parents()[0]._name)
            else:
                par = []
                for p in v.get_parents():
                    par.append(int(p._name))
                nodedesc["parents"] = par

            if v in self.network.semanticsprobs:
                if isinstance(self.network.semanticsprobs[v], str):
                    #print(self.network.semanticsprobs[v])
                    nodedesc["p"] = self.probabilities[self.network.semanticsprobs[v]]
                else:
                    nodedesc["p"] = []
                    for indexprob in self.network.semanticsprobs[v]:
                        nodedesc["p"].append(self.probabilities[indexprob])
            else:
                nodedesc["p"] = []


            if v in self.network.semanticsevidences:
                nodedesc["value"] = True if self.evidences[self.network.semanticsevidences[v]] == "true" else False
            else:
                nodedesc["value"] = []

            tosend.append(nodedesc)
        return tosend


class Opinion():
    """
    Utility function to store a 4-tuple representing a SL opinion and outputting it
    """
    def getBelief(self):
        return mpmath.nstr(self._belief, mpmath.mp.dps)

    def getDisbelief(self):
        return mpmath.nstr(self._disbelief, mpmath.mp.dps)

    def getUncertainty(self):
        return mpmath.nstr(self._uncertainty, mpmath.mp.dps)

    def getBase(self):
        return mpmath.nstr(self._base, mpmath.mp.dps)

    def toList(self):
        return [float(self.getBelief()), float(self.getDisbelief()), float(self.getUncertainty()), float(self.getBase())]

    def __init__(self, b, d, u=None, a=None):
        self._belief = mpmath.mpf(b)
        self._disbelief = mpmath.mpf(d)

        if u is None:
            self._uncertainty = mpmath.mpf(1 - self._belief - self._disbelief)
        else:
            self._uncertainty = mpmath.mpf(u)

        if a is None:
            self._base = mpmath.mpf(0.5)
        else:
            self._base = mpmath.mpf(a)

class DistProbLog:
    """
    Given a ProbProblog object, this class samples the randomly chosen probabilities ntrain times in order to then
    derive beta distributions
    """
    def __init__(self, bnet, ntrain=10):
        self.bn = bnet

        self.samples = {}
        self.opinions = {}
        self.betas = {}
        for p in bnet.getProbabilities():
            self.samples[p] = []

            for i in range(ntrain):
                self.samples[p].append(1 if numpy.random.uniform(0, 1) < bnet.getProbabilities()[p] else 0)

            rcount = sum(self.samples[p])
            scount = ntrain - rcount
            self.opinions[p] = Opinion(rcount / (rcount + scount + 2), scount / (rcount + scount + 2))
            self.betas[p] = from_sl_opinion(self.opinions[p].toList())

    def get_program(self):
        """
        Substitute the various probabilities signposts with SL opinions
        """
        substitutions = {}
        for k in self.opinions:
            substitutions[k] = "w(%s,%s,%s,%s)" % (self.opinions[k].getBelief(),
                                                   self.opinions[k].getDisbelief(),
                                                   self.opinions[k].getUncertainty(),
                                                   self.opinions[k].getBase())
        substitutions.update(self.bn.evidences)
        return Template(self.bn.problogstring).safe_substitute(substitutions)

    def get_program_beta(self):
        """
                Substitute the various probabilities signposts with SL opinions
                """
        substitutions = {}
        for k in self.opinions:
            substitutions[k] = "b(%s,%s)" % (self.betas[k].mean_str(), self.betas[k].variance_str())
        substitutions.update(self.bn.evidences)
        return Template(self.bn.problogstring).safe_substitute(substitutions)

    def _run_matlab_code(self, number):
        """
        Wrapper to run Matlab codes
        """
        octave.addpath(os.getcwd() + MATLABPATH)
        raw = octave.BNfromPY(self.dataToSendToMatlab(), number)

        ret = {}
        idx = 0

        if number == 1:
            for k in sorted(self.bn.network.querynames):
                ret[k] = raw[0][idx]
                idx += 1
        else:
            for k in sorted(self.bn.network.querynames):
                r = raw[idx]
                idx += 1
                ret[k] = [r[0][0], r[0][1], r[0][2], 0.5]

        return ret

    def runSBN(self):
        return self._run_matlab_code(0)

    def runCredal(self):
        return self._run_matlab_code(2)

    def runGBT(self):
        return self._run_matlab_code(3)

    def runVBS(self):
        return self._run_matlab_code(4)

    def dataToSendToMatlab(self):
        tosend = self.bn.dataToSendToMatlab()
        item = 0
        for k, v in sorted(self.bn.network.getNodes().items(), key=lambda x: float(str(x[0]))):
            nodedesc = tosend[item]
            item += 1
            if v in self.bn.network.semanticsprobs:
                if isinstance(self.bn.network.semanticsprobs[v], str):
                    nodedesc["w"] = ((self.opinions[self.bn.network.semanticsprobs[v]]).toList())[:-1]
                else:
                    nodedesc["w"] = []
                    for indexprob in self.bn.network.semanticsprobs[v]:
                        nodedesc["w"].append(((self.opinions[indexprob]).toList())[:-1])
            else:
                nodedesc["w"] = []

        return tosend



class Experiment():
    """
    Class collecting methods for running experiments
    """

    GROUND = "Ground"
    BETAPROBLOG = "BetaProblog"
    HONOLULU = "Honolulu"
    JOSANG = "Josang"
    SBN = "SBN"
    GBT = "GBT"
    CREDAL = "Credal"
    SAMPLESIMP = "MonteCarlo"
    PYBETAPROBLOG = "PyBetaProblog"


    def loadExperiment(picklefile):
        """
        Load from a pickle file
        """
        return pickle.load(open(picklefile, "rb"))

    def __init__(self, picklefile = None):
        """
        Loads from the pickle if passed as parameter, otherwise just creates empty attributes
        """
        if picklefile:
            self = pickle.load(open(picklefile, "rb"))
            print("hello")
        else:
            self._vec_real = []
            self._vec_sl = []
            self._vec_sl_beta = []
            self._vec_sl_beta_cov = []
            self._vec_sbn = []
            self._vec_GBT = []
            self._vec_credal = []
            self.bns = []
            self.net = []
            self._sb = []

            self._results = {}

            self._Nmonte = None
            self._Nnetworks = None
            self._sampleBeta = None
            self._sampleMontecarlo = None

            self._name = None
            self._problogstring = None

            self._is_this_a_bn = None

            self._times = {}

    def _backward_comp_results(self):
        if not hasattr(self, "_results"):
            self._results = None

        if self._results is None:
            self._results = {}

        if hasattr(self, "_vec_real") and self.GROUND not in self._results:
            self._results[self.GROUND] = self._vec_real

        if hasattr(self, "_vec_sl_beta_cov") and self.BETAPROBLOG not in self._results:
            self._results[self.BETAPROBLOG] = self._vec_sl_beta_cov

        if hasattr(self, "_vec_sl_beta") and self.HONOLULU not in self._results:
            self._results[self.HONOLULU] = self._vec_sl_beta

        if hasattr(self, "_vec_sl") and self.JOSANG not in self._results:
            self._results[self.JOSANG] = self._vec_sl

        if hasattr(self, "_vec_sbn") and self.SBN not in self._results:
            self._results[self.SBN] = self._vec_sbn

        if hasattr(self, "_vec_GBT") and self.GBT not in self._results:
            self._results[self.GBT] = self._vec_GBT

        if hasattr(self, "_vec_credal") and self.CREDAL not in self._results:
            self._results[self.CREDAL] = self._vec_credal


    def setup(self, name, content, Nmonte = 10, Nnetworks = 100, sampleBeta = [10], bn=True):
        """
        Storage of attributes
        :param name: name of this experiment: anything
        :param content: file with Bayesian network specified or string otherwise
        :param Nmonte: how often we want to randomly assign probabilities to the template problog string
        :param Nnetworks: how many networks we want to generate
        :param sampleBeta: how many samples to use for create SL opinions
        :param bn: is this a Bayesian network?
        :return:
        """

        self._Nmonte = Nmonte
        self._Nnetworks = Nnetworks
        self._sampleBeta = sampleBeta

        self._name = name
        self._problogstring = None

        self._is_this_a_bn = bn

        self._is_this_a_bn = bn
        if bn:
            self.net = Graph()

            with open(content) as f:
                for line in f:
                    self.net.add_edge(line.rstrip('\n').split(" "))
            self.net.storeFilename(content)
            self._problogstring = self.net.get_problog_string()
        else:
            if not isinstance(content,str):
                raise Exception("Todo")

            self.net = None

            self._problogstring = content

    def _expected_value(self,l):
        """
        Returns the expected value of a SL opinion
        """
        return l[0] + l[2] * l[3]

    def _computer_vector_error(self, one, two):
        res = []
        items = 0
        if len(one) != len(two):
            print(len(one))
            print(len(two))
            raise Exception("wrong length")

        for i in range(len(one)):
            for k, x in one[i].items():
                w = two[i][k]

                p = None
                if isinstance(x, list):
                    p = self._expected_value(x)
                else:
                    p = x

                if isinstance(w, list) or isinstance(w, tuple):
                    res.append((p - float(self._expected_value(w))) ** 2)
                else:
                    res.append((p - float(w)) ** 2)
                items += 1
        return res

    def _compute_error(self, one, two):
        """
        Compute the distance between the two values
        """
        res = 0
        items = 0
        if len(one) != len(two):
            print(len(one))
            print(len(two))
            raise Exception("wrong length")

        for i in range(len(one)):
            for k, x in one[i].items():
                w = two[i][k]

                p = None
                if isinstance(x, list):
                    p = self._expected_value(x)
                else:
                    p = x

                if isinstance(w, list) or isinstance(w, tuple):
                    res += (p - float(self._expected_value(w))) ** 2
                else:
                    res += (p - float(w)) ** 2
                items += 1
        return math.sqrt(float(res) / float(items))

    def _expected_error(self, listw):
        """
        Computed the expected error of a SL opinion
        """
        res = 0
        items = 0
        for w in listw:
            for k, v in w.items():
                items += 1
                res += (float(self._expected_value(v)) * (1.0-float(self._expected_value(v))) * float(v[2]) / (2.0 + float(v[2])))

        return math.sqrt(float(res) / float(items))

    def _counting(self, message, now, max):
        sys.stdout.write("\r%s: %d%%" % (message, int(now / max * 100)))
        sys.stdout.flush()

    def run_choices(self, algs = None):
        if algs == None:
            raise Exception("Choose algorithms")

        self._backward_comp_results()

        if self.bns is None:
            self.bns = []

        if not self.bns or not self._sb: # creating programs
            Nruns = self._Nmonte * self._Nnetworks
            self._times[self.GROUND] = []

            for i in range(Nruns):

                self._counting("Computing the ground truth", i, Nruns)


                b = None
                if b is None or i % self._Nmonte:
                    b = ProbProblog(self._problogstring, self.net)
                    self.bns.append(b)

                # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                begin = time.time()
                ground = b.run()
                end = time.time()
                self._times[self.GROUND].append(end-begin)

                # print("Ground")
                # print(ground)
                self._vec_real.append(ground)
                # print("\n")

                for samples in self._sampleBeta:
                    sb = DistProbLog(b, samples)
                    self._sb.append(sb)

        for alg in algs:
            print("")
            i = 0
            if alg == self.BETAPROBLOG:
                self._results[self.BETAPROBLOG] = []
                self._times[self.BETAPROBLOG] = []
                for sb in self._sb:
                    self._counting("Computing BetaProblog", i, len(self._sb))
                    i+= 1
                    begin = time.time()
                    lance = SLProbLog(sb.get_program_beta(), True).run_beta_cov()
                    end = time.time()
                    self._times[self.BETAPROBLOG].append(end-begin)
                    self._results[self.BETAPROBLOG].append(lance)
            # elif alg == self.PYBETAPROBLOG:
            #     self._results[self.PYBETAPROBLOG] = []
            #     for sb in self._sb:
            #         self._counting("Computing PyBetaProblog", i, len(self._sb))
            #         i+= 1
            #         me = SLProbLog(sb.get_program_beta(), True).myrun()
            #         self._results[self.PYBETAPROBLOG].append(me)
            elif alg == self.HONOLULU:
                self._results[self.HONOLULU] = []
                self._times[self.HONOLULU] = []
                for sb in self._sb:
                    self._counting("Computing Honolulu", i, len(self._sb))
                    i += 1
                    begin = time.time()
                    honolulu = SLProbLog(sb.get_program_beta(), True).run_beta()
                    end = time.time()
                    self._times[self.HONOLULU].append(end-begin)
                    self._results[self.HONOLULU].append(honolulu)

                    # lastid = len(self._results[self.HONOLULU]) - 1
                    # nodiff = True
                    # for k in self._results[self.HONOLULU][lastid]:
                    #     whw = self._results[self.HONOLULU][lastid][k]
                    #     wbw = self._results[self.BETAPROBLOG][lastid][k]
                    #     if abs(whw[0]+0.5*whw[2]-(wbw[0]+0.5*wbw[2])) > 1e-10:
                    #         print(str(k) + ' ' + str(whw) + " " + str(wbw))
                    #         nodiff = False
                    #
                    # if not nodiff:
                    #     print(sb.get_program_beta())

            elif alg == self.SAMPLESIMP:
                self._results[self.SAMPLESIMP] = []
                self._times[self.SAMPLESIMP] = []
                instance = 0
                for sb in self._sb:
                    self._counting("Computing Simple Sampling", i, len(self._sb))
                    i += 1
                    begin = time.time()
                    samplesimp = SLProbLog(sb.get_program(), True).run_sample_cov(samples=100)
                    end = time.time()
                    self._times[self.SAMPLESIMP].append(end-begin)
                    self._results[self.SAMPLESIMP].append(samplesimp)
                    #
                    # lan = self._results[self.BETAPROBLOG][instance]
                    # samp = self._results[self.SAMPLESIMP][instance]
                    #
                    # for k, v in lan.items():
                    #
                    #     meanv = v[0]+0.5*v[2]
                    #     meansamp = samp[k][0] + 0.5*samp[k][2]
                    #
                    #     varv = meanv * (1- meanv) / (2/v[2] + 1)
                    #     varsamp = meansamp * (1-meansamp) / (2/samp[k][2] + 1)
                    #
                    #     errorvar = abs(varv - varsamp)
                    #
                    #     if  errorvar > 0.01:
                    #         print("instance: %d; query: %s; betaproblog: %s; sample: %s; error: %s " % (instance, k, str(varv), str(varsamp), str(errorvar)))

                    instance += 1

            elif alg == self.JOSANG:
                self._results[self.JOSANG] = []
                self._times[self.JOSANG] = []
                for sb in self._sb:
                    self._counting("Computing Josang", i, len(self._sb))
                    i += 1
                    begin = time.time()
                    josang = SLProbLog(sb.get_program(), True).run_SL()
                    end = time.time()
                    self._times[self.JOSANG].append(end-begin)
                    self._results[self.JOSANG].append(josang)
            elif alg == self.SBN:
                if not self._is_this_a_bn:
                    raise Exception("Cannot run SBN as this is not a Bayesian Network")
                self._results[self.SBN] = []
                self._times[self.SBN] = []
                for sb in self._sb:
                    self._counting("Computing SBN", i, len(self._sb))
                    i += 1
                    begin = time.time()
                    subjbaynet = sb.runSBN()
                    end = time.time()
                    self._times[self.SBN].append(end-begin)
                    self._results[self.SBN].append(subjbaynet)
            elif alg == self.GBT:
                if not self._is_this_a_bn:
                    raise Exception("Cannot run GBT as this is not a Bayesian Network")
                self._results[self.GBT] = []
                self._times[self.GBT] = []
                for sb in self._sb:
                    self._counting("Computing GBT", i, len(self._sb))
                    i += 1
                    begin = time.time()
                    gbt = sb.runGBT()
                    end = time.time()
                    self._times[self.GBT].append(end-begin)
                    self._results[self.GBT].append(gbt)
            elif alg == self.CREDAL:
                if not self._is_this_a_bn:
                    raise Exception("Cannot run CREDAL as this is not a Bayesian Network")
                self._results[self.CREDAL] = []
                self._times[self.CREDAL] = []
                for sb in self._sb:
                    self._counting("Computing Credal", i, len(self._sb))
                    i += 1
                    begin = time.time()
                    credal = sb.runCredal()
                    end = time.time()
                    self._times[self.CREDAL].append(end-begin)
                    self._results[self.CREDAL].append(credal)
            else:
                raise Exception("%s unknown" % alg)

            self._store()



    def run(self, oneevidence = False):
        """
        Run the experiment with the given setup
        """
        self._vec_real = []
        self._vec_sl = []
        self._vec_sl_beta = []
        self._vec_sl_beta_cov = []
        self._vec_sbn = []
        self._vec_GBT = []
        self._vec_credal = []
        self.bns = []
        self._sb = []

        Nruns = self._Nmonte * self._Nnetworks

        for i in range(Nruns):

            sys.stdout.write("\r%d%%" % int(i/Nruns*100))
            sys.stdout.flush()

            b = None
            if b is None or i % self._Nmonte:
                b = ProbProblog(self._problogstring, self.net)
                self.bns.append(b)

            #print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            ground = b.run()
            #print("Ground")
            #print(ground)
            self._vec_real.append(ground)
            #print("\n")

            for samples in self._sampleBeta:
                sb = DistProbLog(b, samples)
                self._sb.append(sb)

                # f = open("programs/"+self._name+str(i).zfill(4), "w+")
                # f.write(sb.get_program())
                # f.close()

                self._vec_sl.append(SLProbLog(sb.get_program(), True).run_SL())

                honolulu = SLProbLog(sb.get_program_beta(), True).run_beta()
                #print("Honolulu")
                #print(honolulu)
                self._vec_sl_beta.append(honolulu)
                #print("\n")

                if oneevidence:
                    #print("Lance")
                    lance = SLProbLog(sb.get_program_beta(), True).run_KL(MatlabPath="/../SLProbLog/SLengineMatlab/")
                    #print(lance)
                    #print("\n")
                    self._vec_sl_beta_cov.append(lance)

                if self._is_this_a_bn:
                    #print("SBN")
                    subjbaynet = sb.runSBN()
                    #print(subjbaynet)
                    #print("\n")
                    self._vec_sbn.append(subjbaynet)

                    for k in lance:
                        if abs(self._expected_value(lance[k]) - self._expected_value(subjbaynet[k])) > 0.1:
                            print("error %d" % i)
                            f = open("programs/error" + self._name + str(i).zfill(4), "w+")
                            f.write(sb.get_program())
                            f.close()

                    self._vec_credal.append(sb.runCredal())
                    self._vec_GBT.append(sb.runGBT())

                #print("\n\n")


        print("")
        self._store()


    def _store(self):
        """
        Save the pickle file
        """
        now = datetime.datetime.now()
        self._filename = self._name + "-%s-%s-%s-%s-%s" % (now.year, now.month, now.day, now.hour, now.minute)

        pickle.dump(self, open(self._filename + ".pickle", "wb"))


    def perfbound_beta(self, p, w):

        VERYLARGE = 1e10
        VERYSMALL = 1e-1

        vp = np.array(p)
        vw = np.array(w)

        alpha1 = mpmath.mpf(1) * 2*vw[:,0]/vw[:,2]
        alpha2 = mpmath.mpf(1) * 2*vw[:,1]/vw[:,2]

        todelete = []
        for i in range(len(alpha1)):
            # if np.isinf(alpha1[i]):
            #     alpha1[i] = VERYLARGE
            # if np.isinf(alpha2[i]):
            #     alpha2[i] = VERYLARGE

            if alpha1[i] <= VERYSMALL:
                alpha1[i] == VERYSMALL

            if alpha2[i] <= VERYSMALL:
                alpha2[i] == VERYSMALL


        vp = np.delete(vp, todelete)
        vw = np.delete(vw, todelete, axis=0)
        alpha1 = np.delete(alpha1, todelete)
        alpha2 = np.delete(alpha2, todelete)


        pe = vw[:,0] + vw[:,2] * 0.5



        er = vp - pe
        pp = np.arange(0,1,0.01)
        pa = np.zeros(len(pp))

        mn = alpha1 / (alpha1 + alpha2)

        lb = np.zeros(len(er))
        ub = np.zeros(len(er))

        loc1 = np.where(np.logical_or(np.logical_and(alpha1 > 1, alpha2 > 1),np.logical_and(alpha1 == 1, alpha2 == 1)))[0]
        loc2 = np.where(np.logical_and(alpha1 < 1, alpha2 < 1))[0]
        loc3 = np.where(np.logical_and(alpha1 <= 1, alpha2 > 1))[0]
        lb[loc3] = -1 * mn[loc3]
        loc4 = np.where(np.logical_and(alpha1 > 1, alpha2 <= 1))[0]
        ub[loc4] =  1 - mn[loc4]


        for i in range(len(pp)):
            pp2 = (1 - pp[i])/2 * np.ones(len(loc1))
            lb[loc1] = sc.betaincinv(alpha1[loc1], alpha2[loc1], pp2) - mn[loc1]
            ub[loc1] = sc.betaincinv(alpha1[loc1], alpha2[loc1],  1 - pp2) - mn[loc1]

            if len(loc2) > 0:
                pp2 = pp[i]/2 * np.ones(len(loc2))
                lb[loc2] = sc.betaincinv(alpha1[loc2], alpha2[loc2], pp2) - mn[loc2]
                ub[loc2] = sc.betaincinv(alpha1[loc2], alpha2[loc2], 1- pp2) - mn[loc2]
            if len(loc3) > 0:
                pp2 = pp[i] * np.ones(len(loc3))
                ub[loc3] = sc.betaincinv(alpha1[loc3], alpha2[loc3], pp2) - mn[loc3]
            if len(loc4) > 0:
                pp2 = (1-pp[i]) * np.ones(len(loc4))
                lb[loc4] = sc.betaincinv(alpha1[loc4], alpha2[loc4], pp2) - mn[loc4]

            vin = np.zeros(len(er))
            #notnanidx = np.where(np.logical_and(~np.isnan(lb), ~np.isnan(ub)))[0]
            #vin[np.where(np.logical_and(er[notnanidx] <= ub[notnanidx], er[notnanidx] >= lb[notnanidx]))[0]] = 1
            vin[np.where(np.logical_and(er <= ub, er >= lb))[0]] = 1
            vin[loc2] = 1 - vin[loc2]
            pa[i] = np.sum(vin) / len(er)

        return [pa, pp]


    def _getfilename(self):
        return os.path.basename(self._filename)

    def analise(self, graphs = True, printheader = False):
        """
        Analyse the results and print them to screen
        """

        self._backward_comp_results()

        pspace="\\vphantom{$\\mathcal{S}^{\\beta}\\mathcal{S}_{\\text{SL}}$}"
        labels = {
            self.BETAPROBLOG: "{\\fbox{\\color{cBetaProblog}CPB%s}}" % pspace,
            self.HONOLULU: "{\\color{cHonolulu} \\fbox{$\\mathcal{S}^{\\beta}$%s}}" % pspace,
            self.JOSANG: "{\\color{cJosang}\\fbox{$\\mathcal{S}_{\\text{SL}}$%s}}" %pspace,
            self.SBN: "{\\color{cSBN}\\fbox{SBN%s}}" % pspace,
            self.GBT: "{\\color{cGBT}\\fbox{GBT%s}}" %pspace,
            self.CREDAL: "{\color{cCredal}\\fbox{Credal%s}}" % pspace,
            self.SAMPLESIMP: "{\\color{cMonteCarlo}\\fbox{MC%s}}" %pspace,
            self.PYBETAPROBLOG: "\\fbox{PySLProbLog%s}" %pspace
        }

        coloursdef = {
            self.BETAPROBLOG: "#0000FF",
            self.PYBETAPROBLOG: "#0000FF",
            self.HONOLULU: "#000000",
            self.JOSANG: "#FF0000",
            self.SBN: "#BF00BF",
            self.GBT: "#008000",
            self.CREDAL: "#BFBF00",
            self.SAMPLESIMP: "#B35F00"
        }

        objects = []
        performance = []
        colours = []

        header =  "& &   & "
        #strreal = "& & A & "
        #strpred = "& & P & "

        vreal = []
        vpred = []

        errorbox = []
        timesbox = []
        timesboxlabel = []

        if self._results[self.GROUND] is None or self._results[self.JOSANG] is None or self._results[self.HONOLULU] is None:
            raise Exception("Run first")

        if self.BETAPROBLOG in self._results and len(self._results[self.BETAPROBLOG]) > 0:
            header += " %s & " % labels[self.BETAPROBLOG]
            objects.append(labels[self.BETAPROBLOG])
            colours.append(coloursdef[self.BETAPROBLOG])
            performance.append(self._compute_error(self._results[self.GROUND], self._results[self.BETAPROBLOG]))
            errorbox.append(self._computer_vector_error(self._results[self.GROUND], self._results[self.BETAPROBLOG]))
            timesbox.append(self._times[self.BETAPROBLOG])
            timesboxlabel.append(labels[self.BETAPROBLOG])

            realerr = self._compute_error(self._results[self.GROUND], self._results[self.BETAPROBLOG])
            prederr = self._expected_error(self._results[self.BETAPROBLOG])

            vreal.append(realerr)
            vpred.append(prederr)

            # if self.PYBETAPROBLOG in self._results and len(self._results[self.PYBETAPROBLOG]) > 0:
            #     header += " %s & " % labels[self.PYBETAPROBLOG]
            #     objects.append(labels[self.PYBETAPROBLOG])
            #     colours.append(coloursdef[self.PYBETAPROBLOG])
            #     performance.append(self._compute_error(self._results[self.GROUND], self._results[self.PYBETAPROBLOG]))
            #
            #     realerr = self._compute_error(self._results[self.GROUND], self._results[self.PYBETAPROBLOG])
            #     prederr = self._expected_error(self._results[self.PYBETAPROBLOG])
            #
            #     vreal.append(realerr)
            #     vpred.append(prederr)

            #strreal += "%.4f & " %
            #strpred += "%.4f & " %

        if self.HONOLULU in self._results and len(self._results[self.HONOLULU]) > 0:
            header += "%s & " % labels[self.HONOLULU]
            objects.append(labels[self.HONOLULU])
            colours.append(coloursdef[self.HONOLULU])
            performance.append(self._compute_error(self._results[self.GROUND], self._results[self.HONOLULU]))
            errorbox.append(self._computer_vector_error(self._results[self.GROUND], self._results[self.HONOLULU]))
            #timesbox.append(self._times[self.HONOLULU])
            #timesboxlabel.append(labels[self.HONOLULU])

            realerr = self._compute_error(self._results[self.GROUND], self._results[self.HONOLULU])
            prederr = self._expected_error(self._results[self.HONOLULU])
            vreal.append(realerr)
            vpred.append(prederr)


        if self.JOSANG in self._results and len(self._results[self.JOSANG]) > 0:
            header += "%s & " % labels[self.JOSANG]
            objects.append(labels[self.JOSANG])
            colours.append(coloursdef[self.JOSANG])
            performance.append(self._compute_error(self._results[self.GROUND], self._results[self.JOSANG]))
            errorbox.append(self._computer_vector_error(self._results[self.GROUND], self._results[self.JOSANG]))
            #timesbox.append(self._times[self.JOSANG])
            #timesboxlabel.append(labels[self.JOSANG])

            realerr = self._compute_error(self._results[self.GROUND], self._results[self.JOSANG])
            prederr = self._expected_error(self._results[self.JOSANG])
            vreal.append(realerr)
            vpred.append(prederr)

        if self.SAMPLESIMP in self._results and len(self._results[self.SAMPLESIMP]) > 0:
            header += "%s & " % labels[self.SAMPLESIMP]
            objects.append(labels[self.SAMPLESIMP])
            colours.append(coloursdef[self.SAMPLESIMP])
            performance.append(self._compute_error(self._results[self.GROUND], self._results[self.SAMPLESIMP]))
            errorbox.append(self._computer_vector_error(self._results[self.GROUND], self._results[self.SAMPLESIMP]))
            timesbox.append(self._times[self.SAMPLESIMP])
            timesboxlabel.append(labels[self.SAMPLESIMP])

            realerr = self._compute_error(self._results[self.GROUND], self._results[self.SAMPLESIMP])
            prederr = self._expected_error(self._results[self.SAMPLESIMP])
            vreal.append(realerr)
            vpred.append(prederr)



        if self._is_this_a_bn:
            if self.SBN in self._results and len(self._results[self.SBN]) > 0:
                header += "%s & " % labels[self.SBN]
                objects.append(labels[self.SBN])
                colours.append(coloursdef[self.SBN])
                performance.append(self._compute_error(self._results[self.GROUND], self._results[self.SBN]))
                errorbox.append(
                    self._computer_vector_error(self._results[self.GROUND], self._results[self.SBN]))
                #timesbox.append(self._times[self.SBN])

                realerr = self._compute_error(self._results[self.GROUND], self._results[self.SBN])
                prederr = self._expected_error(self._results[self.SBN])
                vreal.append(realerr)
                vpred.append(prederr)

            if self.GBT in self._results and len(self._results[self.GBT]) > 0:
                header += "%s & " % labels[self.GBT]
                objects.append(labels[self.GBT])
                colours.append(coloursdef[self.GBT])
                performance.append(self._compute_error(self._results[self.GROUND], self._results[self.GBT]))
                errorbox.append(
                    self._computer_vector_error(self._results[self.GROUND], self._results[self.GBT]))
                #timesbox.append(self._times[self.GBT])

                realerr = self._compute_error(self._results[self.GROUND], self._results[self.GBT])
                prederr = self._expected_error(self._results[self.GBT])
                vreal.append(realerr)
                vpred.append(prederr)

            if self.CREDAL in self._results and len(self._results[self.CREDAL]) > 0:
                header += "%s &" % labels[self.CREDAL]
                objects.append(labels[self.CREDAL])
                colours.append(coloursdef[self.CREDAL])
                performance.append(self._compute_error(self._results[self.GROUND], self._results[self.CREDAL]))
                errorbox.append(
                    self._computer_vector_error(self._results[self.GROUND], self._results[self.CREDAL]))
                #timesbox.append(self._times[self.CREDAL])

                realerr = self._compute_error(self._results[self.GROUND], self._results[self.CREDAL])
                prederr = self._expected_error(self._results[self.CREDAL])
                vreal.append(realerr)
                vpred.append(prederr)


        header = header[:-2]
        header += "\n"

        #strreal = strreal[:-2]
        #strreal += "\\\\"

        #strpred = strpred[:-2]
        #strpred += "\\\\"


        if printheader:
            print(header)


        v = self._name.split('-')
        n = "\\%s{}" % v[0]
        s = None

        for i in range(1,len(v)):
            if v[i] == 'Nins':
                s = v[i+1]
                break
            if 'Nins' in v[i]:
                s = v[i][4:]

        maxr = min(vreal)
        posmaxr = [i for i,j in enumerate(vreal) if j == maxr]

        strreal = "%s & %s & A " % (n, s)
        for i in range(len(vreal)):
            if i in posmaxr:
                strreal += "& \\best{%.4f} " % vreal[i]
            else:
                strreal += "& %.4f " % vreal[i]

        strreal += "\\\\"

        strpred = " &  & P "
        for vp in vpred:
            strpred += "& %.4f " %  vp

        strpred += "\\\\"

        print(strreal)
        print(strpred)
        print("\\dashedline")

        y_pos = np.arange(len(objects))

        matplotlib.rcParams['text.usetex'] = True
        params = {'text.latex.preamble': [r'\usepackage{amsmath} \usepackage{xcolor} \definecolor{cBetaProblog}{HTML}{0000FF} \definecolor{cHonolulu}{HTML}{000000} \definecolor{cJosang}{HTML}{FF0000} \definecolor{cSBN}{HTML}{BF00BF} \definecolor{cGBT}{HTML}{008000} \definecolor{cCredal}{HTML}{BFBF00} \definecolor{cMonteCarlo}{HTML}{CC921F}']}
        font = {'family': 'normal',
                'weight': 'normal',
                'size': 10}
        plt.rc('font', **font)
        plt.rcParams.update(params)
        f = plt.figure(figsize=(4,4))
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=True)  # labels along the bottom edge are off
        ax = plt.gca()
        ax.set_ylim([0,12])
        #axes.set_xlim([xmin, xmax])
        #axes.set_ylim([0, 1])
        #plt.ylim(bottom=0.5)
        #ax.set_ylim([0, 0.1])
        b = plt.violinplot(timesbox, [(2 * (x + 1)) for x in range(len(timesbox))],
                   showmeans=False,
                   showmedians=True)

        ax.grid(color="#cdcdcd", linewidth=0.3)
        ax.yaxis.grid(True)
        ax.xaxis.grid(False)
        ax.set_xticks([(2 * (x + 1)) for x in range(len(timesbox))])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)


        # for patch, color in zip(b['boxes'], colours):
        #     patch.set_facecolor(color)

        # def autolabel(rects):
        #     """
        #     Attach a text label above each bar displaying its height
        #     """
        #     for rect in rects:
        #         height = rect.get_height()
        #         ax.text(rect.get_x() + rect.get_width() / 2., 1.005 * height,
        #                 '%.04f' % height,
        #                 ha='center', va='bottom')

        # autolabel(b)
        #
        # h = []
        # xt = []
        # for i in range(len(colours)):
        #     h.append(mpatches.Patch(color=colours[i], label=objects[i]))
        #     xt.append(objects[i])
        #
        # #plt.xticks(np.arange(6), xt)
        # (b,t) = plt.ylim()
        # plt.ylim(top=t+0.1)
        # ax.legend(handles=h, loc='upper right')
        # #plt.xlabel('Desired Confidence')

        plt.setp(ax, xticks=[(2 * (x + 1)) for x in range(len(timesbox))],
                 xticklabels=timesboxlabel)

        plt.ylabel('Execution Times')
        #plt.grid()
        f.savefig(self._getfilename() + "-times.pdf", bbox_inches='tight')

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off

        #plt.xticks(y_pos, objects)

        if graphs:

            octave.addpath(os.getcwd() + MATLABPATH)

            psl_beta_cov = None
            if self.BETAPROBLOG in self._results and len(self._results[self.BETAPROBLOG]) > 0:
                psl_beta_cov = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                                 octave.double(self._toMatlabArray(self._results[self.BETAPROBLOG])))

            psl_beta_cov_py = None
            if self.PYBETAPROBLOG in self._results and len(self._results[self.PYBETAPROBLOG]) > 0:
                psl_beta_cov_py = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                                 octave.double(self._toMatlabArray(self._results[self.PYBETAPROBLOG])))

            psl = None
            if self.JOSANG in self._results and len(self._results[self.JOSANG]) > 0:
                psl = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                        octave.double(self._toMatlabArray(self._results[self.JOSANG])))

            psl_beta = None
            if self.HONOLULU in self._results and len(self._results[self.HONOLULU]) > 0:
                psl_beta = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                          octave.double(self._toMatlabArray(self._results[self.HONOLULU])))

            psl_sample = None
            if self.SAMPLESIMP in self._results and len(self._results[self.SAMPLESIMP]) > 0:
                psl_sample = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                        octave.double(self._toMatlabArray(self._results[self.SAMPLESIMP])))

            psbn = None
            pcredal = None
            pgbt = None
            if self._is_this_a_bn:
                psbn = None
                if self.SBN in self._results and len(self._results[self.SBN]) > 0:
                    psbn = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                         octave.double(self._toMatlabArray(self._results[self.SBN])))

                pcredal = None
                if self.CREDAL in self._results and len(self._results[self.CREDAL]) > 0:
                    pcredal = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                          octave.double(self._toMatlabArray(self._results[self.CREDAL])))

                pgbt = None
                if self.GBT in self._results and len(self._results[self.GBT]) > 0:
                    pgbt = octave.perfbound_beta(octave.double(self._toMatlabArray(self._results[self.GROUND])),
                                             octave.double(self._toMatlabArray(self._results[self.GBT])))


            pp = list(numpy.arange(0,1.0,0.01))

            matplotlib.rcParams['text.usetex'] = True
            f = plt.figure()
            markers_on = 10
            if psl_beta_cov is not None:
                plt.plot(pp, psl_beta_cov[0], ':D', color = '%s' % coloursdef[self.BETAPROBLOG], markevery=markers_on, label=labels[self.BETAPROBLOG])

            if psl_beta_cov_py is not None:
                plt.plot(pp, psl_beta_cov_py[0], ':p', color = '%s' % coloursdef[self.PYBETAPROBLOG], markevery=markers_on, label=labels[self.PYBETAPROBLOG])

            if psl_beta is not None:
                plt.plot(pp, psl_beta[0], ':o', color = '%s' % coloursdef[self.HONOLULU], markevery=markers_on, label=labels[self.HONOLULU])

            if psl is not None:
                plt.plot(pp, psl[0], ':p', color = '%s' % coloursdef[self.JOSANG], markevery=markers_on, label=labels[self.JOSANG])

            if psl_sample is not None:
                plt.plot(pp, psl_sample[0], ':p', color = '%s' % coloursdef[self.SAMPLESIMP], markevery=markers_on, label=labels[self.SAMPLESIMP])

            if self._is_this_a_bn:
                if psbn is not None:
                    plt.plot(pp, psbn[0], ':*', color = '%s' % coloursdef[self.SBN], markevery=markers_on, label=labels[self.SBN])

                if pgbt is not None:
                    plt.plot(pp, pgbt[0], ':v', color = '%s' % coloursdef[self.GBT], markevery=markers_on, label=labels[self.GBT])

                if pcredal is not None:
                    plt.plot(pp, pcredal[0], ':^', color = '%s' % coloursdef[self.CREDAL], markevery=markers_on, label=labels[self.CREDAL])
            plt.legend()
            plt.xlabel('Desired Confidence')
            plt.ylabel('Actual Confidence')
            plt.grid()
            f.savefig(self._getfilename() + ".pdf", bbox_inches='tight')


    def _toMatlabArray(self, l):
        ret = []
        for el in l:
            for k, v in el.items():
                if isinstance(v, list) or isinstance(v, tuple):
                    w = v
                    ret.append([float(w[0]), float(w[1]), float(w[2])])
                else:
                    ret.append([v])
        return ret

