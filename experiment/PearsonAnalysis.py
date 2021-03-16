import sys
sys.path.append('../')
import pickle
import scipy.stats
from SLProbLog.SLProbLog import from_sl_opinion
from SLProbLog.SLProbLog import SLProbLog

def loadfile(f):
    return pickle.load(open(f, "rb"))

class PearsonAnalysis:
    def __init__(self, filename):
        self._filename = filename
        self.simulatedgroundtruth = []
        self.betaproblogres = []
        self.iter_samples = []
        self.pbetaproblog = []
        self.psimulatedground = []
        self.pitermonte = {}
        self.pitermonteground = {}
        self.corrbeta = []
        self.corrmonte = {}

    def run_analysis(self):
        exp = pickle.load(open(self._filename, "rb"))
        self.betaproblogres = exp._results[exp.BETAPROBLOG]
        for idx in range(len(exp._sb)):
            spl = SLProbLog(exp._sb[idx].get_program_beta())
            self.simulatedgroundtruth.append(spl.run_sample(samples=10000))
            self.iter_samples.append(spl.run_sample(samples=200,progression=True))

        for idx in range(len(self.simulatedgroundtruth)):
            for q in exp._results[exp.BETAPROBLOG][0]:
                self.pbetaproblog.append(float(from_sl_opinion(exp._results[exp.BETAPROBLOG][idx][q]).strength()))
                self.psimulatedground.append(float(self.simulatedgroundtruth[idx][q].strength()))

            for iteridx in self.iter_samples[idx]:
                if iteridx not in self.pitermonte:
                    self.pitermonte[iteridx] = []
                if iteridx not in self.pitermonteground:
                    self.pitermonteground[iteridx] = []

                for q in exp._results[exp.BETAPROBLOG][idx]:
                    self.pitermonte[iteridx].append(float(self.iter_samples[idx][iteridx][q].strength()))
                    self.pitermonteground[iteridx].append(float(self.simulatedgroundtruth[idx][q].strength()))

        self.corrbeta = scipy.stats.pearsonr(self.pbetaproblog, self.psimulatedground)

        self.corrmonte = {}
        for iteridx in self.iter_samples[0]:
            self.corrmonte[iteridx] = scipy.stats.pearsonr(self.pitermonte[iteridx], self.pitermonteground[iteridx])

        self.save()


    def save(self):
        pickle.dump(self, open("corr-" + self._filename + ".pickle", "wb"))


if __name__ == '__main__':
    p = PearsonAnalysis(sys.argv[1])
    p.run_analysis()


# #picklefile = "net1-Nins-10-cov-2019-6-11-11-19.pickle"
# picklefile = "net3-Nins-100-cov-2019-7-10-11-43.pickle"
# picklefile = "smoker-Nins100-2019-8-6-13-22.pickle"
# picklefile = "net1-Nins-10-2019-12-11-13-33.pickle"
#
# p = PearsonAnalysis(picklefile)
# p.run_analysis()
# print(p.corrbeta)
# print(p.corrmonte)