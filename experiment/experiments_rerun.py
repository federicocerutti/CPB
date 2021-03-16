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

# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.1507 & 0.1507 & 0.2121 & 0.1501 & 0.1507 & 0.1533 & 0.1610 \\
# & & P & 0.1465 & 0.1895 & 0.1661 & 0.1458 & 0.1466 & 0.0872 & 0.2008 \\
# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.0776 & 0.0776 & 0.1144 & 0.0779 & 0.0776 & 0.0801 & 0.0786 \\
# & & P & 0.0759 & 0.1060 & 0.0820 & 0.0749 & 0.0753 & 0.0378 & 0.1025 \\
# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.0559 & 0.0559 & 0.0851 & 0.0560 & 0.0559 & 0.0610 & 0.0563 \\
# & & P & 0.0565 & 0.0785 & 0.0610 & 0.0558 & 0.0560 & 0.0262 & 0.0760 \\



import sys
sys.path.append('../')

from experiment.experimental_setting import Experiment

files = [
    "net1-Nins-10-2020-4-20-20-11.pickle",
	"net1-Nins-50-2020-4-20-20-41.pickle",
	"net1-Nins-100-2020-4-20-21-11.pickle",
	"net2-Nins-10-2020-4-20-21-51.pickle",
	"net2-Nins-50-2020-4-20-22-31.pickle",
	"net2-Nins-100-2020-4-20-23-10.pickle",
	"net3-Nins-10-2020-4-21-0-11.pickle",
	"net3-Nins-50-2020-4-21-1-12.pickle",
	"net3-Nins-100-2020-4-21-2-14.pickle",
	"smoker-Nins10-2020-4-21-3-4.pickle",
	"smoker-Nins50-2020-4-21-3-53.pickle",
	"smoker-Nins100-2020-4-21-4-43.pickle"
]

clusterrun = int(sys.argv[1])

e = Experiment()
e = Experiment.loadExperiment(files[clusterrun])
e.run_choices(
        [Experiment.SAMPLESIMP])
# if "smorker" in files[clusterrun]:
#     e.run_choices(
#         [Experiment.BETAPROBLOG, Experiment.SAMPLESIMP, Experiment.HONOLULU, Experiment.JOSANG])
# else:
#     e.run_choices(
#         [Experiment.BETAPROBLOG, Experiment.SAMPLESIMP, Experiment.HONOLULU, Experiment.JOSANG, Experiment.SBN,
#          Experiment.GBT, Experiment.CREDAL])

#
# for f in files:
#     e = Experiment.loadExperiment(f)
#     e.run_choices(
#         [Experiment.HONOLULU])
#     #e.analise()
#
# e = Experiment.loadExperiment("net1-Nins-10-cov-2019-6-15-16-30.pickle")
# # e = Experiment.loadExperiment("net1-test-20190620-2019-7-9-23-1.pickle")
# # e._name = "net1-test-20190620"
# # e.run_choices([Experiment.SAMPLESIMP])
# e.analise()

#
# e = Experiment.loadExperiment("net1-Nins-50-cov-2019-6-15-18-7.pickle")
# e.run_choices([Experiment.SAMPLESIMP])
# e.analise()
#
# e = Experiment.loadExperiment("net1-Nins-100-cov-2019-6-15-19-39.pickle")
# e.run_choices([Experiment.SAMPLESIMP])
# e.analise()
#
#
#
# e.setup("net1-Nins-10-cov", fname, 10, 100, [10])
# e.run(oneevidence=True)
# e.analise()
#
# print("")
#
# e.setup("net1-Nins-50-cov", fname, 10, 100, [50])
# e.run(oneevidence=True)
# e.analise()
#
# print("")
#
# e.setup("net1-Nins-100-cov", fname, 10, 100, [100])
# e.run(oneevidence=True)
# #e = Experiment.loadExperiment("net1-Nins-100-cov-2019-6-7-13-59.pickle")
# e.analise()
