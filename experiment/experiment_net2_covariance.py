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
#
#
# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.1424 & 0.1424 & 0.1941 & 0.1421 & 0.1424 & 0.1441 & 0.1509 \\
# & & P & 0.1375 & 0.1686 & 0.1477 & 0.1366 & 0.1382 & 0.1045 & 0.1845 \\
# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.0734 & 0.0734 & 0.1117 & 0.0741 & 0.0734 & 0.0751 & 0.0739 \\
# & & P & 0.0734 & 0.0931 & 0.0789 & 0.0724 & 0.0732 & 0.0485 & 0.0960 \\
# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.0515 & 0.0515 & 0.0807 & 0.0518 & 0.0515 & 0.0538 & 0.0521 \\
# & & P & 0.0529 & 0.0666 & 0.0572 & 0.0522 & 0.0528 & 0.0343 & 0.0694 \\

import sys
sys.path.append('../')

from experiment.experimental_setting import Experiment


fname="networks/net2.file"

e = Experiment()

e = Experiment.loadExperiment("net2-Nins-10-cov-2019-6-15-18-39.pickle")
e.run_choices([Experiment.SAMPLESIMP])
e.analise()

e = Experiment.loadExperiment("net2-Nins-50-cov-2019-6-15-20-25.pickle")
e.run_choices([Experiment.SAMPLESIMP])
e.analise()

e = Experiment.loadExperiment("net2-Nins-100-cov-2019-6-15-21-54.pickle")
e.run_choices([Experiment.SAMPLESIMP])
e.analise()


# e.setup("net2-Nins-10-cov", fname, 10, 100, [10])
# e.run(oneevidence=True)
# e.analise()
#
# print("")
#
# e.setup("net2-Nins-50-cov", fname, 10, 100, [50])
# e.run(oneevidence=True)
# e.analise()
#
# print("")
#
# e.setup("net2-Nins-100-cov", fname, 10, 100, [100])
# e.run(oneevidence=True)
# # e = Experiment.loadExperiment("net1-Nins-100-cov-2019-6-7-13-59.pickle")
# e.analise()
