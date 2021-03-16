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
sys.path.append('../')

from experiment.experimental_setting import Experiment

# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.1524 & 0.1524 & 0.2111 & 0.1524 & 0.1524 & 0.1540 & 0.1671 \\
# & & P & 0.1459 & 0.1619 & 0.1536 & 0.1448 & 0.1460 & 0.0822 & 0.1969 \\
# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.0744 & 0.0744 & 0.1166 & 0.0749 & 0.0744 & 0.0773 & 0.0752 \\
# & & P & 0.0762 & 0.0870 & 0.0761 & 0.0751 & 0.0757 & 0.0354 & 0.0989 \\
# Computing Simple Sampling: 99%& &   &  BetaProbLog & BetaProbLog AAAI19 & SL Operators & Monte Carlo & SBN & GBT & Credal\\
# & & A & 0.0563 & 0.0563 & 0.0890 & 0.0570 & 0.0563 & 0.0598 & 0.0582 \\
# & & P & 0.0571 & 0.0664 & 0.0573 & 0.0564 & 0.0567 & 0.0244 & 0.0739 \\

fname="networks/net3.file"

e = Experiment()

# e = Experiment.loadExperiment("net3-Nins-10-cov-2019-6-15-19-33.pickle")
# e.run_choices([Experiment.SAMPLESIMP])
# e.analise()
#
# e = Experiment.loadExperiment("net3-Nins-50-cov-2019-6-15-21-37.pickle")
# e.run_choices([Experiment.SAMPLESIMP])
# e.analise()

e = Experiment.loadExperiment("net3-Nins-100-cov-2019-6-15-23-47.pickle")
e.run_choices([Experiment.SAMPLESIMP])
e.analise()


# e.setup("net3-Nins-10-cov", fname, 10, 100, [10])
# e.run(oneevidence=True)
# e.analise()
#
# print("")
#
# e.setup("net3-Nins-50-cov", fname, 10, 100, [50])
# e.run(oneevidence=True)
# e.analise()
#
# print("")
#
# e.setup("net3-Nins-100-cov", fname, 10, 100, [100])
# e.run(oneevidence=True)
# # e = Experiment.loadExperiment("net1-Nins-100-cov-2019-6-7-13-59.pickle")
# e.analise()
