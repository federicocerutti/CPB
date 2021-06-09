import mpmath
from flask import Flask
from flask import jsonify
from flask import request
from flask_cors import cross_origin
import logging
logging.basicConfig(level=logging.DEBUG)
from SLProbLog.SLProbLog import SLProbLog


app = Flask(__name__, static_url_path='/static')


@app.route('/run', methods=['POST'])
@cross_origin(origin='localhost', headers=['Content- Type', 'Authorization'])
def run_sl_problog():
    # if request.headers['Content-Type'] == 'text/plain':
    #     return "Text Message: " + request.data
    app.logger.debug("data: " + request.data.decode("utf-8"))
    res = SLProbLog(request.data.decode("utf-8"), True).run_KL()

    jsonres = []
    for k,v in res.items():
        if isinstance(v, list):
            newel = {}
            newel['query'] = k
            newel['belief'] = mpmath.nstr(v[0], mpmath.mp.dps)
            newel['disbelief'] = mpmath.nstr(v[1], mpmath.mp.dps)
            newel['uncertainty'] = mpmath.nstr(v[2], mpmath.mp.dps)
            newel['base'] = mpmath.nstr(v[3], mpmath.mp.dps)
            jsonres.append(newel)

    return jsonify(jsonres)

if __name__ == '__main__':
    app.run(debug=True)
