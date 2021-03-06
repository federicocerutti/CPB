function drawSLOpinion(trianglesvg, betasvg, trianglewidth, posD, betawidth, betaheight, betainitx, betainity, betamaxheight, divbelief, divdisbelief, divuncertainty, divalpha, divbeta, divqualitative, inbel, inunc) {
    /**
     *  Variables to change the size of the canvas/triangle
     */
    var s = Snap("#"+trianglesvg);
    var width = trianglewidth;
    var D = posD;
    var B = [D[0] + width, D[1]];
    var U = [D[0] + width / 2, D[1] - Math.sqrt(Math.pow(width, 2) - Math.pow(width / 2, 2))];

    /**
     * Drawing the triangle
     */
    var lines = s.group(s.line(B[0], B[1], D[0], D[1]), s.line(B[0], B[1], U[0], U[1]), s.line(D[0], D[1], U[0], U[1]));
    lines.attr({
        stroke: "#000",
        strokeWidth: 2
    });

    s.text(D[0], D[1] + 20, "d");
    s.text(B[0], B[1] + 20, "b");
    s.text(U[0], U[1] - 10, "u");

    /**
     * Qualitative representations
     */
    var totalCertainty = 0.05;
    var highCertainty = 0.16;
    var someCertainty = 0.5;

    var unlikely = 0.15;
    var someUnlikely = 0.3;
    var chancesEven = 0.45;
    var someLikely = 0.55;
    var likely = 0.7;
    var vLikely = 0.85;

    var ytotal = D[1] - (D[1] - U[1]) * totalCertainty;
    var yhigh = D[1] - (D[1] - U[1]) * highCertainty;
    var ysome = D[1] - (D[1] - U[1]) * someCertainty;


    var boundaries = s.group(
        s.line(inverseLeftLine(ytotal), ytotal, inverseRightLine(ytotal), ytotal),
        s.line(inverseLeftLine(yhigh), yhigh, inverseRightLine(yhigh), yhigh),
        s.line(inverseLeftLine(ysome), ysome, inverseRightLine(ysome), ysome),
        s.line(D[0] + unlikely * width, D[1], D[0] + unlikely * width, leftLine(D[0] + unlikely * width)),
        s.line(D[0] + someUnlikely * width, D[1], D[0] + someUnlikely * width, leftLine(D[0] + someUnlikely * width)),
        s.line(D[0] + chancesEven * width, D[1], D[0] + chancesEven * width, leftLine(D[0] + chancesEven * width)),
        s.line(D[0] + someLikely * width, D[1], D[0] + someLikely * width, rightLine(D[0] + someLikely * width)),
        s.line(D[0] + likely * width, D[1], D[0] + likely * width, rightLine(D[0] + likely * width)),
        s.line(D[0] + vLikely * width, D[1], D[0] + vLikely * width, rightLine(D[0] + vLikely * width))
    );
    boundaries.attr({
        stroke: "#AAA",
        strokeWidth: 1,
        strokeDasharray: "5 5"
    });

    /**
     * Draggable opinion
     */
    var opinion = s.circle(D[0] + width / 2, D[1] - 0.5 * Math.sqrt(Math.pow(width, 2) - Math.pow(width / 2, 2)), 5).attr({fill: "red"});


    /**
     *  Beta pdf
     */
    var pdf = Snap("#"+betasvg);

    function drawBeta() {
        pdf.clear();

        var width = betawidth, height = betaheight, initx = betainitx, inity = betainity, maxheight = betamaxheight;

        // Get the parameters from the sliders.
        var alpha = parseFloat($('#'+divalpha).text()); //30; //Number($("#alpha-val").text());
        var beta = parseFloat($('#'+divbeta).text());//10; //Number($("#beta-val" ).text());
        var maxY = 7;//Number($("#y-val"    ).text());


        pdf.line(0, height, initx + width, height).attr({stroke: "#000", strokeWidth: 1});
        pdf.line(initx, maxheight - inity + 1, initx, 0).attr({stroke: "#000", strokeWidth: 1});
        pdf.group(
            pdf.line(initx, height - Math.floor(maxY / 3) * height / maxY, initx + width, height - Math.floor(maxY / 3) * height / maxY),
            pdf.line(initx, height - Math.floor(maxY / 3) * 2 * height / maxY, initx + width, height - Math.floor(maxY / 3) * 2 * height / maxY),
            pdf.line(initx, height - Math.floor(maxY / 3) * 3 * height / maxY, initx + width, height - Math.floor(maxY / 3) * 3 * height / maxY),
            pdf.line(initx + unlikely * width, height, initx + unlikely * width, 0),
            pdf.line(initx + someUnlikely * width, height, initx + someUnlikely * width, 0),
            pdf.line(initx + chancesEven * width, height, initx + chancesEven * width, 0),
            pdf.line(initx + someLikely * width, height, initx + someLikely * width, 0),
            pdf.line(initx + likely * width, height, initx + likely * width, 0),
            pdf.line(initx + vLikely * width, height, initx + vLikely * width, 0),
            pdf.line(initx + 1 * width, height, initx + 1 * width, 0)
        ).attr({
            stroke: "#AAA",
            strokeWidth: 1,
            strokeDasharray: "5 5"
        });

        pdf.group(
            pdf.text(initx - 7, maxheight - 5, "0.00"),
            pdf.text(initx + unlikely * width - 7, maxheight - 5, unlikely.toFixed(2)),
            pdf.text(initx + someUnlikely * width - 7, maxheight - 5, someUnlikely.toFixed(2)),
            pdf.text(initx + chancesEven * width - 7, maxheight - 5, chancesEven.toFixed(2)),
            pdf.text(initx + someLikely * width - 7, maxheight - 5, someLikely.toFixed(2)),
            pdf.text(initx + likely * width - 7, maxheight - 5, likely.toFixed(2)),
            pdf.text(initx + vLikely * width - 7, maxheight - 5, vLikely.toFixed(2)),
            pdf.text(initx + 1 * width - 7, maxheight - 5, "1.00"),
            pdf.text(initx - 8, height - Math.floor(maxY / 3) * height / maxY + 2, Math.floor(maxY / 3)),
            pdf.text(initx - 8, height - Math.floor(maxY / 3) * 2 * height / maxY + 2, Math.floor(maxY / 3) * 2),
            pdf.text(initx - 8, height - Math.floor(maxY / 3) * 3 * height / maxY + 2, Math.floor(maxY / 3) * 3)
        ).attr({'font': '400 10px system-ui'});


        /*
        if (alpha > 120)
            alpha = 120;
        if (beta > 120)
            beta = 120;

        if (!isFinite(alpha))
            alpha = 5000;
        if (!isFinite(beta))
            beta = 5000;

        if (isNaN(alpha))
            alpha = 0.01;
        if (isNaN(beta))
            beta = 0.01;

         */

        var points = [];
        var xpx = initx;
        for (var xi = 0; xi <= 1; xi = xi + 1 / width) {
            yi = betaPDF(xi, alpha, beta)
            if (isFinite(yi)) {
                points.push([xpx, height - yi * height / maxY]);
            }
            xpx += 1;
        }

        if (points.length == 0) {
            alpha = parseFloat($('#'+divalpha).text());
            beta = parseFloat($('#'+divbeta).text());
            var mu = alpha / (alpha + beta);


            pdf.line(initx + mu * width, height, initx + mu * width, 0).attr({stroke: "#000", strokeWidth: 4});
        } else {
            pdf.polyline(points).attr({stroke: "#000", strokeWidth: 1.5, fill: "none"});
        }


    }

    drawBeta();


    /**
     * Getting the value of the draggable opinion
     */
    function opinionToText() {
        var x = opinion.getBBox().cx;
        var y = opinion.getBBox().cy;

        var u = Math.abs((D[1] - y) / Math.sqrt(Math.pow(width, 2) - Math.pow(width / 2, 2)));
        var e = Math.abs((x - D[0]) / width);

        if (u < 0.0001)
            u = 0;

        var b = Math.abs(e - u / 2);
        var d = Math.abs(1 - b - u);

        var alpha = 2 / u * b + 1;
        var beta = 2 / u * d + 1;
        if (!isFinite(alpha) || !isFinite(beta) || isNaN(alpha) || isNaN(beta)) {
            if (e <= 0.01) {
                beta = 100000;
                alpha = 0.1;
            }
            if (e > 0.001 && e <= 0.5) {
                beta = 100000;
                alpha = e / (1 - e) * beta;

            }
            if (e > 0.5 && e < 0.999) {
                alpha = 100000;
                beta = (1 - e) / e * alpha;
            }
            if (e >= 0.999) {
                alpha = 100000;
                beta = 0.1;
            }
        }

        /*
        if (alpha < 0)
            alpha = 10000;
        if (beta < 0)
            beta = 10000;

         */

        $('#'+divalpha).html(Math.abs(alpha).toFixed(4));
        $('#'+divbeta).html(Math.abs(beta).toFixed(4));

        $('#'+divbelief).html(b.toFixed(4));
        $('#'+divdisbelief).html(d.toFixed(4));
        $('#'+divuncertainty).html(u.toFixed(4));

        var qualE = "";
        var qualU = "";

        if (e <= 0.001)
            qualE = "Absolutely not";
        if (e > 0.001 && e < unlikely)
            qualE = "Very unlikely";
        if (e >= unlikely && e < someUnlikely)
            qualE = "Unlikely";
        if (e >= someUnlikely && e < chancesEven)
            qualE = "Somewhat unlikely";
        if (e >= chancesEven && e <= someLikely)
            qualE = "Chances about even";
        if (e > someLikely && e <= likely)
            qualE = "Somewhat likely";
        if (e > likely && e <= vLikely)
            qualE = "Likely";
        if (e > vLikely && e < 0.999)
            qualE = "Very likely";
        if (e >= 0.999)
            qualE = "Absolutely";


        if (u <= totalCertainty)
            qualU = "total confidence";
        if (u > totalCertainty && u <= highCertainty)
            qualU = "high confidence";
        if (u > highCertainty && u <= someCertainty)
            qualU = "some confidence";
        if (u > someCertainty && u <= 0.999)
            qualU = "low confidence";
        if (u >= 0.999)
            qualU = "no confidence";

        $('#'+divqualitative).html(qualE + " with " + qualU);

        drawBeta();

    }

    opinionToText();

    function inverseLeftLine(y) {
        var q = D[1] + Math.sqrt(3) * D[0];
        return (-y + q) / Math.sqrt(3);
    }

    function leftLine(x) {
        var q = D[1] + Math.sqrt(3) * D[0];
        return -Math.sqrt(3) * x + q;
    }

    function inverseRightLine(y) {
        var q2 = -Math.sqrt(3) * U[0] + U[1];
        return (y - q2) / Math.sqrt(3)
    }

    function rightLine(x) {
        var q2 = -Math.sqrt(3) * U[0] + U[1];
        return Math.sqrt(3) * x + q2;
    }

    /**
     * What happens when dragging the opinion
     */
    /*var moveTriangle = function (dx, dy, posx, posy) {

        var tdx, tdy;
        var sInvMatrix = this.transform().globalMatrix.invert();
        sInvMatrix.e = sInvMatrix.f = 0;
        tdx = sInvMatrix.x(dx, dy);
        tdy = sInvMatrix.y(dx, dy);

        this.data('x', +this.data('ox') + tdx);
        this.data('y', +this.data('oy') + tdy);
        if (this.data('y') > D[1]) {
            dy = D[1] - this.data('oy');
            this.data('y', D[1])
        }
        if (this.data('y') < U[1]) {
            dy = U[1] - this.data('oy');
            this.data('y', U[1])
        }

        if (this.data('x') < inverseLeftLine(this.data('y'))) {
            dx = inverseLeftLine(this.data('y')) - this.data('ox');
            this.data('x', inverseLeftLine(this.data('y')));
        }

        if (this.data('x') > inverseRightLine(this.data('y'))) {
            dx = inverseRightLine(this.data('y')) - this.data('ox');
            this.data('x', inverseRightLine(this.data('y')));
        }

        this.attr({
            transform: this.data('ot') + (this.data('ot') ? "T" : "t") + [dx, dy]
        });

        opinionToText();

    };
    var start = function () {
        //this.data('origTransform', this.transform().local);
        this.data('ot', this.transform().local);
        this.data('x', this.getBBox().cx);
        this.data('y', this.getBBox().cy);
        this.data('ot', this.transform().local);
        this.data('ox', this.data('x'));
        this.data('oy', this.data('y'));
    };
    var stop = function () {
        console.log('finished dragging');
    };

    opinion.drag(moveTriangle, start, stop);*/

    function drawOpinion(b, u) {
        var e = b + u / 2;
        opinion.remove();
        opinion = s.circle(D[0] + e * width, D[1] - u * (Math.sqrt(Math.pow(width, 2) - Math.pow(width / 2, 2))), 5).attr({fill: "red"});
        //opinion.drag(moveTriangle, start, stop);
        opinionToText();
    }


    function beliefChanged() {
        var newbelief = parseFloat($('#'+divbelief).text());

        if (newbelief < 0)
            newbelief = 0;
        if (newbelief > 1)
            newbelief = 1;

        var olduncertainty = parseFloat($('#'+divuncertainty).text());
        var newuncertainty = olduncertainty;

        if (newbelief + olduncertainty > 1)
            newuncertainty = 1 - newbelief;

        if (newuncertainty < 0)
            newuncertainty = 0;

        drawOpinion(newbelief, newuncertainty);
    }

    function uncertaintyChanged() {
        var newuncertainty = parseFloat($('#'+divuncertainty).text());

        if (newuncertainty < 0)
            newuncertainty = 0;
        if (newuncertainty > 1)
            newuncertainty = 1;

        var oldbelief = parseFloat($('#'+divbelief).text());
        var newbelief = oldbelief;

        if (newuncertainty + oldbelief > 1)
            newbelief = 1 - newuncertainty;

        if (newbelief < 0)
            newbelief = 0;

        drawOpinion(newbelief, newuncertainty);
    }

    function alphabetachanged() {
        var alpha = parseFloat($('#'+divalpha).text()); //30; //Number($("#alpha-val").text());
        var beta = parseFloat($('#'+divbeta).text());//10; //Number($("#beta-val" ).text());

        if (alpha <= 0)
            alpha = 0.1;
        if (beta <= 0)
            beta = 0.1;

        var sx = alpha + beta;
        var mu = alpha / sx;
        var variance = alpha * beta / (Math.pow(alpha + beta, 2) * (alpha + beta + 1));
        var sz = Math.max(mu * (1 - mu) / variance - 1, 1 / mu, 1 / (1 - mu));
        var newbelief = (alpha - 1) / sz;
        var newuncertainty = 2 / sz;

        if (newbelief < 0)
            newbelief = 0;
        if (newbelief > 1)
            newbelief = 1;

        if (newuncertainty < 0)
            newuncertainty = 0;
        if (newuncertainty > 1)
            newuncertainty = 1;

        drawOpinion(newbelief, newuncertainty)
    }

    $('#'+divbelief).change(beliefChanged);
    $('#'+divuncertainty).change(uncertaintyChanged);
    $('#'+divalpha).change(alphabetachanged);
    $('#'+divbeta).change(alphabetachanged);

    /**
     * Clicking on the triangle
     */
    /*s.click(function (event) {
            // Get coords of click minus the offsets due to CSS positioning of SVG such as margins
            var x = event.clientX - $("#svg").offset().left;
            var y = event.clientY - $("#svg").offset().top;
            if (y < D[1] && y > U[1] && x > inverseLeftLine(y) && x < inverseRightLine(y)) {
                opinion.remove();
                opinion = s.circle(x, y, 5).attr({fill: "red"});
                opinion.drag(moveTriangle, start, stop);
                opinionToText();


            }
        }
    );*/

    drawOpinion(inbel, inunc);
}