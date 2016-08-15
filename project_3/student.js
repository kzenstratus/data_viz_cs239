/* Authors: Kevin Zen, and Adam Freymiller
CNET: kzen, afreymiller

Accredidation: 
Kindlmann Algorithm Viz Paper
Kindlmann Lectures for Data Viz Spr 16
Kindlmann Template Code

Collaborated: Amy Sitwala, Daniel Ni, Elliot Levy

// How to create our univariate color gradient
http://stackoverflow.com/questions/37649016/javascript-color-gradient-even-split-blue-to-grey-to-red
--------------------------------------------------------------------------------------
*/

/*
  The basic rules for what needs to be availble from this student.js are:

  dataFinish: will be called once at the end of d3.csv()
  choiceSet: will be called with radioButton changes
  toggleState: will be called with clicks on states or their marks

  Beyond that, you can add and re-structure things as you see fit.
  Most of the code below is based on project 2. Places where its
  especially important to add code are marked with "(your code here)"
*/

// trying to use https://toddmotto.com/mastering-the-module-pattern/
var P3=(function () {

/* variable controlling map geometry; you can reduce this if you think
   it will help your depiction of which states are selected, while not
   creating too distracting a boundary between all the states */
var HexScaling = 1.0; // hexagon scaling (1 == touching)

/* radius of circular marks in bivariate case; change this if you
   think it will make things clearer */
var MarkRadius = 5.0;
/* CmapLegSize and PCALegSize are set in index.html since they
   shouldn't be changed */

/* duration, in milliseconds, of transitions between visualizations */
var TransitionDuration = 500;

/* other variables to track current state of visualization */
var CmapUnivariate = false; // is current colormap univariate?
/* you can add variables more here.  For example, how will you keep
   track of whether a state has been selected in the visualization?
   (your code here) */

/* utility functions that should not need changing */
var lerp = function (w,[a,b]) { return (1.0-w)*a + w*b; }
var unlerp = function (x,[x0,x1]) { return (x-x0)/(x1-x0); }
var minmax = function (arr) {
    var minval=arr[0], maxval=minval;
    arr.map(function (x) {
            minval = Math.min(x, minval);
            maxval = Math.max(x, maxval);
        });
    return [minval, maxval];
}

/* toggleState is called when you click on either a state in the map,
   or its indication in the colormap legend; the passed "state" is the
   two letter state abbreviation.  That means you can select the hexagon
   for the state with d3.select("#" + state + "hex"), and the tickmark
   for the state with d3.select("#" + state + "mark"). How you modify
   the tickmark for the state will probably depend on whether a univariate
   or a bivariate colormap is being used (CmapUnivariate) */
var toggleState = function(state) {
     // feel free to remove this next line (for debugging)
    var stateOutline = d3.select("#mapUS").select("#" + state);
    var mapOutline = d3.select("#cmlMarks").select("#" + state + "mark");
    var currMap = mapOutline.style("stroke");
    var currColor = stateOutline.select("path").style("stroke");

    if(currMap == "rgb(0, 0, 0)"){
        mapOutline.style('stroke', 'white');
        mapOutline.style("stroke-dasharray", ("3", "3"));
        //mapOutline.style("stroke-dasharray", ("3", "3"));
    }
    else{
        mapOutline.style('stroke', 'black');
        mapOutline.style('stroke-dasharray', null);
    }

    if(currColor == "none"){
        stateOutline.select("path").style('stroke', 'white');
        stateOutline.select("path").style("stroke-dasharray", ("3", "3"));
    }
    else{
        stateOutline.select("path").style('stroke', null);
        stateOutline.select("stroke-dasharray", null);
    }

}

/* PCA: computes PCA of given array of arrays.
   uses http://www.numericjs.com for linear algebra */
var PCA = function (dcols) {
    
    if (dcols.length < 3) {
        d3.select("#pcaWarning").html("PCA() needs at least 3 variables (got " + dcols.length+ ")");
        return null;
    }
    /* else got enough variables */
    d3.select("#pcaWarning").html("");
    // dcols: (short) array of (long) data arrays (each element ~ a csv column)
    // drows: (long) array of data vectors (each element ~ a csv row)
    var drows = numeric.transpose(dcols);
    // covar: covariance matrix
    var covar = numeric.dot(dcols,drows);
    /* NOTE: numeric.dot is for matrix multiplication in general,
       which includes matrix-matrix multiply (as above), and
       matrix-vector multiply, as well as
       vector-vector (inner) product, which you might want to use for
       compute coordinates in the basis of PCA eigenvectors */
    // nmeig: numeric.js's eigensystem representation of covar
    var nmeig = numeric.eig(covar);
    /* NOTE: If you see in the javascript console:
       "Uncaught Error: numeric: eigenvalue iteration does not converge -- increase maxiter?"
       then it is likely that one or more values being passed to
       numeric.eig(covar) are not numeric (e.g. "NaN"), which can happen if
       one or more values in dcols are not numeric */
    // evec: array of covariance matrix eigenvectors (unit-length)
    var evec = numeric.transpose(nmeig.E.x);
    // evec: array of corresponding eigenvalues
    var eval = nmeig.lambda.x;
    // esys: zipping up each component of eigensysem into a little object:
    // "l" for eigenvalue, "v" eigenvector, and "mm" for zero-centered range
    // of projections of data into that eigenvector
    var esys = eval.map(function (_,i) {
            var mindot = 0, maxdot = 0;
            drows.map(function (_,j) { // learn range of projections
                    var x = numeric.dot(drows[j],evec[i]);
                    mindot = Math.min(mindot, x);
                    maxdot = Math.max(maxdot, x);
                });
            // center range around zero
            var mmin = Math.min(mindot, -maxdot);
            var mmax = Math.max(-mindot, maxdot);
            // make sure the range itself is non-zero
            if (mmin == mmax) {
                mmin = -1;
                mmax = 1;
            }
            return {"l": eval[i],
                    "v": evec[i],
                    // simplify needlessly precise representation of range
                    "mm": [d3.format(".3f")(mmin), d3.format(".3f")(mmax)]};
        });
    // sort eigensystem in descending eigenvalue order
    esys.sort(function (a,b) {
            var x = a.l; var y = b.l;
            return ((x < y) ? 1 : ((x > y) ? -1 : 0));
        });
    return esys;
}

function getMaxMin(data) {

    // Get max and min from above array
    var min = Math.min.apply(null, data),
        max = Math.max.apply(null, data);

    // Return as a bundle
    return [max, min];
};


/* dataNorm should take an array of scalar data values and return an
   array resulting from two transformations:
   1) subtract out the mean
   2) make the variance 1
   Making the variance 1 means that no data variable will out an outsized
   influence on the PCA just because of a choice of units: multiplying a
   variable by 10 won't change its information content, but it would
   increase that variable's role in a PCA. */

//compute data standard deviation, divide all values by it

var dataNorm = function (arr) { 

    // Manually normalize w.o using packages
    var sum = 0;
    var variance = 0;
    var norm = [];
    // Calculate the average and variance
    var avg = d3.mean(arr);
    var variance = d3.variance(arr);

    var sd = Math.sqrt(variance);

    // normalize all datapoints in array
    for(var k = 0; k < arr.length; k++){
        norm.push((arr[k] - avg)/sd);
    }

    return norm;

}

/* (from Project2 solution) some stuff we can use for each
 * univariate map.  Feel free to ignore/delete this function
 * if you want to structure things differently */
var stuff = function (what, mmGiven) {
    var sel = function(d) {return +d[what]}
    var slc = P3.data.map(sel);
    if (what == "PC0" || what == "PC1" || what == "PC2"){
      var mm = P3.eigX
    } else {
    
        var mm = ((typeof mmGiven === 'undefined')
              ? minmax(slc) // mmGiven not passed, find min,max
              : mmGiven);   // use given mmGiven
    }
    return {"select" : sel,
            "minmax" : mm,
            "cmlscl" : d3.scale.linear().domain(mm).range([0,P3.CmapLegSize-1]),
            };
}

var dataFinish = function (data) {
    /* save data for future reference (for simplicity, from here on
       out P3.data is the only way we'll refer to the data) */
    P3.data = data;

    /* much of the code here is from Project2 reference solution
       http://people.cs.uchicago.edu/~glk/ class/DataVis/p2.js
       but you should feel free to modify/augment/edit it as you
       see fit for your work (your code here) */
    var voteTotMax = 0;
    P3.data.map(function(d) {
            var VT = +d["ObamaVotes"] + +d["RomneyVotes"];
            d["VT"] = VT;
            d["PL"] = +d["ObamaVotes"]/(1.0 + VT);
            voteTotMax = Math.max(voteTotMax, VT);
        });
    P3.data.map(function(d) {
            d["VA"] = 1 - Math.pow(1- d["VT"]/voteTotMax, 3);
        });

    /* learn earnings ranges */
    P3.earnWMinMax = minmax(P3.data.map(function(d) {return +d["WE"]}));
    P3.earnMMinMax = minmax(P3.data.map(function(d) {return +d["ME"]}));

    /* obesity-related things */
    P3.obeseStuff = stuff("OB");
    var _obeseCmap = d3.scale.linear() /* colormap prior to quantization */
        .domain([0,0.4,1])
        .range([d3.rgb(100,200,100), d3.rgb(220,220,210), d3.rgb(130,0,0)]);
    P3.obeseCmap = function(r) {
        var w0 = Math.round(lerp(unlerp(r,P3.obeseStuff["minmax"]), [-0.5, 6.5]));
        return _obeseCmap(unlerp(Math.min(6, w0),[-0.5, 6.5]));
    }

    /* create unemployment colormap */
    P3.unempStuff = stuff("UN");
    P3.unempCmap = d3.scale.linear()
        .domain([0,1/3,2/3,1].map(function(w) {return lerp(w,P3.unempStuff["minmax"]);}))
        .range([d3.rgb(0,0,0), d3.rgb(210,0,0), d3.rgb(255,210,0), d3.rgb(255,255,255)]);

    /* create infant mortality map */
    P3.imortStuff = stuff("IM");
    P3.imortCmap = function(d) {
        var scl = d3.scale.linear().domain(P3.imortStuff["minmax"]);
        return d3.hcl(scl.range([330,-15])(d),
                      25*Math.pow(Math.sin(scl.range([0,3.14159])(d)),2),
                      scl.range([0,100])(d));
    }

    /* create univariate voter maps */
    P3.pleanStuff = stuff("PL", [0,1]);
    var Dhcl = d3.hcl(d3.rgb(0,0,210));
    var Rhcl = d3.hcl(d3.rgb(210,0,0));
    P3.pleanCmap = function(x) {
        return d3.hcl(x < 0.5 ? Rhcl.h : Dhcl.h,
                      (x < 0.5 ? Rhcl.c : Dhcl.c)*
                      (1 - Math.pow(1 - (Math.abs(x-0.5)/0.5),4)),
                      lerp(x,[Rhcl.l,Dhcl.l]));
    }

    /* create bivariate voter map */
    P3.plean2Cmap = function([pl,va]) {
        var col = P3.pleanCmap(pl);
        return d3.hcl(col.h,  lerp(va,[0,col.c]),  lerp(va,[100,col.l]));
    }

    /* create bivariate earnings maps */
    P3.ERcmap = function([mm,ww]) {
        var erw = unlerp(ww,P3.earnWMinMax);
        var erm = unlerp(mm,P3.earnMMinMax);
        return d3.lab(25+40*(erw + erm), 0, 170*(erm - erw));
    }

}

var choiceSet = function (wat,pvars) {
    var max, min;
    // Used to determine non selected feats
    var all_feats = ["OB", "UN", "IM", 'PL', 'VA', 'ME', 'WE', 'GS', 'FB']

    // Filter out non_selected features
    other_feats = all_feats.filter(function(x) {
              return pvars.indexOf(x) == -1;
              });

    if (wat.startsWith("PC")) {

        if (pvars.length < 1) {
            d3.select("#pcaWarning").html("Select at least one variable below for PCA");
            return;
        }

        d3.select("#pcaWarning").html("");
        /* Else we have at least one variable for PCA; so we do that here,
           in the following steps:
           1) make an array (suppose its called "dcols") of the result
           of calling dataNorm() on each data variable (your code here)
           (can be as little as 3 lines) */

          // Normalize all selected variables in all states
          var dcols = pvars.map(function(vars) {
                      var tmp = P3.data.map(function(x){
                          return parseFloat(x[vars]) });
                      return dataNorm(tmp);
                      });

        /* 2) If less than 3 variables were selected for PCA, add to "dcols"
           one or two arrays of zeros, so that PCA() has at least three
           data variables to work on (your code here) (a few lines) */
        if(pvars.length < 3){
            for(var k = pvars.length; k < 3; k++){

                // Create array of zeros to add to dcols if 
                // not enough selected feats
                var zero = Array.apply(null, Array(P3.data.length)).map(Number.prototype.valueOf,0);
                dcols.push(zero);
            }
        }

           /* 3) call PCA(dcols), and add to P3.data the coordinates of each
           datum in the basis of the first three principle components.  Note
           that "var drows = numeric.transpose(dcols)" will get you an array
           of per-state data (row) vectors, and then with
           "P3.data.map(function(d,ii) { })" you can set PCA coordinate
           fields in per-state datum "d" from the dot product between
           drows[ii] and the PCA eigenvectors. Visualizing the PCA
           results should use these PCA coordinates in the same way that
           in the previous project you used the original data variables.
           (your code here) (roughly ~20 lines of code) */

        // Calculate PCA, store the output into P3.eigen
        var output = PCA(dcols);
        // Store this output for future uses
        P3.eigen = output;


        var drows = numeric.transpose(dcols);

          // Derive the eigenvectors we want to use
        var basis = []
        for (var i = 2; i < wat.length; i++){
          basis.push( output[parseInt(wat[i])]['v'] )

          // Store the min and max of the mapped values
          if (i == 2){
              P3.eigX = output[parseInt(wat[i])]["mm"]
            }
          if (i == 3){
              P3.eigY = output[parseInt(wat[i])]["mm"]
          }
        }

        // Set the max and min for the legend
        var temp0 = Math.abs(P3.eigen[0]["mm"][1])
        var temp1 = Math.abs(P3.eigen[1]["mm"][1])
        var temp2 = Math.abs(P3.eigen[2]["mm"][1])
        var tempMax = Math.max(temp0,temp1,temp2)
        P3.eigMM = [-tempMax, tempMax]

        P3.data.map(function(d,ii) {
                  d["xy"]  =  numeric.dot( basis ,drows[ii])
                  // I can't believe I need to manually input
                  // the data as floats because the provided code structure
                  // is so inflexible.
                  d["PC0"] = 0.0
                  d["PC1"] = 0.0
                  d["PC2"] = 0.0
                  switch (wat){
                    // Extract data data from the resulting PCA matrix
                    // Note we re-used "xy" parameters to store this matrix
                    // This shouldn't be an issue since we only use xy to form
                    // the initial hexagons

                    // Case for bivariate
                    case "PC01":
                    case "PC02":
                    case "PC12":
                      d[wat] = d["xy"];
                      break;
                    // Case for univariate
                    case "PC0" :
                    case "PC1" :
                    case "PC2" :
                      d[wat] = d["xy"][0]
                      break;
                  }
        })

        // Store for univariate packaged min max, etc.  
        P3.pc0Stuff = stuff("PC0")
        P3.pc1Stuff = stuff("PC1")
        P3.pc2Stuff = stuff("PC2")


        /* 4) Visualize what the PCA did with the given data variables inside
           the #pcaMarks svg by changing the text element #pcaXX for
           all variables XX (selected via d3.select("#pca" + XX)):
           a) Make the text opaque for the variables actually included in
           the PCA, and transparent for the rest.
           b) For the variables in PCA, move the text to a position that
           indicates how that variable is aligned with the principle
           component(s) shown (one component for PC0, PC1, PC2, and
           two components for PC01, PC02, PC12). Compute this by forming
           a vector of length pvars.length which is all 0s except for 1 at
           the index of XX in pvars, and then using numeric.dot() to get
           the dot product with a principle component eigenvector. Since
           this is the dot product of two unit-length vectors, the result
           should be in [-1,1], which you should map to coordinates
           [30,P3.PCALegSize-30]) in X or [P3.PCALegSize-30,30]) in Y.
           Text is moved by modifying the "transform" attribute to
           "translate(cx,cy)" for position (cx,cy). For variables not
           in the PCA, the text should be moved back to the center at
           (P3.PCALegSize/2,P3.PCALegSize/2).  You can iterate over the
           #pcaXX with "P3.PCAVars.map(function(XX) { })".
           Changes to both opacity and position should also be made via a
           transition of duration TransitionDuration.  (your code here)
           (roughly ~30 lines of code) */



            // Reset all non-used features
            other_feats.map(function(vars){
                    d3.select("#pca" + vars)
                    .transition().duration(TransitionDuration)
                    .attr("opacity", 0)
                    .attr("transform",
                          "translate("+ parseFloat(P3.PCALegSize/2)+
                          ","+ parseFloat(P3.PCALegSize/2)+")")
            })


            pvars.map(function(vars){
                  // Set opacity if in pvars, else set to 0
                  // Default loc params
                  var tmp = Array.apply(null, new Array(dcols.length)).map(Number.prototype.valueOf,0);
                  tmp[pvars.indexOf(vars)] = 1.0

                  // Center everything around 150, 150
                  var loc     = numeric.dot(basis,tmp)
                  var x       = loc[0]
                  var y       = 0.0
                  // Bivariate case
                  if (basis.length > 1){ y = loc[1]}
                  // map the projected data into legend coordinates
                  x = lerp(unlerp(x,[-1,1]),[30, P3.PCALegSize-30])
                  y = lerp(unlerp(y,[-1,1]),[P3.PCALegSize-30, 30])

                  // Change the SVGs
                  d3.select("#pca" + vars)
                    .transition().duration(TransitionDuration)
                    .attr("opacity", 1)
                    .attr("transform", "translate("+ parseFloat(x) + "," + parseFloat(y)+")")
            });

            //propagate univariate PCA map output
            if(wat == "PC0" || wat == "PC1" || wat == "PC2"){
                d3.select("#cmlMarks").selectAll("ellipse")
                .data(P3.data)
                .transition().duration(TransitionDuration)
                .attr("rx", 0.05) // if zero, outline may disappear
                .attr("ry", P3.CmapLegSize/4)
                .attr("cx", function(d){
                    if(wat == "PC0"){
                        max = output[0]['mm'][1];
                        min = output[0]['mm'][0];
                    }
                    else if(wat == "PC1"){
                        max = output[1]['mm'][1];
                        min = output[1]['mm'][0];
                    }
                    else{
                        max = output[2]['mm'][1];
                        min = output[2]['mm'][0];
                    }
                    return lerp(unlerp(d.xy[0], [min, max]), [0, P3.CmapLegSize])})
                .attr("cy", P3.CmapLegSize/2);

                d3.select("#xminlabel").html("<text>" + min + "</text>");
                d3.select("#xmaxlabel").html("<text>" + max + "</text>");
                d3.select("#yminlabel").html("<text></text>");
                d3.select("#ymaxlabel").html("<text></text>");
            }

            else{
                d3.select("#cmlMarks").selectAll("ellipse")
                    .data(P3.data)
                    .transition().duration(TransitionDuration)
                    .attr("rx", MarkRadius).attr("ry", MarkRadius)
                    .attr("cx", function(d) {
                        if(wat == "PC01"){
                            max = output[0]['mm'][1];
                            min = output[0]['mm'][0];
                        }
                        else if(wat == "PC02"){
                            max = output[1]['mm'][1];
                            min = output[1]['mm'][0];
                        }
                        else{
                            max = output[2]['mm'][1];
                            min = output[2]['mm'][0];
                        }
                        return lerp(unlerp(d.xy[0], [min, max]), [0, P3.CmapLegSize])})
                    .attr("cy", function(d) {
                        return lerp(unlerp(d.xy[1], [min, max]), [P3.CmapLegSize, 0])});

                d3.select("#xminlabel").html("<text>" + min + "</text>");
                d3.select("#xmaxlabel").html("<text>" + max + "</text>");
                d3.select("#yminlabel").html("<text>" + min + "</text>");
                d3.select("#ymaxlabel").html("<text>" + max + "</text>");
            }
    }
    else {
        // same as above, reset all non used features
        d3.select("#pcaWarning").html("");
        all_feats.map(function(vars){
                    d3.select("#pca" + vars)
                    .transition().duration(TransitionDuration)
                    .attr("opacity", 0)
                    .attr("transform",
                      "translate("+ parseFloat(P3.PCALegSize/2)+
                          ","+ parseFloat(P3.PCALegSize/2)+")")
            })
        /* else this isn't a PCA visualization, so none of the
           variables are involved in the PCA, so re-center all the PCA
           marks and make them transparent (your code here) (~10 lines) */
    }


    // Colors used for univariate colormap
    botRgb = d3.rgb(127, 0, 255); // Purple
    topRgb = d3.rgb(255, 255, 0); // Yellow
    midRgb = d3.rgb(128,128,128); // Grey

    // Create the univariate colormap
    P3.uniMap = function(x){
        // Color from purple to grey to yellow (complementary colors)
        var color = d3.scale.linear()
          .domain([P3.eigX[0], 0, P3.eigX[1]])
          .range([botRgb, midRgb, topRgb]);
        return color(x)

      }

    // Create the bivariate colormap 2d interpolation using LAB
    P3.biMap = function([ x, y]){
      var mmx = unlerp(x, P3.eigMM);
      var mmy = unlerp(y, P3.eigMM);
      return d3.lab(25+40*(mmx + mmy), 0, 170*(mmy - mmx));
    }
    /* is this a univariate map? */
    CmapUnivariate = (["OB", "UN", "IM", "VU", "PC0", "PC1", "PC2"].indexOf(wat) >= 0);

    /* set the colormapping function */
    var colormap = {"OB" : P3.obeseCmap,
                    "UN" : P3.unempCmap,
                    "IM" : P3.imortCmap,
                    "VU" : P3.pleanCmap,
                    "VB" : P3.plean2Cmap,
                    "ER" : P3.ERcmap,
                    "PC0": P3.uniMap,
                    "PC1": P3.uniMap,
                    "PC2": P3.uniMap,
                    "PC01": P3.biMap,
                    "PC02": P3.biMap,
                    "PC12": P3.biMap,

                    /* anything else? (your code here) */
    }[wat];
    var cml, cmlx, cmly, sel, mmx, mmy;
    if (CmapUnivariate) {
        var stf = {"OB" : P3.obeseStuff,
                   "UN" : P3.unempStuff,
                   "IM" : P3.imortStuff,
                   "VU" : P3.pleanStuff,
                   "PC0": P3.pc0Stuff,
                   "PC1": P3.pc1Stuff,
                   "PC2": P3.pc2Stuff,
                   /* anything else? (your code here) */
        }[wat];
        [cml,mmx,sel] = [stf["cmlscl"], stf["minmax"], stf["select"]];
        mmy = null;
    } else {
        cml = mmx = mmy = sel = null;
    }
    /* handle the bivariate cases */
    switch (wat) {
    case "VB" :
        cmlx = cmly = d3.scale.linear().domain([0, 1]).range([0,P3.CmapLegSize-1]);
        mmx = mmy = [0,1];
        sel = function(d) {return [+d.PL,+d.VA]};
        break;
    case "ER" :
        cmlx = d3.scale.linear().domain(P3.earnMMinMax).range([0,P3.CmapLegSize-1]);
        cmly = d3.scale.linear().domain(P3.earnWMinMax).range([0,P3.CmapLegSize-1]);
        mmx = P3.earnMMinMax;
        mmy = P3.earnWMinMax;
        sel = function(d) {return [+d.ME,+d.WE]};
        break;
    // Hard code in the different cases for PC01, PC02, PC03,
    // They all have the same logic, just different input values
    // because the current method of storing and implementing 
    //this information is inflexible
    case "PC01":
      cmlx = d3.scale.linear().domain(P3.eigMM).range([0,P3.CmapLegSize-1]);
      cmly = d3.scale.linear().domain(P3.eigMM).range([0,P3.CmapLegSize-1]);
      mmx = P3.eigMM;
      mmy = P3.eigMM;
      sel = function(d) {return [d["PC01"][0], d["PC01"][1]]};
      break;
    case "PC02":
      cmlx = d3.scale.linear().domain(P3.eigMM).range([0,P3.CmapLegSize-1]);
      cmly = d3.scale.linear().domain(P3.eigMM).range([0,P3.CmapLegSize-1]);
      mmx = P3.eigMM;
      mmy = P3.eigMM;
      sel = function(d) {return [d["PC02"][0], d["PC02"][1]]};
      break;
    case "PC12":
      cmlx = d3.scale.linear().domain(P3.eigMM).range([0,P3.CmapLegSize-1]);
      cmly = d3.scale.linear().domain(P3.eigMM).range([0,P3.CmapLegSize-1]);
      mmx = P3.eigMM;
      mmy = P3.eigMM;
      sel = function(d) {return [d["PC12"][0], d["PC12"][1]]};
      break;

    }


    /* 1) reapply colorDatum to the "fill" of the states in #mapUS.
       be sure to add a transition that lasts TransitionDuration */
    d3.select("#mapUS").selectAll("path")
        .data(P3.data) /* (your code here) */
        .transition().duration(TransitionDuration) //add the transition
        .style("fill", function(d){ return colormap(sel(d)); });

    /* 2) reset pixels of cmlImage.data, and redisplay it with
       P3.cmlContext.putImageData(P3.cmlImage, 0, 0); */
    if (CmapUnivariate) {
        for (var j=0, k=0, c; j < P3.CmapLegSize; ++j) {
            for (var i=0; i < P3.CmapLegSize; ++i) {
                if (0 == j) {
                    c = d3.rgb(colormap(cml.invert(i)));
                    P3.cmlImage.data[k++] = c.r;
                    P3.cmlImage.data[k++] = c.g;
                    P3.cmlImage.data[k++] = c.b;
                    P3.cmlImage.data[k++] = 255;
                } else {
                    P3.cmlImage.data[k] = P3.cmlImage.data[(k++)-4*P3.CmapLegSize];
                    P3.cmlImage.data[k] = P3.cmlImage.data[(k++)-4*P3.CmapLegSize];
                    P3.cmlImage.data[k] = P3.cmlImage.data[(k++)-4*P3.CmapLegSize];
                    P3.cmlImage.data[k] = 255; k++;
                }
            }
        }
    } else {
        for (var j=0, k=0, c; j < P3.CmapLegSize; ++j) {
            for (var i=0; i < P3.CmapLegSize; ++i) {
                c = d3.rgb(colormap([cmlx.invert(i),
                                     cmly.invert(P3.CmapLegSize-1-j)]));
                P3.cmlImage.data[k++] = c.r;
                P3.cmlImage.data[k++] = c.g;
                P3.cmlImage.data[k++] = c.b;
                P3.cmlImage.data[k++] = 255;
            }
        }
    }
    P3.cmlContext.putImageData(P3.cmlImage, 0, 0);

    /* 3) set d3.select("#xminlabel").html(), and similarly for the other
       three labels, to reflect the range of values that are
       colormapped when displaying "wat".  For univariate maps,
       set xminlabel and yminlabel to show the range, and set
       yminlabel and ymaxlabel to an empty string.  For bivariate
       maps, set all labels to show the X and Y ranges. */
    d3.select("#xminlabel").html("<text>" + parseFloat(mmx[0]).toFixed(2) + "</text>");
    d3.select("#xmaxlabel").html("<text>" + parseFloat(mmx[1]).toFixed(2) + "</text>");
    if (CmapUnivariate) {
        d3.select("#yminlabel").html("<text></text>");
        d3.select("#ymaxlabel").html("<text></text>");
    } else {
        d3.select("#yminlabel").html("<text>" + parseFloat(mmy[0]).toFixed(2) + "</text>");
        d3.select("#ymaxlabel").html("<text>" + parseFloat(mmy[1]).toFixed(2) + "</text>");
    }

    /* 4) update the geometric attributes (rx, ry, cx, cy) of the #cmlMarks
       to indicate the data variables, and any other attributes you want
       to control according to whether the state is selected. Changes should
       happen with a transition of duration TransitionDuration.
       (your code here) (or interspersed below) */
    if (CmapUnivariate) {
        d3.select("#cmlMarks").selectAll("ellipse")
            .data(P3.data)
            .transition().duration(TransitionDuration)
            .attr("rx", 0.05) // if zero, outline may disappear
            .attr("ry", P3.CmapLegSize/4)
            .attr("cx", function(d) { return 0.5+cml(sel(d)); })
            .attr("cy", P3.CmapLegSize/2);
    } else {
        d3.select("#cmlMarks").selectAll("ellipse")
            .data(P3.data)
            .transition().duration(TransitionDuration) 
            .attr("rx", MarkRadius).attr("ry", MarkRadius)
            .attr("cx", function(d) { return 0.5+cmlx(sel(d)[0]); })
            .attr("cy", function(d) { return P3.CmapLegSize-0.5-cmly(sel(d)[1]); });
    }
}

/* shouldn't have to change anything from here on */
return { // the P3 "API"
    HexScaling: HexScaling,
    choiceSet: choiceSet,
    dataFinish: dataFinish,
    toggleState: toggleState,
};

})();






/* Answer questions here. Each should be no more than ~40 words.

#1) Concisely describe and justify your method of indicating, in the map
and in the colormap, whether a state is selected.

In terms of determining whether a state was selected or not, we identified selected states as either having a stroke style 
on their border or not, and then changed the CSS properties of the state to a dashed white border if it had no stroke style,
and removed all CSS border properties if it previously had them. Hypothetically to
solve the problem of mis-selecting a state when you select 6 states around it, as suggested
by Elliot Levy, we can use inward pointing triangle boarder.

In terms of determining whether a state was selected or not, we identified selected states as either having a stroke style 
on their border or not, and then changed the CSS properties of the state to a dashed white border if it had no stroke style,
and removed all CSS border properties if it previously had them. 


#2) In the terminology of "An Algebraic Process for Visualization
Design" (class May 26), what is one "confuser" for PCA as it it used
in this project (i.e. a particular change in the data that will have
no significant effect on the PCA result)?  (hint: think about what
dataNarm() does, and how the covariance matrix is computed).  There
are at least three possible answers.

If we perform a linear transformation (confuser) on the underlying data, normalizing the data
will center the data and produce the same inputs into PCA. Thus we don't see any visual
difference but the underlying data was changed.


#3) Concisely describe and justify your design of colormaps for the
univariate and bivarite cases of PCA visualization.


Univariate: We chose an even split colormap that ranges from purple to grey to yellow.
Since the range of PCA is normalized around 0, we figured this grey center was appropriate
to denote the intrinsic meaning of 0 produced from normalizing the data. 

We know that pca normalized is a ratio variable as all of original feature values are ratio variables
 and due to the intrinsic meaning of zero.
Note we chose the above colors as they are complementary colors, and in the rgb space
can be used to depict linear scales easily, rather than if we used change in hue (eg
went from say purple to blue)

Bivariate: 
Since we're dealing with ratio variables, the continuous changes in brightness employed 
by our color map accurately reflect corresponding differences in the scalar quantities 
we're plotting. 


#4) Based on exploring the data with PCA, what was a result that you found
from PCA of three or four variables?  Describe your result in terms of
which variables had similar vs complementary contributions to the PCA,
or the geographic trend you saw in the map.

Looking at the PCA = 0, OB,UN, VA, ME; OB opposite to UN of the center, showing that
OB and UN have complementary effects on PC0. Geograhically, using the above conditions,
that the east and midwest had similar results and the west had similar results.
 

In the output from performing PCA on OB, UN, IM, and 
PL, each projection onto a two-dimensional plane showed Unemployment
with the greatest distance from the origin (crosshairs), so it contributes the most 
to the variance found within that particular subset.


(extra credit) #5) How did you maximize visual consistency when switching
between PCA and other radio buttons?



*/
