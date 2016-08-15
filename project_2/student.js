/*
Author: Kevin Zen, Adam Freymiller
CNET: kzen, afreymiller

Accredidation: 
Ideas and code structure from Gordan Kindleman and class lectures
https://github.com/mbostock/d3/wiki/Colors
https://bost.ocks.org/mike/join/
https://github.com/mbostock/d3/wiki/Transitions#d3_interpolate
Data: 
  http://www.cdc.gov/nchs/hus/state.htm
  http://www.bls.gov/web/laus/laumstrk.htm
  http://www.cdc.gov/obesity/data/table-adults.html
  http://blog.apps.npr.org/2015/05/11/hex-tile-maps.html
  
Collaboration: Amy Sitwala
*/

var DATADUMP; // Global to store updated data


// Returns the max and min of a target variable
function getMaxMin(data, key) {
    // Get the all the values of target value as a float
    var keyData = data.map(function(x){return parseFloat(x[key])});

    // Get max and min from above array
    var obeseMin = Math.min.apply(null, keyData),
        obeseMax = Math.max.apply(null, keyData);
  
    // Return as a bundle
    return [obeseMax, obeseMin];
};

//propagate tick marks and circles based on univariate vs. bivariate

function generateLinesCircles(category, dataPoints){
    if(category == "U"){
        d3.select("#cmlTicks").selectAll("line").remove();
        d3.select("#cmlCircs").selectAll("circle").remove();
        d3.select("#cmlTicks").selectAll("line")
            .data(dataPoints)
            .enter().append("line")
                .attr("y1", 0)
                .attr("y2", CmapLegSize/2)
                .attr("x1", function(d) {return d * CmapLegSize;} ) // update
                .attr("x2", function(d) {return d * CmapLegSize;} ) // update
                .attr("stroke", "black")
                .attr("stroke-opacity", 0.6)
                .attr("stroke-width", 1.5);

        d3.select("#cmlTicks").attr("display", null);
    }
    else if(category == "B"){
        d3.select("#cmlTicks").selectAll("line").remove();
        d3.select("#cmlCircs").selectAll("circle").remove();
        d3.select("#cmlCircs").selectAll("circle")
            .data(dataPoints)
            .enter().append("circle")
                .attr("r", 4)
                .attr("cx", function(d) { return CmapLegSize * d[0]; }) // update
                .attr("cy", function(d) { return CmapLegSize - (CmapLegSize * d[1]); }) // update
                .attr("stroke", "black")
                .attr("stroke-opacity", 0.6)
                .attr("stroke-width", 1.5)
                .attr("fill", "none");

        d3.select("#cmlCircs").attr("display", null);
    }
}

/* Our interpretation of lerp function
 Credit:http://stackoverflow.com/questions/3080421/javascript-color-gradient
 Answerd by Desau */
function makeInterval(a, b, ratio) {
    return(a + (b-a)*(ratio));
};
  
//helper function for generating labels to the color map
function addLabel(maxX, minX,maxY,minY ) {
    var xMax = d3.selectAll("#cmlLabels").select("#xmaxlabel").select("text");
    var xMin = d3.selectAll("#cmlLabels").select("#xminlabel").select("text");
    var yMax = d3.selectAll("#cmlLabels").select("#ymaxlabel").select("text");
    var yMin = d3.selectAll("#cmlLabels").select("#yminlabel").select("text");

    xMax.text(function(d) { return maxX });

    xMin.text(function(d) { return minX });

    yMax.text(function(d) { return maxY });

    yMin.text(function(d) { return minY });
}


function obesity(data) {

    // Set our color ranges
  
    var lower = d3.rgb(100,200,100);
    var mid   = d3.rgb(220,220,210);
    var upper = d3.rgb(130,0,0);

    // Get the max and min value for our target variable

    var maxMin = getMaxMin(data,'Obesity');
    var max = maxMin[0],  
        min = maxMin[1];

    // Used to map max and min to 0,1 scale
    var scale = d3.scale.linear().domain([min, max]).range([0, 1]);

    // map 0.4, odd partition into 0,1 scale
    var gradient1Scale = d3.scale.linear().domain([0, 0.4]).range([0, 1]);

    var gradient2Scale = d3.scale.linear().domain([0.4, 1]).range([0, 1]);

    addLabel(max, min, ' ', ' ' );


    // Two color gradients to look up colors from based on obesity values
    var gradient1 = d3.interpolateRgb(lower, mid);
    var gradient2 = d3.interpolateRgb(mid, upper);

    // Need to break up interval into two chunks, because linear scaling not consistent
    var temp1 = [1,3,5].map(function(x){return gradient1(x /(14 * 0.4))});
    var temp2 = [7,9,11,13].map(function(x){return gradient2(gradient2Scale(x/14))});
    var obeseGradient = temp1.concat(temp2);


    // Create array of intervals - 7 evenly spaced bars
    var temp = [0,2,4,6,8,10,12,14].map(function(x) {return x * (1/14)});

    // Create an array of thresholds to map obesity values to color
    intervals = temp.map(function(x){return makeInterval(min,max,x)});
  
    // Used to create store the tick marks
    var dataPoints = [];

    // Go through each path, and plot the colors
    d3.select("#mapUS").selectAll("path")
        .style("fill", function(d){
            for(var i = 1; i < intervals.length; i++) {    
                dataPoints.push(scale(d["Obesity"]));
                var start = intervals[i-1];
                var end = intervals[i];

                if ((d["Obesity"] >= start) && (d["Obesity"] <= end)) {
                  return obeseGradient[i-1];
                }
            }
    });
    // Plots the legend for obesity, 7 chunks
    var c   = document.querySelector("canvas"),
        ctx = c.getContext("2d");

    // Go through each 1/7 chunk and fill with color
    for (var i = 0; i < 7; i++ ) {
        ctx.fillStyle = obeseGradient[i];
        ctx.fillRect(i * CmapLegSize/7, 0, CmapLegSize/7, CmapLegSize);
    }
  
    generateLinesCircles("U", dataPoints);
  
};


function unemployment(data) {
    // Make three color gradients 
    var gradient1 = d3.interpolateRgb("#000000", "#e60000");
    var gradient2 = d3.interpolateRgb("#e60000", "#ffe600");
    var gradient3 = d3.interpolateRgb("#ffe600", "#ffffff");

    var maxMin = getMaxMin(data,'Unemployment');
    var max = maxMin[0],  min = maxMin[1];
    var scale = d3.scale.linear().domain([min, max]).range([0, 1]);

    //create the three separate gradients
    var gradient1Scale = d3.scale.linear().domain([0, 1/3]).range([0, 1]);
    var gradient2Scale = d3.scale.linear().domain([1/3, 2/3]).range([0, 1]);
    var gradient3Scale = d3.scale.linear().domain([2/3, 1]).range([0, 1]);

    addLabel(max,min,' ', ' ');

    var dataPoints = [];

    //based on the three separate partitions of [0, 1], propagate the colors from the gradient

    d3.select("#mapUS").selectAll("path")
        .style("fill", function(d){
            dataPoints.push(scale(d["Unemployment"]));
            if(scale(d["Unemployment"]) < (1/3)){
                return gradient1(gradient1Scale(scale(d["Unemployment"])));
            } 
            else if(scale(d["Unemployment"]) < (2/3)) {
                return gradient2(gradient2Scale(scale(d["Unemployment"])));
            }
            else if(scale(d["Unemployment"]) < 1) {
                return gradient3(gradient3Scale(scale(d["Unemployment"])));
            }
         });

    var c   = document.querySelector("canvas"),
        ctx = c.getContext("2d");

    //divide into three regions for gradients
    for(var t = 0; t < 3; t++){
        for(var i = (Math.floor(CmapLegSize/3) * t); i < ((Math.floor(CmapLegSize/3) * t) + Math.floor(CmapLegSize/3)); i++) {
            var tmp = i/CmapLegSize;

            if(i < CmapLegSize*(1/3)){
                ctx.fillStyle = gradient1(gradient1Scale(tmp));
                ctx.fillRect(i, 0, 1, CmapLegSize);
            }
            else if(i < CmapLegSize*(2/3)){
                ctx.fillStyle = gradient2(gradient2Scale(tmp));
                ctx.fillRect(i, 0, 1, CmapLegSize);
            }
            else if(i < CmapLegSize){
                ctx.fillStyle = gradient3(gradient3Scale(tmp));
                ctx.fillRect(i, 0, 1, CmapLegSize);
            }      
        }
    } 
    
    generateLinesCircles("U", dataPoints);

};

function infantMortality(data) {
    //extract max and min of data
    var maxMin = getMaxMin(data,'InfantMortality');
    var min = maxMin[0], max = maxMin[1];

    //create scales
    var angleScale = d3.scale.linear().domain([max, min]).range([330, -15]);
    var generalScale = d3.scale.linear().domain([max, min]).range([0, 1]);
    var mapGenScale = d3.scale.linear().domain([0, CmapLegSize]).range([0, 1]);
    var mapAngleScale = d3.scale.linear().domain([0, CmapLegSize]).range([330, -15]);

    addLabel(min,max, ' ', ' ');

    var dataPoints = [];

    d3.select("#mapUS").selectAll("path")
        .style("fill", function(d){
            dataPoints.push(generalScale(d["InfantMortality"]));
            var hueAngle = (angleScale(d["InfantMortality"]) + 360) % 360;
            var chroma = 25 * Math.pow(Math.sin((generalScale(d["InfantMortality"]) * Math.PI)), 2);
            var luminance = 100 * generalScale(d["InfantMortality"]);

            return d3.hcl(hueAngle, chroma, luminance);
         });

    var c = document.querySelector("canvas");
    var ctx = c.getContext("2d");

    //generate the vertical strips per how the formula evaluates at the given scaled inputs
    for(var i = 0; i < CmapLegSize; i++){
        var hue = (mapAngleScale(i) + 360) % 360;
        var chromaMap = 25 * Math.pow(Math.sin((mapGenScale(i) * Math.PI)), 2);
        var luminance = 100 * mapGenScale(i);
        var color = d3.hcl(hue, chromaMap, luminance);
        ctx.fillStyle = color;
        ctx.fillRect(i, 0, 1, CmapLegSize);
    }

    generateLinesCircles("U", dataPoints);
}

function earningsSymmetric(data) {
    var maxMinMen = getMaxMin(data,'MenEarning'),
        maxMinWomen = getMaxMin(data, 'WomenEarning'),
        totalMax = Math.max(maxMinWomen[0], maxMinMen[0]);

    var dataPoints = [];

    d3.select("#mapUS").selectAll("path")
        .style("fill", function(d) {
            var esm = d["MenEarning"]/totalMax;
            var esw = d["WomenEarning"]/totalMax;
            dataPoints.push([esm, esw]);
            var L = 25 + (40 * (esw + esm));
            var A = 0;
            var B = 170 * (esm - esw);
            var color = d3.lab(L, A, B);
            return color;
         });

    addLabel(totalMax, 0, totalMax, 0)

    var mapXScale = d3.scale.linear().domain([0, CmapLegSize]).range([0, totalMax]);

    var c = document.querySelector("canvas");
    var ctx = c.getContext("2d");

    //propagate to color map based on evaluation based at that scaled pixel
    for(var i = 0; i < CmapLegSize; i++) {
        for(var j = 0; j < CmapLegSize; j++) {
            var esm = mapXScale(i)/totalMax;
            var esw = mapXScale(j)/totalMax;
            var L = 25 + (40 * (esw + esm));
            var A = 0;
            var B = 170 * (esm - esw);
            var color = d3.lab(L, A, B);
            ctx.fillStyle = color;
            ctx.fillRect(i, CmapLegSize-j, 1, 1);
        }
    }

    generateLinesCircles("B", dataPoints);
}

function earningsRecentered(data) {
    var maxMinMen = getMaxMin(data,'MenEarning'),
        maxMinWomen = getMaxMin(data, 'WomenEarning');

    var mapXScale = d3.scale.linear().domain([0, 200]).range([0, maxMinMen[0]]);
    var mapYScale = d3.scale.linear().domain([0, 200]).range([0, maxMinWomen[0]]);

    //get max and min for men and women earnings
    var EMmax = maxMinMen[0], EMmin = maxMinMen[1],
        EWmax = maxMinWomen[0], EWmin = maxMinWomen[1];
    var dataPoints = [];

    addLabel(EMmax, EMmin, EWmax, EWmin);

    d3.select("#mapUS").selectAll("path")
        .style("fill", function(d) {
            var erm = (d["MenEarning"] - EMmin)/(EMmax - EMmin);
            var erw = (d["WomenEarning"] - EWmin)/(EWmax - EWmin);
            dataPoints.push([erm, erw]);
            var L = 25 + (40 * (erw + erm));
            var A = 0;
            var B = 170 * (erm - erw);
            var color = d3.lab(L, A, B);
            return color;
         });

    var c = document.querySelector("canvas");
    var ctx = c.getContext("2d");

    //generate value of formula returned for pixel input
    for(var i = 0; i < 200; i++) {
        for(var j = 0; j < 200; j++) {
            var esm = mapXScale(i)/95;
            var esw = mapYScale(j)/74;
            var L = 25 + (40 * (esw + esm));
            var A = 0;
            var B = 170 * (esm - esw);
            var color = d3.lab(L, A, B);
            ctx.fillStyle = color;
            ctx.fillRect(i, 200-j, 1, 1);
        }
    }

    generateLinesCircles("B", dataPoints);
};

// Calculates HCL values for the political info (Ob/Rm)
function calcPLHCL(pl){

  // Set the range of colors in HCL Space
  var low   = d3.hcl((d3.rgb(210, 0, 0)).toString()),
      upper = d3.hcl((d3.rgb(0, 0, 210)).toString());

  // Perform calculations as specified in the assignment
  if(pl < 0.5){
      var Cscl   = low.c,
          chroma =  Cscl* (1- Math.pow((1- Math.abs(pl-0.5)/0.5),4));
      return d3.hcl(low.h, chroma, upper.l);
  } 
  else {
      var Cscl   = upper.c,
          chroma = Cscl* (1- Math.pow((1- Math.abs(pl-0.5)/0.5),4));
      return d3.hcl(upper.h,chroma,upper.l)
  }
}
// Function for univariate obama romney
function uniObRm(data){
    
    // Map from 0 - 200 to 0-1
    var scale = d3.scale.linear().domain([0, 200]).range([0, 1]);
    var dataPoints = [];
    
    // Add labels to the legend
    addLabel(1, 0, ' ', ' ');

    // returns the HCL value for each path based on PL value
    d3.select("#mapUS").selectAll("path")
        .style("fill", function(d) {
            dataPoints.push(d["PL"]);
            return calcPLHCL(d["PL"]);       
        }); 

    var c = document.querySelector("canvas");
    var ctx = c.getContext("2d");

    // Fill in legend
    for(var i = 0; i < 200; i++) {
        var color = calcPLHCL(scale(i));
        ctx.fillStyle = color;
        ctx.fillRect(i, 0, 1, 200);
    }
    
    generateLinesCircles("U", dataPoints);

}

// Function for calculating bivariate obama romney
function biObRm(data){

  var low = d3.hcl("#ffffff"); // white

  var plScale = d3.scale.linear().domain([0, 200]).range([0, 1]);
  var dataPoints = [];

  addLabel(1,0,1,0);

  d3.select("#mapUS").selectAll("path")
        .style("fill", function(d){
          dataPoints.push([d["PL"], d["VA"]]);

          // upper color based on univariate rv
          var upper = calcPLHCL(d["PL"])

          // makeInterval is our version of lerp
          var chroma = makeInterval(low.c,upper.c,d["VA"]);
          var lumi = makeInterval(low.l,upper.l,d["VA"]);

          return d3.hcl(upper.h,chroma, lumi);
    });

    var c = document.querySelector("canvas");
    var ctx = c.getContext("2d");
    // Fill in legend
    for(var i = 0; i < 200; i++){
        for(var j = 0; j < 200; j++){
            var upper = calcPLHCL(plScale(i));
            var chroma = makeInterval(low.c,upper.c,plScale(j));
            var lumi = makeInterval(low.l,upper.l,plScale(j)); 
            var color = d3.hcl(upper.h, chroma, lumi);
            ctx.fillStyle = color;
            ctx.fillRect(i, 200-j, 1, 1);
        }
    }
    
    generateLinesCircles("B", dataPoints);

};

// Helper function for calculating VA
function calcVA(state,maxVT){
    vn = Math.min(1,state["VT"]/maxVT);
    va = 1- Math.pow((1-vn),3);
    return va;
}

// calculates PL and VT
function calcPL(state){
    ob = parseFloat(state['ObamaVotes']);
    rm = parseFloat(state['RomneyVotes']);
    return [(ob/(1.0+rm+ob)), (ob+rm)]; 
};


function dataFinish(data) {

    // Calculate VT and PL values, stores them in data
    dataPL_VT = data.map(function(x){( x["PL"] = calcPL(x)[0]); 
                                       x["VT"] = calcPL(x)[1];
                                                    return x });
    // Gets the maxVT for use to calculate VA
    var VT = getMaxMin(data,'VT');
    var maxVT = VT[0];
    
    // Calculate the VA values and store them in data
    dataVA = dataPL_VT.map(function(x){( x["VA"] = calcVA(x,maxVT)); return x });

    DATADUMP = dataVA;
};

function choiceSet(wat) {
    var toggle = {}

    toggle["OB"] = obesity;
    toggle["UN"] = unemployment;
    toggle["IM"] = infantMortality;
    toggle["ES"] = earningsSymmetric;
    toggle["ER"] = earningsRecentered;
    toggle["VB"] = biObRm;
    toggle["VU"] = uniObRm;

    toggle[wat](DATADUMP);
    
};
