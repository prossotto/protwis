// * CONSTANTS
const margin = { top: 40, right: 200, bottom: 180, left: 130 };
const w = 1200 - margin.left - margin.right,
  h = 1000 - margin.top - margin.bottom;

// array for data in infobox
let info_data = [];

const signprotmat = {
  // * DATA TRANSFORMATION FUNCTIONS
  data: {
    extractPdbIDs: function(dataset: object) {
      let ret = [];
      Object.keys(dataset).forEach(e => ret.push(e.toUpperCase()));

      return ret;
    },

    objectToArray: function(dataset: object) {
      return Object.keys(dataset).map(key => dataset[key]);
    },

    moveKeyToArray: function(dataset: Array, pdb_ids: Array) {
      for (let i = 0; i < pdb_ids.length; i++) {
        const pdb = pdb_ids[i];
        for (let j = 0; j < dataset[i].length; j++) {
          dataset[i][j].push(pdb);
        }
      }
      return dataset;
    },

    // https://stackoverflow.com/questions/10865025/merge-flatten-an-array-of-arrays-in-javascript/25804569#comment50580016_10865042
    flattenOnce: array => {
      return [].concat(...array);
    },

    labelData: function(data, keys) {
      const data_labeled = data.map(function(e) {
        let obj = {};
        keys.forEach(function(key, i) {
          // comment this out later
          if (key === "sig_gn") {
            obj[key] = Math.floor(e[i] / 10);
            return;
          }
          obj[key] = e[i];
        });
        return obj;
      });

      return data_labeled;
    },

    getInteractionTypes: function(dataset: Array) {
      let int_ty = [];
      for (let i = 0; i < dataset.length; i++) {
        const int_arr = dataset[i].int_ty;
        int_ty.push(int_arr);
      }
      int_ty = signprotmat.data
        .flattenOnce(int_ty)
        .filter((v, i, a) => a.indexOf(v) === i);
      let rm_index = int_ty.indexOf("undefined");
      if (rm_index > -1) {
        int_ty.splice(rm_index, 1);
      }

      return int_ty;
    },

    get_additional_receptors: function(data, xvals, prids) {
      let new_receptor_data = [];

      for (let index = 0; index < data.length; index++) {
        const e1 = data[index]["rec_gn"];
        const e2 = data[index]["rec_id"];
        if (xvals.includes(e1) && prids.includes(e2)) {
          new_receptor_data.push(data[index]);
        }
      }
      return new_receptor_data;
    },

    extractRecSigData: function(data, which_component: string) {
      if (which_component === "rec") {
        return _.uniqBy(data, t => [t.rec_gn, t.pdb_id].join());
      } else if (which_component === "sig") {
        return _.uniqBy(data, t => [t.sig_gn, t.pdb_id].join());
      } else {
        console.log("No component specified...");
      }
    },

    dataTransformationWrapper: function(dataset, keys, pdb_sel) {
      dataset = _.pick(dataset, pdb_sel);
      let pdb_ids = signprotmat.data.extractPdbIDs(dataset);
      let data_t = signprotmat.data.objectToArray(dataset);
      data_t = signprotmat.data.moveKeyToArray(data_t, pdb_ids);
      data_t = signprotmat.data.flattenOnce(data_t);
      data_t = signprotmat.data.labelData(data_t, keys);

      let data_t_rec = signprotmat.data.extractRecSigData(data_t, "rec");
      let data_t_sig = signprotmat.data.extractRecSigData(data_t, "sig");
      let int_ty = signprotmat.data.getInteractionTypes(data_t);

      let return_data = {
        transformed: data_t,
        receptor: data_t_rec,
        signprot: data_t_sig,
        inttypes: int_ty,
        pdbids: pdb_ids
      };

      return return_data;
    }
  },

  // * D3 DRAW FUNCTIONS
  d3: {
    // * SETTING UP SVG FOR OUTPUT
    setup: function() {
      let svg = d3
        .select("body")
        .select("div#content")
        .append("div")
        .classed("svg-container", true) //container class to make it responsive
        .append("svg")
        .attr("preserveAspectRatio", "xMinYMin meet")
        .attr(
          "viewBox",
          "0 0 " +
            (w + margin.left + margin.right) +
            " " +
            (h + 200 + margin.top + margin.bottom)
        )
        .classed("svg-content", true) //class to make it responsive
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      return svg;
    },

    // * SETTING THE X/Y SCALE
    xScale: function(data) {
      let xScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.rec_gn)
            .keys()
            .sort(d3.ascending)
        )
        .range([0, w])
        // .round(true)
        .padding(1);

      return xScale;
    },

    yScale: function(data) {
      let yScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.sig_gn)
            .keys()
            .sort(d3.descending)
        )
        .range([h, 0])
        // .round(true)
        .padding(1);

      return yScale;
    },

    // * SETTING THE PDB/SIG-PROT SCALE
    pdbScale: function(data) {
      let pdbScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.pdb_id)
            .keys()
            .sort(d3.descending)
        )
        .range([300, 0])
        .padding(1);

      return pdbScale;
    },

    sigScale: function(data) {
      let sigScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.pdb_id)
            .keys()
            .sort(d3.descending)
        )
        .range([120, 0])
        .padding(1);

      return sigScale;
    },

    // * SETTING THE COLOR SCALE
    colScale: function(data) {
      let colScale = d3
        .scaleOrdinal()
        .domain(data)
        .range(d3.schemeSet3);

      return colScale;
    },

    // * DEFINING AXIS FOR X/Y AND GRID
    xAxis: function(xScale) {
      let xAxis = d3
        .axisBottom(xScale)
        .tickSize(0)
        .tickPadding(8);

      return xAxis;
    },

    yAxis: function(yScale) {
      let yAxis = d3
        .axisRight(yScale)
        .tickSize(0)
        .tickPadding(8);

      return yAxis;
    },

    xAxisGrid: function(xScale, yScale) {
      let xAxisGrid = d3
        .axisTop(xScale)
        .tickSize(h - yScale.step())
        .tickFormat(d => "");

      return xAxisGrid;
    },

    yAxisGrid: function(xScale, yScale) {
      let yAxisGrid = d3
        .axisRight(yScale)
        .tickSize(w - xScale.step())
        .tickFormat(d => "");

      return yAxisGrid;
    },

    // * ADD TOOLTIP FUNCTIONALITY
    tooltip: function(svg) {
      let tip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function(d) {
          return d.rec_gn + "<br>" + d.sig_gn + "<br>" + d.int_ty;
        });
      svg.call(tip);

      return tip;
    },

    // * RENDER DATA
    renderData: function(
      svg,
      data,
      xScale,
      yScale,
      xAxis,
      yAxis,
      xAxisGrid,
      yAxisGrid,
      colScale,
      pdbScale,
      sigScale,
      tip
    ) {
      let shift_left: number = 7 / 8;
      let shift_top: number = 1 / 8;
      let scale_size: number = shift_left - shift_top;
      let offset: number = 1;
      let each_res;

      svg
        .append("g")
        .attr("id", "interact")
        .selectAll("rects")
        .data(data.transformed)
        .enter()
        .append("rect")
        .attr("x", function(d: any) {
          return xScale(d.rec_gn) - shift_left * xScale.step() + offset;
        })
        .attr("y", function(d: any) {
          return yScale(d.sig_gn) + shift_top * yScale.step() + offset;
        })
        .attr("rx", function() {
          if (data.transformed.length < 15) {
            return 5;
          } else {
            return 3;
          }
        })
        .attr("ry", function() {
          if (data.transformed.length < 15) {
            return 5;
          } else {
            return 3;
          }
        })
        .attr("width", xScale.step() * scale_size)
        .attr("height", yScale.step() * scale_size)
        .attr("fill", function(d: any) {
          if (d.int_ty === undefined) {
            return "none";
          } else {
            // return colScale(d.int_ty[0]);
            return "#0000001A";
          }
        })
        .on("mouseover", function(d) {
          tip.show(d);
        })
        .on("mouseout", function(d) {
          tip.hide();
        })
        .on("click", function(d) {
          let index;
          // let rect_x = d3.event.target.getAttribute('x')
          // let rect_y = d3.event.target.getAttribute('y')
          // console.log(rect_x, rect_y)

          // https://stackoverflow.com/a/20251369/8160230
          // select the rect under cursor
          let curr = d3.select(this);

          // Determine if current rect was clicked before
          let active = d.active ? false : true;

          // Update whether or not the elements are active
          d.active = active;

          // set style in regards to active
          if (d.active) {
            curr.style("stroke", "black").style("stroke-width", 2);
            info_data.push(d);
          } else {
            curr.style("stroke", "none").style("stroke-width", 2);
            index = info_data.indexOf(d);
            info_data.splice(index, 1);
          }
          signprotmat.d3.infoBoxUpdate();
        });

      // * DRAWING AXES
      svg
        .append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
        .call(xAxis)
        .selectAll("text")
        .attr("text-anchor", "end")
        .attr("font-size", "12px")
        .attr("dx", "-5px")
        .attr("dy", "-5px")
        .attr("transform", "rotate(-90)");

      svg
        .append("g")
        .attr("class", "y axis")
        .attr(
          "transform",
          "translate(" + (w - xScale.step()) + "," + yScale.step() / 2 + ")"
        )
        .call(yAxis)
        .selectAll("text")
        .attr("font-size", "12px");

      // * DRAWING GRIDLINES
      svg
        .append("g")
        .attr("class", "x grid")
        .attr("transform", "translate(" + 0 + "," + h + ")")
        .call(xAxisGrid);

      svg
        .append("g")
        .attr("class", "y grid")
        .attr("transform", "translate(" + 0 + "," + yScale.step() + ")")
        .call(yAxisGrid);

      // * ADDITIONAL FIGURE LINES
      // top x line
      svg
        .append("line")
        .style("stroke", "black")
        .attr("x1", 0)
        .attr("y1", yScale.step())
        .attr("x2", w - xScale.step())
        .attr("y2", yScale.step());

      // left y line
      svg
        .append("line")
        .style("stroke", "black")
        .attr("x1", 0)
        .attr("y1", yScale.step())
        .attr("x2", 0)
        .attr("y2", h);

      // * ADD AXIS LABELS
      svg
        .append("text")
        .attr("class", "x axis_label")
        .attr("text-anchor", "end")
        .attr("x", 0)
        .attr("y", h + 15)
        .text("GPCR");

      svg
        .append("text")
        .attr("class", "y axis_label")
        .attr("text-anchor", "begin")
        .attr("x", w - 0.8 * xScale.step())
        .attr("y", 0.8 * yScale.step())
        .text("G-Protein");

      // * ADD INFOBOX ELEMENT
      svg
        .append("g")
        .attr("id", "infobox")
        .attr(
          "transform",
          "translate(-15," + (data.inttypes.length + 2) * 20 + ")"
        );

      // * ADDING COLOR LEGEND
      svg
        .append("g")
        .attr("class", "legendOrdinal")
        .attr("transform", "translate(-30," + yScale.step() + ")");

      let legendOrdinal = d3
        .legendColor()
        .cells(data.inttypes.length)
        .scale(colScale)
        // .cellFilter(function (d) { return d.label !== "undefined" })
        .orient("vertical")
        .labelOffset(-20);

      svg
        .select(".legendOrdinal")
        .call(legendOrdinal)
        .selectAll("rect")
        .attr("rx", 3)
        .attr("ry", 3);
      svg
        .select(".legendOrdinal")
        .selectAll("text")
        .attr("class", "legend")
        .attr("text-anchor", "end");

      // * APPENDING COL TICK ANNOTATION FOR RECEPTOR GNs
      svg
        .append("g")
        .attr("id", "recPDB")
        .attr("transform", "translate(" + 0 + "," + h + ")")
        .selectAll("text")
        .data(data.pdbids)
        .enter()
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10)
        .attr("y", function(d: any) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function(d: any) {
          return d;
        });

      // * APPENDING ROW TICK ANNOTATION FOR SIGPROT GNs
      svg
        .append("g")
        .attr("id", "sigPDB")
        .attr(
          "transform",
          "translate(" + w + "," + yScale.step() + ")rotate(-90)"
        )
        .selectAll("text")
        .data(data.pdbids)
        .enter()
        .append("text")
        .attr("class", "x seq_label")
        .attr("x", function(d: any, i) {
          return 10;
        })
        .attr("y", function(d: any, i) {
          return sigScale.step() * (i + 1);
        })
        .attr("text-anchor", "begin")
        .attr("dy", 65)
        .text(function(d: any) {
          return d;
        });

      // * APPENDING AMINOACID SEQUENCE [RECEPTOR]
      svg
        .append("g")
        .attr("id", "recAA")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
        .append("rect")
        .attr("class", "border-bg")
        .style("fill", "#eaeaea")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());

      each_res = svg
        .select("g#recAA")
        .selectAll("text")
        .data(data.receptor)
        .enter()
        .append("g");

      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", (d: any) => colScale(d.int_ty))
        .attr("x", (d: any) => xScale(d.rec_gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + pdbScale(d.pdb_id) - pdbScale.step())
        .attr("width", xScale.step())
        .attr("height", pdbScale.step());

      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("x", (d: any) => xScale(d.rec_gn))
        .attr("y", (d: any) => pdbScale(d.pdb_id) - pdbScale.step() / 2)
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .text((d: any) => d.rec_aa);

      d3.select("g#recAA")
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());

      // * APPENDING AMINOACID SEQUENCE [SIGPROT]
      svg
        .append("g")
        .attr("id", "sigAA")
        .attr(
          "transform",
          "translate(" +
            (w + (1 / 3) * margin.right) +
            "," +
            yScale.step() / 2 +
            ")"
        )
        .append("rect")
        .style("fill", "#eaeaea")
        .attr("x", 0 + sigScale.step() / 2)
        .attr("y", yScale.step() / 2)
        .attr("width", sigScale.range()[0] - sigScale.step())
        .attr("height", yScale.range()[0] - yScale.step());

      each_res = svg
        .select("g#sigAA")
        .selectAll("text")
        .data(data.signprot)
        .enter()
        .append("g");

      each_res
        .append("rect")
        .style("fill", (d: any) => colScale(d.int_ty))
        .attr("x", (d: any) => sigScale(d.pdb_id) - sigScale.step() / 2)
        .attr("y", (d: any) => yScale(d.sig_gn) - yScale.step() / 2)
        .attr("width", sigScale.step())
        .attr("height", yScale.step());

      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("x", (d: any) => sigScale(d.pdb_id))
        .attr("y", (d: any) => yScale(d.sig_gn))
        .attr("text-anchor", "middle")
        .text((d: any) => d.sig_aa);

      d3.select("g#sigAA")
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", 0 + sigScale.step() / 2)
        .attr("y", yScale.step() / 2)
        .attr("width", sigScale.range()[0] - sigScale.step())
        .attr("height", yScale.range()[0] - yScale.step());

      return svg;
    },

    addReceptor: function(new_data, data, svg) {
      data = _.union(data.transformed, new_data);
      data = signprotmat.data.extractRecSigData(data, "rec");
      let pdb_ids = [];
      _.forEach(_.uniqBy(data, "pdb_id"), value => {
        pdb_ids.push(value["pdb_id"]);
      });

      let pdbScale = signprotmat.d3.pdbScale(data);
      let xScale = signprotmat.d3.xScale(data);

      let selection = svg
        .select("g#recAA")
        .selectAll("text.res_label")
        .data(data);

      let selection_rect = svg
        .select("g#recAA")
        .selectAll("rect.res_rect")
        .data(data);

      let selection_enter = selection.enter().append("g");

      selection_enter
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", "black")
        .attr("x", (d: any) => xScale(d.rec_gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + pdbScale(d.pdb_id) - pdbScale.step())
        .attr("width", xScale.step())
        .attr("height", pdbScale.step())
        .merge(selection_rect)
        .transition()
        .duration(500)
        .attr("x", (d: any) => xScale(d.rec_gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + pdbScale(d.pdb_id) - pdbScale.step())
        .attr("width", xScale.step())
        .attr("height", pdbScale.step());

      selection_enter
        .append("text")
        .attr("class", "res_label")
        .style("fill", "white")
        .attr("x", (d: any) => xScale(d.rec_gn))
        .attr("y", (d: any) => pdbScale(d.pdb_id) - pdbScale.step() / 2)
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .text((d: any) => d.rec_aa)
        .merge(selection)
        .transition()
        .duration(500)
        .attr("x", (d: any) => xScale(d.rec_gn))
        .attr("y", (d: any) => pdbScale(d.pdb_id) - pdbScale.step() / 2);

      selection
        .exit()
        .transition()
        .duration(500)
        .remove();

      selection_rect
        .exit()
        .transition()
        .duration(500)
        .remove();

      selection = svg
        .select("g#recPDB")
        .selectAll("text")
        .data(pdb_ids);

      selection_enter = selection.enter();

      selection_enter
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10)
        .attr("y", function(d: any) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function(d: any) {
          return d;
        })
        .merge(selection)
        .transition()
        .duration(500)
        .attr("x", -10)
        .attr("y", function(d: any) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("dy", 75);

      selection
        .exit()
        .transition()
        .duration(500)
        .remove();

      d3.select("g#recAA")
        .selectAll("rect.border, rect.border-bg")
        .transition()
        .duration(500)
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());
    },

    infoBoxUpdate: function() {
      // create selection and bind data
      let info_box = d3
        .select("g#infobox")
        .selectAll("text")
        .data(info_data);

      // update existing nodes
      info_box
        .attr("y", function(d, i) {
          return i * 15;
        })
        .attr("text-anchor", "end")
        .attr("class", "legend");

      // create nodes for new data
      info_box
        .enter()
        .append("text")
        .attr("y", function(d, i) {
          return i * 15;
        })
        .attr("text-anchor", "end")
        .attr("class", "legend")
        .text(function(d) {
          return d.rec_gn + " : " + d.sig_gn;
        });

      // discard removed nodes
      info_box.exit().remove();

      // print the data again in case it changed
      info_box.text(function(d) {
        return d.rec_gn + " : " + d.sig_gn;
      });
    }
  }
};
