
/**
    * Parses the given PDB string into json
    * @param {String} pdb
    * @returns {Object}
    */
const parsePdb = function parsePdb(pdb) {
    ///// https://github.com/justinmc/parse-pdb ////
    const ATOM_NAME = 'ATOM  ';
    const RESIDUE_NAME = 'SEQRES';
    const pdbLines = pdb.split('\n');
    const atoms = [];
    const seqRes = []; // raw SEQRES entry data
    let residues = []; // individual residue data parsed from SEQRES
    const chains = new Map(); // individual rchaindata parsed from SEQRES

    // Iterate each line looking for atoms
    pdbLines.forEach((pdbLine) => {
        if (pdbLine.substr(0, 6) === ATOM_NAME) {
            // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            atoms.push({
                serial: parseInt(pdbLine.substring(6, 11)),
                name: pdbLine.substring(12, 16).trim(),
                altLoc: pdbLine.substring(16, 17).trim(),
                resName: pdbLine.substring(17, 20).trim(),
                chainID: pdbLine.substring(21, 22).trim(),
                resSeq: parseInt(pdbLine.substring(22, 26)),
                iCode: pdbLine.substring(26, 27).trim(),
                x: parseFloat(pdbLine.substring(30, 38)),
                y: parseFloat(pdbLine.substring(38, 46)),
                z: parseFloat(pdbLine.substring(46, 54)),
                occupancy: parseFloat(pdbLine.substring(54, 60)),
                tempFactor: parseFloat(pdbLine.substring(60, 66)),
                element: pdbLine.substring(76, 78).trim(),
                charge: pdbLine.substring(78, 80).trim(),
            });
        } else if (pdbLine.substr(0, 6) === RESIDUE_NAME) {
            // http://www.wwpdb.org/documentation/file-format-content/format33/sect3.html#SEQRES
            const seqResEntry = {
                serNum: parseInt(pdbLine.substring(7, 10)),
                chainID: pdbLine.substring(11, 12).trim(),
                numRes: parseInt(pdbLine.substring(13, 17)),
                resNames: pdbLine.substring(19, 70).trim().split(' '),
            };
            seqRes.push(seqResEntry);

            residues = residues.concat(seqResEntry.resNames.map(resName => ({
                id: residues.length,
                serNum: seqResEntry.serNum,
                chainID: seqResEntry.chainID,
                resName,
            })));

            if (!chains.get(seqResEntry.chainID)) {
                chains.set(seqResEntry.chainID, {
                    id: chains.size,
                    chainID: seqResEntry.chainID,
                    // No need to save numRes, can just do chain.residues.length
                });
            }
        }
    });

    // Add residues to chains
    chains.forEach((chain) => {
        chain.residues = residues.filter((residue) =>
            residue.chainID === chain.chainID,
        );
    });

    // Add atoms to residues
    residues.forEach((residue) => {
        residue.atoms = atoms.filter((atom) =>
            atom.chainID === residue.chainID && atom.resSeq === residue.serNum,
        );
    });

    return {
        // Raw data from pdb
        atoms,
        seqRes,
        // Derived
        residues, // Array of residue objects
        chains, // Map of chain objects keyed on chainID
    };
}

async function prepare(pdbID) {
    d3.select("#pdbtitle").html(pdbID)

    let url = "	https://files.rcsb.org/download/" + pdbID + ".pdb"

    if (pdbID === '6W9C') {
        url = "https://coronavirus3d.org/structures/6W9C.pdb";
    }

    const response = await fetch(url);

    if (response.ok) { // if HTTP-status is 200-299


        // .then(response => response.text())
        // // .then(data => processStructure(data));
        // .then(data => {
        const pdbString = await response.text();
        const pdbJson = parsePdb(pdbString);

        const chains = pdbJson.chains;

        const chainNames = Array.from(chains.keys());


        // // fill selects with chain names
        const select1 = d3.select("#chain1");
        select1.selectAll('option').remove();

        const select2 = d3.select("#chain2");
        select2.selectAll('option').remove();

        const options1 = select1
            .selectAll('option')
            .data(chainNames).enter()
            .append('option')
            .text(function (d) { return d; });
        const options2 = select2
            .selectAll('option')
            .data(chainNames).enter()
            .append('option')
            .text(function (d) { return d; });

        return chains
    } else {
        alert("HTTP-Error: " + response.status);
        return null;
    }

}

const displayStructureContacts = function (pdbID, chains, divID, chainIDs, maxDistance) {
    pdbID = pdbID || '6W9C';

    const width = 900;
    const height = 900;
    const margin = 100;

    // The radius of the pieplot is half the width or half the height (smallest one). I subtract a bit of margin.
    const radius = Math.min(width, height) / 2 - margin

    const innerRaduis = radius - 40;


    // --------------HELPER FUNCTIONS---------------

    const calcExternalDistances = function (chain1, chain2, chainIDs) {
        let distances = [];

        const l1 = chain1.residues.length;
        const l2 = chain2.residues.length;


        let min = 0, max = 0;

        for (let i = 0; i < l1; i++) {
            let row = new Array(l1);
            row.fill(0);
            distances[i] = row;

            const res1 = chain1.residues[i];

            res1.id = i;
            res1.exlinks = []
            res1.chain = chainIDs[0]

            const caAtoms1 = res1.atoms.filter(a => a.name === CA);
            if (caAtoms1.length < 1) continue;

            for (let j = 0; j < l2; j++) {

                const res2 = chain2.residues[j];

                res2.id = j + l1;
                res2.chain = chainIDs[1]

                res2.exlinks = []

                const caAtoms2 = res2.atoms.filter(a => a.name === CA);
                if (caAtoms2.length < 1) continue;

                const d = getDistanceBetweenAtoms(caAtoms1[0], caAtoms2[0])

                if (i === 0 & j === 0) {
                    min = d; max = d;
                }
                min = Math.min(min, d)
                max = Math.max(max, d)

                row[j] = d;
                res1.exlinks.push({ to: j + l1, distance: d });
            }

            distances[i] = row;

        }
        // console.log(`min ${min}, max:${max}`)


        d3.select("#minmax").html(`Distances: min ${min}, max:${max}, cutoff:${MAX_EX_DISTANCE}`);
    }

    const calcInternalDistances = function (chain) {
        // let distances = [];

        const l = chain.residues.length;

        let min = 100, max = 0;

        for (let i = 0; i < l; i++) {
            let row = new Array(l);
            row.fill(0);
            // distances[i] = row;

            const res1 = chain.residues[i];
            res1.id = i;
            res1.links = []
            const caAtoms1 = res1.atoms.filter(a => a.name === CA);
            if (caAtoms1.length < 1) continue;

            for (let j = i; j < l; j++) {
                const res2 = chain.residues[j];
                const caAtoms2 = res2.atoms.filter(a => a.name === CA);
                if (caAtoms2.length < 1) continue;

                const d = getDistanceBetweenAtoms(caAtoms1[0], caAtoms2[0])

                if (i === 0 & j === 0) {
                    min = d; max = d;
                }
                min = Math.min(min, d)
                max = Math.max(max, d)
                row[j] = d;
                res1.links.push({ to: j, distance: d });
            }

            // distances[i] = row;

        }
        // chain.distances = distances;

        // console.log(`min ${min}, max:${max}`)

        d3.select("#minmax").html(`Distances: min ${min}, max:${max}, cutoff:${MAX_DISTANCE}`);
    }

    const getDistanceBetweenAtoms = function (atom1, atom2) {
        return getDistanceBetweenPoints(atom1.x, atom1.y, atom1.z, atom2.x, atom2.y, atom2.z)
    }

    const getDistanceBetweenPoints = function (x1, y1, z1, x2, y2, z2) {
        return Math.sqrt(
            Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2) + Math.pow((z1 - z2), 2)
        )
    }
    const handleMouseOver = function (d, i) {
        const self = d3.select(this)

        self.classed("highlighted", true)
        // console.log(d)

        if (self.classed("chord")) {

            //find and show ends
            const start = d.source.id;
            const end = d.target.id;
            d3.select("#" + start).classed("highlighted", true);
            d3.select("#" + end).classed("highlighted", true)
        } else {
            const id = d.id;
            if (id.startsWith(RES)) {
                // residue
                // console.log(id)
                let idNum = id.replace(RES, "");

            }
        }
    }
    const handleMouseOut = function (d, i) {
        const self = d3.select(this)

        self.classed("highlighted", false)
        // console.log(d)

        if (self.classed("chord")) {
            //find and show ends
            const start = d.source.id;
            const end = d.target.id;
            d3.select("#" + start).classed("highlighted", false);
            d3.select("#" + end).classed("highlighted", false)
        }
        else {
            //  console.log(d)
        }
    }

    const render = function (chain) {

        const colors = d3.scaleSequential()
            .interpolator(d3.interpolateSpectral)
            .domain([0, chain.residues.length]);


        let residues = chain.residues;

        let data = residues.map((r, i) => {
            return {
                id: RES + r.id,
                len: 1,
                label: r.resName,
                color: colors(i)
            }
        });

        let chords = [];

        residues.forEach(r => {
            r.links.forEach(l => {
                if (l.distance >= MIN_DISTANCE && l.distance <= MAX_DISTANCE) {
                    const from = {
                        id: RES + r.id,
                        start: .5,
                        end: .6
                    };
                    chords.push({
                        value: l.distance,
                        source: from,
                        target: {
                            id: RES + l.to,
                            start: .5,
                            end: .6
                        }
                    })
                }
            });
        })

        // console.log(chords);


        var circos = new Circos({
            container: '#' + divID,
            width: width,
            height: height,
        });
        var configuration = {
            innerRadius: innerRaduis,
            outerRadius: radius,
            class: "residue",
            cornerRadius: 3,
            gap: 0.0001, // in radian
            labels: {
                display: false
            },
            ticks: {
                display: false
            },

        }


        let x = circos
            .layout(data, configuration)
            .text('labels', residues.map(function (d) {
                return {
                    block_id: RES + d.id,
                    position: .5,
                    value: d.resName + ":" + (d.id + 1)
                }
            }), {
                innerRadius: .9,
                outerRadius: 1.3,
                style: {
                    'font-size': 8
                },
                tooltipContent: function (d) {
                    return d.value
                },
                events: {
                    'click.alert': function (datum, index, nodes, event) {
                        window.alert(datum)
                    }
                }
            })
            .chords("l1",
                chords, {
                opacity: 0.4,
                color: 'grey',
                tooltipContent: function (d) {
                    return '<h4>' + d.source.id + ' &gt; ' + d.target.id + ': ' + d.value + '</h4>'
                }
            });

        if (fooMutations.length) {
            x.histogram('his', hist, {
                innerRadius: 1.01,
                outerRadius: 1.4,
                color: 'YlOrRd',
                logScale: true,
                tooltipContent: function (d) {
                    return d.id + ": " + d.value
                }
            })
        }


        circos.render();


        d3.selectAll("path.chord").on("mouseover", handleMouseOver);
        d3.selectAll("path.chord").on("mouseout", handleMouseOut);


        d3.selectAll(".cs-layout path").on("mouseover", handleMouseOver);
        d3.selectAll(".cs-layout path").on("mouseout", handleMouseOut);

    }

    const render2 = function (chainA, chainB) {
        let residues = chainA.residues.concat(chainB.residues);

        const colors = d3.scaleOrdinal(d3.schemeCategory10);
        const alphabet = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"];// dirty trick 

        let data = residues.map((r, i) => {
            return {
                id: RES + r.id,
                chain: r.chain,
                len: 1,
                label: r.resName,
                color: colors(alphabet.indexOf(r.chain.toLowerCase()) + 1)
            }
        });

        let chords = [];

        residues.forEach(r => {
            r.exlinks.forEach(l => {
                if (l.distance > 0 && l.distance <= MAX_EX_DISTANCE) {
                    const from = {
                        id: RES + r.id,
                        start: .5,
                        end: .6
                    };
                    chords.push({
                        value: l.distance,
                        source: from,
                        target: {
                            id: RES + l.to,
                            start: .5,
                            end: .6
                        }
                    })
                }
            });
        })


        let circos = new Circos({
            container: '#' + divID,
            width: width,
            height: height,
        });

        let configuration = {
            innerRadius: innerRaduis,
            outerRadius: radius,
            class: "residue",
            cornerRadius: 3,
            gap: 0.0001, // in radian
            labels: {
                display: false
            },
            ticks: {
                display: false
            }
        }


        let x = circos
            .layout(data, configuration)
            .text('labels', residues.map(function (d) {
                return {
                    block_id: RES + d.id,
                    position: .5,
                    value: d.chain + "|" + d.resName + ":" + (d.id + 1)
                }
            }), {
                innerRadius: .9,
                outerRadius: 1.3,
                style: {
                    'font-size': 8
                },
                tooltipContent: function (d) {
                    return d.value
                },
                events: {
                    'click.alert': function (datum, index, nodes, event) {
                        window.alert(datum)
                    }
                }
            })
            .chords("l1",
                chords, {
                opacity: 0.4,
                color: 'grey',
                tooltipContent: function (d) {
                    return '<h4>' + d.source.id + ' &gt; ' + d.target.id + ': ' + d.value + '</h4>'
                }
            });

        if (fooMutations.length) {
            x.histogram('his', hist, {
                innerRadius: 1.01,
                outerRadius: 1.4,
                color: 'YlOrRd',
                logScale: true,
                tooltipContent: function (d) {
                    return d.id + ": " + d.value
                }
            })
        }


        circos.render();

        d3.selectAll("path.chord").on("mouseover", handleMouseOver);
        d3.selectAll("path.chord").on("mouseout", handleMouseOut);


        d3.selectAll(".cs-layout path").on("mouseover", handleMouseOver);
        d3.selectAll(".cs-layout path").on("mouseout", handleMouseOut);

    }

    // --------------MAIN PART---------------

    d3.select("#" + divID).select("svg").remove();

    let fooMutations = [];
    if (pdbID === '6W9C') {
        fooMutations = [
            {
                "value": 124,
                "id": 249
            },
            {
                "value": 111,
                "id": 240
            },
            {
                "value": 63,
                "id": 221
            },
            {
                "value": 55,
                "id": 206
            },
            {
                "value": 48,
                "id": 207
            },
            {
                "value": 45,
                "id": 36
            },
            {
                "value": 44,
                "id": 259
            },
            {
                "value": 41,
                "id": 63
            },
            {
                "value": 36,
                "id": 77
            },
            {
                "value": 35,
                "id": 277
            },
            {
                "value": 31,
                "id": 19
            },
            {
                "value": 30,
                "id": 171
            },
            {
                "value": 29,
                "id": 119
            },
            {
                "value": 28,
                "id": 92
            },
            {
                "value": 26,
                "id": 288
            },
            {
                "value": 25,
                "id": 291
            },
            {
                "value": 22,
                "id": 255
            },
            {
                "value": 20,
                "id": 77
            },
            {
                "value": 19,
                "id": 84
            },
            {
                "value": 19,
                "id": 157
            },
            {
                "value": 16,
                "id": 289
            },
            {
                "value": 14,
                "id": 74
            },
            {
                "value": 13,
                "id": 137
            },
            {
                "value": 13,
                "id": 223
            },
            {
                "value": 12,
                "id": 4
            },
            {
                "value": 12,
                "id": 195
            },
            {
                "value": 12,
                "id": 301
            },
            {
                "value": 11,
                "id": 69
            },
            {
                "value": 11,
                "id": 129
            },
            {
                "value": 10,
                "id": 261
            },
            {
                "value": 10,
                "id": 263
            },
            {
                "value": 9,
                "id": 34
            },
            {
                "value": 9,
                "id": 75
            },
            {
                "value": 9,
                "id": 116
            },
            {
                "value": 9,
                "id": 254
            },
            {
                "value": 8,
                "id": 44
            },
            {
                "value": 8,
                "id": 86
            },
            {
                "value": 7,
                "id": 120
            },
            {
                "value": 7,
                "id": 191
            },
            {
                "value": 7,
                "id": 197
            },
            {
                "value": 7,
                "id": 200
            },
            {
                "value": 7,
                "id": 223
            },
            {
                "value": 7,
                "id": 247
            },
            {
                "value": 7,
                "id": 292
            },
            {
                "value": 6,
                "id": 18
            },
            {
                "value": 6,
                "id": 33
            },
            {
                "value": 6,
                "id": 61
            },
            {
                "value": 6,
                "id": 153
            },
            {
                "value": 6,
                "id": 198
            },
            {
                "value": 6,
                "id": 231
            },
            {
                "value": 6,
                "id": 299
            },
            {
                "value": 5,
                "id": 3
            },
            {
                "value": 5,
                "id": 26
            },
            {
                "value": 5,
                "id": 62
            },
            {
                "value": 5,
                "id": 68
            },
            {
                "value": 5,
                "id": 119
            },
            {
                "value": 5,
                "id": 123
            },
            {
                "value": 5,
                "id": 200
            },
            {
                "value": 5,
                "id": 203
            },
            {
                "value": 5,
                "id": 275
            },
            {
                "value": 5,
                "id": 290
            },
            {
                "value": 5,
                "id": 311
            },
            {
                "value": 4,
                "id": 42
            },
            {
                "value": 4,
                "id": 65
            },
            {
                "value": 4,
                "id": 66
            },
            {
                "value": 4,
                "id": 80
            },
            {
                "value": 4,
                "id": 129
            },
            {
                "value": 4,
                "id": 211
            },
            {
                "value": 4,
                "id": 222
            },
            {
                "value": 4,
                "id": 257
            },
            {
                "value": 3,
                "id": 20
            },
            {
                "value": 3,
                "id": 28
            },
            {
                "value": 3,
                "id": 39
            },
            {
                "value": 3,
                "id": 47
            },
            {
                "value": 3,
                "id": 54
            },
            {
                "value": 3,
                "id": 76
            },
            {
                "value": 3,
                "id": 86
            },
            {
                "value": 3,
                "id": 89
            },
            {
                "value": 3,
                "id": 108
            },
            {
                "value": 3,
                "id": 109
            },
            {
                "value": 3,
                "id": 130
            },
            {
                "value": 3,
                "id": 146
            },
            {
                "value": 3,
                "id": 168
            },
            {
                "value": 3,
                "id": 187
            },
            {
                "value": 3,
                "id": 196
            },
            {
                "value": 3,
                "id": 222
            },
            {
                "value": 3,
                "id": 246
            },
            {
                "value": 3,
                "id": 247
            },
            {
                "value": 3,
                "id": 250
            },
            {
                "value": 3,
                "id": 295
            },
            {
                "value": 3,
                "id": 296
            },
            {
                "value": 3,
                "id": 297
            },
            {
                "value": 3,
                "id": 308
            },
            {
                "value": 3,
                "id": 309
            },
            {
                "value": 2,
                "id": 9
            },
            {
                "value": 2,
                "id": 13
            },
            {
                "value": 2,
                "id": 16
            },
            {
                "value": 2,
                "id": 20
            },
            {
                "value": 2,
                "id": 22
            },
            {
                "value": 2,
                "id": 22
            },
            {
                "value": 2,
                "id": 24
            },
            {
                "value": 2,
                "id": 25
            },
            {
                "value": 2,
                "id": 28
            },
            {
                "value": 2,
                "id": 29
            },
            {
                "value": 2,
                "id": 33
            },
            {
                "value": 2,
                "id": 37
            },
            {
                "value": 2,
                "id": 37
            },
            {
                "value": 2,
                "id": 47
            },
            {
                "value": 2,
                "id": 49
            },
            {
                "value": 2,
                "id": 56
            },
            {
                "value": 2,
                "id": 57
            },
            {
                "value": 2,
                "id": 71
            },
            {
                "value": 2,
                "id": 74
            },
            {
                "value": 2,
                "id": 75
            },
            {
                "value": 2,
                "id": 82
            },
            {
                "value": 2,
                "id": 90
            },
            {
                "value": 2,
                "id": 92
            },
            {
                "value": 2,
                "id": 96
            },
            {
                "value": 2,
                "id": 100
            },
            {
                "value": 2,
                "id": 104
            },
            {
                "value": 2,
                "id": 108
            },
            {
                "value": 2,
                "id": 126
            },
            {
                "value": 2,
                "id": 138
            },
            {
                "value": 2,
                "id": 149
            },
            {
                "value": 2,
                "id": 156
            },
            {
                "value": 2,
                "id": 161
            },
            {
                "value": 2,
                "id": 169
            },
            {
                "value": 2,
                "id": 170
            },
            {
                "value": 2,
                "id": 177
            },
            {
                "value": 2,
                "id": 196
            },
            {
                "value": 2,
                "id": 201
            },
            {
                "value": 2,
                "id": 206
            },
            {
                "value": 2,
                "id": 208
            },
            {
                "value": 2,
                "id": 225
            },
            {
                "value": 2,
                "id": 235
            },
            {
                "value": 2,
                "id": 281
            },
            {
                "value": 2,
                "id": 289
            },
            {
                "value": 2,
                "id": 293
            },
            {
                "value": 2,
                "id": 313
            },
            {
                "value": 1,
                "id": 10
            },
            {
                "value": 1,
                "id": 11
            },
            {
                "value": 1,
                "id": 12
            },
            {
                "value": 1,
                "id": 19
            },
            {
                "value": 1,
                "id": 25
            },
            {
                "value": 1,
                "id": 29
            },
            {
                "value": 1,
                "id": 49
            },
            {
                "value": 1,
                "id": 50
            },
            {
                "value": 1,
                "id": 51
            },
            {
                "value": 1,
                "id": 54
            },
            {
                "value": 1,
                "id": 55
            },
            {
                "value": 1,
                "id": 55
            },
            {
                "value": 1,
                "id": 65
            },
            {
                "value": 1,
                "id": 66
            },
            {
                "value": 1,
                "id": 73
            },
            {
                "value": 1,
                "id": 90
            },
            {
                "value": 1,
                "id": 94
            },
            {
                "value": 1,
                "id": 98
            },
            {
                "value": 1,
                "id": 102
            },
            {
                "value": 1,
                "id": 107
            },
            {
                "value": 1,
                "id": 110
            },
            {
                "value": 1,
                "id": 112
            },
            {
                "value": 1,
                "id": 115
            },
            {
                "value": 1,
                "id": 124
            },
            {
                "value": 1,
                "id": 126
            },
            {
                "value": 1,
                "id": 128
            },
            {
                "value": 1,
                "id": 135
            },
            {
                "value": 1,
                "id": 144
            },
            {
                "value": 1,
                "id": 144
            },
            {
                "value": 1,
                "id": 151
            },
            {
                "value": 1,
                "id": 151
            },
            {
                "value": 1,
                "id": 152
            },
            {
                "value": 1,
                "id": 159
            },
            {
                "value": 1,
                "id": 160
            },
            {
                "value": 1,
                "id": 175
            },
            {
                "value": 1,
                "id": 202
            },
            {
                "value": 1,
                "id": 204
            },
            {
                "value": 1,
                "id": 205
            },
            {
                "value": 1,
                "id": 209
            },
            {
                "value": 1,
                "id": 219
            },
            {
                "value": 1,
                "id": 220
            },
            {
                "value": 1,
                "id": 222
            },
            {
                "value": 1,
                "id": 228
            },
            {
                "value": 1,
                "id": 230
            },
            {
                "value": 1,
                "id": 235
            },
            {
                "value": 1,
                "id": 238
            },
            {
                "value": 1,
                "id": 241
            },
            {
                "value": 1,
                "id": 242
            },
            {
                "value": 1,
                "id": 248
            },
            {
                "value": 1,
                "id": 254
            },
            {
                "value": 1,
                "id": 255
            },
            {
                "value": 1,
                "id": 263
            },
            {
                "value": 1,
                "id": 267
            },
            {
                "value": 1,
                "id": 276
            },
            {
                "value": 1,
                "id": 278
            },
            {
                "value": 1,
                "id": 285
            },
            {
                "value": 1,
                "id": 291
            },
            {
                "value": 1,
                "id": 298
            },
            {
                "value": 1,
                "id": 299
            },
            {
                "value": 1,
                "id": 300
            },
            {
                "value": 1,
                "id": 303
            },
            {
                "value": 1,
                "id": 306
            },
            {
                "value": 1,
                "id": 307
            },
            {
                "value": 1,
                "id": 313
            }
        ];
    }

    const CA = "CA";

    const MIN_DISTANCE = 0;
    let MAX_DISTANCE = maxDistance || 5;
    let MAX_EX_DISTANCE = maxDistance || 30;

    const RES = "res_"; // html id prefix

    let multiple = chainIDs[0] != chainIDs[1];

    let hist = fooMutations.map(m => {
        m.block_id = RES + m.id;
        m.start = .3;
        m.end = .7;
        return m;
    });

    if (multiple) {
        chainIDs = chainIDs ? chainIDs : ["A", "B"]

        const chA = chains.get(chainIDs[0])
        const chB = chains.get(chainIDs[1])
        calcExternalDistances(chA, chB, chainIDs)
        render2(chA, chB)
    } else {
        const chA = chains.get("A")
        calcInternalDistances(chA);
        render(chA)
    }

}


