OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(0.15396804) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(4.1339388) q[1];
sx q[1];
rz(9.0864656) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6217123) q[0];
sx q[0];
rz(-1.8147239) q[0];
sx q[0];
rz(-1.2432616) q[0];
rz(-pi) q[1];
rz(-2.2489684) q[2];
sx q[2];
rz(-1.8632338) q[2];
sx q[2];
rz(0.48722789) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.675128) q[1];
sx q[1];
rz(-0.29106859) q[1];
sx q[1];
rz(2.9558099) q[1];
x q[2];
rz(0.015720856) q[3];
sx q[3];
rz(-1.058488) q[3];
sx q[3];
rz(-2.6920464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.9677229) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(-0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8006111) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(-1.8992791) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500538) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(-0.48135249) q[0];
x q[1];
rz(-0.1588891) q[2];
sx q[2];
rz(-0.94413589) q[2];
sx q[2];
rz(1.5140669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9638483) q[1];
sx q[1];
rz(-1.4114393) q[1];
sx q[1];
rz(0.53667712) q[1];
rz(-pi) q[2];
rz(2.0598689) q[3];
sx q[3];
rz(-0.92182577) q[3];
sx q[3];
rz(1.058941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3339281) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(2.5578965) q[2];
rz(0.57404533) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(-0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7211001) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(-0.72845355) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(-1.0167936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66660488) q[0];
sx q[0];
rz(-0.089086108) q[0];
sx q[0];
rz(2.7276917) q[0];
rz(-pi) q[1];
rz(-0.74480199) q[2];
sx q[2];
rz(-2.475127) q[2];
sx q[2];
rz(-1.9772066) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9079202) q[1];
sx q[1];
rz(-0.83041149) q[1];
sx q[1];
rz(-0.18632142) q[1];
rz(-pi) q[2];
rz(-1.73909) q[3];
sx q[3];
rz(-2.8561391) q[3];
sx q[3];
rz(1.1807549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-2.0920848) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8124354) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.6695492) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-0.24681117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26537672) q[0];
sx q[0];
rz(-2.2582158) q[0];
sx q[0];
rz(2.9486297) q[0];
rz(-pi) q[1];
rz(2.6631782) q[2];
sx q[2];
rz(-1.3832428) q[2];
sx q[2];
rz(-1.071655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2553195) q[1];
sx q[1];
rz(-1.6233994) q[1];
sx q[1];
rz(0.37087755) q[1];
rz(-pi) q[2];
x q[2];
rz(2.85535) q[3];
sx q[3];
rz(-2.1677368) q[3];
sx q[3];
rz(0.7569353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4776769) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(-1.56303) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(0.24838233) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029593) q[0];
sx q[0];
rz(-0.54134936) q[0];
sx q[0];
rz(3.0754509) q[0];
rz(-pi) q[1];
rz(-2.9551198) q[2];
sx q[2];
rz(-1.4154134) q[2];
sx q[2];
rz(-1.1764256) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0866962) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(-0.30026786) q[1];
x q[2];
rz(0.6487209) q[3];
sx q[3];
rz(-1.3789163) q[3];
sx q[3];
rz(-2.3313525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-0.50271547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0080863) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(2.0603518) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(3.016901) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42449441) q[0];
sx q[0];
rz(-1.4846804) q[0];
sx q[0];
rz(2.9647102) q[0];
rz(1.2776676) q[2];
sx q[2];
rz(-0.98205245) q[2];
sx q[2];
rz(-2.6063906) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9166959) q[1];
sx q[1];
rz(-1.6273013) q[1];
sx q[1];
rz(-2.4042261) q[1];
rz(-pi) q[2];
rz(-0.74495875) q[3];
sx q[3];
rz(-2.8248441) q[3];
sx q[3];
rz(-1.6572286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51320118) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(2.52264) q[2];
rz(-2.0882873) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597647) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(1.0429617) q[0];
rz(-2.6783121) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(-1.0707062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8585513) q[0];
sx q[0];
rz(-1.0867449) q[0];
sx q[0];
rz(-1.5714684) q[0];
rz(-pi) q[1];
rz(-0.3473862) q[2];
sx q[2];
rz(-1.4457448) q[2];
sx q[2];
rz(-1.320968) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5412022) q[1];
sx q[1];
rz(-1.4711079) q[1];
sx q[1];
rz(0.77460918) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0827818) q[3];
sx q[3];
rz(-2.8088514) q[3];
sx q[3];
rz(-2.9174093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-0.94747296) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5652126) q[0];
sx q[0];
rz(-2.1064261) q[0];
sx q[0];
rz(-3.0239848) q[0];
rz(-2.410789) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(-1.3256324) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87654963) q[1];
sx q[1];
rz(-1.3339086) q[1];
sx q[1];
rz(-0.9898647) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.794907) q[3];
sx q[3];
rz(-0.83804916) q[3];
sx q[3];
rz(2.8447062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5802713) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(-2.6718111) q[2];
rz(1.8404768) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-2.9387617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751574) q[0];
sx q[0];
rz(-3.0946819) q[0];
sx q[0];
rz(0.58880083) q[0];
x q[1];
rz(-2.3896396) q[2];
sx q[2];
rz(-2.8402036) q[2];
sx q[2];
rz(-1.6981268) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.087346615) q[1];
sx q[1];
rz(-3.0294703) q[1];
sx q[1];
rz(2.0263158) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8392302) q[3];
sx q[3];
rz(-2.6609169) q[3];
sx q[3];
rz(2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7245076) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(0.99651304) q[2];
rz(2.7881682) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(-2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(-3.0985447) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(-2.8607686) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.330864) q[0];
sx q[0];
rz(-2.4952336) q[0];
sx q[0];
rz(0.75673639) q[0];
rz(-pi) q[1];
rz(0.73239399) q[2];
sx q[2];
rz(-0.28389441) q[2];
sx q[2];
rz(0.64901272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22051375) q[1];
sx q[1];
rz(-1.6069357) q[1];
sx q[1];
rz(-1.6284579) q[1];
rz(-pi) q[2];
rz(-0.19631581) q[3];
sx q[3];
rz(-1.8972978) q[3];
sx q[3];
rz(-1.6895837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-0.62310702) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6476718) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(-0.80550823) q[2];
sx q[2];
rz(-1.1136354) q[2];
sx q[2];
rz(-1.8352933) q[2];
rz(0.86250967) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];