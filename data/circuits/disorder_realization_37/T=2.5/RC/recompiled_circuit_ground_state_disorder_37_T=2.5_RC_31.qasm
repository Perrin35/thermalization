OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3545544) q[0];
sx q[0];
rz(0.7631425) q[0];
sx q[0];
rz(7.2416303) q[0];
rz(0.42674843) q[1];
sx q[1];
rz(-0.42438212) q[1];
sx q[1];
rz(-0.97270614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866731) q[0];
sx q[0];
rz(-1.7717607) q[0];
sx q[0];
rz(1.6660287) q[0];
rz(-pi) q[1];
rz(-1.6390267) q[2];
sx q[2];
rz(-1.7426602) q[2];
sx q[2];
rz(1.9910938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7643775) q[1];
sx q[1];
rz(-0.7684349) q[1];
sx q[1];
rz(-2.7646116) q[1];
x q[2];
rz(-2.8593721) q[3];
sx q[3];
rz(-1.7190386) q[3];
sx q[3];
rz(0.46503286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0613609) q[2];
sx q[2];
rz(-1.9539359) q[2];
sx q[2];
rz(2.3415671) q[2];
rz(-3.0112265) q[3];
sx q[3];
rz(-2.3733807) q[3];
sx q[3];
rz(0.31726328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80113634) q[0];
sx q[0];
rz(-1.9213333) q[0];
sx q[0];
rz(2.1521547) q[0];
rz(2.2438352) q[1];
sx q[1];
rz(-2.0866626) q[1];
sx q[1];
rz(-0.78136939) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5637276) q[0];
sx q[0];
rz(-1.0672309) q[0];
sx q[0];
rz(0.2204075) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23405646) q[2];
sx q[2];
rz(-1.6566725) q[2];
sx q[2];
rz(1.9746922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.20618379) q[1];
sx q[1];
rz(-2.0061135) q[1];
sx q[1];
rz(-2.101442) q[1];
rz(-pi) q[2];
rz(-0.64919169) q[3];
sx q[3];
rz(-2.6731369) q[3];
sx q[3];
rz(2.2476803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.024461688) q[2];
sx q[2];
rz(-1.3776366) q[2];
sx q[2];
rz(0.77829877) q[2];
rz(-0.69027573) q[3];
sx q[3];
rz(-0.92167753) q[3];
sx q[3];
rz(1.5811623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2799601) q[0];
sx q[0];
rz(-2.3065688) q[0];
sx q[0];
rz(-2.7296208) q[0];
rz(-2.1906134) q[1];
sx q[1];
rz(-2.488766) q[1];
sx q[1];
rz(1.9297809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6854343) q[0];
sx q[0];
rz(-0.051101772) q[0];
sx q[0];
rz(0.85275485) q[0];
rz(-1.7936208) q[2];
sx q[2];
rz(-1.6671902) q[2];
sx q[2];
rz(2.4501462) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0101357) q[1];
sx q[1];
rz(-1.5984416) q[1];
sx q[1];
rz(-0.70998876) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36529593) q[3];
sx q[3];
rz(-1.72137) q[3];
sx q[3];
rz(0.07380658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1872306) q[2];
sx q[2];
rz(-0.53202859) q[2];
sx q[2];
rz(1.5032035) q[2];
rz(-3.1014118) q[3];
sx q[3];
rz(-0.92622042) q[3];
sx q[3];
rz(-0.23475501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16315854) q[0];
sx q[0];
rz(-1.6272767) q[0];
sx q[0];
rz(-3.0117595) q[0];
rz(1.2854598) q[1];
sx q[1];
rz(-1.1760271) q[1];
sx q[1];
rz(1.0349549) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1945788) q[0];
sx q[0];
rz(-2.0067375) q[0];
sx q[0];
rz(2.2313124) q[0];
rz(-pi) q[1];
rz(0.90702727) q[2];
sx q[2];
rz(-1.656678) q[2];
sx q[2];
rz(-2.2965388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.538547) q[1];
sx q[1];
rz(-1.6563452) q[1];
sx q[1];
rz(2.892832) q[1];
rz(-pi) q[2];
rz(1.0134936) q[3];
sx q[3];
rz(-0.77269197) q[3];
sx q[3];
rz(-1.1232291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7776362) q[2];
sx q[2];
rz(-1.7549606) q[2];
sx q[2];
rz(0.094495471) q[2];
rz(-2.7206521) q[3];
sx q[3];
rz(-1.7178444) q[3];
sx q[3];
rz(-1.4360992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41973758) q[0];
sx q[0];
rz(-0.57191816) q[0];
sx q[0];
rz(3.1255334) q[0];
rz(-0.84469604) q[1];
sx q[1];
rz(-2.7222996) q[1];
sx q[1];
rz(-0.90763456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3428872) q[0];
sx q[0];
rz(-1.3488975) q[0];
sx q[0];
rz(-2.368244) q[0];
rz(-pi) q[1];
rz(2.0727022) q[2];
sx q[2];
rz(-2.6703983) q[2];
sx q[2];
rz(1.0028749) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.221783) q[1];
sx q[1];
rz(-2.5804829) q[1];
sx q[1];
rz(-2.2235548) q[1];
x q[2];
rz(1.5974371) q[3];
sx q[3];
rz(-1.9147885) q[3];
sx q[3];
rz(0.87061239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57546651) q[2];
sx q[2];
rz(-1.0593654) q[2];
sx q[2];
rz(0.063684138) q[2];
rz(-0.56695402) q[3];
sx q[3];
rz(-2.3072115) q[3];
sx q[3];
rz(-2.5440192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8276234) q[0];
sx q[0];
rz(-0.14179985) q[0];
sx q[0];
rz(1.1578479) q[0];
rz(-1.4843548) q[1];
sx q[1];
rz(-1.4446222) q[1];
sx q[1];
rz(2.8725913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6873574) q[0];
sx q[0];
rz(-1.7116065) q[0];
sx q[0];
rz(1.6160377) q[0];
rz(1.2397175) q[2];
sx q[2];
rz(-0.28713687) q[2];
sx q[2];
rz(-1.3763217) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3693302) q[1];
sx q[1];
rz(-2.9360665) q[1];
sx q[1];
rz(-2.0783246) q[1];
rz(-pi) q[2];
rz(2.6746142) q[3];
sx q[3];
rz(-0.45736936) q[3];
sx q[3];
rz(-2.6757095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77702648) q[2];
sx q[2];
rz(-0.51743162) q[2];
sx q[2];
rz(2.2606134) q[2];
rz(3.0173054) q[3];
sx q[3];
rz(-1.7139939) q[3];
sx q[3];
rz(0.65829268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86150259) q[0];
sx q[0];
rz(-0.26743356) q[0];
sx q[0];
rz(-0.68914831) q[0];
rz(0.46734494) q[1];
sx q[1];
rz(-0.57599774) q[1];
sx q[1];
rz(-2.0435832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91749597) q[0];
sx q[0];
rz(-0.18315975) q[0];
sx q[0];
rz(-0.44395019) q[0];
rz(-pi) q[1];
rz(1.271209) q[2];
sx q[2];
rz(-1.0429842) q[2];
sx q[2];
rz(-0.53003788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6642521) q[1];
sx q[1];
rz(-2.1108529) q[1];
sx q[1];
rz(-1.5332721) q[1];
rz(-pi) q[2];
rz(-2.1995898) q[3];
sx q[3];
rz(-1.1515472) q[3];
sx q[3];
rz(-2.2908831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5757307) q[2];
sx q[2];
rz(-2.2634759) q[2];
sx q[2];
rz(-2.511054) q[2];
rz(-2.5343043) q[3];
sx q[3];
rz(-2.2346965) q[3];
sx q[3];
rz(-2.0489395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8232515) q[0];
sx q[0];
rz(-0.14597758) q[0];
sx q[0];
rz(-1.7952221) q[0];
rz(1.775555) q[1];
sx q[1];
rz(-0.64907688) q[1];
sx q[1];
rz(2.5128561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36489428) q[0];
sx q[0];
rz(-1.8275675) q[0];
sx q[0];
rz(2.5904694) q[0];
rz(0.89968483) q[2];
sx q[2];
rz(-1.9860528) q[2];
sx q[2];
rz(-0.64476162) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6639316) q[1];
sx q[1];
rz(-1.0870502) q[1];
sx q[1];
rz(0.61772924) q[1];
rz(1.3717184) q[3];
sx q[3];
rz(-1.4961096) q[3];
sx q[3];
rz(-1.5916361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8667355) q[2];
sx q[2];
rz(-1.3254415) q[2];
sx q[2];
rz(2.2692915) q[2];
rz(-2.9288779) q[3];
sx q[3];
rz(-1.3958967) q[3];
sx q[3];
rz(-2.5581636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23671959) q[0];
sx q[0];
rz(-1.5308335) q[0];
sx q[0];
rz(-1.8778296) q[0];
rz(2.1243375) q[1];
sx q[1];
rz(-0.89824289) q[1];
sx q[1];
rz(-2.5668626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86704937) q[0];
sx q[0];
rz(-2.6066271) q[0];
sx q[0];
rz(-0.127129) q[0];
rz(-pi) q[1];
x q[1];
rz(1.310964) q[2];
sx q[2];
rz(-0.75922478) q[2];
sx q[2];
rz(0.34789839) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8239261) q[1];
sx q[1];
rz(-1.9374018) q[1];
sx q[1];
rz(1.3701141) q[1];
x q[2];
rz(1.0887126) q[3];
sx q[3];
rz(-1.7974554) q[3];
sx q[3];
rz(0.72666336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42402521) q[2];
sx q[2];
rz(-1.307345) q[2];
sx q[2];
rz(0.032111017) q[2];
rz(1.9854246) q[3];
sx q[3];
rz(-0.38438946) q[3];
sx q[3];
rz(-1.2110075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0403989) q[0];
sx q[0];
rz(-2.5960584) q[0];
sx q[0];
rz(1.1241166) q[0];
rz(2.595937) q[1];
sx q[1];
rz(-2.7898495) q[1];
sx q[1];
rz(-0.55145946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1301981) q[0];
sx q[0];
rz(-1.4470248) q[0];
sx q[0];
rz(-1.3817203) q[0];
rz(1.2966424) q[2];
sx q[2];
rz(-1.621843) q[2];
sx q[2];
rz(0.54948839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.76186) q[1];
sx q[1];
rz(-1.1071536) q[1];
sx q[1];
rz(-2.597888) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48088185) q[3];
sx q[3];
rz(-0.34264287) q[3];
sx q[3];
rz(-0.51183701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55190279) q[2];
sx q[2];
rz(-2.6305113) q[2];
sx q[2];
rz(-1.7299293) q[2];
rz(-2.3949413) q[3];
sx q[3];
rz(-1.6801497) q[3];
sx q[3];
rz(2.1667229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31155561) q[0];
sx q[0];
rz(-1.5604326) q[0];
sx q[0];
rz(-1.5695705) q[0];
rz(-2.742782) q[1];
sx q[1];
rz(-1.3347944) q[1];
sx q[1];
rz(1.4904108) q[1];
rz(-1.5854992) q[2];
sx q[2];
rz(-2.6530545) q[2];
sx q[2];
rz(-1.1966111) q[2];
rz(0.29295425) q[3];
sx q[3];
rz(-2.3200547) q[3];
sx q[3];
rz(-3.0153081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
