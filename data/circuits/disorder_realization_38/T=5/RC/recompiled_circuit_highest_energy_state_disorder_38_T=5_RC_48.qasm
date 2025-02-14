OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6925447) q[0];
sx q[0];
rz(-1.1392051) q[0];
sx q[0];
rz(-0.44266242) q[0];
rz(0.17065419) q[1];
sx q[1];
rz(-2.3499188) q[1];
sx q[1];
rz(0.72905529) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937482) q[0];
sx q[0];
rz(-1.9707075) q[0];
sx q[0];
rz(-2.4104959) q[0];
x q[1];
rz(0.38972028) q[2];
sx q[2];
rz(-1.5699631) q[2];
sx q[2];
rz(0.69609797) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9022791) q[1];
sx q[1];
rz(-2.2444176) q[1];
sx q[1];
rz(-0.31142733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0341845) q[3];
sx q[3];
rz(-0.70176673) q[3];
sx q[3];
rz(-0.90493363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2810104) q[2];
sx q[2];
rz(-0.082110114) q[2];
sx q[2];
rz(2.3561467) q[2];
rz(-0.9663409) q[3];
sx q[3];
rz(-1.0680501) q[3];
sx q[3];
rz(0.6399703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6099089) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(-0.63572836) q[0];
rz(-0.6334148) q[1];
sx q[1];
rz(-2.3737291) q[1];
sx q[1];
rz(2.8107218) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2658871) q[0];
sx q[0];
rz(-0.71015753) q[0];
sx q[0];
rz(0.52010923) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.10373) q[2];
sx q[2];
rz(-1.5469917) q[2];
sx q[2];
rz(-1.3184774) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65407359) q[1];
sx q[1];
rz(-1.6292058) q[1];
sx q[1];
rz(-0.033694555) q[1];
x q[2];
rz(-2.9959566) q[3];
sx q[3];
rz(-3.0560827) q[3];
sx q[3];
rz(0.95216361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50515455) q[2];
sx q[2];
rz(-2.324489) q[2];
sx q[2];
rz(-0.46112296) q[2];
rz(-2.4165706) q[3];
sx q[3];
rz(-2.6889763) q[3];
sx q[3];
rz(-2.4729474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4918936) q[0];
sx q[0];
rz(-1.3839586) q[0];
sx q[0];
rz(-2.6032676) q[0];
rz(-0.0093983924) q[1];
sx q[1];
rz(-2.5098269) q[1];
sx q[1];
rz(1.9422772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.050735) q[0];
sx q[0];
rz(-1.5260625) q[0];
sx q[0];
rz(-1.0629774) q[0];
rz(-pi) q[1];
rz(0.17260562) q[2];
sx q[2];
rz(-0.69059082) q[2];
sx q[2];
rz(-1.1081072) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9925132) q[1];
sx q[1];
rz(-2.1285526) q[1];
sx q[1];
rz(2.7191619) q[1];
x q[2];
rz(-1.2869419) q[3];
sx q[3];
rz(-0.7470658) q[3];
sx q[3];
rz(3.0892854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.37375307) q[2];
sx q[2];
rz(-1.6497296) q[2];
sx q[2];
rz(-2.4191432) q[2];
rz(0.59822285) q[3];
sx q[3];
rz(-3.1066419) q[3];
sx q[3];
rz(1.2178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5825321) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(-3.1225358) q[0];
rz(1.6825153) q[1];
sx q[1];
rz(-0.2225114) q[1];
sx q[1];
rz(0.51024514) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2697578) q[0];
sx q[0];
rz(-1.8352011) q[0];
sx q[0];
rz(-2.1931936) q[0];
rz(2.4195477) q[2];
sx q[2];
rz(-1.8068442) q[2];
sx q[2];
rz(-1.2942734) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4245119) q[1];
sx q[1];
rz(-1.1649017) q[1];
sx q[1];
rz(-2.4688429) q[1];
x q[2];
rz(-1.4278379) q[3];
sx q[3];
rz(-1.1605439) q[3];
sx q[3];
rz(0.15907915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14429188) q[2];
sx q[2];
rz(-1.4520626) q[2];
sx q[2];
rz(-2.9193381) q[2];
rz(3.0620767) q[3];
sx q[3];
rz(-0.54602081) q[3];
sx q[3];
rz(-2.5047746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52614373) q[0];
sx q[0];
rz(-1.0960824) q[0];
sx q[0];
rz(-1.3413606) q[0];
rz(-0.51097956) q[1];
sx q[1];
rz(-1.7781517) q[1];
sx q[1];
rz(2.708639) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47577417) q[0];
sx q[0];
rz(-1.5137714) q[0];
sx q[0];
rz(2.8840469) q[0];
rz(-1.1002925) q[2];
sx q[2];
rz(-0.44657257) q[2];
sx q[2];
rz(-0.64556015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.49952415) q[1];
sx q[1];
rz(-2.7244901) q[1];
sx q[1];
rz(0.92315556) q[1];
x q[2];
rz(-0.28750395) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(2.7903729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0114835) q[2];
sx q[2];
rz(-0.92769647) q[2];
sx q[2];
rz(-0.9978655) q[2];
rz(-2.4326371) q[3];
sx q[3];
rz(-2.7537789) q[3];
sx q[3];
rz(0.97878218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3157432) q[0];
sx q[0];
rz(-0.5873) q[0];
sx q[0];
rz(2.24776) q[0];
rz(1.029344) q[1];
sx q[1];
rz(-2.4004332) q[1];
sx q[1];
rz(0.11480521) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.608716) q[0];
sx q[0];
rz(-2.8136557) q[0];
sx q[0];
rz(-0.84285835) q[0];
rz(0.64115142) q[2];
sx q[2];
rz(-1.8617407) q[2];
sx q[2];
rz(-0.095130446) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56267589) q[1];
sx q[1];
rz(-1.8496184) q[1];
sx q[1];
rz(1.2827966) q[1];
rz(-pi) q[2];
rz(-0.28088681) q[3];
sx q[3];
rz(-1.509519) q[3];
sx q[3];
rz(2.3063502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7041695) q[2];
sx q[2];
rz(-0.88714522) q[2];
sx q[2];
rz(-0.21370055) q[2];
rz(0.63654381) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(-0.73591939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8146424) q[0];
sx q[0];
rz(-0.7433759) q[0];
sx q[0];
rz(-3.0315234) q[0];
rz(-0.07668177) q[1];
sx q[1];
rz(-2.2814467) q[1];
sx q[1];
rz(-2.0032047) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0698893) q[0];
sx q[0];
rz(-1.2217297) q[0];
sx q[0];
rz(-1.8047234) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78779548) q[2];
sx q[2];
rz(-1.5681364) q[2];
sx q[2];
rz(3.015369) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0323022) q[1];
sx q[1];
rz(-1.7545953) q[1];
sx q[1];
rz(1.079783) q[1];
x q[2];
rz(1.7848694) q[3];
sx q[3];
rz(-1.8118637) q[3];
sx q[3];
rz(1.9578707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.25973213) q[2];
sx q[2];
rz(-1.0819819) q[2];
sx q[2];
rz(-0.38121769) q[2];
rz(0.11738736) q[3];
sx q[3];
rz(-0.50506794) q[3];
sx q[3];
rz(-3.05486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36670244) q[0];
sx q[0];
rz(-0.77228868) q[0];
sx q[0];
rz(-2.6430687) q[0];
rz(-2.3122834) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(0.097749762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5978591) q[0];
sx q[0];
rz(-1.2654788) q[0];
sx q[0];
rz(-1.1174563) q[0];
rz(0.4011863) q[2];
sx q[2];
rz(-0.36683914) q[2];
sx q[2];
rz(2.7374379) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5619288) q[1];
sx q[1];
rz(-1.0439928) q[1];
sx q[1];
rz(2.4234523) q[1];
rz(-pi) q[2];
rz(3.1322543) q[3];
sx q[3];
rz(-1.8923602) q[3];
sx q[3];
rz(3.0967979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30663124) q[2];
sx q[2];
rz(-1.9913048) q[2];
sx q[2];
rz(-0.59070307) q[2];
rz(0.090204209) q[3];
sx q[3];
rz(-0.43961757) q[3];
sx q[3];
rz(0.91501594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0331994) q[0];
sx q[0];
rz(-0.84302253) q[0];
sx q[0];
rz(2.3868308) q[0];
rz(-0.33590487) q[1];
sx q[1];
rz(-0.55481189) q[1];
sx q[1];
rz(-2.1051443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91387902) q[0];
sx q[0];
rz(-2.4464985) q[0];
sx q[0];
rz(2.9577425) q[0];
rz(-pi) q[1];
rz(0.48123863) q[2];
sx q[2];
rz(-1.5889152) q[2];
sx q[2];
rz(-2.6157524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7215342) q[1];
sx q[1];
rz(-2.0871401) q[1];
sx q[1];
rz(-1.615585) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.127946) q[3];
sx q[3];
rz(-1.8004775) q[3];
sx q[3];
rz(1.0245293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0466517) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(-0.45647344) q[2];
rz(2.5392695) q[3];
sx q[3];
rz(-0.18776247) q[3];
sx q[3];
rz(-2.203234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8896821) q[0];
sx q[0];
rz(-1.2924117) q[0];
sx q[0];
rz(0.46345261) q[0];
rz(0.015333029) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(-0.32036805) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.244359) q[0];
sx q[0];
rz(-2.9787052) q[0];
sx q[0];
rz(0.49400671) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66530975) q[2];
sx q[2];
rz(-0.39145601) q[2];
sx q[2];
rz(-0.91400439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55904222) q[1];
sx q[1];
rz(-2.1782785) q[1];
sx q[1];
rz(0.59848018) q[1];
rz(-pi) q[2];
rz(1.2804561) q[3];
sx q[3];
rz(-1.2429894) q[3];
sx q[3];
rz(-0.39738712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1128803) q[2];
sx q[2];
rz(-2.4148603) q[2];
sx q[2];
rz(-2.1325364) q[2];
rz(-2.3446337) q[3];
sx q[3];
rz(-1.8877441) q[3];
sx q[3];
rz(0.33826452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0491895) q[0];
sx q[0];
rz(-1.6272463) q[0];
sx q[0];
rz(1.8820681) q[0];
rz(3.0567723) q[1];
sx q[1];
rz(-1.4823352) q[1];
sx q[1];
rz(-1.7099554) q[1];
rz(1.8126429) q[2];
sx q[2];
rz(-1.6326661) q[2];
sx q[2];
rz(2.2018955) q[2];
rz(-2.8390213) q[3];
sx q[3];
rz(-1.745001) q[3];
sx q[3];
rz(0.59636084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
