OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(-2.6383658) q[0];
sx q[0];
rz(0.72416645) q[0];
rz(0.63996285) q[1];
sx q[1];
rz(-0.53007403) q[1];
sx q[1];
rz(-0.78483265) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0137579) q[0];
sx q[0];
rz(-2.7779967) q[0];
sx q[0];
rz(2.512393) q[0];
x q[1];
rz(-2.7230524) q[2];
sx q[2];
rz(-1.4782018) q[2];
sx q[2];
rz(-1.6469524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8514429) q[1];
sx q[1];
rz(-2.5622497) q[1];
sx q[1];
rz(1.9860553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9853398) q[3];
sx q[3];
rz(-1.6016377) q[3];
sx q[3];
rz(0.24658345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.589754) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215866) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(-2.8979229) q[0];
rz(0.63175732) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.7858645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49039098) q[0];
sx q[0];
rz(-0.25679195) q[0];
sx q[0];
rz(1.4635565) q[0];
x q[1];
rz(-0.35661125) q[2];
sx q[2];
rz(-1.9372254) q[2];
sx q[2];
rz(-2.8474142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.3416459) q[1];
sx q[1];
rz(-0.46657545) q[1];
sx q[1];
rz(1.8599618) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59216604) q[3];
sx q[3];
rz(-1.5321931) q[3];
sx q[3];
rz(1.5241227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(0.24965723) q[2];
rz(0.50659531) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(2.8095968) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24519414) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(-2.2431592) q[0];
rz(-1.8067182) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(-1.2737087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0262336) q[0];
sx q[0];
rz(-0.44239487) q[0];
sx q[0];
rz(-0.94985234) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8286335) q[2];
sx q[2];
rz(-0.3096146) q[2];
sx q[2];
rz(1.6329873) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0395567) q[1];
sx q[1];
rz(-1.2219056) q[1];
sx q[1];
rz(-0.7015014) q[1];
rz(-2.418163) q[3];
sx q[3];
rz(-1.9577141) q[3];
sx q[3];
rz(0.5667516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0818103) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(0.95345062) q[2];
rz(0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(-2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2801441) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(3.0674556) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0291245) q[0];
sx q[0];
rz(-1.0974786) q[0];
sx q[0];
rz(2.8549457) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7934302) q[2];
sx q[2];
rz(-2.2194038) q[2];
sx q[2];
rz(1.5319038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84872765) q[1];
sx q[1];
rz(-1.25602) q[1];
sx q[1];
rz(-0.12401144) q[1];
rz(-pi) q[2];
rz(-2.9014663) q[3];
sx q[3];
rz(-1.4014981) q[3];
sx q[3];
rz(-1.5907767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-0.038643535) q[2];
rz(2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53428179) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(-1.779153) q[0];
rz(2.3249987) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(-1.1626676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1122738) q[0];
sx q[0];
rz(-0.11538878) q[0];
sx q[0];
rz(1.832604) q[0];
rz(-pi) q[1];
rz(-1.8936929) q[2];
sx q[2];
rz(-1.731551) q[2];
sx q[2];
rz(-1.9184743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3001544) q[1];
sx q[1];
rz(-1.9342124) q[1];
sx q[1];
rz(0.20519786) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2546222) q[3];
sx q[3];
rz(-0.63347048) q[3];
sx q[3];
rz(-2.5630643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5148619) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-0.14222063) q[2];
rz(-0.90406117) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(-0.18946762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(1.4676771) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(2.8335559) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7361569) q[0];
sx q[0];
rz(-1.6653367) q[0];
sx q[0];
rz(-3.0148538) q[0];
rz(2.1790444) q[2];
sx q[2];
rz(-1.0481917) q[2];
sx q[2];
rz(2.7163497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.088965485) q[1];
sx q[1];
rz(-1.7672156) q[1];
sx q[1];
rz(-0.88167015) q[1];
x q[2];
rz(2.8940053) q[3];
sx q[3];
rz(-2.7831804) q[3];
sx q[3];
rz(-2.4393516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3623111) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(-0.1299783) q[0];
rz(-0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-0.67108363) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6993461) q[0];
sx q[0];
rz(-1.6244349) q[0];
sx q[0];
rz(1.6058812) q[0];
x q[1];
rz(0.030521557) q[2];
sx q[2];
rz(-1.3866716) q[2];
sx q[2];
rz(-2.7437999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1783501) q[1];
sx q[1];
rz(-1.9415783) q[1];
sx q[1];
rz(0.11022719) q[1];
rz(1.5116793) q[3];
sx q[3];
rz(-0.67614188) q[3];
sx q[3];
rz(-2.6638871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(0.028586483) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(-1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639444) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(0.41241616) q[0];
rz(-1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.9746045) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0359356) q[0];
sx q[0];
rz(-1.783839) q[0];
sx q[0];
rz(-0.94353326) q[0];
rz(-2.2799759) q[2];
sx q[2];
rz(-1.3590727) q[2];
sx q[2];
rz(-0.35702969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3762714) q[1];
sx q[1];
rz(-1.5404535) q[1];
sx q[1];
rz(2.9529851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0409045) q[3];
sx q[3];
rz(-1.7003254) q[3];
sx q[3];
rz(-2.0199752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-0.18903014) q[2];
rz(-2.9947301) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(-1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1259595) q[0];
sx q[0];
rz(-1.3110315) q[0];
sx q[0];
rz(-2.2145859) q[0];
rz(1.3828297) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(1.6519201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.888436) q[0];
sx q[0];
rz(-2.2373767) q[0];
sx q[0];
rz(-2.1783834) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0909806) q[2];
sx q[2];
rz(-2.0423371) q[2];
sx q[2];
rz(-0.6363578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5786963) q[1];
sx q[1];
rz(-1.3778731) q[1];
sx q[1];
rz(1.3955411) q[1];
x q[2];
rz(2.6418266) q[3];
sx q[3];
rz(-1.0904113) q[3];
sx q[3];
rz(1.8936611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5513409) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(-1.4302953) q[2];
rz(2.5643505) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(0.28840315) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(-2.9945701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4924016) q[0];
sx q[0];
rz(-2.1426139) q[0];
sx q[0];
rz(2.6834821) q[0];
rz(-pi) q[1];
rz(-1.5006127) q[2];
sx q[2];
rz(-1.8101705) q[2];
sx q[2];
rz(1.0603051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22132561) q[1];
sx q[1];
rz(-2.1888071) q[1];
sx q[1];
rz(3.0926535) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9045194) q[3];
sx q[3];
rz(-0.58963886) q[3];
sx q[3];
rz(1.3414563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(0.75731164) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(-1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1159146) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(0.81746447) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(-0.13851891) q[2];
sx q[2];
rz(-0.45590966) q[2];
sx q[2];
rz(1.6675303) q[2];
rz(-1.0088624) q[3];
sx q[3];
rz(-1.4603793) q[3];
sx q[3];
rz(0.59752656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
