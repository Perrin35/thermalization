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
rz(1.1627816) q[0];
sx q[0];
rz(-0.94036189) q[0];
sx q[0];
rz(-2.9475687) q[0];
rz(2.5972875) q[1];
sx q[1];
rz(-1.5736009) q[1];
sx q[1];
rz(0.1967217) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5618043) q[0];
sx q[0];
rz(-1.4477948) q[0];
sx q[0];
rz(-0.99539006) q[0];
rz(1.5926378) q[2];
sx q[2];
rz(-1.7882203) q[2];
sx q[2];
rz(-0.79628758) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13290545) q[1];
sx q[1];
rz(-0.42294914) q[1];
sx q[1];
rz(-2.4407331) q[1];
rz(-pi) q[2];
rz(1.807735) q[3];
sx q[3];
rz(-0.58962017) q[3];
sx q[3];
rz(0.20649621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0804312) q[2];
sx q[2];
rz(-1.72074) q[2];
sx q[2];
rz(-3.1283992) q[2];
rz(-2.0773928) q[3];
sx q[3];
rz(-2.1049757) q[3];
sx q[3];
rz(-2.9707151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5203633) q[0];
sx q[0];
rz(-2.5907488) q[0];
sx q[0];
rz(0.4826104) q[0];
rz(-0.47166011) q[1];
sx q[1];
rz(-2.1444247) q[1];
sx q[1];
rz(0.68798033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575098) q[0];
sx q[0];
rz(-0.36195746) q[0];
sx q[0];
rz(0.56364001) q[0];
x q[1];
rz(-0.28765388) q[2];
sx q[2];
rz(-0.81865962) q[2];
sx q[2];
rz(-0.56861906) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1222306) q[1];
sx q[1];
rz(-1.5270342) q[1];
sx q[1];
rz(0.47143176) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.08556871) q[3];
sx q[3];
rz(-1.7734756) q[3];
sx q[3];
rz(0.50093716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1241577) q[2];
sx q[2];
rz(-2.0686801) q[2];
sx q[2];
rz(0.2249445) q[2];
rz(-1.0791091) q[3];
sx q[3];
rz(-2.2033117) q[3];
sx q[3];
rz(-2.6911531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9709622) q[0];
sx q[0];
rz(-0.58929515) q[0];
sx q[0];
rz(-1.2901837) q[0];
rz(1.2782485) q[1];
sx q[1];
rz(-1.9854913) q[1];
sx q[1];
rz(-1.8589171) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0541924) q[0];
sx q[0];
rz(-2.4303959) q[0];
sx q[0];
rz(-2.1367226) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.569246) q[2];
sx q[2];
rz(-2.9949014) q[2];
sx q[2];
rz(2.0990505) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45047347) q[1];
sx q[1];
rz(-1.5906003) q[1];
sx q[1];
rz(-3.1234689) q[1];
rz(0.77529737) q[3];
sx q[3];
rz(-0.92758815) q[3];
sx q[3];
rz(2.624529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4625385) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(-0.47492096) q[2];
rz(1.176152) q[3];
sx q[3];
rz(-0.85634309) q[3];
sx q[3];
rz(0.28158751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7684286) q[0];
sx q[0];
rz(-1.3685065) q[0];
sx q[0];
rz(-0.1678634) q[0];
rz(2.8236875) q[1];
sx q[1];
rz(-1.8524421) q[1];
sx q[1];
rz(1.0344523) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7307593) q[0];
sx q[0];
rz(-2.877088) q[0];
sx q[0];
rz(1.704731) q[0];
rz(-pi) q[1];
rz(0.45082966) q[2];
sx q[2];
rz(-0.73720142) q[2];
sx q[2];
rz(-0.99861713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9368091) q[1];
sx q[1];
rz(-2.2310871) q[1];
sx q[1];
rz(-2.5461063) q[1];
x q[2];
rz(-1.6717222) q[3];
sx q[3];
rz(-0.96138326) q[3];
sx q[3];
rz(0.17796365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0340524) q[2];
sx q[2];
rz(-2.5120021) q[2];
sx q[2];
rz(2.8752987) q[2];
rz(0.1693503) q[3];
sx q[3];
rz(-1.063187) q[3];
sx q[3];
rz(-2.9625986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87194815) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(-1.8484775) q[0];
rz(-0.070501892) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(0.806113) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9490813) q[0];
sx q[0];
rz(-1.3638138) q[0];
sx q[0];
rz(3.0144318) q[0];
rz(3.084153) q[2];
sx q[2];
rz(-1.4220793) q[2];
sx q[2];
rz(1.1401389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15547968) q[1];
sx q[1];
rz(-0.9782741) q[1];
sx q[1];
rz(-1.3377124) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53907303) q[3];
sx q[3];
rz(-1.381449) q[3];
sx q[3];
rz(2.2089778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21141323) q[2];
sx q[2];
rz(-1.6665919) q[2];
sx q[2];
rz(0.30964568) q[2];
rz(-2.6304701) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(-2.3206594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822134) q[0];
sx q[0];
rz(-2.5143304) q[0];
sx q[0];
rz(-0.43701592) q[0];
rz(-0.45477319) q[1];
sx q[1];
rz(-2.3450856) q[1];
sx q[1];
rz(2.2589267) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8681141) q[0];
sx q[0];
rz(-1.4820288) q[0];
sx q[0];
rz(-1.8411536) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4910554) q[2];
sx q[2];
rz(-0.3796595) q[2];
sx q[2];
rz(-2.9987631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8375802) q[1];
sx q[1];
rz(-1.092957) q[1];
sx q[1];
rz(-2.7663852) q[1];
rz(-pi) q[2];
rz(2.168247) q[3];
sx q[3];
rz(-2.4472467) q[3];
sx q[3];
rz(-2.2996817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5060045) q[2];
sx q[2];
rz(-1.4281861) q[2];
sx q[2];
rz(0.44090718) q[2];
rz(-1.0718369) q[3];
sx q[3];
rz(-0.25522885) q[3];
sx q[3];
rz(0.60060874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275682) q[0];
sx q[0];
rz(-1.9381645) q[0];
sx q[0];
rz(-1.2822275) q[0];
rz(2.0558689) q[1];
sx q[1];
rz(-1.6174569) q[1];
sx q[1];
rz(1.5132743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3886627) q[0];
sx q[0];
rz(-1.4528028) q[0];
sx q[0];
rz(0.8334882) q[0];
x q[1];
rz(2.7658913) q[2];
sx q[2];
rz(-0.73601572) q[2];
sx q[2];
rz(-1.0003566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6748057) q[1];
sx q[1];
rz(-1.9478096) q[1];
sx q[1];
rz(0.88175718) q[1];
rz(-pi) q[2];
rz(-0.75700883) q[3];
sx q[3];
rz(-1.8534213) q[3];
sx q[3];
rz(-0.69863897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47161237) q[2];
sx q[2];
rz(-2.4358304) q[2];
sx q[2];
rz(2.3213049) q[2];
rz(-2.7465076) q[3];
sx q[3];
rz(-0.79212752) q[3];
sx q[3];
rz(-2.5274966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.26246) q[0];
sx q[0];
rz(-1.6430055) q[0];
sx q[0];
rz(0.60687989) q[0];
rz(-2.795769) q[1];
sx q[1];
rz(-1.1864097) q[1];
sx q[1];
rz(-0.80723673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61582923) q[0];
sx q[0];
rz(-0.35789546) q[0];
sx q[0];
rz(-1.2919904) q[0];
rz(-pi) q[1];
rz(0.99844653) q[2];
sx q[2];
rz(-2.672451) q[2];
sx q[2];
rz(-2.1690705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2507627) q[1];
sx q[1];
rz(-1.8381665) q[1];
sx q[1];
rz(-1.3469338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1764507) q[3];
sx q[3];
rz(-1.0612086) q[3];
sx q[3];
rz(1.4450362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3240933) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(2.6452046) q[2];
rz(3.056622) q[3];
sx q[3];
rz(-0.43046633) q[3];
sx q[3];
rz(-2.5607204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1519415) q[0];
sx q[0];
rz(-1.713151) q[0];
sx q[0];
rz(-2.7570643) q[0];
rz(0.25228581) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(-1.5834454) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2315518) q[0];
sx q[0];
rz(-1.1052026) q[0];
sx q[0];
rz(1.5814073) q[0];
rz(-1.4531288) q[2];
sx q[2];
rz(-1.298438) q[2];
sx q[2];
rz(-0.54335574) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.059824316) q[1];
sx q[1];
rz(-1.5835973) q[1];
sx q[1];
rz(0.98155419) q[1];
rz(-2.2966398) q[3];
sx q[3];
rz(-2.264177) q[3];
sx q[3];
rz(1.4385731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(2.0015049) q[2];
rz(1.358486) q[3];
sx q[3];
rz(-1.2525109) q[3];
sx q[3];
rz(-0.31806773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43050218) q[0];
sx q[0];
rz(-1.4416114) q[0];
sx q[0];
rz(0.39518133) q[0];
rz(1.2218062) q[1];
sx q[1];
rz(-2.2876574) q[1];
sx q[1];
rz(-0.88776678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0096171776) q[0];
sx q[0];
rz(-1.8450301) q[0];
sx q[0];
rz(-2.4188309) q[0];
x q[1];
rz(-1.5639006) q[2];
sx q[2];
rz(-0.47499945) q[2];
sx q[2];
rz(1.4003225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6391082) q[1];
sx q[1];
rz(-0.35668761) q[1];
sx q[1];
rz(-1.6438032) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2492845) q[3];
sx q[3];
rz(-2.8148533) q[3];
sx q[3];
rz(0.48130166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6937574) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(0.57175222) q[2];
rz(1.6311496) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(-1.800764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0157264) q[0];
sx q[0];
rz(-1.3073574) q[0];
sx q[0];
rz(0.14412185) q[0];
rz(1.1852599) q[1];
sx q[1];
rz(-2.2846501) q[1];
sx q[1];
rz(1.2774998) q[1];
rz(-1.324426) q[2];
sx q[2];
rz(-1.194321) q[2];
sx q[2];
rz(1.3512667) q[2];
rz(-0.64438728) q[3];
sx q[3];
rz(-1.3105262) q[3];
sx q[3];
rz(-2.5986828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
