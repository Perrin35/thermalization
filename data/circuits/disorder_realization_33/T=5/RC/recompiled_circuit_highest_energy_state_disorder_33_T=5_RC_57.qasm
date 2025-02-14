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
rz(5.3428234) q[0];
sx q[0];
rz(12.760395) q[0];
rz(2.5972875) q[1];
sx q[1];
rz(-1.5736009) q[1];
sx q[1];
rz(0.1967217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5618043) q[0];
sx q[0];
rz(-1.6937979) q[0];
sx q[0];
rz(-2.1462026) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21747422) q[2];
sx q[2];
rz(-1.5921235) q[2];
sx q[2];
rz(0.77922098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78208343) q[1];
sx q[1];
rz(-1.302914) q[1];
sx q[1];
rz(2.8102576) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9858304) q[3];
sx q[3];
rz(-0.99975306) q[3];
sx q[3];
rz(-3.0653428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0804312) q[2];
sx q[2];
rz(-1.72074) q[2];
sx q[2];
rz(0.013193456) q[2];
rz(-1.0641998) q[3];
sx q[3];
rz(-1.0366169) q[3];
sx q[3];
rz(0.17087759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5203633) q[0];
sx q[0];
rz(-0.55084387) q[0];
sx q[0];
rz(0.4826104) q[0];
rz(-0.47166011) q[1];
sx q[1];
rz(-2.1444247) q[1];
sx q[1];
rz(0.68798033) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9887138) q[0];
sx q[0];
rz(-1.3804624) q[0];
sx q[0];
rz(0.3097663) q[0];
rz(0.28765388) q[2];
sx q[2];
rz(-0.81865962) q[2];
sx q[2];
rz(0.56861906) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6370442) q[1];
sx q[1];
rz(-0.47330644) q[1];
sx q[1];
rz(-3.0454703) q[1];
rz(-1.1766731) q[3];
sx q[3];
rz(-0.21977327) q[3];
sx q[3];
rz(3.0434853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.017435) q[2];
sx q[2];
rz(-1.0729125) q[2];
sx q[2];
rz(0.2249445) q[2];
rz(-2.0624835) q[3];
sx q[3];
rz(-0.93828097) q[3];
sx q[3];
rz(0.45043954) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17063046) q[0];
sx q[0];
rz(-2.5522975) q[0];
sx q[0];
rz(1.851409) q[0];
rz(1.2782485) q[1];
sx q[1];
rz(-1.9854913) q[1];
sx q[1];
rz(-1.8589171) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034866355) q[0];
sx q[0];
rz(-1.9283656) q[0];
sx q[0];
rz(-2.1995972) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.01802) q[2];
sx q[2];
rz(-1.4915492) q[2];
sx q[2];
rz(-1.0956956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45047347) q[1];
sx q[1];
rz(-1.5509924) q[1];
sx q[1];
rz(0.018123716) q[1];
rz(-pi) q[2];
x q[2];
rz(2.380326) q[3];
sx q[3];
rz(-2.1652615) q[3];
sx q[3];
rz(0.52235421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4625385) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(2.6666717) q[2];
rz(-1.9654407) q[3];
sx q[3];
rz(-0.85634309) q[3];
sx q[3];
rz(0.28158751) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7684286) q[0];
sx q[0];
rz(-1.7730862) q[0];
sx q[0];
rz(0.1678634) q[0];
rz(-2.8236875) q[1];
sx q[1];
rz(-1.8524421) q[1];
sx q[1];
rz(2.1071404) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7307593) q[0];
sx q[0];
rz(-2.877088) q[0];
sx q[0];
rz(1.704731) q[0];
x q[1];
rz(-1.9475161) q[2];
sx q[2];
rz(-0.92096201) q[2];
sx q[2];
rz(1.5638994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9368091) q[1];
sx q[1];
rz(-0.91050557) q[1];
sx q[1];
rz(2.5461063) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6717222) q[3];
sx q[3];
rz(-0.96138326) q[3];
sx q[3];
rz(-0.17796365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1075403) q[2];
sx q[2];
rz(-2.5120021) q[2];
sx q[2];
rz(-2.8752987) q[2];
rz(-0.1693503) q[3];
sx q[3];
rz(-1.063187) q[3];
sx q[3];
rz(2.9625986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87194815) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(-1.2931152) q[0];
rz(3.0710908) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(-2.3354796) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4045532) q[0];
sx q[0];
rz(-1.4463639) q[0];
sx q[0];
rz(-1.7794153) q[0];
rz(-pi) q[1];
rz(1.7197554) q[2];
sx q[2];
rz(-1.5139914) q[2];
sx q[2];
rz(2.7194552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2835088) q[1];
sx q[1];
rz(-1.7635937) q[1];
sx q[1];
rz(-2.5363065) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7905037) q[3];
sx q[3];
rz(-1.0423805) q[3];
sx q[3];
rz(2.3912969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21141323) q[2];
sx q[2];
rz(-1.4750007) q[2];
sx q[2];
rz(0.30964568) q[2];
rz(-0.51112255) q[3];
sx q[3];
rz(-0.57356727) q[3];
sx q[3];
rz(-2.3206594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.9822134) q[0];
sx q[0];
rz(-2.5143304) q[0];
sx q[0];
rz(0.43701592) q[0];
rz(0.45477319) q[1];
sx q[1];
rz(-2.3450856) q[1];
sx q[1];
rz(0.88266596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2734786) q[0];
sx q[0];
rz(-1.4820288) q[0];
sx q[0];
rz(1.300439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9493616) q[2];
sx q[2];
rz(-1.6003216) q[2];
sx q[2];
rz(1.6395417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55403472) q[1];
sx q[1];
rz(-1.9022501) q[1];
sx q[1];
rz(-2.0786839) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96787937) q[3];
sx q[3];
rz(-1.2025739) q[3];
sx q[3];
rz(0.24711025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63558811) q[2];
sx q[2];
rz(-1.7134066) q[2];
sx q[2];
rz(2.7006855) q[2];
rz(1.0718369) q[3];
sx q[3];
rz(-2.8863638) q[3];
sx q[3];
rz(-2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275682) q[0];
sx q[0];
rz(-1.9381645) q[0];
sx q[0];
rz(-1.8593651) q[0];
rz(2.0558689) q[1];
sx q[1];
rz(-1.5241357) q[1];
sx q[1];
rz(-1.5132743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31103872) q[0];
sx q[0];
rz(-0.74493248) q[0];
sx q[0];
rz(-5*pi/9) q[0];
x q[1];
rz(1.2499181) q[2];
sx q[2];
rz(-0.89628637) q[2];
sx q[2];
rz(-2.6303076) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3320184) q[1];
sx q[1];
rz(-0.93832131) q[1];
sx q[1];
rz(2.6676086) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75700883) q[3];
sx q[3];
rz(-1.2881713) q[3];
sx q[3];
rz(-0.69863897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6699803) q[2];
sx q[2];
rz(-2.4358304) q[2];
sx q[2];
rz(-0.82028779) q[2];
rz(0.39508501) q[3];
sx q[3];
rz(-2.3494651) q[3];
sx q[3];
rz(-0.61409605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.26246) q[0];
sx q[0];
rz(-1.4985871) q[0];
sx q[0];
rz(-0.60687989) q[0];
rz(-2.795769) q[1];
sx q[1];
rz(-1.1864097) q[1];
sx q[1];
rz(-0.80723673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2169289) q[0];
sx q[0];
rz(-1.6673526) q[0];
sx q[0];
rz(1.9159622) q[0];
x q[1];
rz(-2.1431461) q[2];
sx q[2];
rz(-2.672451) q[2];
sx q[2];
rz(0.97252211) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25995421) q[1];
sx q[1];
rz(-1.3550161) q[1];
sx q[1];
rz(-0.27387932) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96514197) q[3];
sx q[3];
rz(-2.0803841) q[3];
sx q[3];
rz(-1.4450362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81749934) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(-2.6452046) q[2];
rz(0.084970623) q[3];
sx q[3];
rz(-0.43046633) q[3];
sx q[3];
rz(2.5607204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9896511) q[0];
sx q[0];
rz(-1.4284416) q[0];
sx q[0];
rz(-2.7570643) q[0];
rz(-2.8893068) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(1.5581473) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88640942) q[0];
sx q[0];
rz(-2.6758868) q[0];
sx q[0];
rz(0.021115594) q[0];
rz(-pi) q[1];
rz(-2.7436888) q[2];
sx q[2];
rz(-0.29610866) q[2];
sx q[2];
rz(0.95740151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0817683) q[1];
sx q[1];
rz(-1.5579954) q[1];
sx q[1];
rz(0.98155419) q[1];
rz(-0.6739613) q[3];
sx q[3];
rz(-0.95790759) q[3];
sx q[3];
rz(-2.6496668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3462476) q[2];
sx q[2];
rz(-1.8833505) q[2];
sx q[2];
rz(-1.1400878) q[2];
rz(1.358486) q[3];
sx q[3];
rz(-1.2525109) q[3];
sx q[3];
rz(-0.31806773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110905) q[0];
sx q[0];
rz(-1.6999812) q[0];
sx q[0];
rz(-0.39518133) q[0];
rz(-1.9197865) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(0.88776678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0096171776) q[0];
sx q[0];
rz(-1.8450301) q[0];
sx q[0];
rz(-2.4188309) q[0];
x q[1];
rz(0.0035462499) q[2];
sx q[2];
rz(-1.0958091) q[2];
sx q[2];
rz(1.7490243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99988538) q[1];
sx q[1];
rz(-1.5453242) q[1];
sx q[1];
rz(1.2149806) q[1];
rz(-pi) q[2];
rz(1.2597144) q[3];
sx q[3];
rz(-1.4691989) q[3];
sx q[3];
rz(1.3950619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.44783529) q[2];
sx q[2];
rz(-1.8303904) q[2];
sx q[2];
rz(-2.5698404) q[2];
rz(-1.6311496) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(1.800764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0157264) q[0];
sx q[0];
rz(-1.3073574) q[0];
sx q[0];
rz(0.14412185) q[0];
rz(-1.9563328) q[1];
sx q[1];
rz(-2.2846501) q[1];
sx q[1];
rz(1.2774998) q[1];
rz(0.55276362) q[2];
sx q[2];
rz(-2.6949099) q[2];
sx q[2];
rz(0.75134122) q[2];
rz(-2.4972054) q[3];
sx q[3];
rz(-1.8310665) q[3];
sx q[3];
rz(0.54290988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
