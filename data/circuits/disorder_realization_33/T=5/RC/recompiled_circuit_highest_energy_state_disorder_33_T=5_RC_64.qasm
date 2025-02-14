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
rz(-2.9448709) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08840522) q[0];
sx q[0];
rz(-1.0002828) q[0];
sx q[0];
rz(2.9952917) q[0];
rz(-pi) q[1];
rz(3.0430531) q[2];
sx q[2];
rz(-0.218501) q[2];
sx q[2];
rz(0.69536415) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87952033) q[1];
sx q[1];
rz(-1.2517057) q[1];
sx q[1];
rz(1.8533005) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1473925) q[3];
sx q[3];
rz(-1.4399035) q[3];
sx q[3];
rz(1.5623705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0611614) q[2];
sx q[2];
rz(-1.4208527) q[2];
sx q[2];
rz(-3.1283992) q[2];
rz(2.0773928) q[3];
sx q[3];
rz(-1.0366169) q[3];
sx q[3];
rz(0.17087759) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5203633) q[0];
sx q[0];
rz(-0.55084387) q[0];
sx q[0];
rz(2.6589822) q[0];
rz(0.47166011) q[1];
sx q[1];
rz(-2.1444247) q[1];
sx q[1];
rz(-0.68798033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9887138) q[0];
sx q[0];
rz(-1.3804624) q[0];
sx q[0];
rz(0.3097663) q[0];
x q[1];
rz(1.8652164) q[2];
sx q[2];
rz(-2.3465119) q[2];
sx q[2];
rz(0.15995041) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5045484) q[1];
sx q[1];
rz(-0.47330644) q[1];
sx q[1];
rz(0.096122336) q[1];
rz(-pi) q[2];
rz(-1.9649195) q[3];
sx q[3];
rz(-0.21977327) q[3];
sx q[3];
rz(0.098107396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1241577) q[2];
sx q[2];
rz(-1.0729125) q[2];
sx q[2];
rz(-0.2249445) q[2];
rz(2.0624835) q[3];
sx q[3];
rz(-2.2033117) q[3];
sx q[3];
rz(-2.6911531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17063046) q[0];
sx q[0];
rz(-0.58929515) q[0];
sx q[0];
rz(-1.851409) q[0];
rz(1.8633441) q[1];
sx q[1];
rz(-1.1561013) q[1];
sx q[1];
rz(-1.8589171) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7851832) q[0];
sx q[0];
rz(-0.98726596) q[0];
sx q[0];
rz(-2.7088091) q[0];
rz(3.01802) q[2];
sx q[2];
rz(-1.6500435) q[2];
sx q[2];
rz(-1.0956956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29073971) q[1];
sx q[1];
rz(-0.026844414) q[1];
sx q[1];
rz(2.31183) q[1];
rz(2.3219789) q[3];
sx q[3];
rz(-0.96246877) q[3];
sx q[3];
rz(1.5386594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4625385) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(2.6666717) q[2];
rz(-1.176152) q[3];
sx q[3];
rz(-0.85634309) q[3];
sx q[3];
rz(-0.28158751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7684286) q[0];
sx q[0];
rz(-1.3685065) q[0];
sx q[0];
rz(0.1678634) q[0];
rz(2.8236875) q[1];
sx q[1];
rz(-1.8524421) q[1];
sx q[1];
rz(1.0344523) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7307593) q[0];
sx q[0];
rz(-2.877088) q[0];
sx q[0];
rz(-1.704731) q[0];
rz(-pi) q[1];
rz(-1.1940766) q[2];
sx q[2];
rz(-2.2206306) q[2];
sx q[2];
rz(-1.5776933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7718879) q[1];
sx q[1];
rz(-2.2834816) q[1];
sx q[1];
rz(0.94526498) q[1];
rz(-pi) q[2];
rz(2.529783) q[3];
sx q[3];
rz(-1.6535078) q[3];
sx q[3];
rz(-1.3349319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1075403) q[2];
sx q[2];
rz(-2.5120021) q[2];
sx q[2];
rz(-2.8752987) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87194815) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(1.8484775) q[0];
rz(0.070501892) q[1];
sx q[1];
rz(-2.2673456) q[1];
sx q[1];
rz(0.806113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4045532) q[0];
sx q[0];
rz(-1.6952288) q[0];
sx q[0];
rz(-1.3621773) q[0];
x q[1];
rz(-1.9367123) q[2];
sx q[2];
rz(-0.15934773) q[2];
sx q[2];
rz(1.6312576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8580839) q[1];
sx q[1];
rz(-1.3779989) q[1];
sx q[1];
rz(-0.60528614) q[1];
rz(-2.784291) q[3];
sx q[3];
rz(-0.56824486) q[3];
sx q[3];
rz(2.8082591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21141323) q[2];
sx q[2];
rz(-1.4750007) q[2];
sx q[2];
rz(0.30964568) q[2];
rz(2.6304701) q[3];
sx q[3];
rz(-0.57356727) q[3];
sx q[3];
rz(-2.3206594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1593793) q[0];
sx q[0];
rz(-2.5143304) q[0];
sx q[0];
rz(-0.43701592) q[0];
rz(2.6868195) q[1];
sx q[1];
rz(-0.79650703) q[1];
sx q[1];
rz(0.88266596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32188181) q[0];
sx q[0];
rz(-1.30153) q[0];
sx q[0];
rz(-0.092094941) q[0];
rz(-1.9493616) q[2];
sx q[2];
rz(-1.6003216) q[2];
sx q[2];
rz(1.6395417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55403472) q[1];
sx q[1];
rz(-1.9022501) q[1];
sx q[1];
rz(2.0786839) q[1];
rz(2.168247) q[3];
sx q[3];
rz(-2.4472467) q[3];
sx q[3];
rz(-2.2996817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63558811) q[2];
sx q[2];
rz(-1.7134066) q[2];
sx q[2];
rz(2.7006855) q[2];
rz(2.0697557) q[3];
sx q[3];
rz(-2.8863638) q[3];
sx q[3];
rz(2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(3.1275682) q[0];
sx q[0];
rz(-1.9381645) q[0];
sx q[0];
rz(1.8593651) q[0];
rz(1.0857238) q[1];
sx q[1];
rz(-1.5241357) q[1];
sx q[1];
rz(1.5132743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075628932) q[0];
sx q[0];
rz(-0.83978486) q[0];
sx q[0];
rz(-2.9828067) q[0];
rz(-1.8916745) q[2];
sx q[2];
rz(-0.89628637) q[2];
sx q[2];
rz(0.51128507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6172863) q[1];
sx q[1];
rz(-2.3712161) q[1];
sx q[1];
rz(1.0138033) q[1];
rz(0.75700883) q[3];
sx q[3];
rz(-1.2881713) q[3];
sx q[3];
rz(2.4429537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6699803) q[2];
sx q[2];
rz(-0.70576224) q[2];
sx q[2];
rz(-0.82028779) q[2];
rz(-2.7465076) q[3];
sx q[3];
rz(-2.3494651) q[3];
sx q[3];
rz(-0.61409605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.26246) q[0];
sx q[0];
rz(-1.4985871) q[0];
sx q[0];
rz(0.60687989) q[0];
rz(-2.795769) q[1];
sx q[1];
rz(-1.9551829) q[1];
sx q[1];
rz(0.80723673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61582923) q[0];
sx q[0];
rz(-0.35789546) q[0];
sx q[0];
rz(-1.8496022) q[0];
rz(0.99844653) q[2];
sx q[2];
rz(-0.46914161) q[2];
sx q[2];
rz(-0.97252211) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25995421) q[1];
sx q[1];
rz(-1.3550161) q[1];
sx q[1];
rz(0.27387932) q[1];
rz(-pi) q[2];
rz(0.96514197) q[3];
sx q[3];
rz(-1.0612086) q[3];
sx q[3];
rz(-1.4450362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3240933) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(0.49638805) q[2];
rz(-3.056622) q[3];
sx q[3];
rz(-2.7111263) q[3];
sx q[3];
rz(-2.5607204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519415) q[0];
sx q[0];
rz(-1.4284416) q[0];
sx q[0];
rz(-2.7570643) q[0];
rz(0.25228581) q[1];
sx q[1];
rz(-1.7583881) q[1];
sx q[1];
rz(1.5834454) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91004086) q[0];
sx q[0];
rz(-2.0363901) q[0];
sx q[0];
rz(1.5814073) q[0];
rz(-2.7436888) q[2];
sx q[2];
rz(-0.29610866) q[2];
sx q[2];
rz(0.95740151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0817683) q[1];
sx q[1];
rz(-1.5835973) q[1];
sx q[1];
rz(0.98155419) q[1];
rz(-pi) q[2];
rz(2.4676314) q[3];
sx q[3];
rz(-2.1836851) q[3];
sx q[3];
rz(-0.4919258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(-1.1400878) q[2];
rz(-1.358486) q[3];
sx q[3];
rz(-1.8890817) q[3];
sx q[3];
rz(2.8235249) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7110905) q[0];
sx q[0];
rz(-1.4416114) q[0];
sx q[0];
rz(2.7464113) q[0];
rz(1.9197865) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(-0.88776678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3459612) q[0];
sx q[0];
rz(-2.2610616) q[0];
sx q[0];
rz(1.9296586) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0035462499) q[2];
sx q[2];
rz(-2.0457836) q[2];
sx q[2];
rz(1.3925683) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.58037591) q[1];
sx q[1];
rz(-1.2151011) q[1];
sx q[1];
rz(0.027173398) q[1];
x q[2];
rz(-1.8923081) q[3];
sx q[3];
rz(-2.8148533) q[3];
sx q[3];
rz(-2.660291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44783529) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(-2.5698404) q[2];
rz(1.6311496) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(-1.800764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1258662) q[0];
sx q[0];
rz(-1.8342352) q[0];
sx q[0];
rz(-2.9974708) q[0];
rz(-1.9563328) q[1];
sx q[1];
rz(-2.2846501) q[1];
sx q[1];
rz(1.2774998) q[1];
rz(-1.8171667) q[2];
sx q[2];
rz(-1.9472716) q[2];
sx q[2];
rz(-1.7903259) q[2];
rz(0.41729144) q[3];
sx q[3];
rz(-2.4536316) q[3];
sx q[3];
rz(-1.3578547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
