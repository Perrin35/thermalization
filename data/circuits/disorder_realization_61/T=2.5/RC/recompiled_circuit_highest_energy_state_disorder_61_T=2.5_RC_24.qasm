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
rz(-2.6391368) q[0];
sx q[0];
rz(-1.1075736) q[0];
sx q[0];
rz(1.7662319) q[0];
rz(-0.31887588) q[1];
sx q[1];
rz(-1.640929) q[1];
sx q[1];
rz(0.73135102) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1494182) q[0];
sx q[0];
rz(-0.96620027) q[0];
sx q[0];
rz(-0.3913619) q[0];
rz(-pi) q[1];
rz(0.66197239) q[2];
sx q[2];
rz(-1.5740924) q[2];
sx q[2];
rz(0.73288435) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9125376) q[1];
sx q[1];
rz(-2.3542535) q[1];
sx q[1];
rz(0.064219193) q[1];
rz(-0.5238138) q[3];
sx q[3];
rz(-2.9546545) q[3];
sx q[3];
rz(-0.65526774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2607164) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-2.1991849) q[2];
rz(2.2830394) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(-0.71949351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.079916) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(3.0829561) q[0];
rz(-0.05642852) q[1];
sx q[1];
rz(-1.3899048) q[1];
sx q[1];
rz(2.0243534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4898713) q[0];
sx q[0];
rz(-2.3014223) q[0];
sx q[0];
rz(-1.0307167) q[0];
rz(-1.9675206) q[2];
sx q[2];
rz(-2.8182321) q[2];
sx q[2];
rz(0.67229313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8630115) q[1];
sx q[1];
rz(-1.8092691) q[1];
sx q[1];
rz(-3.0568074) q[1];
x q[2];
rz(0.47567077) q[3];
sx q[3];
rz(-1.4632309) q[3];
sx q[3];
rz(-1.2616664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3650018) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(-2.7351725) q[2];
rz(2.8255919) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(0.78280226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2773748) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(1.8055441) q[0];
rz(0.27851963) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(0.96885243) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3144041) q[0];
sx q[0];
rz(-2.1705856) q[0];
sx q[0];
rz(-1.1161854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79040852) q[2];
sx q[2];
rz(-1.1539173) q[2];
sx q[2];
rz(-1.1000772) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67476701) q[1];
sx q[1];
rz(-1.6129426) q[1];
sx q[1];
rz(-1.274363) q[1];
rz(-2.4163867) q[3];
sx q[3];
rz(-1.375633) q[3];
sx q[3];
rz(-2.8233087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1442673) q[2];
sx q[2];
rz(-2.6582025) q[2];
sx q[2];
rz(0.17511314) q[2];
rz(-1.3052321) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(1.2828627) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99821943) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(1.0293707) q[0];
rz(-2.3221305) q[1];
sx q[1];
rz(-1.1623323) q[1];
sx q[1];
rz(0.79684657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5982847) q[0];
sx q[0];
rz(-0.48439041) q[0];
sx q[0];
rz(0.422307) q[0];
rz(1.5605075) q[2];
sx q[2];
rz(-0.52797374) q[2];
sx q[2];
rz(-2.1473403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0681035) q[1];
sx q[1];
rz(-1.3655546) q[1];
sx q[1];
rz(-1.2935335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53189532) q[3];
sx q[3];
rz(-0.73557702) q[3];
sx q[3];
rz(3.0656505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0113819) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(-1.9127362) q[2];
rz(-2.7025488) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6586886) q[0];
sx q[0];
rz(-1.5331601) q[0];
sx q[0];
rz(-2.3332692) q[0];
rz(1.7841548) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(2.435991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7019136) q[0];
sx q[0];
rz(-1.8712303) q[0];
sx q[0];
rz(-1.6070532) q[0];
rz(0.23106261) q[2];
sx q[2];
rz(-0.88132826) q[2];
sx q[2];
rz(1.5514411) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28532923) q[1];
sx q[1];
rz(-1.3176022) q[1];
sx q[1];
rz(-3.1388389) q[1];
rz(-pi) q[2];
rz(0.79885317) q[3];
sx q[3];
rz(-2.0435413) q[3];
sx q[3];
rz(0.96351876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6610873) q[2];
sx q[2];
rz(-1.2837807) q[2];
sx q[2];
rz(2.535848) q[2];
rz(0.12039603) q[3];
sx q[3];
rz(-1.2984637) q[3];
sx q[3];
rz(-2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8814988) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(3.1392198) q[0];
rz(-0.35898769) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(0.40828362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7255032) q[0];
sx q[0];
rz(-1.444284) q[0];
sx q[0];
rz(-3.1055272) q[0];
rz(-pi) q[1];
rz(-1.9469366) q[2];
sx q[2];
rz(-2.500211) q[2];
sx q[2];
rz(-2.1958786) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6039415) q[1];
sx q[1];
rz(-1.8007601) q[1];
sx q[1];
rz(1.5020788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8114641) q[3];
sx q[3];
rz(-2.0272019) q[3];
sx q[3];
rz(1.8965669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6458873) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(0.1296981) q[2];
rz(-2.7221223) q[3];
sx q[3];
rz(-1.2129815) q[3];
sx q[3];
rz(-2.3278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3892589) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(2.4580521) q[0];
rz(2.2031802) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(1.2260812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968564) q[0];
sx q[0];
rz(-0.71030197) q[0];
sx q[0];
rz(-2.5596749) q[0];
rz(-0.60003321) q[2];
sx q[2];
rz(-1.4156315) q[2];
sx q[2];
rz(-0.70916353) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9270806) q[1];
sx q[1];
rz(-2.0190372) q[1];
sx q[1];
rz(-2.8238676) q[1];
rz(-2.5874373) q[3];
sx q[3];
rz(-0.99923493) q[3];
sx q[3];
rz(3.1031648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.73416) q[2];
sx q[2];
rz(-1.4340883) q[2];
sx q[2];
rz(2.4343991) q[2];
rz(1.4736942) q[3];
sx q[3];
rz(-0.67631045) q[3];
sx q[3];
rz(1.428712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4555175) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(2.2177875) q[0];
rz(2.0763981) q[1];
sx q[1];
rz(-1.9520452) q[1];
sx q[1];
rz(1.1346029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77662879) q[0];
sx q[0];
rz(-2.7641649) q[0];
sx q[0];
rz(-1.6799742) q[0];
rz(-pi) q[1];
rz(2.258111) q[2];
sx q[2];
rz(-2.0728014) q[2];
sx q[2];
rz(2.4740681) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90701404) q[1];
sx q[1];
rz(-1.624214) q[1];
sx q[1];
rz(1.0397041) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60875684) q[3];
sx q[3];
rz(-1.8894102) q[3];
sx q[3];
rz(-2.1515338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.137546) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(-2.498632) q[2];
rz(0.17063394) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(1.0954866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74649015) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(-2.3353031) q[0];
rz(0.0099446615) q[1];
sx q[1];
rz(-0.6178304) q[1];
sx q[1];
rz(0.26201216) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952103) q[0];
sx q[0];
rz(-1.4052009) q[0];
sx q[0];
rz(-1.0402334) q[0];
rz(-0.21199439) q[2];
sx q[2];
rz(-2.4851481) q[2];
sx q[2];
rz(0.5613297) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9741874) q[1];
sx q[1];
rz(-1.3648808) q[1];
sx q[1];
rz(-2.0108825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32829653) q[3];
sx q[3];
rz(-2.0496164) q[3];
sx q[3];
rz(-1.4737859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3917196) q[2];
sx q[2];
rz(-0.55730692) q[2];
sx q[2];
rz(-0.14287512) q[2];
rz(2.2376132) q[3];
sx q[3];
rz(-1.6429792) q[3];
sx q[3];
rz(-2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751462) q[0];
sx q[0];
rz(-1.9117993) q[0];
sx q[0];
rz(0.65650702) q[0];
rz(-0.87908602) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(0.62969977) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59081599) q[0];
sx q[0];
rz(-1.970298) q[0];
sx q[0];
rz(2.4062064) q[0];
rz(-pi) q[1];
rz(1.1348502) q[2];
sx q[2];
rz(-2.5647047) q[2];
sx q[2];
rz(-1.0283805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7460664) q[1];
sx q[1];
rz(-1.5256572) q[1];
sx q[1];
rz(-3.1137115) q[1];
rz(-pi) q[2];
rz(-2.8318226) q[3];
sx q[3];
rz(-2.0783721) q[3];
sx q[3];
rz(0.47981706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62275824) q[2];
sx q[2];
rz(-2.3164985) q[2];
sx q[2];
rz(2.4614914) q[2];
rz(2.425219) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(2.0228001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730561) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(2.9987891) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(1.3423781) q[2];
sx q[2];
rz(-2.8932718) q[2];
sx q[2];
rz(0.415033) q[2];
rz(0.96137267) q[3];
sx q[3];
rz(-1.537804) q[3];
sx q[3];
rz(0.15450879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
