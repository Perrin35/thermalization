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
rz(-1.9788111) q[0];
sx q[0];
rz(-2.2012308) q[0];
sx q[0];
rz(-0.19402394) q[0];
rz(2.5972875) q[1];
sx q[1];
rz(-1.5736009) q[1];
sx q[1];
rz(0.1967217) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9636363) q[0];
sx q[0];
rz(-0.58695176) q[0];
sx q[0];
rz(1.3474083) q[0];
rz(-pi) q[1];
rz(0.098539515) q[2];
sx q[2];
rz(-0.218501) q[2];
sx q[2];
rz(-0.69536415) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3595092) q[1];
sx q[1];
rz(-1.302914) q[1];
sx q[1];
rz(0.33133502) q[1];
rz(-0.99420011) q[3];
sx q[3];
rz(-1.4399035) q[3];
sx q[3];
rz(-1.5623705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0611614) q[2];
sx q[2];
rz(-1.72074) q[2];
sx q[2];
rz(-3.1283992) q[2];
rz(-2.0773928) q[3];
sx q[3];
rz(-1.0366169) q[3];
sx q[3];
rz(-0.17087759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-0.99716798) q[1];
sx q[1];
rz(0.68798033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575098) q[0];
sx q[0];
rz(-0.36195746) q[0];
sx q[0];
rz(2.5779526) q[0];
rz(2.8539388) q[2];
sx q[2];
rz(-2.322933) q[2];
sx q[2];
rz(-2.5729736) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6124561) q[1];
sx q[1];
rz(-2.0417401) q[1];
sx q[1];
rz(-1.5216842) q[1];
rz(-pi) q[2];
rz(-1.1766731) q[3];
sx q[3];
rz(-0.21977327) q[3];
sx q[3];
rz(-0.098107396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1241577) q[2];
sx q[2];
rz(-2.0686801) q[2];
sx q[2];
rz(2.9166481) q[2];
rz(1.0791091) q[3];
sx q[3];
rz(-0.93828097) q[3];
sx q[3];
rz(-2.6911531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17063046) q[0];
sx q[0];
rz(-2.5522975) q[0];
sx q[0];
rz(-1.2901837) q[0];
rz(-1.8633441) q[1];
sx q[1];
rz(-1.1561013) q[1];
sx q[1];
rz(1.8589171) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0541924) q[0];
sx q[0];
rz(-2.4303959) q[0];
sx q[0];
rz(2.1367226) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.01802) q[2];
sx q[2];
rz(-1.6500435) q[2];
sx q[2];
rz(-2.045897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29073971) q[1];
sx q[1];
rz(-0.026844414) q[1];
sx q[1];
rz(0.82976262) q[1];
rz(-pi) q[2];
rz(-0.81961378) q[3];
sx q[3];
rz(-2.1791239) q[3];
sx q[3];
rz(1.6029333) q[3];
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
rz(-2.2852496) q[3];
sx q[3];
rz(2.8600051) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3731641) q[0];
sx q[0];
rz(-1.7730862) q[0];
sx q[0];
rz(2.9737293) q[0];
rz(-0.31790512) q[1];
sx q[1];
rz(-1.2891506) q[1];
sx q[1];
rz(-1.0344523) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2721336) q[0];
sx q[0];
rz(-1.3087166) q[0];
sx q[0];
rz(3.1054405) q[0];
x q[1];
rz(-2.690763) q[2];
sx q[2];
rz(-0.73720142) q[2];
sx q[2];
rz(-0.99861713) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2047836) q[1];
sx q[1];
rz(-2.2310871) q[1];
sx q[1];
rz(-2.5461063) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61180964) q[3];
sx q[3];
rz(-1.4880848) q[3];
sx q[3];
rz(1.8066607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0340524) q[2];
sx q[2];
rz(-2.5120021) q[2];
sx q[2];
rz(-0.26629392) q[2];
rz(0.1693503) q[3];
sx q[3];
rz(-1.063187) q[3];
sx q[3];
rz(0.17899409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2696445) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(-1.2931152) q[0];
rz(0.070501892) q[1];
sx q[1];
rz(-2.2673456) q[1];
sx q[1];
rz(-2.3354796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4045532) q[0];
sx q[0];
rz(-1.4463639) q[0];
sx q[0];
rz(1.3621773) q[0];
x q[1];
rz(0.057439645) q[2];
sx q[2];
rz(-1.7195133) q[2];
sx q[2];
rz(1.1401389) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55743417) q[1];
sx q[1];
rz(-0.63156742) q[1];
sx q[1];
rz(0.33051349) q[1];
x q[2];
rz(0.53907303) q[3];
sx q[3];
rz(-1.381449) q[3];
sx q[3];
rz(-2.2089778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21141323) q[2];
sx q[2];
rz(-1.4750007) q[2];
sx q[2];
rz(0.30964568) q[2];
rz(2.6304701) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(-0.82093325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9822134) q[0];
sx q[0];
rz(-2.5143304) q[0];
sx q[0];
rz(0.43701592) q[0];
rz(-0.45477319) q[1];
sx q[1];
rz(-2.3450856) q[1];
sx q[1];
rz(-0.88266596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012271492) q[0];
sx q[0];
rz(-0.28421775) q[0];
sx q[0];
rz(-1.2491262) q[0];
rz(-1.9493616) q[2];
sx q[2];
rz(-1.541271) q[2];
sx q[2];
rz(-1.6395417) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5957634) q[1];
sx q[1];
rz(-2.543151) q[1];
sx q[1];
rz(-0.95495895) q[1];
rz(0.43805505) q[3];
sx q[3];
rz(-2.1282624) q[3];
sx q[3];
rz(-1.5665975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.63558811) q[2];
sx q[2];
rz(-1.4281861) q[2];
sx q[2];
rz(0.44090718) q[2];
rz(1.0718369) q[3];
sx q[3];
rz(-2.8863638) q[3];
sx q[3];
rz(-2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275682) q[0];
sx q[0];
rz(-1.2034282) q[0];
sx q[0];
rz(1.8593651) q[0];
rz(1.0857238) q[1];
sx q[1];
rz(-1.5241357) q[1];
sx q[1];
rz(1.5132743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3886627) q[0];
sx q[0];
rz(-1.4528028) q[0];
sx q[0];
rz(0.8334882) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4413928) q[2];
sx q[2];
rz(-1.8196897) q[2];
sx q[2];
rz(-2.2867416) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6172863) q[1];
sx q[1];
rz(-2.3712161) q[1];
sx q[1];
rz(-2.1277894) q[1];
x q[2];
rz(-1.9508771) q[3];
sx q[3];
rz(-0.8506368) q[3];
sx q[3];
rz(-2.5270568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6699803) q[2];
sx q[2];
rz(-0.70576224) q[2];
sx q[2];
rz(-0.82028779) q[2];
rz(-2.7465076) q[3];
sx q[3];
rz(-0.79212752) q[3];
sx q[3];
rz(-2.5274966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87913269) q[0];
sx q[0];
rz(-1.6430055) q[0];
sx q[0];
rz(0.60687989) q[0];
rz(2.795769) q[1];
sx q[1];
rz(-1.1864097) q[1];
sx q[1];
rz(0.80723673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9246638) q[0];
sx q[0];
rz(-1.47424) q[0];
sx q[0];
rz(-1.9159622) q[0];
rz(2.1431461) q[2];
sx q[2];
rz(-2.672451) q[2];
sx q[2];
rz(-0.97252211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8816384) q[1];
sx q[1];
rz(-1.7865766) q[1];
sx q[1];
rz(-0.27387932) q[1];
rz(-pi) q[2];
rz(0.59698481) q[3];
sx q[3];
rz(-1.0506949) q[3];
sx q[3];
rz(0.4515243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3240933) q[2];
sx q[2];
rz(-1.9976511) q[2];
sx q[2];
rz(2.6452046) q[2];
rz(-3.056622) q[3];
sx q[3];
rz(-0.43046633) q[3];
sx q[3];
rz(2.5607204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896511) q[0];
sx q[0];
rz(-1.713151) q[0];
sx q[0];
rz(0.38452837) q[0];
rz(2.8893068) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(-1.5581473) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2315518) q[0];
sx q[0];
rz(-1.1052026) q[0];
sx q[0];
rz(1.5601853) q[0];
x q[1];
rz(1.4531288) q[2];
sx q[2];
rz(-1.298438) q[2];
sx q[2];
rz(0.54335574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.059824316) q[1];
sx q[1];
rz(-1.5835973) q[1];
sx q[1];
rz(-0.98155419) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4676314) q[3];
sx q[3];
rz(-2.1836851) q[3];
sx q[3];
rz(2.6496668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(1.1400878) q[2];
rz(1.7831066) q[3];
sx q[3];
rz(-1.2525109) q[3];
sx q[3];
rz(-2.8235249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7110905) q[0];
sx q[0];
rz(-1.6999812) q[0];
sx q[0];
rz(2.7464113) q[0];
rz(1.2218062) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(0.88776678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783274) q[0];
sx q[0];
rz(-0.76416956) q[0];
sx q[0];
rz(-0.40212888) q[0];
x q[1];
rz(-1.5776921) q[2];
sx q[2];
rz(-0.47499945) q[2];
sx q[2];
rz(-1.4003225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50248442) q[1];
sx q[1];
rz(-0.35668761) q[1];
sx q[1];
rz(1.4977895) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0349109) q[3];
sx q[3];
rz(-1.8802207) q[3];
sx q[3];
rz(-0.14313652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6937574) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(-2.5698404) q[2];
rz(1.5104431) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(-1.3408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.754517) q[2];
sx q[2];
rz(-1.3420001) q[2];
sx q[2];
rz(-0.31172253) q[2];
rz(2.7243012) q[3];
sx q[3];
rz(-0.68796102) q[3];
sx q[3];
rz(1.7837379) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
