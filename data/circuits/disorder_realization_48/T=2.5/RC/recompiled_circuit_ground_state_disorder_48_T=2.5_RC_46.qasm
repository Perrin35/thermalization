OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1005062) q[0];
sx q[0];
rz(-2.9334928) q[0];
sx q[0];
rz(-2.7383374) q[0];
rz(-0.23303214) q[1];
sx q[1];
rz(-1.4401399) q[1];
sx q[1];
rz(-2.9177102) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96432811) q[0];
sx q[0];
rz(-1.3601662) q[0];
sx q[0];
rz(-0.62646477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7213905) q[2];
sx q[2];
rz(-0.70993844) q[2];
sx q[2];
rz(0.89872724) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76624291) q[1];
sx q[1];
rz(-2.2170984) q[1];
sx q[1];
rz(2.8588041) q[1];
x q[2];
rz(0.68239642) q[3];
sx q[3];
rz(-2.2060437) q[3];
sx q[3];
rz(-1.531383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4441121) q[2];
sx q[2];
rz(-0.31284249) q[2];
sx q[2];
rz(0.035813896) q[2];
rz(2.4114285) q[3];
sx q[3];
rz(-1.9883479) q[3];
sx q[3];
rz(2.9353976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8908454) q[0];
sx q[0];
rz(-0.99531168) q[0];
sx q[0];
rz(-2.9440951) q[0];
rz(-2.6644871) q[1];
sx q[1];
rz(-0.66480607) q[1];
sx q[1];
rz(0.8055996) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8715848) q[0];
sx q[0];
rz(-0.91579899) q[0];
sx q[0];
rz(1.3359469) q[0];
rz(-1.637403) q[2];
sx q[2];
rz(-1.9506467) q[2];
sx q[2];
rz(2.0646281) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2239469) q[1];
sx q[1];
rz(-2.2093281) q[1];
sx q[1];
rz(-2.3238514) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3918341) q[3];
sx q[3];
rz(-0.20515144) q[3];
sx q[3];
rz(1.6501901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8311367) q[2];
sx q[2];
rz(-1.4996935) q[2];
sx q[2];
rz(-1.4031225) q[2];
rz(0.63203114) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(-3.0145751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3056575) q[0];
sx q[0];
rz(-2.5129565) q[0];
sx q[0];
rz(-1.0868616) q[0];
rz(0.19733363) q[1];
sx q[1];
rz(-1.5060164) q[1];
sx q[1];
rz(-2.8089583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5088288) q[0];
sx q[0];
rz(-0.45464215) q[0];
sx q[0];
rz(1.2484545) q[0];
rz(-pi) q[1];
x q[1];
rz(1.640075) q[2];
sx q[2];
rz(-0.73613077) q[2];
sx q[2];
rz(2.2683805) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91326317) q[1];
sx q[1];
rz(-2.8614534) q[1];
sx q[1];
rz(2.254451) q[1];
rz(-pi) q[2];
rz(-0.72686355) q[3];
sx q[3];
rz(-1.3976025) q[3];
sx q[3];
rz(-1.3722591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7670333) q[2];
sx q[2];
rz(-1.550753) q[2];
sx q[2];
rz(1.6732875) q[2];
rz(-3.1324006) q[3];
sx q[3];
rz(-0.46949783) q[3];
sx q[3];
rz(2.3495242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62995768) q[0];
sx q[0];
rz(-2.503643) q[0];
sx q[0];
rz(1.6988423) q[0];
rz(2.1171782) q[1];
sx q[1];
rz(-1.9099648) q[1];
sx q[1];
rz(-0.82694298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8762739) q[0];
sx q[0];
rz(-1.4965222) q[0];
sx q[0];
rz(-1.7552046) q[0];
rz(-pi) q[1];
rz(-2.2587772) q[2];
sx q[2];
rz(-1.6903631) q[2];
sx q[2];
rz(-1.2403229) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0961541) q[1];
sx q[1];
rz(-1.3609582) q[1];
sx q[1];
rz(-2.7010598) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0204861) q[3];
sx q[3];
rz(-0.89205974) q[3];
sx q[3];
rz(0.27287441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3247165) q[2];
sx q[2];
rz(-2.3574895) q[2];
sx q[2];
rz(-1.8640222) q[2];
rz(-0.68814021) q[3];
sx q[3];
rz(-2.1161931) q[3];
sx q[3];
rz(0.96243206) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0012896) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(-0.21743123) q[0];
rz(-2.8561719) q[1];
sx q[1];
rz(-2.5358584) q[1];
sx q[1];
rz(2.6434456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.69553) q[0];
sx q[0];
rz(-2.208606) q[0];
sx q[0];
rz(2.7317156) q[0];
x q[1];
rz(-1.1310546) q[2];
sx q[2];
rz(-0.56098191) q[2];
sx q[2];
rz(-1.9064685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0302311) q[1];
sx q[1];
rz(-1.7788789) q[1];
sx q[1];
rz(2.0672805) q[1];
rz(-pi) q[2];
rz(0.44002779) q[3];
sx q[3];
rz(-1.8002568) q[3];
sx q[3];
rz(-3.0408531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6238326) q[2];
sx q[2];
rz(-2.0443003) q[2];
sx q[2];
rz(2.9916054) q[2];
rz(3.127626) q[3];
sx q[3];
rz(-1.5499127) q[3];
sx q[3];
rz(-0.51378957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6496395) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(0.81277043) q[0];
rz(1.4179519) q[1];
sx q[1];
rz(-1.5128472) q[1];
sx q[1];
rz(1.0636122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9169711) q[0];
sx q[0];
rz(-1.2422891) q[0];
sx q[0];
rz(-2.5909831) q[0];
rz(-1.5539507) q[2];
sx q[2];
rz(-1.7344513) q[2];
sx q[2];
rz(2.9962199) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8955161) q[1];
sx q[1];
rz(-0.59883307) q[1];
sx q[1];
rz(-0.020222874) q[1];
x q[2];
rz(-1.316458) q[3];
sx q[3];
rz(-1.3876378) q[3];
sx q[3];
rz(2.0034127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26760462) q[2];
sx q[2];
rz(-2.5817817) q[2];
sx q[2];
rz(1.416729) q[2];
rz(-3.0806165) q[3];
sx q[3];
rz(-2.2305326) q[3];
sx q[3];
rz(1.4643668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7199719) q[0];
sx q[0];
rz(-2.3201729) q[0];
sx q[0];
rz(-3.022497) q[0];
rz(-2.6630317) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(0.68971577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0499303) q[0];
sx q[0];
rz(-0.99794594) q[0];
sx q[0];
rz(-0.86688231) q[0];
rz(-pi) q[1];
rz(2.5096171) q[2];
sx q[2];
rz(-2.5594829) q[2];
sx q[2];
rz(-1.7012973) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3773681) q[1];
sx q[1];
rz(-2.2768094) q[1];
sx q[1];
rz(-0.059125916) q[1];
x q[2];
rz(-0.14755149) q[3];
sx q[3];
rz(-2.2780212) q[3];
sx q[3];
rz(0.13015166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0488284) q[2];
sx q[2];
rz(-0.060828716) q[2];
sx q[2];
rz(1.0678585) q[2];
rz(-0.16718665) q[3];
sx q[3];
rz(-1.3719631) q[3];
sx q[3];
rz(-3.0021477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.5329413) q[0];
sx q[0];
rz(-2.3280188) q[0];
sx q[0];
rz(0.90840489) q[0];
rz(0.99336973) q[1];
sx q[1];
rz(-2.9790331) q[1];
sx q[1];
rz(-2.2672674) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2539815) q[0];
sx q[0];
rz(-1.7029666) q[0];
sx q[0];
rz(1.553276) q[0];
rz(-pi) q[1];
rz(1.6443917) q[2];
sx q[2];
rz(-2.8831867) q[2];
sx q[2];
rz(-2.0143353) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0093551) q[1];
sx q[1];
rz(-1.435408) q[1];
sx q[1];
rz(-0.05986771) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0215346) q[3];
sx q[3];
rz(-0.62432271) q[3];
sx q[3];
rz(-2.2933427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5270762) q[2];
sx q[2];
rz(-2.9824342) q[2];
sx q[2];
rz(2.2873774) q[2];
rz(0.45520511) q[3];
sx q[3];
rz(-0.56536094) q[3];
sx q[3];
rz(-2.8860886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47827569) q[0];
sx q[0];
rz(-1.9238967) q[0];
sx q[0];
rz(1.664337) q[0];
rz(1.5746337) q[1];
sx q[1];
rz(-0.68395558) q[1];
sx q[1];
rz(2.4050567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4111209) q[0];
sx q[0];
rz(-0.53148848) q[0];
sx q[0];
rz(2.0474252) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8868448) q[2];
sx q[2];
rz(-2.7494135) q[2];
sx q[2];
rz(2.2318411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7334637) q[1];
sx q[1];
rz(-0.40181364) q[1];
sx q[1];
rz(0.6647474) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1062713) q[3];
sx q[3];
rz(-0.97408726) q[3];
sx q[3];
rz(2.2281856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7443202) q[2];
sx q[2];
rz(-0.87817764) q[2];
sx q[2];
rz(1.5273904) q[2];
rz(-2.7496036) q[3];
sx q[3];
rz(-1.9422928) q[3];
sx q[3];
rz(3.0919302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024260661) q[0];
sx q[0];
rz(-1.6680822) q[0];
sx q[0];
rz(0.20183739) q[0];
rz(2.9182538) q[1];
sx q[1];
rz(-0.52450648) q[1];
sx q[1];
rz(-0.39252678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72001624) q[0];
sx q[0];
rz(-0.47552738) q[0];
sx q[0];
rz(-1.8029193) q[0];
x q[1];
rz(1.5129651) q[2];
sx q[2];
rz(-2.4041345) q[2];
sx q[2];
rz(-2.2493169) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0328249) q[1];
sx q[1];
rz(-0.6995753) q[1];
sx q[1];
rz(-0.8731858) q[1];
rz(-2.659117) q[3];
sx q[3];
rz(-1.2996718) q[3];
sx q[3];
rz(1.9737873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7212123) q[2];
sx q[2];
rz(-2.1537697) q[2];
sx q[2];
rz(-2.5610899) q[2];
rz(-0.37603363) q[3];
sx q[3];
rz(-2.1707363) q[3];
sx q[3];
rz(1.6002801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9911983) q[0];
sx q[0];
rz(-1.5924441) q[0];
sx q[0];
rz(1.7572255) q[0];
rz(-1.4906384) q[1];
sx q[1];
rz(-1.8883659) q[1];
sx q[1];
rz(-2.2813003) q[1];
rz(0.80581325) q[2];
sx q[2];
rz(-2.056682) q[2];
sx q[2];
rz(-0.32120612) q[2];
rz(-1.5932455) q[3];
sx q[3];
rz(-1.2590564) q[3];
sx q[3];
rz(1.1127478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
