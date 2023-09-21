OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1169352) q[0];
sx q[0];
rz(-2.5664461) q[0];
sx q[0];
rz(2.2206109) q[0];
rz(-pi) q[1];
rz(-2.0338697) q[2];
sx q[2];
rz(-2.827364) q[2];
sx q[2];
rz(0.83480922) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1957789) q[1];
sx q[1];
rz(-1.073277) q[1];
sx q[1];
rz(1.0832018) q[1];
rz(0.63763036) q[3];
sx q[3];
rz(-2.5507567) q[3];
sx q[3];
rz(1.6319815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2279921) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.4555567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128132) q[0];
sx q[0];
rz(-1.7377186) q[0];
sx q[0];
rz(-2.1524327) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4886841) q[2];
sx q[2];
rz(-1.1973235) q[2];
sx q[2];
rz(0.79370802) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1728954) q[1];
sx q[1];
rz(-2.6186133) q[1];
sx q[1];
rz(-1.5437267) q[1];
x q[2];
rz(1.5442113) q[3];
sx q[3];
rz(-1.9823325) q[3];
sx q[3];
rz(2.6708024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.3519752) q[2];
rz(2.9591566) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(-1.9690537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4370255) q[0];
sx q[0];
rz(-2.4798658) q[0];
sx q[0];
rz(0.5163124) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3861888) q[2];
sx q[2];
rz(-1.3284151) q[2];
sx q[2];
rz(-2.541045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0769656) q[1];
sx q[1];
rz(-0.58276999) q[1];
sx q[1];
rz(2.0240192) q[1];
rz(-pi) q[2];
rz(-3.0828589) q[3];
sx q[3];
rz(-1.3204832) q[3];
sx q[3];
rz(-2.5854923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-0.91119901) q[2];
rz(-2.1905812) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(-3.0134841) q[0];
rz(3.065486) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(-0.52350837) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9273705) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(-0.58332304) q[0];
rz(-2.8115494) q[2];
sx q[2];
rz(-1.651262) q[2];
sx q[2];
rz(0.45804322) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.254863) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(-2.8121594) q[1];
rz(-pi) q[2];
rz(-1.3685437) q[3];
sx q[3];
rz(-2.5407255) q[3];
sx q[3];
rz(2.3893389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(0.564044) q[2];
rz(-0.28856746) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48150912) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447727) q[0];
sx q[0];
rz(-0.20475514) q[0];
sx q[0];
rz(-0.55069189) q[0];
x q[1];
rz(-1.5993824) q[2];
sx q[2];
rz(-2.8386142) q[2];
sx q[2];
rz(-0.44705331) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3494898) q[1];
sx q[1];
rz(-1.2621242) q[1];
sx q[1];
rz(-1.5285138) q[1];
rz(1.7080073) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1296967) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.72702423) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.7927992) q[0];
rz(-0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.7165002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2857367) q[0];
sx q[0];
rz(-0.55186134) q[0];
sx q[0];
rz(-0.43666552) q[0];
rz(3.0503057) q[2];
sx q[2];
rz(-1.8249776) q[2];
sx q[2];
rz(-0.19677481) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9565935) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(-1.1207629) q[1];
x q[2];
rz(-0.96958843) q[3];
sx q[3];
rz(-2.0242656) q[3];
sx q[3];
rz(0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0075334) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(-0.82180506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215295) q[0];
sx q[0];
rz(-0.95525817) q[0];
sx q[0];
rz(-2.2277742) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0242545) q[2];
sx q[2];
rz(-2.3356236) q[2];
sx q[2];
rz(2.1794127) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51470876) q[1];
sx q[1];
rz(-2.2885972) q[1];
sx q[1];
rz(1.2674598) q[1];
x q[2];
rz(0.52845593) q[3];
sx q[3];
rz(-0.65675694) q[3];
sx q[3];
rz(-2.9627851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3372779) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.725175) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.8364505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236429) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(-0.53074093) q[0];
rz(-pi) q[1];
rz(2.4370286) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(1.2003843) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.23755632) q[1];
sx q[1];
rz(-0.95103969) q[1];
sx q[1];
rz(1.1035641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9100902) q[3];
sx q[3];
rz(-1.5486071) q[3];
sx q[3];
rz(2.0458178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(0.70294356) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020333175) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(0.74777491) q[0];
x q[1];
rz(1.0614971) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(-0.73355567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.326509) q[1];
sx q[1];
rz(-2.106296) q[1];
sx q[1];
rz(-2.9498847) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1861107) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(0.51481065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(1.3396324) q[2];
rz(2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(-1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.16383485) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(0.46863619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9194473) q[0];
sx q[0];
rz(-1.7788017) q[0];
sx q[0];
rz(-2.749445) q[0];
rz(-2.8154545) q[2];
sx q[2];
rz(-1.9378127) q[2];
sx q[2];
rz(-2.0740167) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74433078) q[1];
sx q[1];
rz(-1.3052193) q[1];
sx q[1];
rz(3.0708405) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81810276) q[3];
sx q[3];
rz(-2.5460498) q[3];
sx q[3];
rz(-2.9593352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.1432077) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(-1.5293998) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951915) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-0.86482277) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(1.4680396) q[3];
sx q[3];
rz(-2.4874874) q[3];
sx q[3];
rz(-2.7174674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];