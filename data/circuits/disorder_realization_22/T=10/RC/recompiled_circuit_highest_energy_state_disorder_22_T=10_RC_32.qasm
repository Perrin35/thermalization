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
rz(0.67233664) q[0];
sx q[0];
rz(-1.2824751) q[0];
sx q[0];
rz(1.6096492) q[0];
rz(-1.9077644) q[1];
sx q[1];
rz(-1.2169853) q[1];
sx q[1];
rz(1.4438862) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9836228) q[0];
sx q[0];
rz(-0.90216547) q[0];
sx q[0];
rz(-0.10413952) q[0];
rz(-pi) q[1];
rz(-1.9305761) q[2];
sx q[2];
rz(-1.8471247) q[2];
sx q[2];
rz(2.2587549) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0585079) q[1];
sx q[1];
rz(-2.1821686) q[1];
sx q[1];
rz(-1.0669003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.231655) q[3];
sx q[3];
rz(-1.5970486) q[3];
sx q[3];
rz(-1.8363801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7794789) q[2];
sx q[2];
rz(-1.9356091) q[2];
sx q[2];
rz(0.77902478) q[2];
rz(-1.3075167) q[3];
sx q[3];
rz(-2.8179759) q[3];
sx q[3];
rz(1.754508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5883412) q[0];
sx q[0];
rz(-1.3655751) q[0];
sx q[0];
rz(-0.27728444) q[0];
rz(-2.7514027) q[1];
sx q[1];
rz(-2.0398102) q[1];
sx q[1];
rz(1.1276833) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7285293) q[0];
sx q[0];
rz(-1.5908233) q[0];
sx q[0];
rz(-3.0990451) q[0];
rz(2.8088536) q[2];
sx q[2];
rz(-2.840452) q[2];
sx q[2];
rz(2.6469165) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8816205) q[1];
sx q[1];
rz(-1.299159) q[1];
sx q[1];
rz(0.064898811) q[1];
x q[2];
rz(1.1902749) q[3];
sx q[3];
rz(-0.94469686) q[3];
sx q[3];
rz(2.2230119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4358501) q[2];
sx q[2];
rz(-2.5472992) q[2];
sx q[2];
rz(1.1951813) q[2];
rz(2.4863906) q[3];
sx q[3];
rz(-1.4846669) q[3];
sx q[3];
rz(-0.5932194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4707659) q[0];
sx q[0];
rz(-2.4556181) q[0];
sx q[0];
rz(-2.2210448) q[0];
rz(1.8596733) q[1];
sx q[1];
rz(-2.0125407) q[1];
sx q[1];
rz(-1.4528073) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0301967) q[0];
sx q[0];
rz(-1.483484) q[0];
sx q[0];
rz(-3.0285378) q[0];
x q[1];
rz(2.0412943) q[2];
sx q[2];
rz(-0.62242939) q[2];
sx q[2];
rz(-0.4284455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7662168) q[1];
sx q[1];
rz(-1.9060684) q[1];
sx q[1];
rz(-1.6250618) q[1];
rz(-pi) q[2];
rz(3.0759144) q[3];
sx q[3];
rz(-1.6367776) q[3];
sx q[3];
rz(-3.0998041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.526223) q[2];
sx q[2];
rz(-2.0873895) q[2];
sx q[2];
rz(0.39895454) q[2];
rz(-0.51359549) q[3];
sx q[3];
rz(-0.31702888) q[3];
sx q[3];
rz(-1.1494466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.32617635) q[0];
sx q[0];
rz(-0.57891095) q[0];
sx q[0];
rz(1.2166566) q[0];
rz(0.51219621) q[1];
sx q[1];
rz(-2.0307816) q[1];
sx q[1];
rz(2.0373352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062631) q[0];
sx q[0];
rz(-1.2699513) q[0];
sx q[0];
rz(1.8604688) q[0];
x q[1];
rz(1.3894559) q[2];
sx q[2];
rz(-1.7415819) q[2];
sx q[2];
rz(-2.6362399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3901375) q[1];
sx q[1];
rz(-1.631076) q[1];
sx q[1];
rz(-1.7675318) q[1];
rz(-pi) q[2];
rz(-2.2858972) q[3];
sx q[3];
rz(-0.33791204) q[3];
sx q[3];
rz(-2.5355946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5451374) q[2];
sx q[2];
rz(-1.1918273) q[2];
sx q[2];
rz(0.54984251) q[2];
rz(0.95692974) q[3];
sx q[3];
rz(-2.3423829) q[3];
sx q[3];
rz(-1.4889312) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6560219) q[0];
sx q[0];
rz(-0.85000426) q[0];
sx q[0];
rz(-0.75918424) q[0];
rz(0.94866577) q[1];
sx q[1];
rz(-1.1439088) q[1];
sx q[1];
rz(-2.8099828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0873957) q[0];
sx q[0];
rz(-0.83354649) q[0];
sx q[0];
rz(1.9760608) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.070721432) q[2];
sx q[2];
rz(-2.6658733) q[2];
sx q[2];
rz(2.4782319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6575106) q[1];
sx q[1];
rz(-0.23536029) q[1];
sx q[1];
rz(-0.85296191) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4344352) q[3];
sx q[3];
rz(-1.9788187) q[3];
sx q[3];
rz(1.3566164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9552976) q[2];
sx q[2];
rz(-2.1897327) q[2];
sx q[2];
rz(-0.9430421) q[2];
rz(-0.17213639) q[3];
sx q[3];
rz(-2.19682) q[3];
sx q[3];
rz(2.9183563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.405769) q[0];
sx q[0];
rz(-0.74546927) q[0];
sx q[0];
rz(-1.5752342) q[0];
rz(2.9834421) q[1];
sx q[1];
rz(-0.76845303) q[1];
sx q[1];
rz(0.92811531) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1386618) q[0];
sx q[0];
rz(-1.2671906) q[0];
sx q[0];
rz(-1.4360484) q[0];
rz(-pi) q[1];
rz(-2.1238223) q[2];
sx q[2];
rz(-1.206351) q[2];
sx q[2];
rz(1.5096957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5852768) q[1];
sx q[1];
rz(-0.92570526) q[1];
sx q[1];
rz(2.8200802) q[1];
x q[2];
rz(-2.8148303) q[3];
sx q[3];
rz(-2.770677) q[3];
sx q[3];
rz(2.4075631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1930262) q[2];
sx q[2];
rz(-0.28634772) q[2];
sx q[2];
rz(2.7890653) q[2];
rz(-0.36892095) q[3];
sx q[3];
rz(-0.58266321) q[3];
sx q[3];
rz(-2.2841456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7848659) q[0];
sx q[0];
rz(-3.1251188) q[0];
sx q[0];
rz(0.67100588) q[0];
rz(-3.0352133) q[1];
sx q[1];
rz(-0.89768592) q[1];
sx q[1];
rz(-0.62972128) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24983588) q[0];
sx q[0];
rz(-1.5895004) q[0];
sx q[0];
rz(-2.9500701) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1415954) q[2];
sx q[2];
rz(-1.1941807) q[2];
sx q[2];
rz(-1.9887599) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2196673) q[1];
sx q[1];
rz(-1.6451854) q[1];
sx q[1];
rz(-0.32385918) q[1];
rz(-1.8811781) q[3];
sx q[3];
rz(-0.93573007) q[3];
sx q[3];
rz(0.47582754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0715535) q[2];
sx q[2];
rz(-2.2420501) q[2];
sx q[2];
rz(1.0413991) q[2];
rz(-1.533016) q[3];
sx q[3];
rz(-0.93419111) q[3];
sx q[3];
rz(-1.9937203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1958985) q[0];
sx q[0];
rz(-0.65106374) q[0];
sx q[0];
rz(0.4796637) q[0];
rz(0.02462968) q[1];
sx q[1];
rz(-1.004091) q[1];
sx q[1];
rz(3.0020795) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937353) q[0];
sx q[0];
rz(-1.4850886) q[0];
sx q[0];
rz(1.6802963) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.191338) q[2];
sx q[2];
rz(-2.5436311) q[2];
sx q[2];
rz(2.4520252) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7558328) q[1];
sx q[1];
rz(-0.91948286) q[1];
sx q[1];
rz(-1.3612222) q[1];
rz(-pi) q[2];
rz(1.6253396) q[3];
sx q[3];
rz(-2.2404376) q[3];
sx q[3];
rz(-3.0564133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2554243) q[2];
sx q[2];
rz(-1.5443799) q[2];
sx q[2];
rz(1.9891948) q[2];
rz(1.7216916) q[3];
sx q[3];
rz(-2.4072188) q[3];
sx q[3];
rz(1.3091807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060147978) q[0];
sx q[0];
rz(-1.6763433) q[0];
sx q[0];
rz(-1.6434705) q[0];
rz(-0.76088798) q[1];
sx q[1];
rz(-2.1983169) q[1];
sx q[1];
rz(1.4563742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0022572) q[0];
sx q[0];
rz(-2.5482448) q[0];
sx q[0];
rz(-0.32440455) q[0];
rz(2.7414906) q[2];
sx q[2];
rz(-1.7662373) q[2];
sx q[2];
rz(0.2646499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7356811) q[1];
sx q[1];
rz(-1.0764437) q[1];
sx q[1];
rz(2.6895061) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7237503) q[3];
sx q[3];
rz(-1.0438232) q[3];
sx q[3];
rz(2.9067519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8196408) q[2];
sx q[2];
rz(-0.4343304) q[2];
sx q[2];
rz(-2.1257909) q[2];
rz(-2.2470233) q[3];
sx q[3];
rz(-1.2457448) q[3];
sx q[3];
rz(2.3999124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1557409) q[0];
sx q[0];
rz(-0.72991252) q[0];
sx q[0];
rz(-2.1515382) q[0];
rz(0.92297018) q[1];
sx q[1];
rz(-1.762941) q[1];
sx q[1];
rz(-2.1754481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7434563) q[0];
sx q[0];
rz(-1.4706637) q[0];
sx q[0];
rz(-3.0707573) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82426333) q[2];
sx q[2];
rz(-2.5474862) q[2];
sx q[2];
rz(0.73005331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50033874) q[1];
sx q[1];
rz(-0.49420824) q[1];
sx q[1];
rz(-0.92730913) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0279492) q[3];
sx q[3];
rz(-1.7688278) q[3];
sx q[3];
rz(-2.6012492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82603377) q[2];
sx q[2];
rz(-2.3580599) q[2];
sx q[2];
rz(2.3890736) q[2];
rz(0.13713947) q[3];
sx q[3];
rz(-2.5004041) q[3];
sx q[3];
rz(-0.34998676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.242545) q[0];
sx q[0];
rz(-0.62472961) q[0];
sx q[0];
rz(-0.20150264) q[0];
rz(1.7789727) q[1];
sx q[1];
rz(-0.87281223) q[1];
sx q[1];
rz(2.9760117) q[1];
rz(1.4819862) q[2];
sx q[2];
rz(-1.6203151) q[2];
sx q[2];
rz(2.0697894) q[2];
rz(-0.84861075) q[3];
sx q[3];
rz(-0.52042421) q[3];
sx q[3];
rz(2.0757794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
