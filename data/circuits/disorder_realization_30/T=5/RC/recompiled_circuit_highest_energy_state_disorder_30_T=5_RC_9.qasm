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
rz(-2.8415866) q[0];
sx q[0];
rz(-2.0553148) q[0];
sx q[0];
rz(0.80817428) q[0];
rz(-1.1440682) q[1];
sx q[1];
rz(-0.77163982) q[1];
sx q[1];
rz(-1.5789403) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15641744) q[0];
sx q[0];
rz(-2.1083207) q[0];
sx q[0];
rz(0.82348768) q[0];
x q[1];
rz(-2.0879395) q[2];
sx q[2];
rz(-2.8364193) q[2];
sx q[2];
rz(0.15811731) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0275845) q[1];
sx q[1];
rz(-1.6817087) q[1];
sx q[1];
rz(-1.0618889) q[1];
rz(-pi) q[2];
rz(-1.6049625) q[3];
sx q[3];
rz(-1.5251659) q[3];
sx q[3];
rz(0.51866097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.943104) q[2];
sx q[2];
rz(-0.21185943) q[2];
sx q[2];
rz(2.1073821) q[2];
rz(2.4487623) q[3];
sx q[3];
rz(-2.0878017) q[3];
sx q[3];
rz(-2.0809765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8387872) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(-2.5651108) q[0];
rz(0.020596404) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(-1.0284665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45555025) q[0];
sx q[0];
rz(-2.8032113) q[0];
sx q[0];
rz(-0.69940059) q[0];
rz(-pi) q[1];
rz(1.1972381) q[2];
sx q[2];
rz(-1.8578863) q[2];
sx q[2];
rz(-2.8127828) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0513269) q[1];
sx q[1];
rz(-1.169186) q[1];
sx q[1];
rz(0.36220596) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1787492) q[3];
sx q[3];
rz(-0.96091753) q[3];
sx q[3];
rz(1.4489685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91745201) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(-1.2646328) q[2];
rz(-2.1680016) q[3];
sx q[3];
rz(-1.6907938) q[3];
sx q[3];
rz(2.8784331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50315404) q[0];
sx q[0];
rz(-1.3042903) q[0];
sx q[0];
rz(0.85357443) q[0];
rz(2.0897934) q[1];
sx q[1];
rz(-2.167326) q[1];
sx q[1];
rz(3.1265756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5511709) q[0];
sx q[0];
rz(-2.430924) q[0];
sx q[0];
rz(1.5573614) q[0];
rz(-0.096108561) q[2];
sx q[2];
rz(-1.5920361) q[2];
sx q[2];
rz(0.48890314) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7370809) q[1];
sx q[1];
rz(-0.85857262) q[1];
sx q[1];
rz(1.4732857) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.851753) q[3];
sx q[3];
rz(-1.6320758) q[3];
sx q[3];
rz(-2.0743528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0217648) q[2];
sx q[2];
rz(-1.4889577) q[2];
sx q[2];
rz(-0.28142288) q[2];
rz(-1.4540539) q[3];
sx q[3];
rz(-1.2097996) q[3];
sx q[3];
rz(-2.37229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5467095) q[0];
sx q[0];
rz(-1.2458206) q[0];
sx q[0];
rz(-0.17387867) q[0];
rz(-1.8100544) q[1];
sx q[1];
rz(-1.7200108) q[1];
sx q[1];
rz(1.7132267) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9857432) q[0];
sx q[0];
rz(-1.5850889) q[0];
sx q[0];
rz(-1.3784268) q[0];
x q[1];
rz(-1.6971484) q[2];
sx q[2];
rz(-1.9280806) q[2];
sx q[2];
rz(2.7972691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5133068) q[1];
sx q[1];
rz(-0.9137972) q[1];
sx q[1];
rz(-1.6194067) q[1];
rz(-pi) q[2];
rz(1.9589728) q[3];
sx q[3];
rz(-1.2283162) q[3];
sx q[3];
rz(1.6714718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4398769) q[2];
sx q[2];
rz(-2.3150847) q[2];
sx q[2];
rz(-2.7244549) q[2];
rz(-0.16962984) q[3];
sx q[3];
rz(-0.093234213) q[3];
sx q[3];
rz(-2.4197742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193264) q[0];
sx q[0];
rz(-2.7665783) q[0];
sx q[0];
rz(2.8322423) q[0];
rz(-0.7695235) q[1];
sx q[1];
rz(-2.0517495) q[1];
sx q[1];
rz(-1.5187029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8474591) q[0];
sx q[0];
rz(-1.8222229) q[0];
sx q[0];
rz(-0.88788106) q[0];
x q[1];
rz(-0.74414545) q[2];
sx q[2];
rz(-2.2501906) q[2];
sx q[2];
rz(-1.3199922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.57598439) q[1];
sx q[1];
rz(-0.19153015) q[1];
sx q[1];
rz(1.7267141) q[1];
x q[2];
rz(-2.2045787) q[3];
sx q[3];
rz(-2.975836) q[3];
sx q[3];
rz(-2.9637333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3994483) q[2];
sx q[2];
rz(-2.1592906) q[2];
sx q[2];
rz(0.99545109) q[2];
rz(0.71850145) q[3];
sx q[3];
rz(-1.813443) q[3];
sx q[3];
rz(2.3464581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8566078) q[0];
sx q[0];
rz(-1.4549078) q[0];
sx q[0];
rz(1.6108151) q[0];
rz(-1.6948304) q[1];
sx q[1];
rz(-1.4434573) q[1];
sx q[1];
rz(1.7036499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7890678) q[0];
sx q[0];
rz(-1.1630327) q[0];
sx q[0];
rz(3.1192794) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9447359) q[2];
sx q[2];
rz(-2.1293921) q[2];
sx q[2];
rz(2.2064759) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9077545) q[1];
sx q[1];
rz(-1.4871039) q[1];
sx q[1];
rz(-0.5096883) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5737278) q[3];
sx q[3];
rz(-0.81122196) q[3];
sx q[3];
rz(-2.7377759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97856727) q[2];
sx q[2];
rz(-1.359553) q[2];
sx q[2];
rz(-1.2320409) q[2];
rz(1.9949404) q[3];
sx q[3];
rz(-1.8012643) q[3];
sx q[3];
rz(0.31057096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50800407) q[0];
sx q[0];
rz(-0.81729832) q[0];
sx q[0];
rz(-2.001413) q[0];
rz(-2.2445402) q[1];
sx q[1];
rz(-1.8042118) q[1];
sx q[1];
rz(-3.1111029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59647467) q[0];
sx q[0];
rz(-1.6870572) q[0];
sx q[0];
rz(1.7932939) q[0];
x q[1];
rz(-2.4989481) q[2];
sx q[2];
rz(-1.8360999) q[2];
sx q[2];
rz(-2.93769) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5329689) q[1];
sx q[1];
rz(-0.86809671) q[1];
sx q[1];
rz(-2.7547902) q[1];
x q[2];
rz(-1.072149) q[3];
sx q[3];
rz(-0.79550084) q[3];
sx q[3];
rz(1.7444057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9474779) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(1.6537559) q[3];
sx q[3];
rz(-1.0349118) q[3];
sx q[3];
rz(1.8208767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79539913) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(1.7561703) q[0];
rz(-1.9253383) q[1];
sx q[1];
rz(-1.6938035) q[1];
sx q[1];
rz(2.6522955) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70883553) q[0];
sx q[0];
rz(-2.2422195) q[0];
sx q[0];
rz(-1.5455724) q[0];
rz(-0.036542372) q[2];
sx q[2];
rz(-1.9068235) q[2];
sx q[2];
rz(2.450409) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7304436) q[1];
sx q[1];
rz(-0.48161067) q[1];
sx q[1];
rz(-1.7754619) q[1];
rz(1.8200666) q[3];
sx q[3];
rz(-1.6695654) q[3];
sx q[3];
rz(0.45757142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3936002) q[2];
sx q[2];
rz(-1.9990498) q[2];
sx q[2];
rz(-0.46621123) q[2];
rz(1.5318058) q[3];
sx q[3];
rz(-1.6377056) q[3];
sx q[3];
rz(2.8209749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521249) q[0];
sx q[0];
rz(-2.2418699) q[0];
sx q[0];
rz(3.0925282) q[0];
rz(2.9009254) q[1];
sx q[1];
rz(-2.1235695) q[1];
sx q[1];
rz(0.022445591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72180333) q[0];
sx q[0];
rz(-1.5947486) q[0];
sx q[0];
rz(0.99138685) q[0];
rz(2.4149553) q[2];
sx q[2];
rz(-0.49984806) q[2];
sx q[2];
rz(2.8195087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3901523) q[1];
sx q[1];
rz(-0.94822394) q[1];
sx q[1];
rz(2.1639813) q[1];
x q[2];
rz(-2.4565036) q[3];
sx q[3];
rz(-1.2108608) q[3];
sx q[3];
rz(-2.0956831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13751328) q[2];
sx q[2];
rz(-1.9065607) q[2];
sx q[2];
rz(2.1732886) q[2];
rz(-0.83640313) q[3];
sx q[3];
rz(-2.8439549) q[3];
sx q[3];
rz(-1.3692726) q[3];
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
rz(1.4150998) q[0];
sx q[0];
rz(-1.8510171) q[0];
sx q[0];
rz(2.7628164) q[0];
rz(-3.1355766) q[1];
sx q[1];
rz(-0.52979398) q[1];
sx q[1];
rz(2.5078497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0519052) q[0];
sx q[0];
rz(-1.9080592) q[0];
sx q[0];
rz(-2.428672) q[0];
x q[1];
rz(-0.1600645) q[2];
sx q[2];
rz(-1.0032005) q[2];
sx q[2];
rz(2.4307323) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4937456) q[1];
sx q[1];
rz(-2.4253107) q[1];
sx q[1];
rz(2.358308) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2489088) q[3];
sx q[3];
rz(-2.0555946) q[3];
sx q[3];
rz(-2.3546653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4358431) q[2];
sx q[2];
rz(-1.6195932) q[2];
sx q[2];
rz(-0.61326927) q[2];
rz(-0.67522007) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(0.81775445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9153862) q[0];
sx q[0];
rz(-1.362726) q[0];
sx q[0];
rz(1.2363731) q[0];
rz(-0.18648237) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(-1.5507422) q[2];
sx q[2];
rz(-1.1112083) q[2];
sx q[2];
rz(0.61136897) q[2];
rz(2.1929838) q[3];
sx q[3];
rz(-1.4962248) q[3];
sx q[3];
rz(0.20635508) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
