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
rz(-0.041115887) q[0];
sx q[0];
rz(-1.9408722) q[0];
sx q[0];
rz(1.6706985) q[0];
rz(-0.38493758) q[1];
sx q[1];
rz(-0.5783143) q[1];
sx q[1];
rz(-0.90775031) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2157057) q[0];
sx q[0];
rz(-0.74630794) q[0];
sx q[0];
rz(-2.4624636) q[0];
rz(-2.3458751) q[2];
sx q[2];
rz(-2.7331686) q[2];
sx q[2];
rz(-1.0701831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8239535) q[1];
sx q[1];
rz(-1.7574218) q[1];
sx q[1];
rz(-0.49368787) q[1];
x q[2];
rz(3.0201245) q[3];
sx q[3];
rz(-1.5163912) q[3];
sx q[3];
rz(-1.7732946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3140728) q[2];
sx q[2];
rz(-1.4853518) q[2];
sx q[2];
rz(-2.8208288) q[2];
rz(-3.1406506) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(0.48085406) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28474057) q[0];
sx q[0];
rz(-2.4212403) q[0];
sx q[0];
rz(2.5208933) q[0];
rz(-1.2868098) q[1];
sx q[1];
rz(-1.6544673) q[1];
sx q[1];
rz(0.99498814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10886562) q[0];
sx q[0];
rz(-1.0720729) q[0];
sx q[0];
rz(2.9176955) q[0];
rz(2.0109977) q[2];
sx q[2];
rz(-1.2792247) q[2];
sx q[2];
rz(2.7879813) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80985132) q[1];
sx q[1];
rz(-2.6552422) q[1];
sx q[1];
rz(1.9736109) q[1];
rz(-pi) q[2];
rz(-0.50813679) q[3];
sx q[3];
rz(-2.2656815) q[3];
sx q[3];
rz(0.14665237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93028012) q[2];
sx q[2];
rz(-1.899753) q[2];
sx q[2];
rz(2.9105183) q[2];
rz(1.5288345) q[3];
sx q[3];
rz(-1.4569747) q[3];
sx q[3];
rz(0.55135623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7174317) q[0];
sx q[0];
rz(-1.9743974) q[0];
sx q[0];
rz(-2.7447847) q[0];
rz(1.1948168) q[1];
sx q[1];
rz(-1.2620986) q[1];
sx q[1];
rz(-1.6669115) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2541324) q[0];
sx q[0];
rz(-1.6586705) q[0];
sx q[0];
rz(-0.011482419) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72309317) q[2];
sx q[2];
rz(-2.4567219) q[2];
sx q[2];
rz(-1.9346646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35016216) q[1];
sx q[1];
rz(-1.2315948) q[1];
sx q[1];
rz(-0.8703402) q[1];
rz(1.9649404) q[3];
sx q[3];
rz(-1.5817465) q[3];
sx q[3];
rz(-1.8979069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2954703) q[2];
sx q[2];
rz(-0.84540558) q[2];
sx q[2];
rz(1.4836614) q[2];
rz(2.1814003) q[3];
sx q[3];
rz(-0.63513297) q[3];
sx q[3];
rz(0.6944164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8089777) q[0];
sx q[0];
rz(-1.3354477) q[0];
sx q[0];
rz(2.4131925) q[0];
rz(-1.2737466) q[1];
sx q[1];
rz(-1.961901) q[1];
sx q[1];
rz(1.5625928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0141143) q[0];
sx q[0];
rz(-1.593361) q[0];
sx q[0];
rz(-0.0033408611) q[0];
x q[1];
rz(-0.043508675) q[2];
sx q[2];
rz(-2.1183421) q[2];
sx q[2];
rz(2.7491153) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82484805) q[1];
sx q[1];
rz(-1.4364752) q[1];
sx q[1];
rz(-0.21726852) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1636435) q[3];
sx q[3];
rz(-2.4302357) q[3];
sx q[3];
rz(2.9748981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1100715) q[2];
sx q[2];
rz(-1.7379802) q[2];
sx q[2];
rz(-0.20315476) q[2];
rz(-1.7848232) q[3];
sx q[3];
rz(-1.9316831) q[3];
sx q[3];
rz(1.320896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2585231) q[0];
sx q[0];
rz(-1.8061545) q[0];
sx q[0];
rz(-1.6254599) q[0];
rz(-0.37172231) q[1];
sx q[1];
rz(-2.5860791) q[1];
sx q[1];
rz(-0.98977596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76565874) q[0];
sx q[0];
rz(-1.759937) q[0];
sx q[0];
rz(-2.2098757) q[0];
rz(-1.4577652) q[2];
sx q[2];
rz(-2.340132) q[2];
sx q[2];
rz(-2.46704) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.02640662) q[1];
sx q[1];
rz(-0.37749981) q[1];
sx q[1];
rz(2.0707692) q[1];
rz(-1.2618213) q[3];
sx q[3];
rz(-2.0403962) q[3];
sx q[3];
rz(2.0275786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61701361) q[2];
sx q[2];
rz(-0.22713529) q[2];
sx q[2];
rz(-3.0830834) q[2];
rz(1.2837563) q[3];
sx q[3];
rz(-2.3144898) q[3];
sx q[3];
rz(-1.9730998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7002895) q[0];
sx q[0];
rz(-2.0952201) q[0];
sx q[0];
rz(2.6562279) q[0];
rz(-0.46733388) q[1];
sx q[1];
rz(-0.58715564) q[1];
sx q[1];
rz(0.68624085) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29083911) q[0];
sx q[0];
rz(-0.79054442) q[0];
sx q[0];
rz(-0.45773069) q[0];
x q[1];
rz(0.94450343) q[2];
sx q[2];
rz(-1.7472225) q[2];
sx q[2];
rz(-3.021332) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31364775) q[1];
sx q[1];
rz(-2.0454496) q[1];
sx q[1];
rz(0.56874911) q[1];
x q[2];
rz(2.521311) q[3];
sx q[3];
rz(-0.93343319) q[3];
sx q[3];
rz(2.8220774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.301959) q[2];
sx q[2];
rz(-1.0884103) q[2];
sx q[2];
rz(0.23183091) q[2];
rz(0.92136446) q[3];
sx q[3];
rz(-1.6145584) q[3];
sx q[3];
rz(3.11619) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63474083) q[0];
sx q[0];
rz(-2.1379819) q[0];
sx q[0];
rz(0.64742175) q[0];
rz(-2.1000775) q[1];
sx q[1];
rz(-1.2728649) q[1];
sx q[1];
rz(2.0288859) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9346973) q[0];
sx q[0];
rz(-2.0152038) q[0];
sx q[0];
rz(-1.0133176) q[0];
x q[1];
rz(1.4278863) q[2];
sx q[2];
rz(-1.630703) q[2];
sx q[2];
rz(0.026345677) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7156421) q[1];
sx q[1];
rz(-2.1248769) q[1];
sx q[1];
rz(1.7765846) q[1];
rz(-1.5912368) q[3];
sx q[3];
rz(-2.6764565) q[3];
sx q[3];
rz(2.6488369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.57704321) q[2];
sx q[2];
rz(-2.3304522) q[2];
sx q[2];
rz(2.2570293) q[2];
rz(1.0386284) q[3];
sx q[3];
rz(-1.1218772) q[3];
sx q[3];
rz(-1.60892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7811979) q[0];
sx q[0];
rz(-2.9304657) q[0];
sx q[0];
rz(-1.1398844) q[0];
rz(0.79849157) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(-3.0671157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1482753) q[0];
sx q[0];
rz(-1.8145942) q[0];
sx q[0];
rz(-0.93945844) q[0];
x q[1];
rz(-1.1604606) q[2];
sx q[2];
rz(-0.85121819) q[2];
sx q[2];
rz(-1.5026523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3466693) q[1];
sx q[1];
rz(-2.2058869) q[1];
sx q[1];
rz(-1.1785202) q[1];
rz(-pi) q[2];
rz(-2.8387855) q[3];
sx q[3];
rz(-1.9085371) q[3];
sx q[3];
rz(0.95888058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8188339) q[2];
sx q[2];
rz(-2.5311311) q[2];
sx q[2];
rz(-2.9739213) q[2];
rz(-2.3882315) q[3];
sx q[3];
rz(-1.2246917) q[3];
sx q[3];
rz(1.0838449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3566256) q[0];
sx q[0];
rz(-1.4545119) q[0];
sx q[0];
rz(-1.0474569) q[0];
rz(2.6197703) q[1];
sx q[1];
rz(-1.1818886) q[1];
sx q[1];
rz(-2.3375915) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3814731) q[0];
sx q[0];
rz(-1.9895041) q[0];
sx q[0];
rz(0.19180723) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1754403) q[2];
sx q[2];
rz(-2.1198303) q[2];
sx q[2];
rz(2.6556478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2968353) q[1];
sx q[1];
rz(-1.8342819) q[1];
sx q[1];
rz(0.46864639) q[1];
rz(-pi) q[2];
rz(-0.70254962) q[3];
sx q[3];
rz(-1.8773729) q[3];
sx q[3];
rz(-2.2114913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95253402) q[2];
sx q[2];
rz(-1.6510094) q[2];
sx q[2];
rz(2.789759) q[2];
rz(-3.0686038) q[3];
sx q[3];
rz(-1.1806386) q[3];
sx q[3];
rz(-2.4019057) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9214582) q[0];
sx q[0];
rz(-0.84311068) q[0];
sx q[0];
rz(-1.5716918) q[0];
rz(2.4193071) q[1];
sx q[1];
rz(-1.1415569) q[1];
sx q[1];
rz(-1.7438181) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1565052) q[0];
sx q[0];
rz(-0.48195515) q[0];
sx q[0];
rz(-2.5768649) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7977925) q[2];
sx q[2];
rz(-0.8778866) q[2];
sx q[2];
rz(2.5347112) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.367041) q[1];
sx q[1];
rz(-1.4907717) q[1];
sx q[1];
rz(-2.5888279) q[1];
rz(-1.151284) q[3];
sx q[3];
rz(-1.4503917) q[3];
sx q[3];
rz(-1.6199978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.34652823) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(1.2947003) q[2];
rz(1.0558225) q[3];
sx q[3];
rz(-2.081223) q[3];
sx q[3];
rz(2.136371) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9697363) q[0];
sx q[0];
rz(-2.2315401) q[0];
sx q[0];
rz(1.6730614) q[0];
rz(2.3079474) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(-0.8236105) q[2];
sx q[2];
rz(-1.8737741) q[2];
sx q[2];
rz(-0.57536716) q[2];
rz(2.3232719) q[3];
sx q[3];
rz(-2.5723895) q[3];
sx q[3];
rz(1.7583917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
