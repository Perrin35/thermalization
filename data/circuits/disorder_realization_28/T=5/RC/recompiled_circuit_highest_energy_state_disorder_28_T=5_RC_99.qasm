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
rz(-2.2313843) q[0];
sx q[0];
rz(3.892133) q[0];
sx q[0];
rz(12.005796) q[0];
rz(1.9380467) q[1];
sx q[1];
rz(0.37625852) q[1];
sx q[1];
rz(9.6179554) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7734186) q[0];
sx q[0];
rz(-0.76245284) q[0];
sx q[0];
rz(0.66526757) q[0];
x q[1];
rz(1.7142263) q[2];
sx q[2];
rz(-2.4312115) q[2];
sx q[2];
rz(1.8153035) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79598017) q[1];
sx q[1];
rz(-2.4761845) q[1];
sx q[1];
rz(1.7980952) q[1];
rz(-pi) q[2];
x q[2];
rz(2.812593) q[3];
sx q[3];
rz(-1.7339249) q[3];
sx q[3];
rz(-1.8150869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91118497) q[2];
sx q[2];
rz(-0.62778968) q[2];
sx q[2];
rz(2.5837303) q[2];
rz(1.2281536) q[3];
sx q[3];
rz(-1.4598673) q[3];
sx q[3];
rz(-0.094999878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2896344) q[0];
sx q[0];
rz(-1.8486706) q[0];
sx q[0];
rz(1.5767545) q[0];
rz(1.9948888) q[1];
sx q[1];
rz(-1.4253989) q[1];
sx q[1];
rz(2.1099405) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0165591) q[0];
sx q[0];
rz(-1.8047338) q[0];
sx q[0];
rz(2.1159322) q[0];
rz(-pi) q[1];
rz(2.7286058) q[2];
sx q[2];
rz(-1.3527762) q[2];
sx q[2];
rz(0.08139164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7647142) q[1];
sx q[1];
rz(-0.8974613) q[1];
sx q[1];
rz(-2.7002579) q[1];
rz(-pi) q[2];
rz(-2.7835763) q[3];
sx q[3];
rz(-2.6081134) q[3];
sx q[3];
rz(-1.8915576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8138294) q[2];
sx q[2];
rz(-0.29080614) q[2];
sx q[2];
rz(1.7247058) q[2];
rz(-0.82301569) q[3];
sx q[3];
rz(-1.6133285) q[3];
sx q[3];
rz(-0.64457646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7773892) q[0];
sx q[0];
rz(-2.0222029) q[0];
sx q[0];
rz(0.35225824) q[0];
rz(0.38905713) q[1];
sx q[1];
rz(-1.16951) q[1];
sx q[1];
rz(-2.9152117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8181747) q[0];
sx q[0];
rz(-1.3238088) q[0];
sx q[0];
rz(3.0564708) q[0];
rz(-0.34778123) q[2];
sx q[2];
rz(-2.116733) q[2];
sx q[2];
rz(0.40218654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6701874) q[1];
sx q[1];
rz(-0.67882631) q[1];
sx q[1];
rz(-0.77969867) q[1];
x q[2];
rz(-2.705615) q[3];
sx q[3];
rz(-1.6395901) q[3];
sx q[3];
rz(1.9668369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9809197) q[2];
sx q[2];
rz(-1.0272762) q[2];
sx q[2];
rz(-2.7799907) q[2];
rz(-0.30038878) q[3];
sx q[3];
rz(-1.8301679) q[3];
sx q[3];
rz(2.1786407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.921973) q[0];
sx q[0];
rz(-1.0924871) q[0];
sx q[0];
rz(-0.12119448) q[0];
rz(-1.9425862) q[1];
sx q[1];
rz(-2.0620748) q[1];
sx q[1];
rz(1.2746864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.072678) q[0];
sx q[0];
rz(-1.6175748) q[0];
sx q[0];
rz(-0.22225265) q[0];
x q[1];
rz(-1.5546459) q[2];
sx q[2];
rz(-1.3653737) q[2];
sx q[2];
rz(-1.1753163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1772985) q[1];
sx q[1];
rz(-1.7907447) q[1];
sx q[1];
rz(2.5581058) q[1];
x q[2];
rz(-0.14583662) q[3];
sx q[3];
rz(-1.2088299) q[3];
sx q[3];
rz(1.2545619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.092209665) q[2];
sx q[2];
rz(-1.4468687) q[2];
sx q[2];
rz(1.764074) q[2];
rz(-1.1285909) q[3];
sx q[3];
rz(-0.94213525) q[3];
sx q[3];
rz(0.82730627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89151299) q[0];
sx q[0];
rz(-0.87123195) q[0];
sx q[0];
rz(-2.0569892) q[0];
rz(2.4705823) q[1];
sx q[1];
rz(-1.6295461) q[1];
sx q[1];
rz(-3.0674518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0468354) q[0];
sx q[0];
rz(-1.0124542) q[0];
sx q[0];
rz(1.2461016) q[0];
x q[1];
rz(0.29080342) q[2];
sx q[2];
rz(-1.8755696) q[2];
sx q[2];
rz(-2.1220611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74846327) q[1];
sx q[1];
rz(-0.65915758) q[1];
sx q[1];
rz(-1.1200957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1187333) q[3];
sx q[3];
rz(-0.30395711) q[3];
sx q[3];
rz(0.077465103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.025617754) q[2];
sx q[2];
rz(-2.3927549) q[2];
sx q[2];
rz(-0.96308127) q[2];
rz(1.3747619) q[3];
sx q[3];
rz(-1.129496) q[3];
sx q[3];
rz(0.53205427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.9328203) q[0];
sx q[0];
rz(-2.8498579) q[0];
sx q[0];
rz(1.2579086) q[0];
rz(2.3014297) q[1];
sx q[1];
rz(-1.9636619) q[1];
sx q[1];
rz(2.449583) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35273472) q[0];
sx q[0];
rz(-1.8957355) q[0];
sx q[0];
rz(-2.8558225) q[0];
rz(-pi) q[1];
rz(-2.3756723) q[2];
sx q[2];
rz(-1.4408852) q[2];
sx q[2];
rz(-2.9875362) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4547448) q[1];
sx q[1];
rz(-2.1917289) q[1];
sx q[1];
rz(-2.5140155) q[1];
rz(-pi) q[2];
rz(-0.70051381) q[3];
sx q[3];
rz(-2.3729602) q[3];
sx q[3];
rz(-0.045846102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17795263) q[2];
sx q[2];
rz(-1.0963564) q[2];
sx q[2];
rz(0.5160416) q[2];
rz(-1.3456723) q[3];
sx q[3];
rz(-0.44107744) q[3];
sx q[3];
rz(-1.0015944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26821414) q[0];
sx q[0];
rz(-0.28609797) q[0];
sx q[0];
rz(-1.0299261) q[0];
rz(0.65451199) q[1];
sx q[1];
rz(-1.3348568) q[1];
sx q[1];
rz(0.17394224) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2692341) q[0];
sx q[0];
rz(-2.1790811) q[0];
sx q[0];
rz(-2.48003) q[0];
rz(-2.6268466) q[2];
sx q[2];
rz(-0.87640773) q[2];
sx q[2];
rz(2.1026762) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6817045) q[1];
sx q[1];
rz(-0.68716421) q[1];
sx q[1];
rz(-2.0089441) q[1];
rz(-pi) q[2];
rz(0.44486041) q[3];
sx q[3];
rz(-2.5813817) q[3];
sx q[3];
rz(-2.1261393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30764636) q[2];
sx q[2];
rz(-0.53457326) q[2];
sx q[2];
rz(1.8024811) q[2];
rz(-0.80859679) q[3];
sx q[3];
rz(-1.1491038) q[3];
sx q[3];
rz(1.2861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22063743) q[0];
sx q[0];
rz(-1.0458825) q[0];
sx q[0];
rz(1.9761696) q[0];
rz(2.9479345) q[1];
sx q[1];
rz(-2.3083189) q[1];
sx q[1];
rz(1.9906893) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347558) q[0];
sx q[0];
rz(-2.2118705) q[0];
sx q[0];
rz(1.8888372) q[0];
x q[1];
rz(-1.3940349) q[2];
sx q[2];
rz(-0.77406787) q[2];
sx q[2];
rz(-1.34684) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6713359) q[1];
sx q[1];
rz(-1.0733979) q[1];
sx q[1];
rz(-0.48723826) q[1];
rz(-0.64784785) q[3];
sx q[3];
rz(-1.0937249) q[3];
sx q[3];
rz(-1.939524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2151486) q[2];
sx q[2];
rz(-1.0215267) q[2];
sx q[2];
rz(1.273217) q[2];
rz(0.5145973) q[3];
sx q[3];
rz(-2.8223346) q[3];
sx q[3];
rz(-2.4014421) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1505245) q[0];
sx q[0];
rz(-3.0758698) q[0];
sx q[0];
rz(-2.7142628) q[0];
rz(-1.4409675) q[1];
sx q[1];
rz(-1.7040375) q[1];
sx q[1];
rz(2.2429121) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3220427) q[0];
sx q[0];
rz(-1.8801831) q[0];
sx q[0];
rz(-1.1197107) q[0];
rz(-2.1601474) q[2];
sx q[2];
rz(-1.169983) q[2];
sx q[2];
rz(-0.0055731853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8391621) q[1];
sx q[1];
rz(-1.3962455) q[1];
sx q[1];
rz(-1.8086834) q[1];
x q[2];
rz(2.6014448) q[3];
sx q[3];
rz(-1.9753878) q[3];
sx q[3];
rz(0.57309421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.906189) q[2];
sx q[2];
rz(-2.0960505) q[2];
sx q[2];
rz(-1.40353) q[2];
rz(-0.67052001) q[3];
sx q[3];
rz(-1.5454005) q[3];
sx q[3];
rz(-1.8286573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32187605) q[0];
sx q[0];
rz(-2.1724367) q[0];
sx q[0];
rz(-2.4153391) q[0];
rz(-2.4655474) q[1];
sx q[1];
rz(-1.36422) q[1];
sx q[1];
rz(0.77686754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60706106) q[0];
sx q[0];
rz(-2.2103469) q[0];
sx q[0];
rz(3.0846918) q[0];
rz(1.6157584) q[2];
sx q[2];
rz(-1.7167712) q[2];
sx q[2];
rz(3.1342497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.59275) q[1];
sx q[1];
rz(-1.2075066) q[1];
sx q[1];
rz(-0.32925683) q[1];
x q[2];
rz(-0.75001287) q[3];
sx q[3];
rz(-0.74392825) q[3];
sx q[3];
rz(-0.14694005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0016510222) q[2];
sx q[2];
rz(-0.17639128) q[2];
sx q[2];
rz(1.6528992) q[2];
rz(2.0148924) q[3];
sx q[3];
rz(-1.1688346) q[3];
sx q[3];
rz(0.53762976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61459944) q[0];
sx q[0];
rz(-2.499883) q[0];
sx q[0];
rz(-1.5747621) q[0];
rz(2.3552786) q[1];
sx q[1];
rz(-0.17335261) q[1];
sx q[1];
rz(-0.39549624) q[1];
rz(-0.69063557) q[2];
sx q[2];
rz(-0.80123676) q[2];
sx q[2];
rz(-1.5581808) q[2];
rz(1.3287104) q[3];
sx q[3];
rz(-0.81012204) q[3];
sx q[3];
rz(1.3835081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
