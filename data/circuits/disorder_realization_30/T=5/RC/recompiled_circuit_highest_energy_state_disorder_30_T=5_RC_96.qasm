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
rz(1.5626524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.232036) q[0];
sx q[0];
rz(-2.252451) q[0];
sx q[0];
rz(0.85083346) q[0];
rz(-pi) q[1];
rz(0.15450124) q[2];
sx q[2];
rz(-1.3065632) q[2];
sx q[2];
rz(-2.7619948) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7464802) q[1];
sx q[1];
rz(-1.0653138) q[1];
sx q[1];
rz(-3.014747) q[1];
rz(-pi) q[2];
rz(-0.64227958) q[3];
sx q[3];
rz(-3.0845957) q[3];
sx q[3];
rz(1.1617203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19848862) q[2];
sx q[2];
rz(-0.21185943) q[2];
sx q[2];
rz(-2.1073821) q[2];
rz(2.4487623) q[3];
sx q[3];
rz(-2.0878017) q[3];
sx q[3];
rz(-2.0809765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8387872) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(2.5651108) q[0];
rz(0.020596404) q[1];
sx q[1];
rz(-0.3978022) q[1];
sx q[1];
rz(-2.1131262) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45555025) q[0];
sx q[0];
rz(-0.33838135) q[0];
sx q[0];
rz(-2.4421921) q[0];
rz(-pi) q[1];
rz(-0.307085) q[2];
sx q[2];
rz(-1.9283617) q[2];
sx q[2];
rz(-1.7890695) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2814863) q[1];
sx q[1];
rz(-2.6074643) q[1];
sx q[1];
rz(2.2661282) q[1];
x q[2];
rz(1.1787492) q[3];
sx q[3];
rz(-0.96091753) q[3];
sx q[3];
rz(-1.4489685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91745201) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(1.8769598) q[2];
rz(-0.97359109) q[3];
sx q[3];
rz(-1.4507989) q[3];
sx q[3];
rz(-0.26315954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50315404) q[0];
sx q[0];
rz(-1.8373024) q[0];
sx q[0];
rz(-0.85357443) q[0];
rz(2.0897934) q[1];
sx q[1];
rz(-2.167326) q[1];
sx q[1];
rz(3.1265756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1117843) q[0];
sx q[0];
rz(-1.5795603) q[0];
sx q[0];
rz(2.2814204) q[0];
rz(-3.0454841) q[2];
sx q[2];
rz(-1.5920361) q[2];
sx q[2];
rz(-0.48890314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2440363) q[1];
sx q[1];
rz(-1.4970395) q[1];
sx q[1];
rz(-0.71457992) q[1];
x q[2];
rz(2.851753) q[3];
sx q[3];
rz(-1.6320758) q[3];
sx q[3];
rz(-1.0672399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1198279) q[2];
sx q[2];
rz(-1.652635) q[2];
sx q[2];
rz(0.28142288) q[2];
rz(-1.6875387) q[3];
sx q[3];
rz(-1.2097996) q[3];
sx q[3];
rz(-0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467095) q[0];
sx q[0];
rz(-1.895772) q[0];
sx q[0];
rz(-2.967714) q[0];
rz(1.8100544) q[1];
sx q[1];
rz(-1.7200108) q[1];
sx q[1];
rz(-1.7132267) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15584942) q[0];
sx q[0];
rz(-1.5850889) q[0];
sx q[0];
rz(-1.3784268) q[0];
rz(2.7816781) q[2];
sx q[2];
rz(-1.4524621) q[2];
sx q[2];
rz(-1.9595166) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5927939) q[1];
sx q[1];
rz(-2.4830635) q[1];
sx q[1];
rz(-0.062915398) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7739206) q[3];
sx q[3];
rz(-1.9353364) q[3];
sx q[3];
rz(-3.1057949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7017158) q[2];
sx q[2];
rz(-0.82650799) q[2];
sx q[2];
rz(2.7244549) q[2];
rz(0.16962984) q[3];
sx q[3];
rz(-0.093234213) q[3];
sx q[3];
rz(-0.72181845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7222662) q[0];
sx q[0];
rz(-0.37501431) q[0];
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
rz(1.5681388) q[0];
sx q[0];
rz(-0.72069695) q[0];
sx q[0];
rz(-1.9573523) q[0];
rz(2.4028087) q[2];
sx q[2];
rz(-2.1257943) q[2];
sx q[2];
rz(-0.273663) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1479234) q[1];
sx q[1];
rz(-1.5412314) q[1];
sx q[1];
rz(-1.7600585) q[1];
x q[2];
rz(1.7047911) q[3];
sx q[3];
rz(-1.6686642) q[3];
sx q[3];
rz(2.0201473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3994483) q[2];
sx q[2];
rz(-0.98230201) q[2];
sx q[2];
rz(-2.1461416) q[2];
rz(-0.71850145) q[3];
sx q[3];
rz(-1.3281497) q[3];
sx q[3];
rz(-0.79513454) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8566078) q[0];
sx q[0];
rz(-1.4549078) q[0];
sx q[0];
rz(-1.6108151) q[0];
rz(1.4467622) q[1];
sx q[1];
rz(-1.6981354) q[1];
sx q[1];
rz(1.4379427) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3525248) q[0];
sx q[0];
rz(-1.1630327) q[0];
sx q[0];
rz(-0.022313281) q[0];
rz(-pi) q[1];
rz(2.6126768) q[2];
sx q[2];
rz(-0.66096604) q[2];
sx q[2];
rz(-1.5724374) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48541486) q[1];
sx q[1];
rz(-2.6256769) q[1];
sx q[1];
rz(0.1702711) q[1];
rz(3.1385057) q[3];
sx q[3];
rz(-0.75957889) q[3];
sx q[3];
rz(0.39955968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1630254) q[2];
sx q[2];
rz(-1.359553) q[2];
sx q[2];
rz(-1.9095518) q[2];
rz(-1.1466522) q[3];
sx q[3];
rz(-1.3403284) q[3];
sx q[3];
rz(-0.31057096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6935007) q[0];
sx q[0];
rz(-0.25059487) q[0];
sx q[0];
rz(1.084024) q[0];
rz(-pi) q[1];
rz(0.6426446) q[2];
sx q[2];
rz(-1.3054928) q[2];
sx q[2];
rz(2.93769) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6086238) q[1];
sx q[1];
rz(-0.86809671) q[1];
sx q[1];
rz(-2.7547902) q[1];
rz(-0.45400158) q[3];
sx q[3];
rz(-0.89278614) q[3];
sx q[3];
rz(-2.0583722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19411479) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(-1.4878368) q[3];
sx q[3];
rz(-1.0349118) q[3];
sx q[3];
rz(1.8208767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3461935) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(1.7561703) q[0];
rz(1.2162544) q[1];
sx q[1];
rz(-1.6938035) q[1];
sx q[1];
rz(2.6522955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84626694) q[0];
sx q[0];
rz(-1.5510484) q[0];
sx q[0];
rz(-2.4700145) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4665746) q[2];
sx q[2];
rz(-0.33793338) q[2];
sx q[2];
rz(0.80160415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8000187) q[1];
sx q[1];
rz(-1.6650781) q[1];
sx q[1];
rz(-1.0977919) q[1];
x q[2];
rz(1.1888483) q[3];
sx q[3];
rz(-0.26774613) q[3];
sx q[3];
rz(-1.4827888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.74799246) q[2];
sx q[2];
rz(-1.9990498) q[2];
sx q[2];
rz(2.6753814) q[2];
rz(-1.6097869) q[3];
sx q[3];
rz(-1.5038871) q[3];
sx q[3];
rz(-2.8209749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4894678) q[0];
sx q[0];
rz(-2.2418699) q[0];
sx q[0];
rz(0.049064431) q[0];
rz(0.24066726) q[1];
sx q[1];
rz(-2.1235695) q[1];
sx q[1];
rz(3.1191471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81240679) q[0];
sx q[0];
rz(-0.57984771) q[0];
sx q[0];
rz(-1.5270698) q[0];
rz(-pi) q[1];
rz(-1.2227549) q[2];
sx q[2];
rz(-1.204426) q[2];
sx q[2];
rz(-2.6719246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53246236) q[1];
sx q[1];
rz(-2.3098619) q[1];
sx q[1];
rz(2.4799281) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1184228) q[3];
sx q[3];
rz(-2.2045772) q[3];
sx q[3];
rz(2.3364802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13751328) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(-2.1732886) q[2];
rz(2.3051895) q[3];
sx q[3];
rz(-0.29763779) q[3];
sx q[3];
rz(1.3692726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7264929) q[0];
sx q[0];
rz(-1.8510171) q[0];
sx q[0];
rz(-2.7628164) q[0];
rz(3.1355766) q[1];
sx q[1];
rz(-2.6117987) q[1];
sx q[1];
rz(-0.63374296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0519052) q[0];
sx q[0];
rz(-1.9080592) q[0];
sx q[0];
rz(2.428672) q[0];
x q[1];
rz(0.1600645) q[2];
sx q[2];
rz(-1.0032005) q[2];
sx q[2];
rz(0.71086037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6478471) q[1];
sx q[1];
rz(-0.71628191) q[1];
sx q[1];
rz(-0.78328461) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87234906) q[3];
sx q[3];
rz(-2.3309338) q[3];
sx q[3];
rz(1.8332421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70574957) q[2];
sx q[2];
rz(-1.5219995) q[2];
sx q[2];
rz(-0.61326927) q[2];
rz(2.4663726) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(0.81775445) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9153862) q[0];
sx q[0];
rz(-1.7788667) q[0];
sx q[0];
rz(-1.9052196) q[0];
rz(-0.18648237) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(0.45966799) q[2];
sx q[2];
rz(-1.5887693) q[2];
sx q[2];
rz(2.1732687) q[2];
rz(-2.1929838) q[3];
sx q[3];
rz(-1.6453679) q[3];
sx q[3];
rz(-2.9352376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
