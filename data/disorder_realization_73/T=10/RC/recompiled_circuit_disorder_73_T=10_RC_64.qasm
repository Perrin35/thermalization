OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(3.6448195) q[0];
sx q[0];
rz(10.148944) q[0];
rz(6.9231482) q[1];
sx q[1];
rz(5.7531113) q[1];
sx q[1];
rz(2.35676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0137579) q[0];
sx q[0];
rz(-0.36359596) q[0];
sx q[0];
rz(-2.512393) q[0];
rz(-0.22462331) q[2];
sx q[2];
rz(-2.7135239) q[2];
sx q[2];
rz(-3.012804) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9468294) q[1];
sx q[1];
rz(-2.0954872) q[1];
sx q[1];
rz(0.25804934) q[1];
rz(1.1562528) q[3];
sx q[3];
rz(-1.6016377) q[3];
sx q[3];
rz(2.8950092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(-3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(-1.7547866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215866) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(-0.24366972) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.7858645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60123721) q[0];
sx q[0];
rz(-1.8260801) q[0];
sx q[0];
rz(-3.113494) q[0];
rz(1.9594876) q[2];
sx q[2];
rz(-1.2388065) q[2];
sx q[2];
rz(1.4093083) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6527378) q[1];
sx q[1];
rz(-1.4421717) q[1];
sx q[1];
rz(2.0205523) q[1];
rz(1.5242819) q[3];
sx q[3];
rz(-2.1624613) q[3];
sx q[3];
rz(0.072629645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0791066) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(-0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24519414) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(2.2431592) q[0];
rz(1.8067182) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(-1.2737087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1184517) q[0];
sx q[0];
rz(-1.8225192) q[0];
sx q[0];
rz(-1.203042) q[0];
x q[1];
rz(1.4726228) q[2];
sx q[2];
rz(-1.2766826) q[2];
sx q[2];
rz(1.3054747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85325235) q[1];
sx q[1];
rz(-2.3715092) q[1];
sx q[1];
rz(0.5132765) q[1];
rz(-0.7234296) q[3];
sx q[3];
rz(-1.1838786) q[3];
sx q[3];
rz(-2.5748411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0597824) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(-2.188142) q[2];
rz(0.034514286) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(-2.8934073) q[0];
rz(2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(0.074137069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.603133) q[0];
sx q[0];
rz(-0.54766253) q[0];
sx q[0];
rz(2.0752226) q[0];
rz(-pi) q[1];
rz(-0.89216994) q[2];
sx q[2];
rz(-1.2954419) q[2];
sx q[2];
rz(-2.8868669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84872765) q[1];
sx q[1];
rz(-1.25602) q[1];
sx q[1];
rz(0.12401144) q[1];
rz(-pi) q[2];
rz(-0.24012633) q[3];
sx q[3];
rz(-1.4014981) q[3];
sx q[3];
rz(1.5907767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-0.40428287) q[2];
sx q[2];
rz(-3.1029491) q[2];
rz(0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(2.3249987) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(-1.978925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8602596) q[0];
sx q[0];
rz(-1.5409924) q[0];
sx q[0];
rz(1.6822862) q[0];
x q[1];
rz(1.0983724) q[2];
sx q[2];
rz(-0.35944164) q[2];
sx q[2];
rz(3.0430832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3711277) q[1];
sx q[1];
rz(-2.7264997) q[1];
sx q[1];
rz(2.0626555) q[1];
rz(0.43443067) q[3];
sx q[3];
rz(-2.0475004) q[3];
sx q[3];
rz(-1.3694976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5148619) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(0.90406117) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(-2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11480039) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.6739155) q[0];
rz(-2.5698075) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(2.8335559) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.98826) q[0];
sx q[0];
rz(-1.6969661) q[0];
sx q[0];
rz(1.4754962) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61200895) q[2];
sx q[2];
rz(-2.0888121) q[2];
sx q[2];
rz(-0.81105622) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0526272) q[1];
sx q[1];
rz(-1.7672156) q[1];
sx q[1];
rz(2.2599225) q[1];
x q[2];
rz(-1.4792535) q[3];
sx q[3];
rz(-1.9178101) q[3];
sx q[3];
rz(-2.7029944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77928153) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(-0.71427304) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66184735) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(-0.1299783) q[0];
rz(-3.1107483) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(2.470509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2790047) q[0];
sx q[0];
rz(-3.0775078) q[0];
sx q[0];
rz(-0.57871731) q[0];
x q[1];
rz(-1.7332156) q[2];
sx q[2];
rz(-2.9549837) q[2];
sx q[2];
rz(2.9090372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1783501) q[1];
sx q[1];
rz(-1.9415783) q[1];
sx q[1];
rz(3.0313655) q[1];
rz(-pi) q[2];
rz(3.0942261) q[3];
sx q[3];
rz(-0.8960552) q[3];
sx q[3];
rz(2.588152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.122763) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(-2.5637131) q[2];
rz(3.1130062) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776483) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(-2.7291765) q[0];
rz(-1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.9746045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8186504) q[0];
sx q[0];
rz(-2.4837821) q[0];
sx q[0];
rz(1.9239182) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86161676) q[2];
sx q[2];
rz(-1.78252) q[2];
sx q[2];
rz(-0.35702969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.036901722) q[1];
sx q[1];
rz(-2.9505886) q[1];
sx q[1];
rz(2.9810993) q[1];
rz(-pi) q[2];
rz(-0.1498296) q[3];
sx q[3];
rz(-2.095788) q[3];
sx q[3];
rz(-2.6168952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-0.18903014) q[2];
rz(-0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-2.2145859) q[0];
rz(1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(-1.4896726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72367523) q[0];
sx q[0];
rz(-2.0360887) q[0];
sx q[0];
rz(-2.377541) q[0];
x q[1];
rz(0.050612014) q[2];
sx q[2];
rz(-1.0992556) q[2];
sx q[2];
rz(2.5052349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3085732) q[1];
sx q[1];
rz(-0.25990572) q[1];
sx q[1];
rz(-2.412917) q[1];
rz(-pi) q[2];
rz(-2.1065815) q[3];
sx q[3];
rz(-1.1318558) q[3];
sx q[3];
rz(2.5715695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5902517) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(1.7112973) q[2];
rz(-2.5643505) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(2.9945701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64919103) q[0];
sx q[0];
rz(-0.99897879) q[0];
sx q[0];
rz(0.4581106) q[0];
x q[1];
rz(-1.5006127) q[2];
sx q[2];
rz(-1.3314221) q[2];
sx q[2];
rz(-1.0603051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8204931) q[1];
sx q[1];
rz(-1.610678) q[1];
sx q[1];
rz(0.95221968) q[1];
rz(-0.57659984) q[3];
sx q[3];
rz(-1.4398265) q[3];
sx q[3];
rz(3.1104345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0599351) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(0.74404136) q[2];
rz(0.75731164) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(1.5031917) q[2];
sx q[2];
rz(-2.0220145) q[2];
sx q[2];
rz(1.8215712) q[2];
rz(-0.13027262) q[3];
sx q[3];
rz(-1.0126922) q[3];
sx q[3];
rz(-1.0425413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
