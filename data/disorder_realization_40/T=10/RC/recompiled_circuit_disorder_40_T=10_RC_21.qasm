OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(-3.0298046) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(0.15375528) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.789334) q[0];
sx q[0];
rz(-2.6129122) q[0];
sx q[0];
rz(-0.60469158) q[0];
rz(-pi) q[1];
rz(-1.329374) q[2];
sx q[2];
rz(-1.4960939) q[2];
sx q[2];
rz(-0.15969294) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0220713) q[1];
sx q[1];
rz(-1.6011366) q[1];
sx q[1];
rz(1.8896709) q[1];
rz(0.33711707) q[3];
sx q[3];
rz(-1.1532591) q[3];
sx q[3];
rz(-1.2085714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6979606) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.704818) q[2];
rz(-2.4076961) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-2.2170128) q[0];
rz(-2.1444767) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(1.3234214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.130587) q[0];
sx q[0];
rz(-0.6479833) q[0];
sx q[0];
rz(-2.381071) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5263444) q[2];
sx q[2];
rz(-1.2215081) q[2];
sx q[2];
rz(-0.312422) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70222774) q[1];
sx q[1];
rz(-0.80184466) q[1];
sx q[1];
rz(2.761809) q[1];
rz(2.7534361) q[3];
sx q[3];
rz(-2.1106488) q[3];
sx q[3];
rz(-2.848958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.20415846) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(2.5126863) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(-1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73308289) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(0.68840233) q[0];
rz(-0.06772659) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(0.53007954) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2358658) q[0];
sx q[0];
rz(-1.6775963) q[0];
sx q[0];
rz(-1.8624767) q[0];
rz(-2.2377551) q[2];
sx q[2];
rz(-2.505216) q[2];
sx q[2];
rz(-2.9449376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0878151) q[1];
sx q[1];
rz(-2.3936845) q[1];
sx q[1];
rz(1.0653711) q[1];
x q[2];
rz(-1.3313053) q[3];
sx q[3];
rz(-1.665984) q[3];
sx q[3];
rz(2.5483607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80785859) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(0.90144908) q[2];
rz(2.3060913) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19105844) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(2.3572671) q[0];
rz(0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(-3.004946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2211321) q[0];
sx q[0];
rz(-1.7014628) q[0];
sx q[0];
rz(2.0544102) q[0];
rz(-pi) q[1];
rz(-1.253445) q[2];
sx q[2];
rz(-1.1270521) q[2];
sx q[2];
rz(-1.9217984) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1130484) q[1];
sx q[1];
rz(-1.5469157) q[1];
sx q[1];
rz(1.6391812) q[1];
rz(-pi) q[2];
rz(0.99478787) q[3];
sx q[3];
rz(-1.6487062) q[3];
sx q[3];
rz(0.60914492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1293929) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(0.56048918) q[2];
rz(-0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43301582) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.7339773) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6038937) q[0];
sx q[0];
rz(-1.1802117) q[0];
sx q[0];
rz(1.9498528) q[0];
rz(-2.0166964) q[2];
sx q[2];
rz(-0.6302399) q[2];
sx q[2];
rz(1.557204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.447532) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(-1.3259757) q[1];
x q[2];
rz(-1.0293343) q[3];
sx q[3];
rz(-2.7280305) q[3];
sx q[3];
rz(3.001861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6901107) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(0.042479854) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-2.5094331) q[3];
sx q[3];
rz(-2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38934389) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(-0.4831627) q[0];
rz(2.0893611) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(0.56484708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0676346) q[0];
sx q[0];
rz(-1.7059776) q[0];
sx q[0];
rz(-1.0413175) q[0];
x q[1];
rz(0.075750307) q[2];
sx q[2];
rz(-1.9905914) q[2];
sx q[2];
rz(2.9043353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55826742) q[1];
sx q[1];
rz(-2.7100483) q[1];
sx q[1];
rz(1.7130909) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5552093) q[3];
sx q[3];
rz(-2.8660503) q[3];
sx q[3];
rz(1.1875718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57006449) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.338039) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(-2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.9682482) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(0.54779732) q[0];
rz(0.785218) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-2.8731667) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82314202) q[0];
sx q[0];
rz(-0.25932352) q[0];
sx q[0];
rz(0.32487049) q[0];
rz(-pi) q[1];
rz(2.3586876) q[2];
sx q[2];
rz(-2.1856538) q[2];
sx q[2];
rz(-2.8564786) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0722326) q[1];
sx q[1];
rz(-1.6688804) q[1];
sx q[1];
rz(2.1588438) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0915756) q[3];
sx q[3];
rz(-0.5939393) q[3];
sx q[3];
rz(-1.493243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(2.7653149) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(-3.0686839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5557264) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-1.0193753) q[0];
rz(0.85340071) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(0.44874915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4067626) q[0];
sx q[0];
rz(-1.8039628) q[0];
sx q[0];
rz(1.083311) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.87254) q[2];
sx q[2];
rz(-1.3240959) q[2];
sx q[2];
rz(1.6517284) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4510348) q[1];
sx q[1];
rz(-2.096855) q[1];
sx q[1];
rz(-1.8655538) q[1];
x q[2];
rz(2.2711146) q[3];
sx q[3];
rz(-0.92910367) q[3];
sx q[3];
rz(-2.2023647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.86924187) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(-2.8273919) q[2];
rz(-0.82434404) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(-0.79469386) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0712414) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.2605793) q[0];
rz(2.4977327) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(-3.1226645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32301329) q[0];
sx q[0];
rz(-1.5886663) q[0];
sx q[0];
rz(-2.8779526) q[0];
rz(-0.55982121) q[2];
sx q[2];
rz(-2.9580742) q[2];
sx q[2];
rz(-0.57932094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85266528) q[1];
sx q[1];
rz(-1.7551433) q[1];
sx q[1];
rz(-1.5194555) q[1];
x q[2];
rz(-1.8065679) q[3];
sx q[3];
rz(-1.2887495) q[3];
sx q[3];
rz(-2.3896133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5921322) q[2];
sx q[2];
rz(-2.7533054) q[2];
sx q[2];
rz(-2.4712759) q[2];
rz(-2.629225) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5877514) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(-0.49945369) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(-2.0589028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854541) q[0];
sx q[0];
rz(-0.98422613) q[0];
sx q[0];
rz(-1.6295208) q[0];
rz(-pi) q[1];
x q[1];
rz(0.033109025) q[2];
sx q[2];
rz(-2.3220255) q[2];
sx q[2];
rz(2.5493252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.078552695) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(0.54363721) q[1];
x q[2];
rz(1.8977676) q[3];
sx q[3];
rz(-1.2083112) q[3];
sx q[3];
rz(-0.31839759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3988951) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(0.51188525) q[2];
rz(-2.7486457) q[3];
sx q[3];
rz(-1.4018551) q[3];
sx q[3];
rz(1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.95505161) q[0];
sx q[0];
rz(-1.825009) q[0];
sx q[0];
rz(0.64074989) q[0];
rz(-2.4004249) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(2.1096061) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(0.74647222) q[3];
sx q[3];
rz(-1.4440047) q[3];
sx q[3];
rz(0.11228893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
