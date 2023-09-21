OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(-2.764954) q[0];
sx q[0];
rz(-0.11178804) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(0.15375528) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68054799) q[0];
sx q[0];
rz(-1.8616315) q[0];
sx q[0];
rz(0.44797795) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2674238) q[2];
sx q[2];
rz(-0.25250013) q[2];
sx q[2];
rz(1.7054103) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54290463) q[1];
sx q[1];
rz(-0.3202657) q[1];
sx q[1];
rz(-1.4742875) q[1];
rz(-2.0102242) q[3];
sx q[3];
rz(-1.2636375) q[3];
sx q[3];
rz(0.22104056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.443632) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(-1.4367746) q[2];
rz(-2.4076961) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(-2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2519418) q[0];
sx q[0];
rz(-1.9152859) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(-0.997116) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(-1.8181713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2572718) q[0];
sx q[0];
rz(-2.0233676) q[0];
sx q[0];
rz(1.0898468) q[0];
rz(-pi) q[1];
rz(-1.6152482) q[2];
sx q[2];
rz(-1.9200846) q[2];
sx q[2];
rz(2.8291707) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1393226) q[1];
sx q[1];
rz(-1.3011258) q[1];
sx q[1];
rz(-2.3766999) q[1];
rz(2.1453342) q[3];
sx q[3];
rz(-1.2401476) q[3];
sx q[3];
rz(1.4853256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9374342) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(2.3834174) q[2];
rz(0.6289064) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(-1.988407) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73308289) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(2.4531903) q[0];
rz(3.0738661) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-2.6115131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2358658) q[0];
sx q[0];
rz(-1.6775963) q[0];
sx q[0];
rz(1.8624767) q[0];
rz(-1.0447787) q[2];
sx q[2];
rz(-1.194343) q[2];
sx q[2];
rz(-1.9386171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59280076) q[1];
sx q[1];
rz(-2.2081516) q[1];
sx q[1];
rz(0.42216502) q[1];
rz(-pi) q[2];
rz(0.097966627) q[3];
sx q[3];
rz(-1.3324105) q[3];
sx q[3];
rz(1.0007678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3337341) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(2.3060913) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(1.3114595) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(2.3572671) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(-0.13664666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9204606) q[0];
sx q[0];
rz(-1.7014628) q[0];
sx q[0];
rz(-2.0544102) q[0];
x q[1];
rz(1.253445) q[2];
sx q[2];
rz(-1.1270521) q[2];
sx q[2];
rz(-1.2197942) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5977051) q[1];
sx q[1];
rz(-1.6391616) q[1];
sx q[1];
rz(3.1176561) q[1];
rz(-pi) q[2];
rz(0.09282077) q[3];
sx q[3];
rz(-2.1448359) q[3];
sx q[3];
rz(2.1294347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(2.5811035) q[2];
rz(3.1292606) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43301582) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-2.3204455) q[0];
rz(-0.87617809) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.4076153) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6038937) q[0];
sx q[0];
rz(-1.1802117) q[0];
sx q[0];
rz(1.9498528) q[0];
x q[1];
rz(-2.0166964) q[2];
sx q[2];
rz(-0.6302399) q[2];
sx q[2];
rz(1.557204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40286139) q[1];
sx q[1];
rz(-1.0051454) q[1];
sx q[1];
rz(-0.15927844) q[1];
x q[2];
rz(1.9305265) q[3];
sx q[3];
rz(-1.779428) q[3];
sx q[3];
rz(1.9344575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6901107) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(-3.0991128) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7522488) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(-0.4831627) q[0];
rz(1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(2.5767456) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42445499) q[0];
sx q[0];
rz(-2.0949445) q[0];
sx q[0];
rz(-2.9852887) q[0];
x q[1];
rz(-1.9916612) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(1.7771306) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1419303) q[1];
sx q[1];
rz(-1.6301486) q[1];
sx q[1];
rz(1.998494) q[1];
x q[2];
rz(0.23128831) q[3];
sx q[3];
rz(-1.4196718) q[3];
sx q[3];
rz(2.1895727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.57006449) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.8035536) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(-0.6750955) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9682482) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(2.5937953) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(2.8731667) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82314202) q[0];
sx q[0];
rz(-2.8822691) q[0];
sx q[0];
rz(-0.32487049) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78290501) q[2];
sx q[2];
rz(-0.9559388) q[2];
sx q[2];
rz(-0.2851141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7859898) q[1];
sx q[1];
rz(-0.59521788) q[1];
sx q[1];
rz(-1.3952414) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3018474) q[3];
sx q[3];
rz(-1.0511304) q[3];
sx q[3];
rz(0.93320751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5238374) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(-2.3274373) q[2];
rz(0.37627775) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(3.0686839) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58586621) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-1.0193753) q[0];
rz(-0.85340071) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(-0.44874915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73483) q[0];
sx q[0];
rz(-1.8039628) q[0];
sx q[0];
rz(-1.083311) q[0];
rz(-0.86782311) q[2];
sx q[2];
rz(-0.38735577) q[2];
sx q[2];
rz(-2.5572436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6905578) q[1];
sx q[1];
rz(-2.096855) q[1];
sx q[1];
rz(1.8655538) q[1];
rz(-pi) q[2];
rz(-2.2711146) q[3];
sx q[3];
rz(-0.92910367) q[3];
sx q[3];
rz(2.2023647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.86924187) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(0.82434404) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(2.3468988) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.8810133) q[0];
rz(-0.64385995) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(3.1226645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2429598) q[0];
sx q[0];
rz(-1.3071994) q[0];
sx q[0];
rz(1.5893057) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15599613) q[2];
sx q[2];
rz(-1.6678572) q[2];
sx q[2];
rz(2.7023466) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85266528) q[1];
sx q[1];
rz(-1.7551433) q[1];
sx q[1];
rz(-1.6221371) q[1];
rz(-0.28963611) q[3];
sx q[3];
rz(-1.3445065) q[3];
sx q[3];
rz(-2.3895404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54946047) q[2];
sx q[2];
rz(-2.7533054) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(2.629225) q[3];
sx q[3];
rz(-1.3918326) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5877514) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-0.49945369) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(-1.0826899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64718819) q[0];
sx q[0];
rz(-1.6196961) q[0];
sx q[0];
rz(-2.5542269) q[0];
rz(1.6062276) q[2];
sx q[2];
rz(-2.3897768) q[2];
sx q[2];
rz(0.54377901) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4808828) q[1];
sx q[1];
rz(-1.0272659) q[1];
sx q[1];
rz(1.5488312) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38090221) q[3];
sx q[3];
rz(-1.2657832) q[3];
sx q[3];
rz(2.0088793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(0.51188525) q[2];
rz(2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.186541) q[0];
sx q[0];
rz(-1.825009) q[0];
sx q[0];
rz(0.64074989) q[0];
rz(-2.4004249) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(2.1869833) q[2];
sx q[2];
rz(-1.90345) q[2];
sx q[2];
rz(-0.30751139) q[2];
rz(-1.742733) q[3];
sx q[3];
rz(-0.83172432) q[3];
sx q[3];
rz(-1.5749501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];