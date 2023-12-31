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
rz(0.11178804) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(0.15375528) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.114404) q[0];
sx q[0];
rz(-1.9986885) q[0];
sx q[0];
rz(-1.2501636) q[0];
rz(-1.8122187) q[2];
sx q[2];
rz(-1.6454988) q[2];
sx q[2];
rz(-0.15969294) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54290463) q[1];
sx q[1];
rz(-0.3202657) q[1];
sx q[1];
rz(-1.6673052) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0102242) q[3];
sx q[3];
rz(-1.8779552) q[3];
sx q[3];
rz(2.9205521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6979606) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(-1.4367746) q[2];
rz(2.4076961) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-2.2170128) q[0];
rz(2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(-1.8181713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01100563) q[0];
sx q[0];
rz(-2.4936094) q[0];
sx q[0];
rz(2.381071) q[0];
rz(-0.3496062) q[2];
sx q[2];
rz(-1.6125624) q[2];
sx q[2];
rz(-1.2431527) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9604608) q[1];
sx q[1];
rz(-2.301553) q[1];
sx q[1];
rz(-1.9366656) q[1];
rz(2.1453342) q[3];
sx q[3];
rz(-1.901445) q[3];
sx q[3];
rz(1.656267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(0.75817529) q[2];
rz(-0.6289064) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(-1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.4085098) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(-0.68840233) q[0];
rz(0.06772659) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-0.53007954) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63307525) q[0];
sx q[0];
rz(-1.2808262) q[0];
sx q[0];
rz(-0.11147186) q[0];
rz(-pi) q[1];
rz(2.0968139) q[2];
sx q[2];
rz(-1.194343) q[2];
sx q[2];
rz(-1.9386171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5487919) q[1];
sx q[1];
rz(-0.93344102) q[1];
sx q[1];
rz(-2.7194276) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9534555) q[3];
sx q[3];
rz(-0.2573765) q[3];
sx q[3];
rz(2.5352258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3337341) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(2.2401436) q[2];
rz(0.83550134) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(-1.3114595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(-0.78432551) q[0];
rz(0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(-3.004946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2480859) q[0];
sx q[0];
rz(-0.49960217) q[0];
sx q[0];
rz(1.2953555) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8881477) q[2];
sx q[2];
rz(-1.1270521) q[2];
sx q[2];
rz(-1.2197942) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1130484) q[1];
sx q[1];
rz(-1.5946769) q[1];
sx q[1];
rz(1.5024115) q[1];
rz(-1.4284381) q[3];
sx q[3];
rz(-0.58066237) q[3];
sx q[3];
rz(-0.84238392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1293929) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(-2.5811035) q[2];
rz(-0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-2.0509317) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43301582) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(-2.2654146) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.7339773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1836023) q[0];
sx q[0];
rz(-1.9200268) q[0];
sx q[0];
rz(0.41718418) q[0];
rz(-pi) q[1];
rz(-2.1528835) q[2];
sx q[2];
rz(-1.313813) q[2];
sx q[2];
rz(-2.7866521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0820497) q[1];
sx q[1];
rz(-1.4364916) q[1];
sx q[1];
rz(2.1422269) q[1];
x q[2];
rz(2.9191454) q[3];
sx q[3];
rz(-1.9223833) q[3];
sx q[3];
rz(0.44140154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6901107) q[2];
sx q[2];
rz(-1.9268945) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38934389) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(0.4831627) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(0.56484708) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7171377) q[0];
sx q[0];
rz(-2.0949445) q[0];
sx q[0];
rz(-2.9852887) q[0];
rz(-pi) q[1];
rz(1.9916612) q[2];
sx q[2];
rz(-1.6399584) q[2];
sx q[2];
rz(1.7771306) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7397592) q[1];
sx q[1];
rz(-1.9976915) q[1];
sx q[1];
rz(3.0763807) q[1];
rz(-pi) q[2];
rz(-0.23128831) q[3];
sx q[3];
rz(-1.7219208) q[3];
sx q[3];
rz(-0.95201991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57006449) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(-1.338039) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(-0.6750955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682482) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(-2.5937953) q[0];
rz(-0.785218) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(-2.8731667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82314202) q[0];
sx q[0];
rz(-0.25932352) q[0];
sx q[0];
rz(0.32487049) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7873017) q[2];
sx q[2];
rz(-0.95677081) q[2];
sx q[2];
rz(-1.3348483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56663471) q[1];
sx q[1];
rz(-2.1556427) q[1];
sx q[1];
rz(-0.11771867) q[1];
rz(2.0500171) q[3];
sx q[3];
rz(-0.5939393) q[3];
sx q[3];
rz(1.493243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5238374) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-2.3274373) q[2];
rz(0.37627775) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(1.0193753) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(2.6928435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9578581) q[0];
sx q[0];
rz(-1.0976037) q[0];
sx q[0];
rz(-0.26259043) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8837187) q[2];
sx q[2];
rz(-1.2784625) q[2];
sx q[2];
rz(2.984798) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.031520695) q[1];
sx q[1];
rz(-1.8247461) q[1];
sx q[1];
rz(-2.5961848) q[1];
rz(-0.87047808) q[3];
sx q[3];
rz(-2.212489) q[3];
sx q[3];
rz(2.2023647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2723508) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0712414) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.8810133) q[0];
rz(-2.4977327) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(0.018928122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32301329) q[0];
sx q[0];
rz(-1.5529263) q[0];
sx q[0];
rz(0.26364003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.47255) q[2];
sx q[2];
rz(-1.7260523) q[2];
sx q[2];
rz(1.9948024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85266528) q[1];
sx q[1];
rz(-1.7551433) q[1];
sx q[1];
rz(-1.6221371) q[1];
rz(2.8519565) q[3];
sx q[3];
rz(-1.7970861) q[3];
sx q[3];
rz(2.3895404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(0.67031676) q[2];
rz(0.51236764) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(1.5669426) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(-1.0826899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854541) q[0];
sx q[0];
rz(-2.1573665) q[0];
sx q[0];
rz(-1.5120718) q[0];
x q[1];
rz(-0.81929368) q[2];
sx q[2];
rz(-1.594992) q[2];
sx q[2];
rz(-1.0011315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6607099) q[1];
sx q[1];
rz(-1.0272659) q[1];
sx q[1];
rz(-1.5488312) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38090221) q[3];
sx q[3];
rz(-1.8758095) q[3];
sx q[3];
rz(-1.1327133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(2.7486457) q[3];
sx q[3];
rz(-1.4018551) q[3];
sx q[3];
rz(2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.186541) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(2.4004249) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(0.95460931) q[2];
sx q[2];
rz(-1.2381427) q[2];
sx q[2];
rz(2.8340813) q[2];
rz(-1.3988597) q[3];
sx q[3];
rz(-2.3098683) q[3];
sx q[3];
rz(1.5666425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
