OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1928756) q[0];
sx q[0];
rz(-2.2964381) q[0];
sx q[0];
rz(-1.2851508) q[0];
x q[1];
rz(1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(2.8710499) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.018651389) q[1];
sx q[1];
rz(-0.39995799) q[1];
sx q[1];
rz(-0.33756983) q[1];
x q[2];
rz(-2.5472766) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(0.017410226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801385) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-2.8161312) q[0];
rz(1.356396) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.9869841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017529537) q[0];
sx q[0];
rz(-1.5520099) q[0];
sx q[0];
rz(-1.5879052) q[0];
x q[1];
rz(1.0500533) q[2];
sx q[2];
rz(-2.4467231) q[2];
sx q[2];
rz(2.3827202) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6007538) q[1];
sx q[1];
rz(-2.3327017) q[1];
sx q[1];
rz(-1.4536588) q[1];
rz(-0.18493821) q[3];
sx q[3];
rz(-1.8293081) q[3];
sx q[3];
rz(-2.1188494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31323355) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-1.0916969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900345) q[0];
sx q[0];
rz(-1.0150195) q[0];
sx q[0];
rz(-0.65727289) q[0];
rz(-pi) q[1];
rz(1.0513564) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(1.6348334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25635438) q[1];
sx q[1];
rz(-1.2043722) q[1];
sx q[1];
rz(-0.79856915) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(0.88095218) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-2.8312347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46582857) q[0];
sx q[0];
rz(-2.7390263) q[0];
sx q[0];
rz(-0.34253828) q[0];
x q[1];
rz(0.68508673) q[2];
sx q[2];
rz(-1.6685467) q[2];
sx q[2];
rz(0.098066559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1102317) q[1];
sx q[1];
rz(-2.7831315) q[1];
sx q[1];
rz(-1.6237753) q[1];
rz(1.849732) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(-1.0517373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0115396) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(0.87019428) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2244959) q[0];
sx q[0];
rz(-1.774569) q[0];
sx q[0];
rz(-2.8208371) q[0];
x q[1];
rz(-2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(0.86415926) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89228499) q[1];
sx q[1];
rz(-1.2586437) q[1];
sx q[1];
rz(-2.409163) q[1];
x q[2];
rz(-1.3051885) q[3];
sx q[3];
rz(-2.5634273) q[3];
sx q[3];
rz(2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(0.67374054) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34981397) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(2.8836024) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.4917096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0598037) q[0];
sx q[0];
rz(-0.44133082) q[0];
sx q[0];
rz(2.9102737) q[0];
rz(-pi) q[1];
rz(-2.6368124) q[2];
sx q[2];
rz(-2.3961888) q[2];
sx q[2];
rz(0.2573075) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7548435) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(-0.90653231) q[1];
rz(1.9353988) q[3];
sx q[3];
rz(-1.6657262) q[3];
sx q[3];
rz(1.3473542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-0.091726124) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(-2.2475524) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3180852) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-3.1076028) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720374) q[0];
sx q[0];
rz(-0.79622686) q[0];
sx q[0];
rz(-1.3433334) q[0];
rz(2.6007973) q[2];
sx q[2];
rz(-1.0058837) q[2];
sx q[2];
rz(0.62513798) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1696724) q[1];
sx q[1];
rz(-1.5202513) q[1];
sx q[1];
rz(-2.2488942) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7185228) q[3];
sx q[3];
rz(-0.96102321) q[3];
sx q[3];
rz(-0.33965286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6035446) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2440764) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(-2.774033) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.6961018) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.53795) q[0];
sx q[0];
rz(-2.0041487) q[0];
sx q[0];
rz(3.135878) q[0];
x q[1];
rz(-2.3169575) q[2];
sx q[2];
rz(-2.3687009) q[2];
sx q[2];
rz(-0.075721272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3280914) q[1];
sx q[1];
rz(-2.5655167) q[1];
sx q[1];
rz(-0.2920132) q[1];
rz(-pi) q[2];
rz(0.66283488) q[3];
sx q[3];
rz(-1.0532866) q[3];
sx q[3];
rz(2.585632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(0.30977419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4944086) q[0];
sx q[0];
rz(-1.697288) q[0];
sx q[0];
rz(2.6148318) q[0];
rz(-pi) q[1];
rz(1.5179713) q[2];
sx q[2];
rz(-2.9644358) q[2];
sx q[2];
rz(-1.0747386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27948353) q[1];
sx q[1];
rz(-0.21394193) q[1];
sx q[1];
rz(-1.2941542) q[1];
rz(-1.314332) q[3];
sx q[3];
rz(-2.1520352) q[3];
sx q[3];
rz(-0.47282156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(-1.9343728) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(1.6419798) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6584872) q[0];
sx q[0];
rz(-2.8327836) q[0];
sx q[0];
rz(-1.2312066) q[0];
x q[1];
rz(-0.070812289) q[2];
sx q[2];
rz(-0.76528463) q[2];
sx q[2];
rz(-0.27809696) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75913945) q[1];
sx q[1];
rz(-1.1083974) q[1];
sx q[1];
rz(0.23666246) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5564735) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(-2.1567791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.88400921) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-3.070667) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(-0.82540853) q[2];
sx q[2];
rz(-2.5054629) q[2];
sx q[2];
rz(2.7873743) q[2];
rz(-0.72017097) q[3];
sx q[3];
rz(-1.1838893) q[3];
sx q[3];
rz(0.67223687) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
