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
rz(0.69005203) q[0];
sx q[0];
rz(-0.85059387) q[0];
sx q[0];
rz(0.2001065) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(1.8097635) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1855186) q[0];
sx q[0];
rz(-1.3393434) q[0];
sx q[0];
rz(2.0218019) q[0];
rz(1.0806141) q[2];
sx q[2];
rz(-2.1128383) q[2];
sx q[2];
rz(1.3203743) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7847474) q[1];
sx q[1];
rz(-0.6797528) q[1];
sx q[1];
rz(-3.0650782) q[1];
rz(-pi) q[2];
rz(1.41134) q[3];
sx q[3];
rz(-0.48122367) q[3];
sx q[3];
rz(3.0414097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.072210463) q[2];
sx q[2];
rz(-2.78077) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(-0.34590736) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2717993) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(2.2858009) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(2.3537297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44875328) q[0];
sx q[0];
rz(-1.8215382) q[0];
sx q[0];
rz(3.1186799) q[0];
rz(-pi) q[1];
rz(-1.5283952) q[2];
sx q[2];
rz(-1.2268492) q[2];
sx q[2];
rz(-2.6478772) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7447506) q[1];
sx q[1];
rz(-1.9991268) q[1];
sx q[1];
rz(2.7672) q[1];
x q[2];
rz(-0.26420748) q[3];
sx q[3];
rz(-2.4572861) q[3];
sx q[3];
rz(-2.5591889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0756125) q[2];
sx q[2];
rz(-1.2300666) q[2];
sx q[2];
rz(-0.39609972) q[2];
rz(1.570545) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.9539255) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(3.0809825) q[0];
rz(0.68471471) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(0.33904591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44797541) q[0];
sx q[0];
rz(-2.1537499) q[0];
sx q[0];
rz(-1.5459803) q[0];
x q[1];
rz(1.4475559) q[2];
sx q[2];
rz(-1.3831105) q[2];
sx q[2];
rz(-2.0579684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5506997) q[1];
sx q[1];
rz(-1.8807481) q[1];
sx q[1];
rz(2.8551686) q[1];
rz(-pi) q[2];
rz(2.364758) q[3];
sx q[3];
rz(-1.562444) q[3];
sx q[3];
rz(-1.2119105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.8637916) q[2];
sx q[2];
rz(-0.47046146) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(-0.40801868) q[0];
rz(1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(1.8203576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-0.59998673) q[0];
sx q[0];
rz(2.5866051) q[0];
rz(-pi) q[1];
rz(3.0555326) q[2];
sx q[2];
rz(-0.91740184) q[2];
sx q[2];
rz(1.6376405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8712052) q[1];
sx q[1];
rz(-3.0000038) q[1];
sx q[1];
rz(-2.6006727) q[1];
rz(0.77797555) q[3];
sx q[3];
rz(-2.5730592) q[3];
sx q[3];
rz(-1.3028627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9185751) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(0.72285405) q[2];
rz(-2.2919848) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1866813) q[0];
sx q[0];
rz(-2.6246922) q[0];
sx q[0];
rz(-2.5816259) q[0];
rz(-2.9365183) q[1];
sx q[1];
rz(-1.8849025) q[1];
sx q[1];
rz(-2.8505039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4093035) q[0];
sx q[0];
rz(-1.3025874) q[0];
sx q[0];
rz(-1.3951673) q[0];
x q[1];
rz(-1.6339602) q[2];
sx q[2];
rz(-0.77253714) q[2];
sx q[2];
rz(0.26320266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2771137) q[1];
sx q[1];
rz(-0.8833589) q[1];
sx q[1];
rz(2.8836125) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3897393) q[3];
sx q[3];
rz(-1.6883779) q[3];
sx q[3];
rz(2.8428809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9186972) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(0.29829868) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(-1.6303308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30763141) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(0.6024012) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(-2.1098302) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4435972) q[0];
sx q[0];
rz(-1.5549608) q[0];
sx q[0];
rz(-1.3965194) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6844031) q[2];
sx q[2];
rz(-2.6962099) q[2];
sx q[2];
rz(-0.066748652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25615869) q[1];
sx q[1];
rz(-0.51894655) q[1];
sx q[1];
rz(2.7536254) q[1];
x q[2];
rz(-2.5661181) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(-0.3994715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80456698) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.6229013) q[2];
rz(-0.8484146) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0254211) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(2.1436932) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(1.3708699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8557963) q[0];
sx q[0];
rz(-1.5842452) q[0];
sx q[0];
rz(2.5105259) q[0];
rz(-pi) q[1];
rz(0.15372686) q[2];
sx q[2];
rz(-1.2187374) q[2];
sx q[2];
rz(-1.5444666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9129907) q[1];
sx q[1];
rz(-2.1739042) q[1];
sx q[1];
rz(-1.2772389) q[1];
x q[2];
rz(0.73565817) q[3];
sx q[3];
rz(-2.1950985) q[3];
sx q[3];
rz(-1.5746547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-3.1119463) q[2];
rz(-0.61521411) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(-0.7830559) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(-1.7436183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0692978) q[0];
sx q[0];
rz(-2.0136478) q[0];
sx q[0];
rz(2.686727) q[0];
x q[1];
rz(0.88770788) q[2];
sx q[2];
rz(-0.36044932) q[2];
sx q[2];
rz(0.53149022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82247558) q[1];
sx q[1];
rz(-0.90293316) q[1];
sx q[1];
rz(0.0112337) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21381883) q[3];
sx q[3];
rz(-0.79584661) q[3];
sx q[3];
rz(-2.4406432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29584259) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.3207377) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(-2.1871908) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(2.0571713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99710871) q[0];
sx q[0];
rz(-2.2786744) q[0];
sx q[0];
rz(3.0860712) q[0];
rz(-pi) q[1];
rz(-2.3060572) q[2];
sx q[2];
rz(-0.37097574) q[2];
sx q[2];
rz(-2.6717466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7560278) q[1];
sx q[1];
rz(-1.4984908) q[1];
sx q[1];
rz(1.0977618) q[1];
x q[2];
rz(2.3929651) q[3];
sx q[3];
rz(-1.7065062) q[3];
sx q[3];
rz(2.3627797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0584917) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(-1.3433749) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042353543) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.7310671) q[0];
rz(-2.9208185) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-0.9332307) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8932874) q[0];
sx q[0];
rz(-0.74549216) q[0];
sx q[0];
rz(0.020123578) q[0];
rz(-pi) q[1];
rz(1.2104697) q[2];
sx q[2];
rz(-0.98305741) q[2];
sx q[2];
rz(0.34229842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5367914) q[1];
sx q[1];
rz(-1.7030506) q[1];
sx q[1];
rz(2.3459646) q[1];
rz(-1.7228863) q[3];
sx q[3];
rz(-2.503431) q[3];
sx q[3];
rz(-1.9877951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.71020469) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(0.63146511) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5030293) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-2.6999264) q[1];
sx q[1];
rz(-1.8001945) q[1];
sx q[1];
rz(1.1460907) q[1];
rz(-1.491811) q[2];
sx q[2];
rz(-0.38182237) q[2];
sx q[2];
rz(0.67425722) q[2];
rz(-0.56959116) q[3];
sx q[3];
rz(-1.6574331) q[3];
sx q[3];
rz(2.1059753) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
