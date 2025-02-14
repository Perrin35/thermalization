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
rz(-0.96002785) q[0];
sx q[0];
rz(-2.5224944) q[0];
sx q[0];
rz(1.7382789) q[0];
rz(2.4461441) q[1];
sx q[1];
rz(-0.56013501) q[1];
sx q[1];
rz(-2.4897895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8389002) q[0];
sx q[0];
rz(-2.3639285) q[0];
sx q[0];
rz(2.0582786) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.02044631) q[2];
sx q[2];
rz(-0.281535) q[2];
sx q[2];
rz(-0.26007465) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5320325) q[1];
sx q[1];
rz(-1.3390307) q[1];
sx q[1];
rz(1.5941117) q[1];
rz(-pi) q[2];
rz(0.27068287) q[3];
sx q[3];
rz(-2.3029165) q[3];
sx q[3];
rz(1.0974778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92242509) q[2];
sx q[2];
rz(-1.3895915) q[2];
sx q[2];
rz(2.8233042) q[2];
rz(2.2394032) q[3];
sx q[3];
rz(-2.2212494) q[3];
sx q[3];
rz(2.2778146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1282463) q[0];
sx q[0];
rz(-0.25307578) q[0];
sx q[0];
rz(2.6213562) q[0];
rz(-0.086961374) q[1];
sx q[1];
rz(-1.0528916) q[1];
sx q[1];
rz(2.3435977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28198081) q[0];
sx q[0];
rz(-0.3726633) q[0];
sx q[0];
rz(-0.57639846) q[0];
rz(-pi) q[1];
rz(2.3981061) q[2];
sx q[2];
rz(-1.6891589) q[2];
sx q[2];
rz(-1.3918882) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.89435951) q[1];
sx q[1];
rz(-1.1825359) q[1];
sx q[1];
rz(-1.4187993) q[1];
rz(-pi) q[2];
rz(1.1677443) q[3];
sx q[3];
rz(-1.155793) q[3];
sx q[3];
rz(2.5248506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8331376) q[2];
sx q[2];
rz(-1.986958) q[2];
sx q[2];
rz(2.7412565) q[2];
rz(0.18276246) q[3];
sx q[3];
rz(-0.57923135) q[3];
sx q[3];
rz(-1.8496752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0188145) q[0];
sx q[0];
rz(-3.0570539) q[0];
sx q[0];
rz(-1.9345181) q[0];
rz(1.5030376) q[1];
sx q[1];
rz(-1.8820347) q[1];
sx q[1];
rz(-0.17301339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4543145) q[0];
sx q[0];
rz(-1.955324) q[0];
sx q[0];
rz(-3.049535) q[0];
rz(-pi) q[1];
rz(1.7373213) q[2];
sx q[2];
rz(-2.5534667) q[2];
sx q[2];
rz(-3.0808133) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26721482) q[1];
sx q[1];
rz(-0.80433955) q[1];
sx q[1];
rz(-0.72093236) q[1];
rz(-pi) q[2];
rz(-0.2263308) q[3];
sx q[3];
rz(-0.76068288) q[3];
sx q[3];
rz(-2.5193391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59306899) q[2];
sx q[2];
rz(-0.58716455) q[2];
sx q[2];
rz(1.0347838) q[2];
rz(0.89094025) q[3];
sx q[3];
rz(-2.2726629) q[3];
sx q[3];
rz(-1.9152036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0767176) q[0];
sx q[0];
rz(-0.95219505) q[0];
sx q[0];
rz(-3.0127443) q[0];
rz(0.42875641) q[1];
sx q[1];
rz(-1.7196451) q[1];
sx q[1];
rz(-1.6645924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1554028) q[0];
sx q[0];
rz(-2.7895067) q[0];
sx q[0];
rz(0.19908631) q[0];
rz(-pi) q[1];
rz(-1.1294133) q[2];
sx q[2];
rz(-1.980482) q[2];
sx q[2];
rz(0.022155174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13969888) q[1];
sx q[1];
rz(-2.668758) q[1];
sx q[1];
rz(1.0919149) q[1];
rz(-pi) q[2];
rz(0.97954025) q[3];
sx q[3];
rz(-1.5698804) q[3];
sx q[3];
rz(2.9360848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19026549) q[2];
sx q[2];
rz(-1.4172047) q[2];
sx q[2];
rz(1.3235271) q[2];
rz(1.4642814) q[3];
sx q[3];
rz(-2.5530294) q[3];
sx q[3];
rz(2.467449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4071963) q[0];
sx q[0];
rz(-0.36235991) q[0];
sx q[0];
rz(1.4417484) q[0];
rz(0.4231407) q[1];
sx q[1];
rz(-2.1839881) q[1];
sx q[1];
rz(0.29022455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5231756) q[0];
sx q[0];
rz(-1.2783861) q[0];
sx q[0];
rz(-0.66177701) q[0];
rz(-pi) q[1];
rz(2.1802203) q[2];
sx q[2];
rz(-1.0842241) q[2];
sx q[2];
rz(-0.79604641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6685953) q[1];
sx q[1];
rz(-1.5221704) q[1];
sx q[1];
rz(-0.12811382) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30778389) q[3];
sx q[3];
rz(-1.3213385) q[3];
sx q[3];
rz(-1.4741999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4238802) q[2];
sx q[2];
rz(-0.78705698) q[2];
sx q[2];
rz(-0.94669739) q[2];
rz(-0.26149073) q[3];
sx q[3];
rz(-1.3638834) q[3];
sx q[3];
rz(0.3869032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.43712) q[0];
sx q[0];
rz(-1.8081212) q[0];
sx q[0];
rz(-2.3719846) q[0];
rz(2.9656124) q[1];
sx q[1];
rz(-1.8006005) q[1];
sx q[1];
rz(-1.5699068) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0933857) q[0];
sx q[0];
rz(-0.95952672) q[0];
sx q[0];
rz(0.78637357) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4909582) q[2];
sx q[2];
rz(-2.4208899) q[2];
sx q[2];
rz(-1.0716764) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45794931) q[1];
sx q[1];
rz(-1.7548326) q[1];
sx q[1];
rz(-0.37269784) q[1];
rz(-pi) q[2];
rz(1.0553618) q[3];
sx q[3];
rz(-1.1222708) q[3];
sx q[3];
rz(0.12628638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0104597) q[2];
sx q[2];
rz(-0.95644462) q[2];
sx q[2];
rz(1.825911) q[2];
rz(-1.9192421) q[3];
sx q[3];
rz(-2.0171916) q[3];
sx q[3];
rz(0.30092064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2216126) q[0];
sx q[0];
rz(-1.6918809) q[0];
sx q[0];
rz(-1.3395039) q[0];
rz(1.1969396) q[1];
sx q[1];
rz(-2.0868128) q[1];
sx q[1];
rz(-0.96492499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1227568) q[0];
sx q[0];
rz(-2.0792476) q[0];
sx q[0];
rz(1.6792137) q[0];
x q[1];
rz(1.9222505) q[2];
sx q[2];
rz(-1.2932475) q[2];
sx q[2];
rz(0.48257839) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33789148) q[1];
sx q[1];
rz(-1.6543904) q[1];
sx q[1];
rz(-0.22050942) q[1];
rz(-1.3960725) q[3];
sx q[3];
rz(-2.0479322) q[3];
sx q[3];
rz(-1.9486547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40184608) q[2];
sx q[2];
rz(-1.2082938) q[2];
sx q[2];
rz(1.9604663) q[2];
rz(-0.014569672) q[3];
sx q[3];
rz(-0.40112344) q[3];
sx q[3];
rz(1.1723664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75083098) q[0];
sx q[0];
rz(-1.2068692) q[0];
sx q[0];
rz(-0.75585756) q[0];
rz(-0.66468704) q[1];
sx q[1];
rz(-1.1988147) q[1];
sx q[1];
rz(1.2896779) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.792218) q[0];
sx q[0];
rz(-1.5070033) q[0];
sx q[0];
rz(0.79453461) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6634438) q[2];
sx q[2];
rz(-0.59403803) q[2];
sx q[2];
rz(2.4817634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7482704) q[1];
sx q[1];
rz(-2.7549208) q[1];
sx q[1];
rz(2.3063117) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3068376) q[3];
sx q[3];
rz(-2.8182013) q[3];
sx q[3];
rz(0.33228126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0083302) q[2];
sx q[2];
rz(-1.8066758) q[2];
sx q[2];
rz(1.5745715) q[2];
rz(1.3957006) q[3];
sx q[3];
rz(-1.402366) q[3];
sx q[3];
rz(0.76656669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181327) q[0];
sx q[0];
rz(-2.8657931) q[0];
sx q[0];
rz(-2.7511399) q[0];
rz(-2.0402724) q[1];
sx q[1];
rz(-1.7951671) q[1];
sx q[1];
rz(-2.2094545) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9207805) q[0];
sx q[0];
rz(-1.5340065) q[0];
sx q[0];
rz(-0.027618577) q[0];
x q[1];
rz(-0.53932346) q[2];
sx q[2];
rz(-1.612609) q[2];
sx q[2];
rz(1.5635179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80492067) q[1];
sx q[1];
rz(-1.3778566) q[1];
sx q[1];
rz(2.8590917) q[1];
rz(-pi) q[2];
rz(1.8161536) q[3];
sx q[3];
rz(-1.1096769) q[3];
sx q[3];
rz(-0.86129118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.091888) q[2];
sx q[2];
rz(-0.69685093) q[2];
sx q[2];
rz(1.4616802) q[2];
rz(2.0563431) q[3];
sx q[3];
rz(-2.0094252) q[3];
sx q[3];
rz(0.6404883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853448) q[0];
sx q[0];
rz(-0.65971056) q[0];
sx q[0];
rz(-1.1261384) q[0];
rz(1.0477192) q[1];
sx q[1];
rz(-0.62647096) q[1];
sx q[1];
rz(-1.5315717) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5384993) q[0];
sx q[0];
rz(-1.5694008) q[0];
sx q[0];
rz(-3.1183232) q[0];
x q[1];
rz(-1.80349) q[2];
sx q[2];
rz(-1.1232716) q[2];
sx q[2];
rz(1.9767424) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23352392) q[1];
sx q[1];
rz(-1.6485813) q[1];
sx q[1];
rz(-3.0766355) q[1];
rz(-0.13981847) q[3];
sx q[3];
rz(-0.56666683) q[3];
sx q[3];
rz(-1.293928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5213726) q[2];
sx q[2];
rz(-2.3282101) q[2];
sx q[2];
rz(-0.20624557) q[2];
rz(0.5591875) q[3];
sx q[3];
rz(-1.5234448) q[3];
sx q[3];
rz(-1.1243813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54022057) q[0];
sx q[0];
rz(-1.843353) q[0];
sx q[0];
rz(-2.699615) q[0];
rz(-1.1454918) q[1];
sx q[1];
rz(-1.835123) q[1];
sx q[1];
rz(1.8225972) q[1];
rz(-0.27362846) q[2];
sx q[2];
rz(-1.3009334) q[2];
sx q[2];
rz(1.0197659) q[2];
rz(-0.19552874) q[3];
sx q[3];
rz(-0.73368213) q[3];
sx q[3];
rz(2.1140221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
