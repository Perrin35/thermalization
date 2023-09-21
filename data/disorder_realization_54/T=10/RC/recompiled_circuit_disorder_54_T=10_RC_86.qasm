OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(-0.33049345) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(2.4324774) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.37492) q[0];
sx q[0];
rz(-1.8639038) q[0];
sx q[0];
rz(2.2105182) q[0];
rz(-pi) q[1];
rz(1.2572631) q[2];
sx q[2];
rz(-1.5402113) q[2];
sx q[2];
rz(2.7262296) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.658537) q[1];
sx q[1];
rz(-0.72209789) q[1];
sx q[1];
rz(-0.87941054) q[1];
rz(2.764774) q[3];
sx q[3];
rz(-1.0816649) q[3];
sx q[3];
rz(-2.7917142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7314529) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(1.6072134) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(-0.62227917) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(-0.91631779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371671) q[0];
sx q[0];
rz(-0.93351782) q[0];
sx q[0];
rz(0.62563719) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83444886) q[2];
sx q[2];
rz(-1.6685899) q[2];
sx q[2];
rz(2.9546839) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.037109) q[1];
sx q[1];
rz(-2.6625405) q[1];
sx q[1];
rz(-0.80161174) q[1];
rz(-0.022577062) q[3];
sx q[3];
rz(-1.9687679) q[3];
sx q[3];
rz(-2.857739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24094412) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(0.82143482) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18773742) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(-1.0082555) q[0];
rz(-0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(-2.6170513) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0073111) q[0];
sx q[0];
rz(-1.6005922) q[0];
sx q[0];
rz(-1.7344463) q[0];
x q[1];
rz(0.67851615) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(3.0845272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9959065) q[1];
sx q[1];
rz(-1.5090764) q[1];
sx q[1];
rz(-1.5476336) q[1];
rz(-pi) q[2];
rz(2.9933661) q[3];
sx q[3];
rz(-1.1550511) q[3];
sx q[3];
rz(1.9849329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(1.3211936) q[3];
sx q[3];
rz(-1.7318055) q[3];
sx q[3];
rz(-2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33048531) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(3.0526429) q[0];
rz(2.6308909) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8331497) q[0];
sx q[0];
rz(-1.2763378) q[0];
sx q[0];
rz(1.4008646) q[0];
x q[1];
rz(-0.34822779) q[2];
sx q[2];
rz(-0.59154445) q[2];
sx q[2];
rz(-2.0915742) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.888962) q[1];
sx q[1];
rz(-1.2446212) q[1];
sx q[1];
rz(0.77146448) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48456405) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(2.518667) q[2];
rz(1.1359435) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(2.6749271) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(2.5581397) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(1.6437644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72901112) q[0];
sx q[0];
rz(-2.3238365) q[0];
sx q[0];
rz(1.9304995) q[0];
rz(-pi) q[1];
rz(1.4719109) q[2];
sx q[2];
rz(-2.3660198) q[2];
sx q[2];
rz(1.0727739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.076482) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-0.70365023) q[1];
rz(-pi) q[2];
rz(0.28108092) q[3];
sx q[3];
rz(-2.6087458) q[3];
sx q[3];
rz(-1.7588774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65537611) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(1.5779457) q[2];
rz(0.90562138) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.6435796) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.7061589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2611321) q[0];
sx q[0];
rz(-1.7186223) q[0];
sx q[0];
rz(0.79974215) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4749182) q[2];
sx q[2];
rz(-0.55170689) q[2];
sx q[2];
rz(-0.72223896) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1286436) q[1];
sx q[1];
rz(-1.9777021) q[1];
sx q[1];
rz(0.13601555) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4581497) q[3];
sx q[3];
rz(-0.23922353) q[3];
sx q[3];
rz(2.3000172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(-3.0701239) q[2];
rz(-1.4525157) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(2.563971) q[0];
rz(1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-1.0095899) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44430915) q[0];
sx q[0];
rz(-2.5464006) q[0];
sx q[0];
rz(-0.62023456) q[0];
rz(-pi) q[1];
rz(0.52092123) q[2];
sx q[2];
rz(-2.705057) q[2];
sx q[2];
rz(-1.6233363) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0040782) q[1];
sx q[1];
rz(-2.5234748) q[1];
sx q[1];
rz(2.532258) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76241242) q[3];
sx q[3];
rz(-2.0250642) q[3];
sx q[3];
rz(-1.7319958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(-1.1996777) q[2];
rz(0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49884477) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(2.902466) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(2.696864) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0895773) q[0];
sx q[0];
rz(-1.9964295) q[0];
sx q[0];
rz(0.94469597) q[0];
x q[1];
rz(-0.76079255) q[2];
sx q[2];
rz(-0.82423254) q[2];
sx q[2];
rz(-0.21812083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4152894) q[1];
sx q[1];
rz(-2.2518034) q[1];
sx q[1];
rz(-2.9773832) q[1];
rz(1.7041676) q[3];
sx q[3];
rz(-1.1629472) q[3];
sx q[3];
rz(-2.9123902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(2.3262809) q[2];
rz(1.0347962) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(-1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(0.28717336) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-0.33219355) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0567386) q[0];
sx q[0];
rz(-0.85002725) q[0];
sx q[0];
rz(-0.42215729) q[0];
rz(2.2696482) q[2];
sx q[2];
rz(-1.6064062) q[2];
sx q[2];
rz(-1.2138838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48391446) q[1];
sx q[1];
rz(-0.49424833) q[1];
sx q[1];
rz(2.5792522) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6094749) q[3];
sx q[3];
rz(-1.2174774) q[3];
sx q[3];
rz(-1.4766828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93280783) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(-1.2906637) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(0.65657842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-2.8870781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018325018) q[0];
sx q[0];
rz(-1.8013445) q[0];
sx q[0];
rz(0.075320764) q[0];
rz(1.8149257) q[2];
sx q[2];
rz(-1.0955398) q[2];
sx q[2];
rz(0.82820669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.49452457) q[1];
sx q[1];
rz(-1.2459323) q[1];
sx q[1];
rz(-1.4855794) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7032095) q[3];
sx q[3];
rz(-1.4808726) q[3];
sx q[3];
rz(-0.59059483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(1.7808328) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(1.9706456) q[2];
sx q[2];
rz(-0.56213899) q[2];
sx q[2];
rz(3.0664372) q[2];
rz(1.249282) q[3];
sx q[3];
rz(-0.63890639) q[3];
sx q[3];
rz(-1.3756868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
