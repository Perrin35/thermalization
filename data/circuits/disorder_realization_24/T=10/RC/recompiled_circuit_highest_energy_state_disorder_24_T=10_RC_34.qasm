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
rz(-0.10676323) q[0];
sx q[0];
rz(4.4730555) q[0];
sx q[0];
rz(8.0743481) q[0];
rz(-2.8332233) q[1];
sx q[1];
rz(-2.8928533) q[1];
sx q[1];
rz(-2.1623936) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0686183) q[0];
sx q[0];
rz(-2.3280297) q[0];
sx q[0];
rz(-1.4692151) q[0];
rz(-pi) q[1];
rz(2.2297284) q[2];
sx q[2];
rz(-1.2871337) q[2];
sx q[2];
rz(-2.5549169) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82047082) q[1];
sx q[1];
rz(-0.36013569) q[1];
sx q[1];
rz(-2.9481357) q[1];
rz(-pi) q[2];
rz(-2.7704846) q[3];
sx q[3];
rz(-1.733193) q[3];
sx q[3];
rz(-2.028167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.2110445) q[2];
sx q[2];
rz(2.5085874) q[2];
rz(2.4999319) q[3];
sx q[3];
rz(-0.46052614) q[3];
sx q[3];
rz(2.3708926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68384701) q[0];
sx q[0];
rz(-2.1774543) q[0];
sx q[0];
rz(0.82474166) q[0];
rz(-1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(-2.8214084) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8267248) q[0];
sx q[0];
rz(-1.1629675) q[0];
sx q[0];
rz(-2.639222) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0616706) q[2];
sx q[2];
rz(-1.3855644) q[2];
sx q[2];
rz(2.9502466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9588543) q[1];
sx q[1];
rz(-3.0575648) q[1];
sx q[1];
rz(-2.8699204) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3489487) q[3];
sx q[3];
rz(-1.6288405) q[3];
sx q[3];
rz(-1.0674509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1227293) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(-1.2790206) q[2];
rz(-3.0458798) q[3];
sx q[3];
rz(-1.0749823) q[3];
sx q[3];
rz(-1.3191282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1529481) q[0];
sx q[0];
rz(-0.64866346) q[0];
sx q[0];
rz(2.1107819) q[0];
rz(-0.55054322) q[1];
sx q[1];
rz(-0.85854733) q[1];
sx q[1];
rz(1.8969089) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36280945) q[0];
sx q[0];
rz(-2.9735357) q[0];
sx q[0];
rz(-0.3248949) q[0];
x q[1];
rz(-2.9739506) q[2];
sx q[2];
rz(-2.1812005) q[2];
sx q[2];
rz(-2.8053225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0708187) q[1];
sx q[1];
rz(-2.4386775) q[1];
sx q[1];
rz(1.2513112) q[1];
rz(1.1641509) q[3];
sx q[3];
rz(-0.79504025) q[3];
sx q[3];
rz(-2.1926978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9024258) q[2];
sx q[2];
rz(-2.5776165) q[2];
sx q[2];
rz(-1.172056) q[2];
rz(-0.82550448) q[3];
sx q[3];
rz(-1.8768616) q[3];
sx q[3];
rz(-2.4765305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0498407) q[0];
sx q[0];
rz(-1.4565775) q[0];
sx q[0];
rz(0.41393429) q[0];
rz(0.44231689) q[1];
sx q[1];
rz(-2.1374173) q[1];
sx q[1];
rz(0.54534674) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6028311) q[0];
sx q[0];
rz(-1.520938) q[0];
sx q[0];
rz(2.8661348) q[0];
x q[1];
rz(-2.3661883) q[2];
sx q[2];
rz(-2.0606962) q[2];
sx q[2];
rz(-2.6012914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0132349) q[1];
sx q[1];
rz(-1.8473133) q[1];
sx q[1];
rz(-2.1091631) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14800565) q[3];
sx q[3];
rz(-2.0750351) q[3];
sx q[3];
rz(-1.809169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7004194) q[2];
sx q[2];
rz(-1.2607231) q[2];
sx q[2];
rz(-0.1612266) q[2];
rz(-1.5876596) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(0.59066311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0597543) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(2.4181714) q[0];
rz(1.8266034) q[1];
sx q[1];
rz(-1.2774717) q[1];
sx q[1];
rz(2.1037197) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6110566) q[0];
sx q[0];
rz(-2.4145472) q[0];
sx q[0];
rz(-3.0226991) q[0];
rz(2.1234799) q[2];
sx q[2];
rz(-0.92490126) q[2];
sx q[2];
rz(-1.2864662) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2813869) q[1];
sx q[1];
rz(-2.0067999) q[1];
sx q[1];
rz(-1.2748177) q[1];
rz(-pi) q[2];
rz(0.24017374) q[3];
sx q[3];
rz(-2.5824912) q[3];
sx q[3];
rz(-0.19776519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3967241) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(-0.57207668) q[2];
rz(-0.62503302) q[3];
sx q[3];
rz(-1.0038989) q[3];
sx q[3];
rz(-1.2264138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6201651) q[0];
sx q[0];
rz(-3.1075952) q[0];
sx q[0];
rz(0.52325621) q[0];
rz(-1.322586) q[1];
sx q[1];
rz(-1.6736284) q[1];
sx q[1];
rz(-0.62017131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5145743) q[0];
sx q[0];
rz(-2.4051599) q[0];
sx q[0];
rz(-0.66177841) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46862009) q[2];
sx q[2];
rz(-1.5053362) q[2];
sx q[2];
rz(2.5602788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7249696) q[1];
sx q[1];
rz(-0.94601224) q[1];
sx q[1];
rz(-0.10819541) q[1];
rz(-pi) q[2];
rz(2.4661651) q[3];
sx q[3];
rz(-0.27598652) q[3];
sx q[3];
rz(2.7819862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72914499) q[2];
sx q[2];
rz(-2.3690988) q[2];
sx q[2];
rz(1.8577417) q[2];
rz(-0.9681975) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(1.5865954) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17484434) q[0];
sx q[0];
rz(-1.1356249) q[0];
sx q[0];
rz(-1.211776) q[0];
rz(2.1719596) q[1];
sx q[1];
rz(-1.3215093) q[1];
sx q[1];
rz(1.2744354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2454555) q[0];
sx q[0];
rz(-1.7217759) q[0];
sx q[0];
rz(0.87489382) q[0];
rz(-pi) q[1];
rz(-1.6658804) q[2];
sx q[2];
rz(-0.78721744) q[2];
sx q[2];
rz(1.9958391) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.74531534) q[1];
sx q[1];
rz(-2.4419129) q[1];
sx q[1];
rz(-2.0513644) q[1];
rz(1.2479374) q[3];
sx q[3];
rz(-1.1311232) q[3];
sx q[3];
rz(-0.46931258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4659287) q[2];
sx q[2];
rz(-1.6589386) q[2];
sx q[2];
rz(-2.8181804) q[2];
rz(0.0023500738) q[3];
sx q[3];
rz(-1.4582783) q[3];
sx q[3];
rz(0.6215483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49331409) q[0];
sx q[0];
rz(-1.6596153) q[0];
sx q[0];
rz(1.7907273) q[0];
rz(-1.3510652) q[1];
sx q[1];
rz(-1.2850782) q[1];
sx q[1];
rz(1.2877119) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13437102) q[0];
sx q[0];
rz(-2.8343081) q[0];
sx q[0];
rz(-1.1728806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85629927) q[2];
sx q[2];
rz(-1.7921471) q[2];
sx q[2];
rz(2.3846087) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1786036) q[1];
sx q[1];
rz(-1.4006536) q[1];
sx q[1];
rz(2.8135706) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.051750465) q[3];
sx q[3];
rz(-2.7922575) q[3];
sx q[3];
rz(2.3037095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6172341) q[2];
sx q[2];
rz(-1.1523767) q[2];
sx q[2];
rz(0.99346811) q[2];
rz(-1.0348381) q[3];
sx q[3];
rz(-0.60641685) q[3];
sx q[3];
rz(-0.16916999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.33191037) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(1.6140953) q[0];
rz(-0.88036674) q[1];
sx q[1];
rz(-0.4207193) q[1];
sx q[1];
rz(-2.8299433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0904734) q[0];
sx q[0];
rz(-1.3731806) q[0];
sx q[0];
rz(1.2725426) q[0];
rz(-pi) q[1];
rz(-3.1252091) q[2];
sx q[2];
rz(-1.3394622) q[2];
sx q[2];
rz(-2.6838059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.27585718) q[1];
sx q[1];
rz(-1.2051264) q[1];
sx q[1];
rz(-0.72876282) q[1];
x q[2];
rz(1.1470471) q[3];
sx q[3];
rz(-2.4422788) q[3];
sx q[3];
rz(-0.55845234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46939048) q[2];
sx q[2];
rz(-0.92350525) q[2];
sx q[2];
rz(1.3908609) q[2];
rz(-1.5270799) q[3];
sx q[3];
rz(-1.047784) q[3];
sx q[3];
rz(2.0670149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63956082) q[0];
sx q[0];
rz(-1.9968411) q[0];
sx q[0];
rz(-0.61258739) q[0];
rz(0.9901498) q[1];
sx q[1];
rz(-1.9691111) q[1];
sx q[1];
rz(-1.6506857) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51524893) q[0];
sx q[0];
rz(-2.307461) q[0];
sx q[0];
rz(-2.5596928) q[0];
rz(-pi) q[1];
rz(2.2305626) q[2];
sx q[2];
rz(-1.7740031) q[2];
sx q[2];
rz(0.045762941) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4348169) q[1];
sx q[1];
rz(-1.5770327) q[1];
sx q[1];
rz(1.732279) q[1];
rz(-pi) q[2];
rz(-0.69802876) q[3];
sx q[3];
rz(-1.3778316) q[3];
sx q[3];
rz(0.86856996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7154197) q[2];
sx q[2];
rz(-2.3121068) q[2];
sx q[2];
rz(-0.034991525) q[2];
rz(1.70111) q[3];
sx q[3];
rz(-2.248843) q[3];
sx q[3];
rz(-2.6006202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(0.69915019) q[0];
sx q[0];
rz(-0.89542605) q[0];
sx q[0];
rz(-0.17740346) q[0];
rz(1.198248) q[1];
sx q[1];
rz(-1.5180963) q[1];
sx q[1];
rz(0.95388283) q[1];
rz(-0.69474557) q[2];
sx q[2];
rz(-0.80949819) q[2];
sx q[2];
rz(-0.037485952) q[2];
rz(-0.51552897) q[3];
sx q[3];
rz(-0.78918189) q[3];
sx q[3];
rz(0.73317151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
