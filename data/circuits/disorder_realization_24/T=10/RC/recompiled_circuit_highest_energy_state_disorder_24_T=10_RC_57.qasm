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
rz(3.0348294) q[0];
sx q[0];
rz(-1.3314629) q[0];
sx q[0];
rz(-1.7911628) q[0];
rz(0.30836937) q[1];
sx q[1];
rz(6.0344459) q[1];
sx q[1];
rz(11.587172) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0672507) q[0];
sx q[0];
rz(-0.76267159) q[0];
sx q[0];
rz(-3.0347155) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0151695) q[2];
sx q[2];
rz(-2.4326486) q[2];
sx q[2];
rz(2.5042748) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1147656) q[1];
sx q[1];
rz(-1.9239167) q[1];
sx q[1];
rz(-1.4985282) q[1];
rz(-pi) q[2];
rz(-2.7704846) q[3];
sx q[3];
rz(-1.733193) q[3];
sx q[3];
rz(1.1134256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.9305482) q[2];
sx q[2];
rz(0.63300526) q[2];
rz(-0.64166075) q[3];
sx q[3];
rz(-2.6810665) q[3];
sx q[3];
rz(-2.3708926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4577456) q[0];
sx q[0];
rz(-0.96413833) q[0];
sx q[0];
rz(0.82474166) q[0];
rz(1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(2.8214084) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8267248) q[0];
sx q[0];
rz(-1.1629675) q[0];
sx q[0];
rz(2.639222) q[0];
rz(-pi) q[1];
rz(1.2037781) q[2];
sx q[2];
rz(-2.6026155) q[2];
sx q[2];
rz(2.0808329) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0517438) q[1];
sx q[1];
rz(-1.6517354) q[1];
sx q[1];
rz(1.593394) q[1];
x q[2];
rz(-1.8289964) q[3];
sx q[3];
rz(-2.9123983) q[3];
sx q[3];
rz(-2.3865139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1227293) q[2];
sx q[2];
rz(-1.2168695) q[2];
sx q[2];
rz(-1.2790206) q[2];
rz(-3.0458798) q[3];
sx q[3];
rz(-2.0666104) q[3];
sx q[3];
rz(-1.8224645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1529481) q[0];
sx q[0];
rz(-0.64866346) q[0];
sx q[0];
rz(-2.1107819) q[0];
rz(0.55054322) q[1];
sx q[1];
rz(-0.85854733) q[1];
sx q[1];
rz(1.2446838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6129788) q[0];
sx q[0];
rz(-1.6242149) q[0];
sx q[0];
rz(-0.15941711) q[0];
rz(-2.1878476) q[2];
sx q[2];
rz(-1.4336406) q[2];
sx q[2];
rz(1.1378261) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0708187) q[1];
sx q[1];
rz(-2.4386775) q[1];
sx q[1];
rz(1.8902814) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38328992) q[3];
sx q[3];
rz(-2.2858738) q[3];
sx q[3];
rz(-1.6413309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9024258) q[2];
sx q[2];
rz(-2.5776165) q[2];
sx q[2];
rz(-1.172056) q[2];
rz(2.3160882) q[3];
sx q[3];
rz(-1.2647311) q[3];
sx q[3];
rz(-0.66506213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0498407) q[0];
sx q[0];
rz(-1.4565775) q[0];
sx q[0];
rz(2.7276584) q[0];
rz(-2.6992758) q[1];
sx q[1];
rz(-2.1374173) q[1];
sx q[1];
rz(-2.5962459) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53876153) q[0];
sx q[0];
rz(-1.520938) q[0];
sx q[0];
rz(-2.8661348) q[0];
rz(-2.4905861) q[2];
sx q[2];
rz(-2.2525313) q[2];
sx q[2];
rz(1.663547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0132349) q[1];
sx q[1];
rz(-1.8473133) q[1];
sx q[1];
rz(2.1091631) q[1];
rz(-pi) q[2];
rz(2.0796989) q[3];
sx q[3];
rz(-1.4413222) q[3];
sx q[3];
rz(2.9751301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7004194) q[2];
sx q[2];
rz(-1.2607231) q[2];
sx q[2];
rz(0.1612266) q[2];
rz(1.553933) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(0.59066311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081838354) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(-2.4181714) q[0];
rz(1.8266034) q[1];
sx q[1];
rz(-1.864121) q[1];
sx q[1];
rz(1.0378729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6110566) q[0];
sx q[0];
rz(-0.7270455) q[0];
sx q[0];
rz(3.0226991) q[0];
x q[1];
rz(-2.1234799) q[2];
sx q[2];
rz(-2.2166914) q[2];
sx q[2];
rz(1.8551265) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41748504) q[1];
sx q[1];
rz(-1.8383664) q[1];
sx q[1];
rz(2.6883283) q[1];
rz(-0.54609046) q[3];
sx q[3];
rz(-1.4442863) q[3];
sx q[3];
rz(1.5777335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7448685) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(2.569516) q[2];
rz(-0.62503302) q[3];
sx q[3];
rz(-1.0038989) q[3];
sx q[3];
rz(1.9151789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142757) q[0];
sx q[0];
rz(-0.033997424) q[0];
sx q[0];
rz(2.6183364) q[0];
rz(-1.322586) q[1];
sx q[1];
rz(-1.4679642) q[1];
sx q[1];
rz(0.62017131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4374728) q[0];
sx q[0];
rz(-1.0123555) q[0];
sx q[0];
rz(-2.0790786) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46862009) q[2];
sx q[2];
rz(-1.6362564) q[2];
sx q[2];
rz(0.58131389) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7249696) q[1];
sx q[1];
rz(-2.1955804) q[1];
sx q[1];
rz(0.10819541) q[1];
rz(-pi) q[2];
rz(2.9240588) q[3];
sx q[3];
rz(-1.7420046) q[3];
sx q[3];
rz(1.273716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4124477) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(-1.8577417) q[2];
rz(-2.1733952) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(-1.5865954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17484434) q[0];
sx q[0];
rz(-2.0059678) q[0];
sx q[0];
rz(1.9298166) q[0];
rz(-2.1719596) q[1];
sx q[1];
rz(-1.3215093) q[1];
sx q[1];
rz(-1.2744354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14718854) q[0];
sx q[0];
rz(-2.4321803) q[0];
sx q[0];
rz(1.8038007) q[0];
rz(-1.6658804) q[2];
sx q[2];
rz(-0.78721744) q[2];
sx q[2];
rz(-1.1457535) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44597065) q[1];
sx q[1];
rz(-1.8730764) q[1];
sx q[1];
rz(2.2119889) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5481765) q[3];
sx q[3];
rz(-2.6024266) q[3];
sx q[3];
rz(2.9447458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6756639) q[2];
sx q[2];
rz(-1.482654) q[2];
sx q[2];
rz(-0.3234123) q[2];
rz(3.1392426) q[3];
sx q[3];
rz(-1.4582783) q[3];
sx q[3];
rz(-0.6215483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49331409) q[0];
sx q[0];
rz(-1.4819773) q[0];
sx q[0];
rz(-1.7907273) q[0];
rz(1.7905275) q[1];
sx q[1];
rz(-1.8565145) q[1];
sx q[1];
rz(-1.2877119) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5918811) q[0];
sx q[0];
rz(-1.8533819) q[0];
sx q[0];
rz(-0.12235276) q[0];
rz(-pi) q[1];
rz(-1.9016198) q[2];
sx q[2];
rz(-0.74218732) q[2];
sx q[2];
rz(2.0798707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85387522) q[1];
sx q[1];
rz(-2.7734904) q[1];
sx q[1];
rz(2.6516857) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34890449) q[3];
sx q[3];
rz(-1.5530905) q[3];
sx q[3];
rz(-0.68428333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6172341) q[2];
sx q[2];
rz(-1.1523767) q[2];
sx q[2];
rz(0.99346811) q[2];
rz(-2.1067545) q[3];
sx q[3];
rz(-0.60641685) q[3];
sx q[3];
rz(-2.9724227) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8096823) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(-1.6140953) q[0];
rz(-2.2612259) q[1];
sx q[1];
rz(-0.4207193) q[1];
sx q[1];
rz(-0.31164935) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051119251) q[0];
sx q[0];
rz(-1.3731806) q[0];
sx q[0];
rz(-1.86905) q[0];
rz(-pi) q[1];
rz(0.016383532) q[2];
sx q[2];
rz(-1.3394622) q[2];
sx q[2];
rz(-2.6838059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5376725) q[1];
sx q[1];
rz(-2.2419856) q[1];
sx q[1];
rz(-1.0965986) q[1];
rz(-1.9945456) q[3];
sx q[3];
rz(-2.4422788) q[3];
sx q[3];
rz(-0.55845234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46939048) q[2];
sx q[2];
rz(-0.92350525) q[2];
sx q[2];
rz(-1.3908609) q[2];
rz(-1.5270799) q[3];
sx q[3];
rz(-2.0938087) q[3];
sx q[3];
rz(-2.0670149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020318) q[0];
sx q[0];
rz(-1.1447516) q[0];
sx q[0];
rz(-2.5290053) q[0];
rz(-2.1514429) q[1];
sx q[1];
rz(-1.1724816) q[1];
sx q[1];
rz(1.6506857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6263437) q[0];
sx q[0];
rz(-2.307461) q[0];
sx q[0];
rz(-2.5596928) q[0];
rz(-pi) q[1];
rz(-0.91103001) q[2];
sx q[2];
rz(-1.3675895) q[2];
sx q[2];
rz(-0.045762941) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13496357) q[1];
sx q[1];
rz(-1.7322759) q[1];
sx q[1];
rz(3.1352741) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4435639) q[3];
sx q[3];
rz(-1.3778316) q[3];
sx q[3];
rz(-2.2730227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7154197) q[2];
sx q[2];
rz(-0.82948589) q[2];
sx q[2];
rz(0.034991525) q[2];
rz(1.70111) q[3];
sx q[3];
rz(-2.248843) q[3];
sx q[3];
rz(0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424425) q[0];
sx q[0];
rz(-0.89542605) q[0];
sx q[0];
rz(-0.17740346) q[0];
rz(-1.9433446) q[1];
sx q[1];
rz(-1.5180963) q[1];
sx q[1];
rz(0.95388283) q[1];
rz(-2.4468471) q[2];
sx q[2];
rz(-2.3320945) q[2];
sx q[2];
rz(3.1041067) q[2];
rz(0.71975868) q[3];
sx q[3];
rz(-1.9282785) q[3];
sx q[3];
rz(-0.45784605) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
