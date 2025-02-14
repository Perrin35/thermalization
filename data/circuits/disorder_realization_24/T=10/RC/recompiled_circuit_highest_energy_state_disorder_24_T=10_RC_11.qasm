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
rz(0.97919908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0672507) q[0];
sx q[0];
rz(-0.76267159) q[0];
sx q[0];
rz(-3.0347155) q[0];
x q[1];
rz(-2.2297284) q[2];
sx q[2];
rz(-1.8544589) q[2];
sx q[2];
rz(-2.5549169) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5725928) q[1];
sx q[1];
rz(-1.6385983) q[1];
sx q[1];
rz(0.35396934) q[1];
rz(-pi) q[2];
rz(2.7172543) q[3];
sx q[3];
rz(-0.40357129) q[3];
sx q[3];
rz(2.2904325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.9305482) q[2];
sx q[2];
rz(-2.5085874) q[2];
rz(2.4999319) q[3];
sx q[3];
rz(-0.46052614) q[3];
sx q[3];
rz(2.3708926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4577456) q[0];
sx q[0];
rz(-0.96413833) q[0];
sx q[0];
rz(-0.82474166) q[0];
rz(-1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(0.32018426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6711222) q[0];
sx q[0];
rz(-2.0286386) q[0];
sx q[0];
rz(1.1127959) q[0];
x q[1];
rz(-1.0616706) q[2];
sx q[2];
rz(-1.3855644) q[2];
sx q[2];
rz(-2.9502466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.089848891) q[1];
sx q[1];
rz(-1.4898572) q[1];
sx q[1];
rz(1.5481987) q[1];
x q[2];
rz(3.0820936) q[3];
sx q[3];
rz(-1.7922641) q[3];
sx q[3];
rz(-0.49026107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1227293) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(1.8625721) q[2];
rz(3.0458798) q[3];
sx q[3];
rz(-2.0666104) q[3];
sx q[3];
rz(-1.3191282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.9886446) q[0];
sx q[0];
rz(-0.64866346) q[0];
sx q[0];
rz(2.1107819) q[0];
rz(-0.55054322) q[1];
sx q[1];
rz(-0.85854733) q[1];
sx q[1];
rz(1.8969089) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7787832) q[0];
sx q[0];
rz(-0.16805695) q[0];
sx q[0];
rz(-2.8166978) q[0];
rz(2.9739506) q[2];
sx q[2];
rz(-0.9603921) q[2];
sx q[2];
rz(0.33627015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.479949) q[1];
sx q[1];
rz(-2.2315761) q[1];
sx q[1];
rz(-0.26008545) q[1];
x q[2];
rz(-0.81826415) q[3];
sx q[3];
rz(-1.8570559) q[3];
sx q[3];
rz(2.8125826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2391669) q[2];
sx q[2];
rz(-2.5776165) q[2];
sx q[2];
rz(-1.9695367) q[2];
rz(-0.82550448) q[3];
sx q[3];
rz(-1.8768616) q[3];
sx q[3];
rz(-2.4765305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.091752) q[0];
sx q[0];
rz(-1.6850152) q[0];
sx q[0];
rz(0.41393429) q[0];
rz(0.44231689) q[1];
sx q[1];
rz(-2.1374173) q[1];
sx q[1];
rz(0.54534674) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2065679) q[0];
sx q[0];
rz(-0.27982084) q[0];
sx q[0];
rz(-0.18144515) q[0];
rz(0.6510066) q[2];
sx q[2];
rz(-0.88906139) q[2];
sx q[2];
rz(-1.663547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1283577) q[1];
sx q[1];
rz(-1.8473133) q[1];
sx q[1];
rz(-2.1091631) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0618938) q[3];
sx q[3];
rz(-1.7002704) q[3];
sx q[3];
rz(2.9751301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7004194) q[2];
sx q[2];
rz(-1.2607231) q[2];
sx q[2];
rz(0.1612266) q[2];
rz(-1.5876596) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(-2.5509295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081838354) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(-0.72342122) q[0];
rz(1.8266034) q[1];
sx q[1];
rz(-1.2774717) q[1];
sx q[1];
rz(-1.0378729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903477) q[0];
sx q[0];
rz(-1.6497166) q[0];
sx q[0];
rz(0.7235288) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60836253) q[2];
sx q[2];
rz(-2.3181097) q[2];
sx q[2];
rz(-1.057511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7241076) q[1];
sx q[1];
rz(-1.8383664) q[1];
sx q[1];
rz(2.6883283) q[1];
rz(-2.5955022) q[3];
sx q[3];
rz(-1.4442863) q[3];
sx q[3];
rz(1.5638592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7448685) q[2];
sx q[2];
rz(-1.1149656) q[2];
sx q[2];
rz(0.57207668) q[2];
rz(-2.5165596) q[3];
sx q[3];
rz(-2.1376938) q[3];
sx q[3];
rz(1.9151789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142757) q[0];
sx q[0];
rz(-0.033997424) q[0];
sx q[0];
rz(-0.52325621) q[0];
rz(1.8190067) q[1];
sx q[1];
rz(-1.4679642) q[1];
sx q[1];
rz(-2.5214213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70411982) q[0];
sx q[0];
rz(-1.0123555) q[0];
sx q[0];
rz(-2.0790786) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1441346) q[2];
sx q[2];
rz(-2.6687593) q[2];
sx q[2];
rz(2.0236058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.090724548) q[1];
sx q[1];
rz(-1.6584937) q[1];
sx q[1];
rz(0.94322738) q[1];
x q[2];
rz(-1.7460538) q[3];
sx q[3];
rz(-1.7850998) q[3];
sx q[3];
rz(-2.8068723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4124477) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(-1.8577417) q[2];
rz(-0.9681975) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(-1.5549972) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9667483) q[0];
sx q[0];
rz(-1.1356249) q[0];
sx q[0];
rz(-1.9298166) q[0];
rz(0.9696331) q[1];
sx q[1];
rz(-1.3215093) q[1];
sx q[1];
rz(-1.2744354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14718854) q[0];
sx q[0];
rz(-0.70941234) q[0];
sx q[0];
rz(-1.3377919) q[0];
x q[1];
rz(-0.095000141) q[2];
sx q[2];
rz(-2.3534905) q[2];
sx q[2];
rz(2.1301477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3962773) q[1];
sx q[1];
rz(-0.69967979) q[1];
sx q[1];
rz(-1.0902283) q[1];
rz(1.8936552) q[3];
sx q[3];
rz(-2.0104694) q[3];
sx q[3];
rz(-0.46931258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6756639) q[2];
sx q[2];
rz(-1.482654) q[2];
sx q[2];
rz(2.8181804) q[2];
rz(-3.1392426) q[3];
sx q[3];
rz(-1.4582783) q[3];
sx q[3];
rz(0.6215483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-1.3510652) q[1];
sx q[1];
rz(-1.8565145) q[1];
sx q[1];
rz(-1.2877119) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0072216) q[0];
sx q[0];
rz(-0.30728455) q[0];
sx q[0];
rz(1.1728806) q[0];
rz(-pi) q[1];
rz(-1.2399729) q[2];
sx q[2];
rz(-0.74218732) q[2];
sx q[2];
rz(-2.0798707) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8069597) q[1];
sx q[1];
rz(-1.2476842) q[1];
sx q[1];
rz(1.7503225) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.051750465) q[3];
sx q[3];
rz(-0.34933511) q[3];
sx q[3];
rz(0.83788315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52435851) q[2];
sx q[2];
rz(-1.989216) q[2];
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
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8096823) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(-1.5274973) q[0];
rz(-2.2612259) q[1];
sx q[1];
rz(-0.4207193) q[1];
sx q[1];
rz(-0.31164935) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0904734) q[0];
sx q[0];
rz(-1.3731806) q[0];
sx q[0];
rz(1.86905) q[0];
rz(-pi) q[1];
rz(-1.6402354) q[2];
sx q[2];
rz(-2.9096894) q[2];
sx q[2];
rz(2.7551485) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2276409) q[1];
sx q[1];
rz(-2.3415509) q[1];
sx q[1];
rz(0.52180078) q[1];
rz(-pi) q[2];
rz(1.1470471) q[3];
sx q[3];
rz(-2.4422788) q[3];
sx q[3];
rz(-0.55845234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46939048) q[2];
sx q[2];
rz(-2.2180874) q[2];
sx q[2];
rz(1.7507318) q[2];
rz(1.6145128) q[3];
sx q[3];
rz(-1.047784) q[3];
sx q[3];
rz(-1.0745777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020318) q[0];
sx q[0];
rz(-1.9968411) q[0];
sx q[0];
rz(2.5290053) q[0];
rz(2.1514429) q[1];
sx q[1];
rz(-1.9691111) q[1];
sx q[1];
rz(1.6506857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6698786) q[0];
sx q[0];
rz(-1.151514) q[0];
sx q[0];
rz(0.74434481) q[0];
x q[1];
rz(-1.8950997) q[2];
sx q[2];
rz(-0.68585472) q[2];
sx q[2];
rz(-1.7795479) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4348169) q[1];
sx q[1];
rz(-1.56456) q[1];
sx q[1];
rz(-1.732279) q[1];
rz(-pi) q[2];
rz(0.69802876) q[3];
sx q[3];
rz(-1.7637611) q[3];
sx q[3];
rz(-2.2730227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.426173) q[2];
sx q[2];
rz(-0.82948589) q[2];
sx q[2];
rz(-3.1066011) q[2];
rz(1.70111) q[3];
sx q[3];
rz(-2.248843) q[3];
sx q[3];
rz(0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424425) q[0];
sx q[0];
rz(-2.2461666) q[0];
sx q[0];
rz(2.9641892) q[0];
rz(1.198248) q[1];
sx q[1];
rz(-1.5180963) q[1];
sx q[1];
rz(0.95388283) q[1];
rz(2.4468471) q[2];
sx q[2];
rz(-0.80949819) q[2];
sx q[2];
rz(-0.037485952) q[2];
rz(2.421834) q[3];
sx q[3];
rz(-1.2133141) q[3];
sx q[3];
rz(2.6837466) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
