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
rz(-0.35204044) q[0];
sx q[0];
rz(-0.8249324) q[0];
sx q[0];
rz(0.53476778) q[0];
rz(0.98400247) q[1];
sx q[1];
rz(-2.6360631) q[1];
sx q[1];
rz(1.2619789) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0540491) q[0];
sx q[0];
rz(-0.66960483) q[0];
sx q[0];
rz(2.7101507) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2625669) q[2];
sx q[2];
rz(-0.75415416) q[2];
sx q[2];
rz(-1.1471495) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83916503) q[1];
sx q[1];
rz(-1.622597) q[1];
sx q[1];
rz(0.54755819) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7746968) q[3];
sx q[3];
rz(-1.9528021) q[3];
sx q[3];
rz(-2.6873858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24903211) q[2];
sx q[2];
rz(-2.6292215) q[2];
sx q[2];
rz(1.6585635) q[2];
rz(0.28462166) q[3];
sx q[3];
rz(-0.62323815) q[3];
sx q[3];
rz(1.0641789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738709) q[0];
sx q[0];
rz(-0.72672788) q[0];
sx q[0];
rz(3.100585) q[0];
rz(1.2076591) q[1];
sx q[1];
rz(-2.9137847) q[1];
sx q[1];
rz(-1.4606732) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2940833) q[0];
sx q[0];
rz(-2.9929711) q[0];
sx q[0];
rz(2.9411208) q[0];
x q[1];
rz(2.5778077) q[2];
sx q[2];
rz(-0.85104686) q[2];
sx q[2];
rz(-1.189718) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8833958) q[1];
sx q[1];
rz(-2.8821695) q[1];
sx q[1];
rz(-1.4197465) q[1];
rz(-pi) q[2];
rz(-1.3364661) q[3];
sx q[3];
rz(-0.70584345) q[3];
sx q[3];
rz(2.2587551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4760806) q[2];
sx q[2];
rz(-1.4613232) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(-0.0811854) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(1.7056874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3944893) q[0];
sx q[0];
rz(-3.1293226) q[0];
sx q[0];
rz(0.61007208) q[0];
rz(-1.8286145) q[1];
sx q[1];
rz(-2.4459631) q[1];
sx q[1];
rz(-2.5909766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75567882) q[0];
sx q[0];
rz(-2.3649594) q[0];
sx q[0];
rz(2.617111) q[0];
x q[1];
rz(-1.2989013) q[2];
sx q[2];
rz(-1.8823132) q[2];
sx q[2];
rz(-2.2443983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0455835) q[1];
sx q[1];
rz(-0.35638816) q[1];
sx q[1];
rz(2.136904) q[1];
rz(2.315963) q[3];
sx q[3];
rz(-1.6761002) q[3];
sx q[3];
rz(-1.8838175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6637471) q[2];
sx q[2];
rz(-2.9691594) q[2];
sx q[2];
rz(1.0663859) q[2];
rz(1.1586698) q[3];
sx q[3];
rz(-1.3830769) q[3];
sx q[3];
rz(-3.0995479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8313507) q[0];
sx q[0];
rz(-2.5284335) q[0];
sx q[0];
rz(-0.9170652) q[0];
rz(0.65545583) q[1];
sx q[1];
rz(-2.2323445) q[1];
sx q[1];
rz(2.7520032) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.675466) q[0];
sx q[0];
rz(-1.6524757) q[0];
sx q[0];
rz(1.3838751) q[0];
rz(2.5823103) q[2];
sx q[2];
rz(-2.1765095) q[2];
sx q[2];
rz(-2.9756211) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0181737) q[1];
sx q[1];
rz(-1.5416073) q[1];
sx q[1];
rz(1.1694714) q[1];
rz(-pi) q[2];
rz(0.97352685) q[3];
sx q[3];
rz(-1.105666) q[3];
sx q[3];
rz(-1.2965508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83170825) q[2];
sx q[2];
rz(-0.94253057) q[2];
sx q[2];
rz(1.7892828) q[2];
rz(0.74812198) q[3];
sx q[3];
rz(-1.3040521) q[3];
sx q[3];
rz(1.6149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378976) q[0];
sx q[0];
rz(-2.5717773) q[0];
sx q[0];
rz(1.3414398) q[0];
rz(-0.98126423) q[1];
sx q[1];
rz(-1.2813247) q[1];
sx q[1];
rz(0.70995465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34559743) q[0];
sx q[0];
rz(-1.5516722) q[0];
sx q[0];
rz(-0.32581331) q[0];
x q[1];
rz(-0.99468987) q[2];
sx q[2];
rz(-1.7719977) q[2];
sx q[2];
rz(-2.3735769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4932351) q[1];
sx q[1];
rz(-1.7410189) q[1];
sx q[1];
rz(1.3551718) q[1];
rz(-pi) q[2];
rz(-0.59779915) q[3];
sx q[3];
rz(-1.4922659) q[3];
sx q[3];
rz(0.30818916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95825163) q[2];
sx q[2];
rz(-0.7603344) q[2];
sx q[2];
rz(2.23488) q[2];
rz(0.18236154) q[3];
sx q[3];
rz(-1.7014818) q[3];
sx q[3];
rz(0.85842925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3458503) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(2.9398651) q[0];
rz(1.7851104) q[1];
sx q[1];
rz(-0.67142612) q[1];
sx q[1];
rz(1.1511525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0990382) q[0];
sx q[0];
rz(-1.7793613) q[0];
sx q[0];
rz(-1.1369571) q[0];
x q[1];
rz(-0.62080748) q[2];
sx q[2];
rz(-2.4352031) q[2];
sx q[2];
rz(2.8548129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56111873) q[1];
sx q[1];
rz(-0.72787933) q[1];
sx q[1];
rz(-2.5125458) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51064193) q[3];
sx q[3];
rz(-1.770442) q[3];
sx q[3];
rz(-2.8996403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98728937) q[2];
sx q[2];
rz(-0.58649784) q[2];
sx q[2];
rz(-2.7941373) q[2];
rz(-0.82184982) q[3];
sx q[3];
rz(-1.776266) q[3];
sx q[3];
rz(-1.52389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8007941) q[0];
sx q[0];
rz(-0.64685416) q[0];
sx q[0];
rz(-1.1232173) q[0];
rz(-0.42876354) q[1];
sx q[1];
rz(-0.88974297) q[1];
sx q[1];
rz(-2.9926328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4377047) q[0];
sx q[0];
rz(-2.0876679) q[0];
sx q[0];
rz(-0.025512841) q[0];
rz(-pi) q[1];
rz(-1.3690022) q[2];
sx q[2];
rz(-1.861146) q[2];
sx q[2];
rz(1.5767911) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.625714) q[1];
sx q[1];
rz(-0.82059723) q[1];
sx q[1];
rz(1.3681812) q[1];
rz(1.4895053) q[3];
sx q[3];
rz(-0.40697655) q[3];
sx q[3];
rz(1.6124042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.011977) q[2];
sx q[2];
rz(-0.35795438) q[2];
sx q[2];
rz(0.32148662) q[2];
rz(1.1164411) q[3];
sx q[3];
rz(-0.8588841) q[3];
sx q[3];
rz(1.727625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9880992) q[0];
sx q[0];
rz(-1.3445925) q[0];
sx q[0];
rz(0.022911428) q[0];
rz(-1.5494391) q[1];
sx q[1];
rz(-2.0982845) q[1];
sx q[1];
rz(-1.9291838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1003636) q[0];
sx q[0];
rz(-2.1080906) q[0];
sx q[0];
rz(0.30502086) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7801082) q[2];
sx q[2];
rz(-0.87280849) q[2];
sx q[2];
rz(0.43025508) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0689736) q[1];
sx q[1];
rz(-0.5104465) q[1];
sx q[1];
rz(-2.4859395) q[1];
rz(-pi) q[2];
rz(-1.6818265) q[3];
sx q[3];
rz(-0.71655203) q[3];
sx q[3];
rz(-0.46886231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.93078485) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(-2.1442294) q[2];
rz(-1.8831683) q[3];
sx q[3];
rz(-1.1910028) q[3];
sx q[3];
rz(-1.9112126) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9835984) q[0];
sx q[0];
rz(-2.4857434) q[0];
sx q[0];
rz(-2.6672145) q[0];
rz(2.5899218) q[1];
sx q[1];
rz(-2.3730979) q[1];
sx q[1];
rz(2.4840568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6303355) q[0];
sx q[0];
rz(-1.5310107) q[0];
sx q[0];
rz(1.5471094) q[0];
rz(1.1908866) q[2];
sx q[2];
rz(-0.43336855) q[2];
sx q[2];
rz(-2.4182408) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2567946) q[1];
sx q[1];
rz(-1.9391283) q[1];
sx q[1];
rz(-0.017466768) q[1];
rz(-0.94352888) q[3];
sx q[3];
rz(-2.9637058) q[3];
sx q[3];
rz(-1.8454144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9115596) q[2];
sx q[2];
rz(-2.7342789) q[2];
sx q[2];
rz(1.8302906) q[2];
rz(0.71410549) q[3];
sx q[3];
rz(-1.2823391) q[3];
sx q[3];
rz(-0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0062200935) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(2.5332992) q[0];
rz(2.3237806) q[1];
sx q[1];
rz(-1.9655971) q[1];
sx q[1];
rz(-2.3086595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9271256) q[0];
sx q[0];
rz(-1.1499078) q[0];
sx q[0];
rz(1.6273999) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67939679) q[2];
sx q[2];
rz(-2.7459807) q[2];
sx q[2];
rz(-2.6406956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74216182) q[1];
sx q[1];
rz(-0.99282904) q[1];
sx q[1];
rz(-1.7294243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6021483) q[3];
sx q[3];
rz(-2.1636181) q[3];
sx q[3];
rz(-1.624799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1828764) q[2];
sx q[2];
rz(-0.35280886) q[2];
sx q[2];
rz(2.555441) q[2];
rz(-0.90306774) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(3.0557475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9027973) q[0];
sx q[0];
rz(-1.1302523) q[0];
sx q[0];
rz(-0.94373066) q[0];
rz(2.0441652) q[1];
sx q[1];
rz(-2.8891017) q[1];
sx q[1];
rz(-0.15695922) q[1];
rz(-1.9465994) q[2];
sx q[2];
rz(-2.7766418) q[2];
sx q[2];
rz(3.1212213) q[2];
rz(-2.1395352) q[3];
sx q[3];
rz(-2.3593223) q[3];
sx q[3];
rz(0.13025688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
