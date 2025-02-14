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
rz(-1.2404233) q[0];
sx q[0];
rz(1.9440396) q[0];
sx q[0];
rz(9.64111) q[0];
rz(0.95739111) q[1];
sx q[1];
rz(-2.4654145) q[1];
sx q[1];
rz(-1.6300936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0974329) q[0];
sx q[0];
rz(-0.55632797) q[0];
sx q[0];
rz(-3.1233643) q[0];
x q[1];
rz(0.58403973) q[2];
sx q[2];
rz(-0.55070832) q[2];
sx q[2];
rz(-0.60930291) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2701157) q[1];
sx q[1];
rz(-1.2387215) q[1];
sx q[1];
rz(2.6144652) q[1];
x q[2];
rz(-0.86801784) q[3];
sx q[3];
rz(-1.2642164) q[3];
sx q[3];
rz(2.0748695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1067074) q[2];
sx q[2];
rz(-2.7896176) q[2];
sx q[2];
rz(-2.0931639) q[2];
rz(-0.18167051) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(2.488193) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.492391) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(-2.6948068) q[0];
rz(-0.88042879) q[1];
sx q[1];
rz(-1.3648405) q[1];
sx q[1];
rz(-0.78278881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370011) q[0];
sx q[0];
rz(-0.25213045) q[0];
sx q[0];
rz(1.3229516) q[0];
rz(-pi) q[1];
rz(1.9982463) q[2];
sx q[2];
rz(-0.99787092) q[2];
sx q[2];
rz(-1.7597511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4383653) q[1];
sx q[1];
rz(-1.2709054) q[1];
sx q[1];
rz(-0.48734003) q[1];
rz(-2.6683025) q[3];
sx q[3];
rz(-2.1968699) q[3];
sx q[3];
rz(-1.2155217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.11144) q[2];
sx q[2];
rz(-1.6609265) q[2];
sx q[2];
rz(2.1622369) q[2];
rz(-2.7566946) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(0.16429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6556743) q[0];
sx q[0];
rz(-0.05412183) q[0];
sx q[0];
rz(-0.78980494) q[0];
rz(-2.9549331) q[1];
sx q[1];
rz(-1.4013545) q[1];
sx q[1];
rz(-2.1479215) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06834767) q[0];
sx q[0];
rz(-0.6093502) q[0];
sx q[0];
rz(-2.6960899) q[0];
x q[1];
rz(-0.31618677) q[2];
sx q[2];
rz(-0.63717604) q[2];
sx q[2];
rz(0.030046163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6911654) q[1];
sx q[1];
rz(-1.9948927) q[1];
sx q[1];
rz(1.0662088) q[1];
x q[2];
rz(0.30672725) q[3];
sx q[3];
rz(-2.3969458) q[3];
sx q[3];
rz(0.70017231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32187244) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(-2.7703088) q[2];
rz(0.34879455) q[3];
sx q[3];
rz(-2.0472417) q[3];
sx q[3];
rz(-0.5955407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45469859) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(-1.7373079) q[0];
rz(2.4644201) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(-1.4926532) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66144511) q[0];
sx q[0];
rz(-0.3269402) q[0];
sx q[0];
rz(2.4133555) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8437496) q[2];
sx q[2];
rz(-1.1561511) q[2];
sx q[2];
rz(-2.9306102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75912913) q[1];
sx q[1];
rz(-1.8378403) q[1];
sx q[1];
rz(1.3401396) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31217137) q[3];
sx q[3];
rz(-2.2709284) q[3];
sx q[3];
rz(2.2362102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6877785) q[2];
sx q[2];
rz(-2.1731589) q[2];
sx q[2];
rz(-2.7933534) q[2];
rz(1.6866775) q[3];
sx q[3];
rz(-1.7125407) q[3];
sx q[3];
rz(-1.7961563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763879) q[0];
sx q[0];
rz(-0.56469733) q[0];
sx q[0];
rz(-2.2494466) q[0];
rz(2.6773894) q[1];
sx q[1];
rz(-1.8966388) q[1];
sx q[1];
rz(-1.2947327) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2703122) q[0];
sx q[0];
rz(-2.0593606) q[0];
sx q[0];
rz(-2.110092) q[0];
x q[1];
rz(1.6991529) q[2];
sx q[2];
rz(-2.1048628) q[2];
sx q[2];
rz(-2.4160224) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.059493493) q[1];
sx q[1];
rz(-1.1319185) q[1];
sx q[1];
rz(-2.782269) q[1];
x q[2];
rz(2.2512359) q[3];
sx q[3];
rz(-2.0687752) q[3];
sx q[3];
rz(0.80163664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(2.5757705) q[2];
rz(-3.0602509) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.464798) q[0];
sx q[0];
rz(-0.066635266) q[0];
sx q[0];
rz(-1.5555405) q[0];
rz(-1.0606891) q[1];
sx q[1];
rz(-1.5763177) q[1];
sx q[1];
rz(2.5097844) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2151129) q[0];
sx q[0];
rz(-2.8936912) q[0];
sx q[0];
rz(-0.32899022) q[0];
rz(-pi) q[1];
rz(-0.49922322) q[2];
sx q[2];
rz(-0.92485917) q[2];
sx q[2];
rz(1.8725841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0630707) q[1];
sx q[1];
rz(-2.1025189) q[1];
sx q[1];
rz(-0.39549455) q[1];
rz(-pi) q[2];
rz(1.7177714) q[3];
sx q[3];
rz(-2.2490361) q[3];
sx q[3];
rz(3.105046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1598728) q[2];
sx q[2];
rz(-2.0231415) q[2];
sx q[2];
rz(0.49883207) q[2];
rz(-1.8286797) q[3];
sx q[3];
rz(-0.73176089) q[3];
sx q[3];
rz(-1.7074728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071534) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(-0.69951192) q[0];
rz(2.6761159) q[1];
sx q[1];
rz(-2.270348) q[1];
sx q[1];
rz(2.2043998) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4995183) q[0];
sx q[0];
rz(-2.8694911) q[0];
sx q[0];
rz(-2.4689552) q[0];
x q[1];
rz(1.6752536) q[2];
sx q[2];
rz(-1.3323931) q[2];
sx q[2];
rz(-1.6676211) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22058567) q[1];
sx q[1];
rz(-1.3594419) q[1];
sx q[1];
rz(1.2799954) q[1];
rz(0.77009691) q[3];
sx q[3];
rz(-2.7316964) q[3];
sx q[3];
rz(0.65576762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.297544) q[2];
sx q[2];
rz(-1.8397477) q[2];
sx q[2];
rz(-2.76827) q[2];
rz(-1.1897872) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(-2.6313307) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74476403) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(-2.3936791) q[0];
rz(0.76639908) q[1];
sx q[1];
rz(-2.8728569) q[1];
sx q[1];
rz(-3.1386197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2073313) q[0];
sx q[0];
rz(-2.2483279) q[0];
sx q[0];
rz(0.24768655) q[0];
rz(-pi) q[1];
rz(0.019549088) q[2];
sx q[2];
rz(-1.4998933) q[2];
sx q[2];
rz(0.30724684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.06719477) q[1];
sx q[1];
rz(-1.41681) q[1];
sx q[1];
rz(0.043937307) q[1];
rz(0.19488867) q[3];
sx q[3];
rz(-1.3759383) q[3];
sx q[3];
rz(3.0611567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.030674) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(2.0745011) q[2];
rz(-3.057632) q[3];
sx q[3];
rz(-2.6559918) q[3];
sx q[3];
rz(0.77897227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344051) q[0];
sx q[0];
rz(-2.2028956) q[0];
sx q[0];
rz(3.0352266) q[0];
rz(-2.1616409) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(-2.3756557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04127114) q[0];
sx q[0];
rz(-1.0491956) q[0];
sx q[0];
rz(1.0688416) q[0];
rz(0.19975234) q[2];
sx q[2];
rz(-1.6589266) q[2];
sx q[2];
rz(-0.55921184) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83726604) q[1];
sx q[1];
rz(-2.2464753) q[1];
sx q[1];
rz(1.5262239) q[1];
rz(-pi) q[2];
rz(-3.1274904) q[3];
sx q[3];
rz(-1.5745134) q[3];
sx q[3];
rz(-0.93059082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6672259) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(0.70029798) q[2];
rz(-2.4750366) q[3];
sx q[3];
rz(-1.8066112) q[3];
sx q[3];
rz(-2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.67548442) q[0];
sx q[0];
rz(-1.0270783) q[0];
sx q[0];
rz(-1.3960557) q[0];
rz(-1.667977) q[1];
sx q[1];
rz(-1.2584078) q[1];
sx q[1];
rz(2.4748763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5294801) q[0];
sx q[0];
rz(-2.1938938) q[0];
sx q[0];
rz(-2.3701131) q[0];
rz(-1.2856977) q[2];
sx q[2];
rz(-2.210493) q[2];
sx q[2];
rz(2.6486047) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4426431) q[1];
sx q[1];
rz(-1.8869699) q[1];
sx q[1];
rz(1.7085646) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75717775) q[3];
sx q[3];
rz(-0.72523967) q[3];
sx q[3];
rz(-2.4717769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.066976808) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(-2.2508049) q[2];
rz(2.7095419) q[3];
sx q[3];
rz(-2.0712349) q[3];
sx q[3];
rz(2.3945358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4074832) q[0];
sx q[0];
rz(-1.643184) q[0];
sx q[0];
rz(-1.2947422) q[0];
rz(1.2474077) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(1.1449849) q[2];
sx q[2];
rz(-2.7663284) q[2];
sx q[2];
rz(-0.13079499) q[2];
rz(-2.6525146) q[3];
sx q[3];
rz(-1.5413956) q[3];
sx q[3];
rz(1.2914381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
