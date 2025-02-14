OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.85668808) q[0];
sx q[0];
rz(4.9520725) q[0];
sx q[0];
rz(7.0926275) q[0];
rz(0.59960214) q[1];
sx q[1];
rz(-2.0166346) q[1];
sx q[1];
rz(-0.52176276) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4057465) q[0];
sx q[0];
rz(-3.0015115) q[0];
sx q[0];
rz(-1.5625885) q[0];
rz(-pi) q[1];
rz(-0.56057616) q[2];
sx q[2];
rz(-1.8235221) q[2];
sx q[2];
rz(-1.8552903) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1612137) q[1];
sx q[1];
rz(-1.8424826) q[1];
sx q[1];
rz(1.1030688) q[1];
rz(-pi) q[2];
rz(0.14484804) q[3];
sx q[3];
rz(-2.600935) q[3];
sx q[3];
rz(0.40504211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1521505) q[2];
sx q[2];
rz(-2.4336954) q[2];
sx q[2];
rz(-1.36261) q[2];
rz(-1.4710434) q[3];
sx q[3];
rz(-1.5444642) q[3];
sx q[3];
rz(-0.41489261) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7489557) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(1.0697399) q[0];
rz(2.3906129) q[1];
sx q[1];
rz(-2.2396478) q[1];
sx q[1];
rz(-0.78838563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1845368) q[0];
sx q[0];
rz(-0.95648492) q[0];
sx q[0];
rz(1.3411103) q[0];
rz(0.090615409) q[2];
sx q[2];
rz(-0.66747016) q[2];
sx q[2];
rz(-0.27809503) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0505817) q[1];
sx q[1];
rz(-1.7643223) q[1];
sx q[1];
rz(-2.0975608) q[1];
rz(1.1343003) q[3];
sx q[3];
rz(-2.0406764) q[3];
sx q[3];
rz(-1.3677471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.02701935) q[2];
sx q[2];
rz(-0.35832778) q[2];
sx q[2];
rz(2.2349854) q[2];
rz(-2.1410227) q[3];
sx q[3];
rz(-2.0228736) q[3];
sx q[3];
rz(0.51893836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8949316) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(-1.3038127) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.3273393) q[1];
sx q[1];
rz(-1.7283641) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6400063) q[0];
sx q[0];
rz(-1.4297856) q[0];
sx q[0];
rz(1.6580576) q[0];
x q[1];
rz(-0.5858461) q[2];
sx q[2];
rz(-0.72524643) q[2];
sx q[2];
rz(0.21373978) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6709137) q[1];
sx q[1];
rz(-1.7470198) q[1];
sx q[1];
rz(2.6754968) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.71536) q[3];
sx q[3];
rz(-2.4213397) q[3];
sx q[3];
rz(1.4831869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7580737) q[2];
sx q[2];
rz(-1.588409) q[2];
sx q[2];
rz(-3.0124532) q[2];
rz(-2.6376851) q[3];
sx q[3];
rz(-1.4174856) q[3];
sx q[3];
rz(2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1103766) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(2.2953798) q[0];
rz(-2.5695678) q[1];
sx q[1];
rz(-1.5049479) q[1];
sx q[1];
rz(-1.9073073) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4847149) q[0];
sx q[0];
rz(-1.6678255) q[0];
sx q[0];
rz(-0.33330864) q[0];
x q[1];
rz(-2.0300806) q[2];
sx q[2];
rz(-1.3142576) q[2];
sx q[2];
rz(0.81748325) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9920791) q[1];
sx q[1];
rz(-1.7798335) q[1];
sx q[1];
rz(1.4884218) q[1];
rz(-2.0999072) q[3];
sx q[3];
rz(-2.2941781) q[3];
sx q[3];
rz(0.46284562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7531551) q[2];
sx q[2];
rz(-0.93418241) q[2];
sx q[2];
rz(-1.6969121) q[2];
rz(-1.8709024) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-1.1893893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.5020849) q[0];
sx q[0];
rz(-1.4690228) q[0];
sx q[0];
rz(2.0462659) q[0];
rz(-3.0785839) q[1];
sx q[1];
rz(-1.4728225) q[1];
sx q[1];
rz(-0.94620401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.584884) q[0];
sx q[0];
rz(-2.1153304) q[0];
sx q[0];
rz(-0.03277112) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49503742) q[2];
sx q[2];
rz(-2.4433141) q[2];
sx q[2];
rz(-2.3586522) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19095382) q[1];
sx q[1];
rz(-0.51114782) q[1];
sx q[1];
rz(0.010015476) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0466411) q[3];
sx q[3];
rz(-0.44455179) q[3];
sx q[3];
rz(1.9421645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5064064) q[2];
sx q[2];
rz(-2.0621641) q[2];
sx q[2];
rz(1.5080473) q[2];
rz(-2.9888198) q[3];
sx q[3];
rz(-1.2341876) q[3];
sx q[3];
rz(2.2085786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372811) q[0];
sx q[0];
rz(-2.5560684) q[0];
sx q[0];
rz(3.0329419) q[0];
rz(3.0142504) q[1];
sx q[1];
rz(-1.8904949) q[1];
sx q[1];
rz(-2.2551575) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21289794) q[0];
sx q[0];
rz(-1.5542287) q[0];
sx q[0];
rz(2.6221656) q[0];
rz(-pi) q[1];
rz(-0.059659307) q[2];
sx q[2];
rz(-1.3666412) q[2];
sx q[2];
rz(-2.8722389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8005714) q[1];
sx q[1];
rz(-1.6280975) q[1];
sx q[1];
rz(-1.7898466) q[1];
rz(0.57678595) q[3];
sx q[3];
rz(-0.20142074) q[3];
sx q[3];
rz(0.52514986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.844937) q[2];
sx q[2];
rz(-2.3462494) q[2];
sx q[2];
rz(-2.8080158) q[2];
rz(0.9345471) q[3];
sx q[3];
rz(-2.3881674) q[3];
sx q[3];
rz(0.88753382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6880671) q[0];
sx q[0];
rz(-1.4997361) q[0];
sx q[0];
rz(0.30853477) q[0];
rz(-0.60641369) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(0.44953129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21643695) q[0];
sx q[0];
rz(-1.0569658) q[0];
sx q[0];
rz(-1.4975182) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3541127) q[2];
sx q[2];
rz(-1.9132987) q[2];
sx q[2];
rz(-2.1975348) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3492185) q[1];
sx q[1];
rz(-1.4687531) q[1];
sx q[1];
rz(0.95178638) q[1];
rz(-pi) q[2];
rz(-0.87089296) q[3];
sx q[3];
rz(-1.1364391) q[3];
sx q[3];
rz(0.30444452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.061607925) q[2];
sx q[2];
rz(-1.4915165) q[2];
sx q[2];
rz(1.1770581) q[2];
rz(0.70147771) q[3];
sx q[3];
rz(-2.8886815) q[3];
sx q[3];
rz(0.85566068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45706448) q[0];
sx q[0];
rz(-2.6236911) q[0];
sx q[0];
rz(-1.49217) q[0];
rz(1.0555142) q[1];
sx q[1];
rz(-1.5857453) q[1];
sx q[1];
rz(2.7141056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.939291) q[0];
sx q[0];
rz(-2.2965676) q[0];
sx q[0];
rz(0.85406749) q[0];
rz(-pi) q[1];
rz(1.7321109) q[2];
sx q[2];
rz(-1.9967134) q[2];
sx q[2];
rz(-2.6420455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52054469) q[1];
sx q[1];
rz(-1.8836023) q[1];
sx q[1];
rz(-1.5473844) q[1];
rz(-2.4343723) q[3];
sx q[3];
rz(-2.7364991) q[3];
sx q[3];
rz(1.9735379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18275729) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(1.2048362) q[2];
rz(0.13293535) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(-1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28250113) q[0];
sx q[0];
rz(-1.4097255) q[0];
sx q[0];
rz(2.6891563) q[0];
rz(2.2826507) q[1];
sx q[1];
rz(-0.57933885) q[1];
sx q[1];
rz(-1.6890242) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7410788) q[0];
sx q[0];
rz(-2.3483334) q[0];
sx q[0];
rz(2.443072) q[0];
x q[1];
rz(2.8088687) q[2];
sx q[2];
rz(-1.0411658) q[2];
sx q[2];
rz(0.67654787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13551417) q[1];
sx q[1];
rz(-2.0767414) q[1];
sx q[1];
rz(-3.1385026) q[1];
rz(-pi) q[2];
rz(-0.36859103) q[3];
sx q[3];
rz(-1.2342216) q[3];
sx q[3];
rz(-2.9934538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9026044) q[2];
sx q[2];
rz(-1.3024104) q[2];
sx q[2];
rz(-2.7461309) q[2];
rz(-2.0579193) q[3];
sx q[3];
rz(-0.62267059) q[3];
sx q[3];
rz(-1.7284988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.285242) q[0];
sx q[0];
rz(-0.1150035) q[0];
sx q[0];
rz(-1.2913936) q[0];
rz(2.3639288) q[1];
sx q[1];
rz(-1.5823369) q[1];
sx q[1];
rz(1.2819598) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1384821) q[0];
sx q[0];
rz(-1.7958475) q[0];
sx q[0];
rz(2.9182994) q[0];
rz(-pi) q[1];
rz(0.28996946) q[2];
sx q[2];
rz(-2.4735138) q[2];
sx q[2];
rz(-1.5886605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.229515) q[1];
sx q[1];
rz(-1.34206) q[1];
sx q[1];
rz(-2.5994615) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7601717) q[3];
sx q[3];
rz(-2.6145589) q[3];
sx q[3];
rz(2.1552212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.164244) q[2];
sx q[2];
rz(-0.35173309) q[2];
sx q[2];
rz(-2.742761) q[2];
rz(2.2073958) q[3];
sx q[3];
rz(-0.7312921) q[3];
sx q[3];
rz(-3.0493693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6964523) q[0];
sx q[0];
rz(-2.2284989) q[0];
sx q[0];
rz(-3.1365119) q[0];
rz(-2.483881) q[1];
sx q[1];
rz(-1.4124159) q[1];
sx q[1];
rz(1.1884069) q[1];
rz(-0.80398139) q[2];
sx q[2];
rz(-1.3406499) q[2];
sx q[2];
rz(-1.7743877) q[2];
rz(-3.0791238) q[3];
sx q[3];
rz(-2.1556839) q[3];
sx q[3];
rz(1.5869303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
