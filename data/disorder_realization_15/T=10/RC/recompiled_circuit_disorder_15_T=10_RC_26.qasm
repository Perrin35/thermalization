OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29785922) q[0];
sx q[0];
rz(-2.5279186) q[0];
sx q[0];
rz(2.4224129) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52553015) q[0];
sx q[0];
rz(-1.5126192) q[0];
sx q[0];
rz(1.9921897) q[0];
rz(-2.9814331) q[2];
sx q[2];
rz(-0.9848435) q[2];
sx q[2];
rz(1.9507267) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72585427) q[1];
sx q[1];
rz(-1.1032747) q[1];
sx q[1];
rz(-2.095925) q[1];
x q[2];
rz(-3.1036166) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(1.9248885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44895479) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(0.56837481) q[2];
rz(-1.9521936) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73873591) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(0.96827132) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(1.2526858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68080901) q[0];
sx q[0];
rz(-1.660581) q[0];
sx q[0];
rz(0.16016527) q[0];
rz(-2.5622796) q[2];
sx q[2];
rz(-2.1462153) q[2];
sx q[2];
rz(-0.79798165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3167047) q[1];
sx q[1];
rz(-1.0893634) q[1];
sx q[1];
rz(-0.96932051) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1329123) q[3];
sx q[3];
rz(-1.7403733) q[3];
sx q[3];
rz(1.358043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.12053717) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(-0.34376124) q[2];
rz(-0.3343285) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(-3.1345471) q[0];
rz(-2.7650611) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(-0.71281707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5755771) q[0];
sx q[0];
rz(-0.30656719) q[0];
sx q[0];
rz(2.3335639) q[0];
x q[1];
rz(0.17194925) q[2];
sx q[2];
rz(-1.9438582) q[2];
sx q[2];
rz(-2.3200777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30046001) q[1];
sx q[1];
rz(-1.3679879) q[1];
sx q[1];
rz(2.1122785) q[1];
rz(1.4180693) q[3];
sx q[3];
rz(-2.6543791) q[3];
sx q[3];
rz(-1.0917851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(-0.66398579) q[2];
rz(-1.239423) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-0.57408875) q[0];
rz(0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(-2.2132197) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5260122) q[0];
sx q[0];
rz(-2.2768339) q[0];
sx q[0];
rz(-3.1403149) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33889126) q[2];
sx q[2];
rz(-2.3232984) q[2];
sx q[2];
rz(1.3202867) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1410364) q[1];
sx q[1];
rz(-0.49266854) q[1];
sx q[1];
rz(-1.0286742) q[1];
x q[2];
rz(2.1888234) q[3];
sx q[3];
rz(-2.3828155) q[3];
sx q[3];
rz(-3.0892059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0488247) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.9862004) q[2];
rz(0.14285764) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(-0.75240451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.083754152) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(0.074247867) q[0];
rz(1.4986562) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-0.9517076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041391011) q[0];
sx q[0];
rz(-0.75267422) q[0];
sx q[0];
rz(1.676883) q[0];
rz(-pi) q[1];
rz(-2.1553467) q[2];
sx q[2];
rz(-0.84976053) q[2];
sx q[2];
rz(0.69794387) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4068027) q[1];
sx q[1];
rz(-2.4996335) q[1];
sx q[1];
rz(1.8491247) q[1];
rz(-2.1577155) q[3];
sx q[3];
rz(-0.65423274) q[3];
sx q[3];
rz(-2.4901842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64530659) q[2];
sx q[2];
rz(-0.22660613) q[2];
sx q[2];
rz(-1.1523694) q[2];
rz(1.0460098) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(-0.0013105198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-0.9933269) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(-2.6519725) q[0];
rz(-1.024225) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(-1.6061868) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09774694) q[0];
sx q[0];
rz(-1.6698208) q[0];
sx q[0];
rz(1.6188341) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68131955) q[2];
sx q[2];
rz(-2.1207003) q[2];
sx q[2];
rz(-0.51598179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9163497) q[1];
sx q[1];
rz(-0.32685977) q[1];
sx q[1];
rz(-2.0845695) q[1];
rz(-pi) q[2];
rz(-0.225004) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(-1.0596421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17807047) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(0.17573389) q[2];
rz(-1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068129383) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(-1.2699132) q[0];
rz(2.9636256) q[1];
sx q[1];
rz(-2.3033419) q[1];
sx q[1];
rz(1.3659182) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458347) q[0];
sx q[0];
rz(-2.2582044) q[0];
sx q[0];
rz(1.13152) q[0];
x q[1];
rz(2.5125011) q[2];
sx q[2];
rz(-2.7436939) q[2];
sx q[2];
rz(1.2049904) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7682225) q[1];
sx q[1];
rz(-1.3942413) q[1];
sx q[1];
rz(-2.8819487) q[1];
rz(2.6122983) q[3];
sx q[3];
rz(-1.4636453) q[3];
sx q[3];
rz(-1.2308434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0697249) q[2];
sx q[2];
rz(-0.29399997) q[2];
sx q[2];
rz(-2.243637) q[2];
rz(0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(0.34415054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9550069) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(2.8198077) q[0];
rz(2.2161662) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(-0.61703533) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0641545) q[0];
sx q[0];
rz(-1.7700717) q[0];
sx q[0];
rz(-1.827924) q[0];
rz(-1.7865137) q[2];
sx q[2];
rz(-1.5383913) q[2];
sx q[2];
rz(1.0499133) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.9279234) q[1];
sx q[1];
rz(-2.546715) q[1];
sx q[1];
rz(2.2713714) q[1];
rz(-pi) q[2];
rz(1.6156322) q[3];
sx q[3];
rz(-0.82190824) q[3];
sx q[3];
rz(-0.84079784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59163219) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(-0.11403306) q[2];
rz(-0.36241254) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(-1.9053649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(1.3569008) q[0];
rz(0.82018954) q[1];
sx q[1];
rz(-2.7892022) q[1];
sx q[1];
rz(1.6581416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8394422) q[0];
sx q[0];
rz(-1.654002) q[0];
sx q[0];
rz(0.010362576) q[0];
rz(-0.16295095) q[2];
sx q[2];
rz(-0.59646791) q[2];
sx q[2];
rz(-0.60566723) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7328651) q[1];
sx q[1];
rz(-0.36204007) q[1];
sx q[1];
rz(1.15508) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5591194) q[3];
sx q[3];
rz(-2.5081722) q[3];
sx q[3];
rz(2.7137091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1411529) q[2];
sx q[2];
rz(-0.65594643) q[2];
sx q[2];
rz(1.8971987) q[2];
rz(2.7164298) q[3];
sx q[3];
rz(-0.77912283) q[3];
sx q[3];
rz(-1.6311133) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4193831) q[0];
sx q[0];
rz(-0.49383759) q[0];
sx q[0];
rz(0.37049946) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(1.3409021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1271034) q[0];
sx q[0];
rz(-1.2047486) q[0];
sx q[0];
rz(-2.8949379) q[0];
rz(-1.4341899) q[2];
sx q[2];
rz(-2.280526) q[2];
sx q[2];
rz(-3.1106069) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51532981) q[1];
sx q[1];
rz(-1.9326107) q[1];
sx q[1];
rz(0.26730178) q[1];
x q[2];
rz(2.9756536) q[3];
sx q[3];
rz(-0.12291848) q[3];
sx q[3];
rz(2.3636706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6293388) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(1.5226927) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.4168134) q[3];
sx q[3];
rz(-0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387977) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(2.8181656) q[1];
sx q[1];
rz(-1.2558116) q[1];
sx q[1];
rz(-1.5423923) q[1];
rz(-1.9771489) q[2];
sx q[2];
rz(-1.5487557) q[2];
sx q[2];
rz(1.1974481) q[2];
rz(2.2446185) q[3];
sx q[3];
rz(-0.57861181) q[3];
sx q[3];
rz(1.2275916) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
