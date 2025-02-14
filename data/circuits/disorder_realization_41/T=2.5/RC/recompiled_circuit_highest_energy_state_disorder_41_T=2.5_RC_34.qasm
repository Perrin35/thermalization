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
rz(-0.37201878) q[0];
sx q[0];
rz(-2.7899185) q[0];
sx q[0];
rz(-0.055211842) q[0];
rz(1.4456324) q[1];
sx q[1];
rz(-0.90336019) q[1];
sx q[1];
rz(-0.13394314) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3066912) q[0];
sx q[0];
rz(-1.4708232) q[0];
sx q[0];
rz(-2.9213219) q[0];
x q[1];
rz(-0.58007095) q[2];
sx q[2];
rz(-1.7385387) q[2];
sx q[2];
rz(1.3684966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81649924) q[1];
sx q[1];
rz(-1.9032005) q[1];
sx q[1];
rz(1.2219882) q[1];
x q[2];
rz(2.0961742) q[3];
sx q[3];
rz(-1.6046383) q[3];
sx q[3];
rz(-2.6325838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.967531) q[2];
sx q[2];
rz(-1.2207737) q[2];
sx q[2];
rz(2.3270712) q[2];
rz(-0.19541611) q[3];
sx q[3];
rz(-2.312909) q[3];
sx q[3];
rz(-1.0403847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8622417) q[0];
sx q[0];
rz(-2.8793654) q[0];
sx q[0];
rz(-0.97214118) q[0];
rz(0.60130087) q[1];
sx q[1];
rz(-2.1763132) q[1];
sx q[1];
rz(1.7960637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0079673) q[0];
sx q[0];
rz(-1.04038) q[0];
sx q[0];
rz(-0.47904695) q[0];
x q[1];
rz(-0.43429476) q[2];
sx q[2];
rz(-0.58048297) q[2];
sx q[2];
rz(0.24809294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3516772) q[1];
sx q[1];
rz(-0.43444217) q[1];
sx q[1];
rz(-1.2277568) q[1];
rz(-pi) q[2];
rz(0.034288716) q[3];
sx q[3];
rz(-0.18119745) q[3];
sx q[3];
rz(-0.56668355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1874275) q[2];
sx q[2];
rz(-1.2792055) q[2];
sx q[2];
rz(0.97174755) q[2];
rz(2.1238972) q[3];
sx q[3];
rz(-1.965799) q[3];
sx q[3];
rz(0.051232256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63075066) q[0];
sx q[0];
rz(-0.96931163) q[0];
sx q[0];
rz(2.4208659) q[0];
rz(3.0156056) q[1];
sx q[1];
rz(-2.4469913) q[1];
sx q[1];
rz(-0.40649498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6744248) q[0];
sx q[0];
rz(-0.95646383) q[0];
sx q[0];
rz(-2.7697639) q[0];
rz(-2.1265246) q[2];
sx q[2];
rz(-0.90733084) q[2];
sx q[2];
rz(0.86554722) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.49660044) q[1];
sx q[1];
rz(-0.82372249) q[1];
sx q[1];
rz(-0.83622593) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2307736) q[3];
sx q[3];
rz(-2.6916457) q[3];
sx q[3];
rz(-1.6641973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6005818) q[2];
sx q[2];
rz(-1.7966248) q[2];
sx q[2];
rz(2.7715032) q[2];
rz(-3.0626512) q[3];
sx q[3];
rz(-0.97349662) q[3];
sx q[3];
rz(2.6551042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.9321891) q[0];
sx q[0];
rz(-1.3511193) q[0];
sx q[0];
rz(2.6655647) q[0];
rz(1.1141106) q[1];
sx q[1];
rz(-0.94534355) q[1];
sx q[1];
rz(1.6260737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6081656) q[0];
sx q[0];
rz(-1.2100056) q[0];
sx q[0];
rz(-1.3460623) q[0];
rz(0.59882382) q[2];
sx q[2];
rz(-2.235417) q[2];
sx q[2];
rz(-1.7395626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5596931) q[1];
sx q[1];
rz(-1.6554313) q[1];
sx q[1];
rz(0.35178784) q[1];
x q[2];
rz(0.8533303) q[3];
sx q[3];
rz(-1.8195517) q[3];
sx q[3];
rz(2.0847818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49426207) q[2];
sx q[2];
rz(-1.5212955) q[2];
sx q[2];
rz(-2.1935513) q[2];
rz(0.4387795) q[3];
sx q[3];
rz(-2.572757) q[3];
sx q[3];
rz(-0.89890629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5206443) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(-2.7137252) q[0];
rz(-0.95871344) q[1];
sx q[1];
rz(-1.2370647) q[1];
sx q[1];
rz(-2.0895035) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3593345) q[0];
sx q[0];
rz(-2.3365031) q[0];
sx q[0];
rz(2.2710618) q[0];
rz(0.97276997) q[2];
sx q[2];
rz(-2.5199119) q[2];
sx q[2];
rz(-1.5369111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31378779) q[1];
sx q[1];
rz(-0.99936411) q[1];
sx q[1];
rz(-3.0302583) q[1];
rz(-pi) q[2];
rz(0.73486272) q[3];
sx q[3];
rz(-0.70847337) q[3];
sx q[3];
rz(2.3264309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9713356) q[2];
sx q[2];
rz(-1.3054138) q[2];
sx q[2];
rz(-2.7830284) q[2];
rz(1.1809008) q[3];
sx q[3];
rz(-2.2660393) q[3];
sx q[3];
rz(2.6023279) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3233258) q[0];
sx q[0];
rz(-1.8805255) q[0];
sx q[0];
rz(0.64875025) q[0];
rz(-2.5541041) q[1];
sx q[1];
rz(-1.8193388) q[1];
sx q[1];
rz(0.23744753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1797111) q[0];
sx q[0];
rz(-1.474232) q[0];
sx q[0];
rz(-1.5386816) q[0];
rz(-1.6726137) q[2];
sx q[2];
rz(-1.8378496) q[2];
sx q[2];
rz(1.7783527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12585576) q[1];
sx q[1];
rz(-2.3185711) q[1];
sx q[1];
rz(-1.5344844) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6471938) q[3];
sx q[3];
rz(-2.2174617) q[3];
sx q[3];
rz(1.4495415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6610403) q[2];
sx q[2];
rz(-2.48017) q[2];
sx q[2];
rz(0.95575571) q[2];
rz(1.7620979) q[3];
sx q[3];
rz(-2.5199514) q[3];
sx q[3];
rz(-0.1563589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.5697923) q[0];
sx q[0];
rz(-3.1389696) q[0];
sx q[0];
rz(3.0963335) q[0];
rz(-2.3472002) q[1];
sx q[1];
rz(-2.2519799) q[1];
sx q[1];
rz(0.20283595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9211333) q[0];
sx q[0];
rz(-2.5510802) q[0];
sx q[0];
rz(-2.8648225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7603821) q[2];
sx q[2];
rz(-0.84824296) q[2];
sx q[2];
rz(-2.0675532) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.07405532) q[1];
sx q[1];
rz(-0.54880868) q[1];
sx q[1];
rz(-2.9850053) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.037153745) q[3];
sx q[3];
rz(-0.84913466) q[3];
sx q[3];
rz(-0.29216684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63991919) q[2];
sx q[2];
rz(-1.4036274) q[2];
sx q[2];
rz(-2.5301834) q[2];
rz(-1.1008788) q[3];
sx q[3];
rz(-2.5067063) q[3];
sx q[3];
rz(2.8618405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.47073498) q[0];
sx q[0];
rz(-1.1962471) q[0];
sx q[0];
rz(1.083495) q[0];
rz(0.96393839) q[1];
sx q[1];
rz(-1.0716535) q[1];
sx q[1];
rz(-0.49577698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7115323) q[0];
sx q[0];
rz(-1.4976131) q[0];
sx q[0];
rz(-0.021223193) q[0];
rz(-1.8476358) q[2];
sx q[2];
rz(-1.7466892) q[2];
sx q[2];
rz(-1.7664282) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.816905) q[1];
sx q[1];
rz(-2.4084615) q[1];
sx q[1];
rz(0.92589652) q[1];
x q[2];
rz(0.98306174) q[3];
sx q[3];
rz(-1.5458298) q[3];
sx q[3];
rz(-0.72265305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.701391) q[2];
sx q[2];
rz(-0.76875606) q[2];
sx q[2];
rz(-2.4526147) q[2];
rz(1.5135328) q[3];
sx q[3];
rz(-2.3115034) q[3];
sx q[3];
rz(1.0015063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5772407) q[0];
sx q[0];
rz(-1.5654726) q[0];
sx q[0];
rz(-2.054457) q[0];
rz(-2.3785036) q[1];
sx q[1];
rz(-1.3975846) q[1];
sx q[1];
rz(2.9147002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54880035) q[0];
sx q[0];
rz(-0.88946402) q[0];
sx q[0];
rz(2.8517904) q[0];
rz(-pi) q[1];
rz(-2.8434038) q[2];
sx q[2];
rz(-1.2623566) q[2];
sx q[2];
rz(1.0210263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5530738) q[1];
sx q[1];
rz(-1.5344193) q[1];
sx q[1];
rz(-2.9168771) q[1];
rz(-0.97051981) q[3];
sx q[3];
rz(-1.1423938) q[3];
sx q[3];
rz(0.71479931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6202966) q[2];
sx q[2];
rz(-2.2038348) q[2];
sx q[2];
rz(2.1916981) q[2];
rz(1.0319483) q[3];
sx q[3];
rz(-2.2622006) q[3];
sx q[3];
rz(0.43752813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501605) q[0];
sx q[0];
rz(-0.10364769) q[0];
sx q[0];
rz(1.9895122) q[0];
rz(-3.0488455) q[1];
sx q[1];
rz(-0.68789613) q[1];
sx q[1];
rz(2.6377717) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4277305) q[0];
sx q[0];
rz(-2.6127338) q[0];
sx q[0];
rz(1.4725757) q[0];
rz(-pi) q[1];
rz(-2.4359792) q[2];
sx q[2];
rz(-2.793987) q[2];
sx q[2];
rz(-3.0726031) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31695932) q[1];
sx q[1];
rz(-1.0976296) q[1];
sx q[1];
rz(1.6184357) q[1];
rz(-pi) q[2];
rz(-2.2600365) q[3];
sx q[3];
rz(-2.6668352) q[3];
sx q[3];
rz(0.17533824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6759701) q[2];
sx q[2];
rz(-2.0268107) q[2];
sx q[2];
rz(-2.5176804) q[2];
rz(-0.55656773) q[3];
sx q[3];
rz(-0.68816319) q[3];
sx q[3];
rz(-2.9437959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.9482166) q[0];
sx q[0];
rz(-0.93127903) q[0];
sx q[0];
rz(-2.2122526) q[0];
rz(2.9931862) q[1];
sx q[1];
rz(-1.8971309) q[1];
sx q[1];
rz(2.94577) q[1];
rz(0.23870809) q[2];
sx q[2];
rz(-2.6663824) q[2];
sx q[2];
rz(2.0486365) q[2];
rz(-1.8529057) q[3];
sx q[3];
rz(-2.8587881) q[3];
sx q[3];
rz(2.7360873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
