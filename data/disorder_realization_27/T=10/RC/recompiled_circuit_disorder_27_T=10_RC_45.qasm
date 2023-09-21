OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1448016) q[0];
sx q[0];
rz(0.15455833) q[0];
sx q[0];
rz(6.9757087) q[0];
rz(5.0737557) q[1];
sx q[1];
rz(4.3901246) q[1];
sx q[1];
rz(7.6683383) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9557578) q[0];
sx q[0];
rz(-1.8896566) q[0];
sx q[0];
rz(-2.3715109) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7117549) q[2];
sx q[2];
rz(-2.5463383) q[2];
sx q[2];
rz(2.0855479) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7117118) q[1];
sx q[1];
rz(-1.3379339) q[1];
sx q[1];
rz(-1.7377322) q[1];
rz(-pi) q[2];
rz(1.0783644) q[3];
sx q[3];
rz(-2.1636117) q[3];
sx q[3];
rz(-1.4584695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2549071) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(2.0390959) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7339864) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(-2.6876887) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(-1.227238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7176712) q[0];
sx q[0];
rz(-0.14980355) q[0];
sx q[0];
rz(-2.1013837) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4279459) q[2];
sx q[2];
rz(-0.64672856) q[2];
sx q[2];
rz(1.5734067) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2482359) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(-0.30171079) q[1];
rz(-pi) q[2];
rz(-2.482588) q[3];
sx q[3];
rz(-0.25203029) q[3];
sx q[3];
rz(-0.27740955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(0.56742898) q[2];
rz(2.7764017) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(-0.96810961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.658618) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(2.2429402) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-2.8083037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3086739) q[0];
sx q[0];
rz(-1.1090288) q[0];
sx q[0];
rz(-0.69899107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3938445) q[2];
sx q[2];
rz(-1.6672009) q[2];
sx q[2];
rz(2.5071438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2644314) q[1];
sx q[1];
rz(-0.96211551) q[1];
sx q[1];
rz(-0.18552893) q[1];
rz(-1.0201449) q[3];
sx q[3];
rz(-1.9198951) q[3];
sx q[3];
rz(-2.246644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4553392) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(1.1509482) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.6999917) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(-2.4568795) q[0];
rz(-2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.9365786) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8367774) q[0];
sx q[0];
rz(-0.69201058) q[0];
sx q[0];
rz(-1.6230323) q[0];
rz(-pi) q[1];
rz(1.8393458) q[2];
sx q[2];
rz(-0.66187243) q[2];
sx q[2];
rz(0.72999398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4971784) q[1];
sx q[1];
rz(-2.2481611) q[1];
sx q[1];
rz(-0.35269423) q[1];
rz(1.4865925) q[3];
sx q[3];
rz(-0.70662543) q[3];
sx q[3];
rz(-0.018761793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(-0.37115804) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(-0.13312419) q[0];
rz(-0.99331028) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(-0.55508074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017437) q[0];
sx q[0];
rz(-1.886133) q[0];
sx q[0];
rz(-0.01339162) q[0];
rz(-pi) q[1];
rz(-1.416989) q[2];
sx q[2];
rz(-0.4193192) q[2];
sx q[2];
rz(-0.16659444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3913369) q[1];
sx q[1];
rz(-1.5515944) q[1];
sx q[1];
rz(2.0160497) q[1];
rz(-pi) q[2];
rz(-1.6229779) q[3];
sx q[3];
rz(-1.6515886) q[3];
sx q[3];
rz(0.83524708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30620265) q[2];
sx q[2];
rz(-2.1357048) q[2];
sx q[2];
rz(-0.13892697) q[2];
rz(2.1991918) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(0.55148235) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(2.561835) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(-1.5396083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0149536) q[0];
sx q[0];
rz(-0.74288988) q[0];
sx q[0];
rz(-1.8735621) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6164262) q[2];
sx q[2];
rz(-1.4403733) q[2];
sx q[2];
rz(-1.6861196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0582038) q[1];
sx q[1];
rz(-1.1440047) q[1];
sx q[1];
rz(2.8298488) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88605373) q[3];
sx q[3];
rz(-1.3430809) q[3];
sx q[3];
rz(-0.39623228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55398983) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(2.8721151) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(-3.0814734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(2.457298) q[0];
rz(-3.0220095) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(-2.6228242) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2001901) q[0];
sx q[0];
rz(-1.4369643) q[0];
sx q[0];
rz(-0.83394136) q[0];
rz(-2.3709626) q[2];
sx q[2];
rz(-1.4118886) q[2];
sx q[2];
rz(-2.4385914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3541975) q[1];
sx q[1];
rz(-1.1095611) q[1];
sx q[1];
rz(-2.4128777) q[1];
rz(-pi) q[2];
rz(0.026168907) q[3];
sx q[3];
rz(-1.9713638) q[3];
sx q[3];
rz(-0.27163423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(-2.5781412) q[2];
rz(3.0900132) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1812487) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(-1.3990336) q[0];
rz(2.3545806) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(-2.3972437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9104011) q[0];
sx q[0];
rz(-1.7185128) q[0];
sx q[0];
rz(-1.5157248) q[0];
x q[1];
rz(-1.3049576) q[2];
sx q[2];
rz(-2.611534) q[2];
sx q[2];
rz(-2.838138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8213615) q[1];
sx q[1];
rz(-1.7944971) q[1];
sx q[1];
rz(0.51364586) q[1];
x q[2];
rz(0.048150496) q[3];
sx q[3];
rz(-2.6865494) q[3];
sx q[3];
rz(-2.9330848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7156334) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(0.60950935) q[2];
rz(0.65731796) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1881926) q[0];
sx q[0];
rz(-3.0472026) q[0];
sx q[0];
rz(1.5040065) q[0];
rz(1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(2.3666568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79561728) q[0];
sx q[0];
rz(-0.91751639) q[0];
sx q[0];
rz(-0.4972636) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9853285) q[2];
sx q[2];
rz(-2.0093577) q[2];
sx q[2];
rz(1.4807448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.51674023) q[1];
sx q[1];
rz(-1.7588741) q[1];
sx q[1];
rz(2.6810758) q[1];
x q[2];
rz(-1.7899412) q[3];
sx q[3];
rz(-1.7100167) q[3];
sx q[3];
rz(-0.62300357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.187414) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(-0.15787086) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0697486) q[0];
sx q[0];
rz(-0.97244111) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(-0.65746039) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(1.0459895) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2803672) q[0];
sx q[0];
rz(-1.6770289) q[0];
sx q[0];
rz(-0.36061128) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59098737) q[2];
sx q[2];
rz(-1.5467484) q[2];
sx q[2];
rz(1.1514593) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5221094) q[1];
sx q[1];
rz(-1.2088747) q[1];
sx q[1];
rz(-0.98720179) q[1];
x q[2];
rz(1.5973741) q[3];
sx q[3];
rz(-1.9760895) q[3];
sx q[3];
rz(-1.6534896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41632286) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(-1.0894758) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(-0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2789223) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(2.7325148) q[2];
sx q[2];
rz(-1.8553875) q[2];
sx q[2];
rz(2.6864048) q[2];
rz(3.0388721) q[3];
sx q[3];
rz(-2.5892047) q[3];
sx q[3];
rz(1.8451167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
