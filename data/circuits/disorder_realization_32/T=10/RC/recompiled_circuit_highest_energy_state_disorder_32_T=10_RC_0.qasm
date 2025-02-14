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
rz(1.2122756) q[0];
sx q[0];
rz(-2.0097998) q[0];
sx q[0];
rz(1.3728859) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(4.5166587) q[1];
sx q[1];
rz(5.4969129) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776913) q[0];
sx q[0];
rz(-2.7161975) q[0];
sx q[0];
rz(2.0869654) q[0];
rz(-pi) q[1];
rz(-2.588618) q[2];
sx q[2];
rz(-2.1510501) q[2];
sx q[2];
rz(1.0667104) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.3226124) q[1];
sx q[1];
rz(-0.29085813) q[1];
sx q[1];
rz(1.4664654) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1699675) q[3];
sx q[3];
rz(-1.9135981) q[3];
sx q[3];
rz(-0.179571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11975153) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(-1.2500259) q[2];
rz(2.405808) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(-1.638394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982872) q[0];
sx q[0];
rz(-1.1730288) q[0];
sx q[0];
rz(0.071320891) q[0];
rz(1.9885063) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(0.52871314) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77069717) q[0];
sx q[0];
rz(-2.5636854) q[0];
sx q[0];
rz(0.96630567) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9527093) q[2];
sx q[2];
rz(-1.3391277) q[2];
sx q[2];
rz(0.29555368) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1518883) q[1];
sx q[1];
rz(-2.3373649) q[1];
sx q[1];
rz(-0.81664919) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1558481) q[3];
sx q[3];
rz(-1.7768716) q[3];
sx q[3];
rz(1.8725439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(0.028623494) q[2];
rz(-2.7875767) q[3];
sx q[3];
rz(-0.75494868) q[3];
sx q[3];
rz(-0.55907512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(-0.89135528) q[0];
rz(2.640653) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(-3.0827789) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0863773) q[0];
sx q[0];
rz(-0.098345938) q[0];
sx q[0];
rz(2.6340061) q[0];
rz(-pi) q[1];
rz(2.8562653) q[2];
sx q[2];
rz(-1.0614191) q[2];
sx q[2];
rz(-2.2338108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8938287) q[1];
sx q[1];
rz(-0.86834748) q[1];
sx q[1];
rz(0.6599627) q[1];
x q[2];
rz(1.5427049) q[3];
sx q[3];
rz(-1.7511611) q[3];
sx q[3];
rz(-1.9244058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3567051) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(2.9050262) q[2];
rz(0.92075721) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22223602) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(-0.87523571) q[0];
rz(1.596176) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(0.045305591) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0939498) q[0];
sx q[0];
rz(-1.6730621) q[0];
sx q[0];
rz(3.0759252) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95618177) q[2];
sx q[2];
rz(-0.76757694) q[2];
sx q[2];
rz(-1.1078579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.873535) q[1];
sx q[1];
rz(-0.90374871) q[1];
sx q[1];
rz(-0.11063834) q[1];
x q[2];
rz(1.8604467) q[3];
sx q[3];
rz(-2.5389606) q[3];
sx q[3];
rz(1.2887736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8594325) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(-1.1939987) q[2];
rz(-2.2512839) q[3];
sx q[3];
rz(-1.9186391) q[3];
sx q[3];
rz(-0.94849006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0328261) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(-2.4591675) q[0];
rz(-1.2884864) q[1];
sx q[1];
rz(-1.5623743) q[1];
sx q[1];
rz(-1.8103745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3159861) q[0];
sx q[0];
rz(-0.84418143) q[0];
sx q[0];
rz(0.021922317) q[0];
x q[1];
rz(2.140215) q[2];
sx q[2];
rz(-0.45854202) q[2];
sx q[2];
rz(-1.4024782) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4466616) q[1];
sx q[1];
rz(-2.7547421) q[1];
sx q[1];
rz(-2.7617161) q[1];
x q[2];
rz(-1.0430452) q[3];
sx q[3];
rz(-1.966779) q[3];
sx q[3];
rz(2.2137303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1149301) q[2];
sx q[2];
rz(-0.59938359) q[2];
sx q[2];
rz(-2.4746573) q[2];
rz(-2.5866348) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(-2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894593) q[0];
sx q[0];
rz(-1.2586559) q[0];
sx q[0];
rz(2.4378648) q[0];
rz(-0.16920371) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(-0.15883787) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.279351) q[0];
sx q[0];
rz(-0.18687525) q[0];
sx q[0];
rz(-2.274155) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9877404) q[2];
sx q[2];
rz(-1.6922608) q[2];
sx q[2];
rz(2.6328994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.006268) q[1];
sx q[1];
rz(-0.49990955) q[1];
sx q[1];
rz(-0.27854021) q[1];
x q[2];
rz(1.7991285) q[3];
sx q[3];
rz(-0.94551802) q[3];
sx q[3];
rz(-2.5061553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3682897) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(-0.85757315) q[2];
rz(1.9722021) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17539772) q[0];
sx q[0];
rz(-2.3747787) q[0];
sx q[0];
rz(2.4229557) q[0];
rz(2.8596558) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(2.4729572) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40067264) q[0];
sx q[0];
rz(-2.1869279) q[0];
sx q[0];
rz(-0.4192062) q[0];
rz(-pi) q[1];
rz(0.94271548) q[2];
sx q[2];
rz(-1.8657902) q[2];
sx q[2];
rz(-1.4843954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.201828) q[1];
sx q[1];
rz(-1.9376905) q[1];
sx q[1];
rz(0.35122996) q[1];
rz(-1.578978) q[3];
sx q[3];
rz(-0.93957179) q[3];
sx q[3];
rz(-0.62731987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4947027) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(1.0364214) q[2];
rz(2.5070665) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(2.3488267) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6805639) q[0];
sx q[0];
rz(-0.11756086) q[0];
sx q[0];
rz(0.72599894) q[0];
rz(1.9793824) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(-1.1964218) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40563075) q[0];
sx q[0];
rz(-2.9655632) q[0];
sx q[0];
rz(1.4077296) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4216719) q[2];
sx q[2];
rz(-0.31055488) q[2];
sx q[2];
rz(-2.9322185) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4359522) q[1];
sx q[1];
rz(-2.0355823) q[1];
sx q[1];
rz(0.56805261) q[1];
rz(-0.4287339) q[3];
sx q[3];
rz(-0.7168684) q[3];
sx q[3];
rz(-0.69566876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2339345) q[2];
sx q[2];
rz(-1.5479167) q[2];
sx q[2];
rz(0.93351239) q[2];
rz(-0.61698169) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(-0.84701076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9762978) q[0];
sx q[0];
rz(-0.10460654) q[0];
sx q[0];
rz(3.0883375) q[0];
rz(2.8540197) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(2.8823749) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0236146) q[0];
sx q[0];
rz(-0.058246944) q[0];
sx q[0];
rz(-1.144676) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0854112) q[2];
sx q[2];
rz(-0.26517235) q[2];
sx q[2];
rz(-2.3725703) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99154918) q[1];
sx q[1];
rz(-0.75356149) q[1];
sx q[1];
rz(1.6607453) q[1];
rz(0.15786981) q[3];
sx q[3];
rz(-0.77220687) q[3];
sx q[3];
rz(1.7534353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58664924) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(-2.9602236) q[2];
rz(2.7465316) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(-2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7419389) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(0.3821061) q[0];
rz(0.29640472) q[1];
sx q[1];
rz(-2.1370685) q[1];
sx q[1];
rz(1.2688676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73055501) q[0];
sx q[0];
rz(-2.0417622) q[0];
sx q[0];
rz(-2.6593326) q[0];
rz(-pi) q[1];
rz(-0.13258719) q[2];
sx q[2];
rz(-0.90182038) q[2];
sx q[2];
rz(-0.35945177) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.813432) q[1];
sx q[1];
rz(-2.8012271) q[1];
sx q[1];
rz(1.634738) q[1];
rz(-pi) q[2];
rz(-0.45370729) q[3];
sx q[3];
rz(-1.6804916) q[3];
sx q[3];
rz(3.0806993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2962013) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(-0.38983795) q[2];
rz(-1.2975533) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(-2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525477) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(1.7851495) q[1];
sx q[1];
rz(-1.0379797) q[1];
sx q[1];
rz(-1.5300068) q[1];
rz(2.155921) q[2];
sx q[2];
rz(-0.73230848) q[2];
sx q[2];
rz(-2.6072235) q[2];
rz(2.3572902) q[3];
sx q[3];
rz(-1.3275066) q[3];
sx q[3];
rz(2.2307997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
