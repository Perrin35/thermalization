OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(4.2545118) q[1];
sx q[1];
rz(1.7634044) q[1];
sx q[1];
rz(7.4982285) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62984798) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(-1.758979) q[0];
rz(-1.2826074) q[2];
sx q[2];
rz(-0.93027861) q[2];
sx q[2];
rz(-3.1079907) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4814264) q[1];
sx q[1];
rz(-2.5170442) q[1];
sx q[1];
rz(0.98510965) q[1];
rz(-pi) q[2];
rz(2.0753161) q[3];
sx q[3];
rz(-0.30116943) q[3];
sx q[3];
rz(-1.0217713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.5391301) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8388226) q[0];
sx q[0];
rz(-1.4937703) q[0];
sx q[0];
rz(2.1588615) q[0];
x q[1];
rz(-1.7878754) q[2];
sx q[2];
rz(-2.2976544) q[2];
sx q[2];
rz(-0.94490563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(-2.944988) q[1];
x q[2];
rz(0.42373557) q[3];
sx q[3];
rz(-1.7859922) q[3];
sx q[3];
rz(-0.37950294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.845528) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(2.5862397) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1486737) q[0];
sx q[0];
rz(-0.86704555) q[0];
sx q[0];
rz(-0.91627319) q[0];
rz(-0.31366445) q[2];
sx q[2];
rz(-2.1137538) q[2];
sx q[2];
rz(-1.107159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4761423) q[1];
sx q[1];
rz(-2.1995771) q[1];
sx q[1];
rz(-0.20035845) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2461353) q[3];
sx q[3];
rz(-2.5723296) q[3];
sx q[3];
rz(-0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811532) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(-2.326791) q[0];
rz(1.3793777) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(0.25517685) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164453) q[0];
sx q[0];
rz(-1.2108004) q[0];
sx q[0];
rz(-1.8871334) q[0];
x q[1];
rz(0.63919477) q[2];
sx q[2];
rz(-1.8393469) q[2];
sx q[2];
rz(-0.07721363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.081101) q[1];
sx q[1];
rz(-1.2730518) q[1];
sx q[1];
rz(2.9213419) q[1];
rz(1.049794) q[3];
sx q[3];
rz(-1.61506) q[3];
sx q[3];
rz(1.9998159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2531551) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.859905) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(0.056093562) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1144975) q[0];
sx q[0];
rz(-1.9341015) q[0];
sx q[0];
rz(2.5278805) q[0];
rz(0.9390097) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(-2.3914571) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6519421) q[1];
sx q[1];
rz(-2.838476) q[1];
sx q[1];
rz(-0.048708212) q[1];
x q[2];
rz(-1.7124743) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(-0.39839881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(-0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(1.8796857) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.944825) q[0];
sx q[0];
rz(-1.7260572) q[0];
sx q[0];
rz(1.0374271) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2013821) q[2];
sx q[2];
rz(-1.5988837) q[2];
sx q[2];
rz(-2.0036151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8050025) q[1];
sx q[1];
rz(-2.046236) q[1];
sx q[1];
rz(-0.57979433) q[1];
x q[2];
rz(1.8618705) q[3];
sx q[3];
rz(-1.4143412) q[3];
sx q[3];
rz(0.94369704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(1.0423638) q[2];
rz(0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(-1.8235122) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689432) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(0.48604301) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91673135) q[2];
sx q[2];
rz(-1.6628633) q[2];
sx q[2];
rz(0.23840657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2513189) q[1];
sx q[1];
rz(-2.2379025) q[1];
sx q[1];
rz(3.1097079) q[1];
rz(-pi) q[2];
rz(-0.20603541) q[3];
sx q[3];
rz(-1.9427951) q[3];
sx q[3];
rz(2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(0.73658529) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37305957) q[0];
sx q[0];
rz(-0.031950843) q[0];
sx q[0];
rz(1.91747) q[0];
rz(-1.2484776) q[2];
sx q[2];
rz(-1.467448) q[2];
sx q[2];
rz(-1.581574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8229586) q[1];
sx q[1];
rz(-1.4689323) q[1];
sx q[1];
rz(0.46551367) q[1];
x q[2];
rz(2.8183476) q[3];
sx q[3];
rz(-1.1439699) q[3];
sx q[3];
rz(-2.1904898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(2.4576808) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443611) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(2.1445403) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(-0.7448147) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9841524) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(1.1293344) q[0];
rz(-pi) q[1];
rz(-2.7460329) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(-1.2566483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77764952) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(-1.2121483) q[1];
x q[2];
rz(-0.39448491) q[3];
sx q[3];
rz(-2.3186765) q[3];
sx q[3];
rz(-0.2713954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(1.194681) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.9706479) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02689657) q[0];
sx q[0];
rz(-2.7141502) q[0];
sx q[0];
rz(-0.052274152) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3659336) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(-1.0722216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7217692) q[1];
sx q[1];
rz(-0.024200736) q[1];
sx q[1];
rz(1.2060549) q[1];
rz(-pi) q[2];
rz(-2.714614) q[3];
sx q[3];
rz(-1.8378165) q[3];
sx q[3];
rz(0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(0.56636089) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(0.099427632) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-2.2291017) q[2];
sx q[2];
rz(-0.9463263) q[2];
sx q[2];
rz(1.8778388) q[2];
rz(1.2373274) q[3];
sx q[3];
rz(-1.5860535) q[3];
sx q[3];
rz(0.64762583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
