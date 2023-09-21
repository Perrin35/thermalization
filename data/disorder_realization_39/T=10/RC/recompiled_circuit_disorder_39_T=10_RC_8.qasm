OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(-0.46407035) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(1.8571412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0425134) q[0];
sx q[0];
rz(-2.6455542) q[0];
sx q[0];
rz(2.2201559) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.651865) q[2];
sx q[2];
rz(-1.5999319) q[2];
sx q[2];
rz(-0.65087989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7317036) q[1];
sx q[1];
rz(-1.1176425) q[1];
sx q[1];
rz(1.975592) q[1];
rz(-2.9908882) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(2.5022262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-2.822067) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(2.615036) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(0.79663509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5027673) q[0];
sx q[0];
rz(-0.728038) q[0];
sx q[0];
rz(1.0685705) q[0];
x q[1];
rz(-0.30226207) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(2.9994534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14256515) q[1];
sx q[1];
rz(-1.5836645) q[1];
sx q[1];
rz(2.5343115) q[1];
x q[2];
rz(-0.91005743) q[3];
sx q[3];
rz(-2.4168192) q[3];
sx q[3];
rz(1.4409325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(-2.3550418) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-2.6942159) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-0.70297855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709135) q[0];
sx q[0];
rz(-1.850607) q[0];
sx q[0];
rz(-0.12165102) q[0];
rz(-pi) q[1];
rz(-2.5419652) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(-1.3404913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30333334) q[1];
sx q[1];
rz(-1.5240182) q[1];
sx q[1];
rz(-0.10951885) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94743518) q[3];
sx q[3];
rz(-1.4294335) q[3];
sx q[3];
rz(-0.75368222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(0.91153574) q[2];
rz(0.91397816) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(0.56328303) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8877836) q[0];
sx q[0];
rz(-2.1953708) q[0];
sx q[0];
rz(2.8028691) q[0];
rz(-1.9070542) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(-1.3627571) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.024151) q[1];
sx q[1];
rz(-1.5643331) q[1];
sx q[1];
rz(-0.27177377) q[1];
rz(-pi) q[2];
rz(1.2761649) q[3];
sx q[3];
rz(-2.2664824) q[3];
sx q[3];
rz(-0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48878601) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(0.164786) q[2];
rz(-2.9131043) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(2.1951108) q[3];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(-0.35119855) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(-1.9794827) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0162109) q[0];
sx q[0];
rz(-2.1124766) q[0];
sx q[0];
rz(1.1075695) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94050546) q[2];
sx q[2];
rz(-2.5224707) q[2];
sx q[2];
rz(1.3319912) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8309161) q[1];
sx q[1];
rz(-1.1569996) q[1];
sx q[1];
rz(-0.78891854) q[1];
x q[2];
rz(1.5557628) q[3];
sx q[3];
rz(-2.2923277) q[3];
sx q[3];
rz(-0.3258457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6953485) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(1.8943141) q[2];
rz(0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-0.22918992) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86876774) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(1.8355339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7040946) q[0];
sx q[0];
rz(-1.1878345) q[0];
sx q[0];
rz(-2.1279446) q[0];
rz(0.46361228) q[2];
sx q[2];
rz(-1.7525275) q[2];
sx q[2];
rz(0.52464991) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2893147) q[1];
sx q[1];
rz(-1.5783974) q[1];
sx q[1];
rz(0.72035933) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58879866) q[3];
sx q[3];
rz(-1.6949123) q[3];
sx q[3];
rz(2.0037946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(-0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.612161) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(-2.7244862) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313275) q[0];
sx q[0];
rz(-1.350704) q[0];
sx q[0];
rz(-3.1165645) q[0];
rz(-pi) q[1];
x q[1];
rz(1.276937) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(0.35640946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6228094) q[1];
sx q[1];
rz(-0.86474027) q[1];
sx q[1];
rz(2.519636) q[1];
x q[2];
rz(1.2715289) q[3];
sx q[3];
rz(-2.2021658) q[3];
sx q[3];
rz(2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32507867) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(0.11238012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053514078) q[0];
sx q[0];
rz(-1.3980165) q[0];
sx q[0];
rz(2.0515576) q[0];
x q[1];
rz(2.8939395) q[2];
sx q[2];
rz(-2.7981749) q[2];
sx q[2];
rz(-2.0695956) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97265128) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(-2.6886743) q[1];
x q[2];
rz(1.3120679) q[3];
sx q[3];
rz(-1.6885763) q[3];
sx q[3];
rz(3.0974914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0743951) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(1.2049234) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444721) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(0.034974139) q[0];
rz(-0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(2.2299178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.002279) q[0];
sx q[0];
rz(-0.3542491) q[0];
sx q[0];
rz(-2.2646963) q[0];
rz(-1.6976835) q[2];
sx q[2];
rz(-2.1107026) q[2];
sx q[2];
rz(-1.9829696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1940143) q[1];
sx q[1];
rz(-1.7742426) q[1];
sx q[1];
rz(0.46781637) q[1];
x q[2];
rz(-1.0802286) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(-2.3232943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(0.24924499) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(0.58139873) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(1.4153597) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0534131) q[0];
sx q[0];
rz(-1.3706494) q[0];
sx q[0];
rz(-2.8739909) q[0];
rz(-pi) q[1];
rz(-2.084311) q[2];
sx q[2];
rz(-2.1272749) q[2];
sx q[2];
rz(-1.8884115) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49299875) q[1];
sx q[1];
rz(-0.10222888) q[1];
sx q[1];
rz(-0.74536721) q[1];
rz(-pi) q[2];
rz(-0.38378999) q[3];
sx q[3];
rz(-1.0597611) q[3];
sx q[3];
rz(0.25110652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(1.7278016) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(-1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5851371) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(1.8150868) q[2];
sx q[2];
rz(-0.46637022) q[2];
sx q[2];
rz(-2.7844219) q[2];
rz(-2.1429569) q[3];
sx q[3];
rz(-1.640366) q[3];
sx q[3];
rz(0.58983005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];