OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(4.2098213) q[0];
sx q[0];
rz(9.8888483) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099079236) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(2.2201559) q[0];
rz(3.0797144) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(-2.2762736) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1661108) q[1];
sx q[1];
rz(-1.9327285) q[1];
sx q[1];
rz(2.6544177) q[1];
rz(-pi) q[2];
rz(0.10847096) q[3];
sx q[3];
rz(-2.9900108) q[3];
sx q[3];
rz(0.82419318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39711943) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-2.822067) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(0.52655667) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(-2.3449576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731199) q[0];
sx q[0];
rz(-2.1935049) q[0];
sx q[0];
rz(0.40533439) q[0];
x q[1];
rz(-0.57467069) q[2];
sx q[2];
rz(-1.7386912) q[2];
sx q[2];
rz(-1.6811973) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7318774) q[1];
sx q[1];
rz(-2.5341923) q[1];
sx q[1];
rz(0.022547988) q[1];
rz(-2.2315352) q[3];
sx q[3];
rz(-0.72477341) q[3];
sx q[3];
rz(1.4409325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(-2.3550418) q[2];
rz(2.6484047) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-2.6942159) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-0.70297855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706791) q[0];
sx q[0];
rz(-1.2909856) q[0];
sx q[0];
rz(-3.0199416) q[0];
rz(2.643232) q[2];
sx q[2];
rz(-1.8866072) q[2];
sx q[2];
rz(-0.74893803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30333334) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(3.0320738) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.968077) q[3];
sx q[3];
rz(-2.1870038) q[3];
sx q[3];
rz(0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-2.3655868) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(0.56328303) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3455428) q[0];
sx q[0];
rz(-0.69957083) q[0];
sx q[0];
rz(-2.0027341) q[0];
x q[1];
rz(-0.40043719) q[2];
sx q[2];
rz(-0.73095989) q[2];
sx q[2];
rz(-1.8797344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.024151) q[1];
sx q[1];
rz(-1.5643331) q[1];
sx q[1];
rz(0.27177377) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2761649) q[3];
sx q[3];
rz(-2.2664824) q[3];
sx q[3];
rz(0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(-2.9768067) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8673458) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-0.85246032) q[0];
rz(-0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.9794827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4974385) q[0];
sx q[0];
rz(-2.4442406) q[0];
sx q[0];
rz(-2.5028412) q[0];
rz(-pi) q[1];
x q[1];
rz(1.048462) q[2];
sx q[2];
rz(-1.9198717) q[2];
sx q[2];
rz(2.366684) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1197549) q[1];
sx q[1];
rz(-0.86958414) q[1];
sx q[1];
rz(0.55418684) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5858298) q[3];
sx q[3];
rz(-2.2923277) q[3];
sx q[3];
rz(0.3258457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6953485) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(-0.22918992) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(0.75063467) q[0];
rz(2.0320832) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40697843) q[0];
sx q[0];
rz(-0.66440551) q[0];
sx q[0];
rz(-2.2218496) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6779804) q[2];
sx q[2];
rz(-1.3890651) q[2];
sx q[2];
rz(-0.52464991) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29014153) q[1];
sx q[1];
rz(-2.4212004) q[1];
sx q[1];
rz(-0.011522567) q[1];
rz(0.2209729) q[3];
sx q[3];
rz(-2.5413725) q[3];
sx q[3];
rz(2.5252987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7541472) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(2.5409017) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(2.7887662) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(-2.0948998) q[0];
rz(-1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(-2.7244862) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313275) q[0];
sx q[0];
rz(-1.7908887) q[0];
sx q[0];
rz(-3.1165645) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8344526) q[2];
sx q[2];
rz(-1.2899613) q[2];
sx q[2];
rz(1.839523) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6542146) q[1];
sx q[1];
rz(-1.1114792) q[1];
sx q[1];
rz(-0.76141255) q[1];
x q[2];
rz(-1.2715289) q[3];
sx q[3];
rz(-2.2021658) q[3];
sx q[3];
rz(-2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(0.55316365) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.9518071) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.534879) q[0];
sx q[0];
rz(-1.0977854) q[0];
sx q[0];
rz(0.19434778) q[0];
rz(-pi) q[1];
rz(-2.8078812) q[2];
sx q[2];
rz(-1.4881655) q[2];
sx q[2];
rz(2.4090648) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29490678) q[1];
sx q[1];
rz(-1.2290188) q[1];
sx q[1];
rz(-2.3218367) q[1];
x q[2];
rz(1.1376082) q[3];
sx q[3];
rz(-0.28372753) q[3];
sx q[3];
rz(-1.9445436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0743951) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(0.91167489) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23110403) q[0];
sx q[0];
rz(-1.7945053) q[0];
sx q[0];
rz(1.847812) q[0];
rz(2.5981204) q[2];
sx q[2];
rz(-1.6795571) q[2];
sx q[2];
rz(0.47765884) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1940143) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(2.6737763) q[1];
rz(-pi) q[2];
rz(-2.6058042) q[3];
sx q[3];
rz(-2.0013323) q[3];
sx q[3];
rz(2.146194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-2.8923477) q[2];
rz(-0.76672673) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(-0.37208474) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(1.7262329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0534131) q[0];
sx q[0];
rz(-1.3706494) q[0];
sx q[0];
rz(-2.8739909) q[0];
rz(-2.4731589) q[2];
sx q[2];
rz(-2.4032776) q[2];
sx q[2];
rz(0.43503209) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.806554) q[1];
sx q[1];
rz(-1.5015258) q[1];
sx q[1];
rz(-3.0663504) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1595702) q[3];
sx q[3];
rz(-2.5128799) q[3];
sx q[3];
rz(-2.2002937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(-0.20467219) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(-3.0204308) q[2];
sx q[2];
rz(-2.0222752) q[2];
sx q[2];
rz(-3.0565699) q[2];
rz(-0.99863573) q[3];
sx q[3];
rz(-1.5012267) q[3];
sx q[3];
rz(-2.5517626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
