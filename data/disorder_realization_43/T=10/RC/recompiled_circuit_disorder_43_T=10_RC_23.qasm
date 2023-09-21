OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(-1.3347081) q[0];
rz(2.788738) q[1];
sx q[1];
rz(3.3021441) q[1];
sx q[1];
rz(8.4488206) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2655576) q[0];
sx q[0];
rz(-1.8119438) q[0];
sx q[0];
rz(0.59224706) q[0];
rz(-pi) q[1];
rz(2.1337778) q[2];
sx q[2];
rz(-1.8554167) q[2];
sx q[2];
rz(-2.9993338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.576697) q[1];
sx q[1];
rz(-1.7045867) q[1];
sx q[1];
rz(-0.96058515) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7482412) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(2.2497183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(-2.3577918) q[2];
rz(-2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(0.17856199) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(-3.1352502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3865249) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-3.0088739) q[0];
x q[1];
rz(0.96887529) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.7768163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.43286846) q[1];
sx q[1];
rz(-1.8808865) q[1];
sx q[1];
rz(-1.6354145) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5793521) q[3];
sx q[3];
rz(-1.7842818) q[3];
sx q[3];
rz(-1.436304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94770849) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(2.2634899) q[2];
rz(2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(2.7745461) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(0.095741622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75671065) q[0];
sx q[0];
rz(-1.4178935) q[0];
sx q[0];
rz(-1.3923313) q[0];
x q[1];
rz(-0.50600608) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(0.72088748) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0248191) q[1];
sx q[1];
rz(-2.1110592) q[1];
sx q[1];
rz(0.54622548) q[1];
x q[2];
rz(0.057283244) q[3];
sx q[3];
rz(-2.5896642) q[3];
sx q[3];
rz(-2.3416167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(0.21162027) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(-0.37477469) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(-2.7950177) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(2.1077164) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047065145) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(-2.4497776) q[0];
rz(-pi) q[1];
rz(-0.83861645) q[2];
sx q[2];
rz(-0.71574434) q[2];
sx q[2];
rz(-1.0775623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4775866) q[1];
sx q[1];
rz(-0.2020745) q[1];
sx q[1];
rz(-1.9007773) q[1];
x q[2];
rz(-1.8978118) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(1.3467849) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-2.4436061) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(-1.746009) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(-2.5674852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2643471) q[0];
sx q[0];
rz(-0.97752042) q[0];
sx q[0];
rz(-2.9869153) q[0];
x q[1];
rz(0.026002361) q[2];
sx q[2];
rz(-2.4261195) q[2];
sx q[2];
rz(-2.0040087) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1787604) q[1];
sx q[1];
rz(-1.8492336) q[1];
sx q[1];
rz(2.7774485) q[1];
rz(-0.19458171) q[3];
sx q[3];
rz(-2.2268725) q[3];
sx q[3];
rz(1.3692828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(-2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-0.73227698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2039316) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(2.926814) q[0];
x q[1];
rz(-2.3789669) q[2];
sx q[2];
rz(-1.0566933) q[2];
sx q[2];
rz(-2.5800173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9975035) q[1];
sx q[1];
rz(-2.2170076) q[1];
sx q[1];
rz(-1.287582) q[1];
x q[2];
rz(-0.332571) q[3];
sx q[3];
rz(-2.5453574) q[3];
sx q[3];
rz(0.81070825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(-0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(0.77004534) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5891588) q[0];
sx q[0];
rz(-1.7434412) q[0];
sx q[0];
rz(0.60683672) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40600834) q[2];
sx q[2];
rz(-2.1210665) q[2];
sx q[2];
rz(-1.5232435) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29994446) q[1];
sx q[1];
rz(-2.3023459) q[1];
sx q[1];
rz(-0.70920918) q[1];
rz(-pi) q[2];
rz(-0.55925925) q[3];
sx q[3];
rz(-0.57563215) q[3];
sx q[3];
rz(-1.5881133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4975171) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(2.1264123) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5148233) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-2.7752005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.217427) q[0];
sx q[0];
rz(-0.96091849) q[0];
sx q[0];
rz(-2.2141371) q[0];
x q[1];
rz(1.1364469) q[2];
sx q[2];
rz(-0.98698101) q[2];
sx q[2];
rz(2.7455612) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.061325039) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(-2.8894043) q[1];
rz(-pi) q[2];
rz(1.6493158) q[3];
sx q[3];
rz(-0.55320569) q[3];
sx q[3];
rz(-1.5403403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0176795) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(2.8477342) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-0.88395399) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-0.79137897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521097) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(0.86274685) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4653761) q[2];
sx q[2];
rz(-1.1306136) q[2];
sx q[2];
rz(2.344775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3190805) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(0.49088571) q[1];
rz(1.295624) q[3];
sx q[3];
rz(-2.5857946) q[3];
sx q[3];
rz(2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.7609319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35683435) q[0];
sx q[0];
rz(-2.1958302) q[0];
sx q[0];
rz(2.9805095) q[0];
rz(1.5153377) q[2];
sx q[2];
rz(-0.43606731) q[2];
sx q[2];
rz(-0.01576327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32523649) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(0.34646323) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7396183) q[3];
sx q[3];
rz(-1.7572548) q[3];
sx q[3];
rz(-2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-2.805368) q[2];
rz(-2.0119038) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(-1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(-0.54429383) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(-0.82472807) q[2];
sx q[2];
rz(-1.8246973) q[2];
sx q[2];
rz(1.4551103) q[2];
rz(-0.23562283) q[3];
sx q[3];
rz(-2.1274673) q[3];
sx q[3];
rz(-2.6475788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];