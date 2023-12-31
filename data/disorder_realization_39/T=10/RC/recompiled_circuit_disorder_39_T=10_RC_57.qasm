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
rz(2.6775223) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0811631) q[0];
sx q[0];
rz(-1.2788749) q[0];
sx q[0];
rz(1.1638327) q[0];
x q[1];
rz(-2.651865) q[2];
sx q[2];
rz(-1.5416607) q[2];
sx q[2];
rz(-2.4907128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97548188) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(-2.6544177) q[1];
x q[2];
rz(0.15070446) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(2.5022262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.39711943) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(2.822067) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-2.615036) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(0.79663509) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731199) q[0];
sx q[0];
rz(-2.1935049) q[0];
sx q[0];
rz(-2.7362583) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57467069) q[2];
sx q[2];
rz(-1.7386912) q[2];
sx q[2];
rz(-1.4603953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7318774) q[1];
sx q[1];
rz(-2.5341923) q[1];
sx q[1];
rz(-0.022547988) q[1];
x q[2];
rz(-2.2315352) q[3];
sx q[3];
rz(-2.4168192) q[3];
sx q[3];
rz(-1.4409325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.440381) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-0.70297855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5876578) q[0];
sx q[0];
rz(-2.8371187) q[0];
sx q[0];
rz(-1.1712043) q[0];
rz(-pi) q[1];
rz(-0.59962745) q[2];
sx q[2];
rz(-0.58279524) q[2];
sx q[2];
rz(1.8011013) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86537251) q[1];
sx q[1];
rz(-0.11905383) q[1];
sx q[1];
rz(2.7369376) q[1];
x q[2];
rz(2.968077) q[3];
sx q[3];
rz(-2.1870038) q[3];
sx q[3];
rz(-0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26677033) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(-0.91397816) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-0.77600586) q[0];
rz(-1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(0.56328303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8877836) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(-0.33872351) q[0];
x q[1];
rz(-0.69022501) q[2];
sx q[2];
rz(-1.8340655) q[2];
sx q[2];
rz(-3.1379679) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5717585) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(0.024072577) q[1];
x q[2];
rz(-2.8068845) q[3];
sx q[3];
rz(-0.74581205) q[3];
sx q[3];
rz(1.0583744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48878601) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8673458) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(-2.7903941) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441372) q[0];
sx q[0];
rz(-1.1778957) q[0];
sx q[0];
rz(2.5494954) q[0];
x q[1];
rz(-1.048462) q[2];
sx q[2];
rz(-1.221721) q[2];
sx q[2];
rz(-0.77490865) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1197549) q[1];
sx q[1];
rz(-0.86958414) q[1];
sx q[1];
rz(-0.55418684) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.01708548) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(-2.7929896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(1.2472786) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-2.9124027) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(0.75063467) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.437498) q[0];
sx q[0];
rz(-1.1878345) q[0];
sx q[0];
rz(1.0136481) q[0];
x q[1];
rz(0.46361228) q[2];
sx q[2];
rz(-1.3890651) q[2];
sx q[2];
rz(2.6169427) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8667824) q[1];
sx q[1];
rz(-0.85046235) q[1];
sx q[1];
rz(-1.5809098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.552794) q[3];
sx q[3];
rz(-1.4466803) q[3];
sx q[3];
rz(-2.0037946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951185) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(-2.0948998) q[0];
rz(-1.5294317) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(2.7244862) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55506599) q[0];
sx q[0];
rz(-1.546372) q[0];
sx q[0];
rz(1.3506372) q[0];
rz(-pi) q[1];
rz(0.30714005) q[2];
sx q[2];
rz(-1.2899613) q[2];
sx q[2];
rz(-1.839523) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6542146) q[1];
sx q[1];
rz(-2.0301135) q[1];
sx q[1];
rz(-0.76141255) q[1];
x q[2];
rz(-0.38325558) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(-2.3243429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1356915) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(0.55316365) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(0.52545351) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32507867) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(1.4272383) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(3.0292125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.534879) q[0];
sx q[0];
rz(-2.0438072) q[0];
sx q[0];
rz(2.9472449) q[0];
x q[1];
rz(-0.33371146) q[2];
sx q[2];
rz(-1.4881655) q[2];
sx q[2];
rz(0.73252788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1689414) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(0.45291839) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3120679) q[3];
sx q[3];
rz(-1.4530164) q[3];
sx q[3];
rz(-0.044101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(1.3486264) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(2.2299178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1393136) q[0];
sx q[0];
rz(-0.3542491) q[0];
sx q[0];
rz(2.2646963) q[0];
x q[1];
rz(-1.6976835) q[2];
sx q[2];
rz(-1.03089) q[2];
sx q[2];
rz(1.9829696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94757838) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(-0.46781637) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0613641) q[3];
sx q[3];
rz(-2.053223) q[3];
sx q[3];
rz(-2.3232943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(2.7695079) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6045195) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(1.3634691) q[0];
rz(-0.66843372) q[2];
sx q[2];
rz(-2.4032776) q[2];
sx q[2];
rz(2.7065606) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.806554) q[1];
sx q[1];
rz(-1.6400669) q[1];
sx q[1];
rz(-0.075242234) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98202242) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(-2.2002937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(1.0958825) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(0.12116184) q[2];
sx q[2];
rz(-2.0222752) q[2];
sx q[2];
rz(-3.0565699) q[2];
rz(2.1429569) q[3];
sx q[3];
rz(-1.5012267) q[3];
sx q[3];
rz(-2.5517626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
