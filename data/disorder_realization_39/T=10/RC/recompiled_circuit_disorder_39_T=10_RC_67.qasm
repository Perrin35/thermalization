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
rz(1.9595454) q[1];
sx q[1];
rz(-0.067117604) q[1];
sx q[1];
rz(1.2844515) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61304898) q[0];
sx q[0];
rz(-1.9595946) q[0];
sx q[0];
rz(0.31624985) q[0];
rz(-pi) q[1];
x q[1];
rz(1.60381) q[2];
sx q[2];
rz(-1.0812949) q[2];
sx q[2];
rz(0.93544338) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1661108) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(-2.6544177) q[1];
rz(-pi) q[2];
rz(-1.554261) q[3];
sx q[3];
rz(-1.7214805) q[3];
sx q[3];
rz(-0.93391234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(0.31952566) q[2];
rz(2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(2.4676676) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.3449576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8684727) q[0];
sx q[0];
rz(-0.94808775) q[0];
sx q[0];
rz(2.7362583) q[0];
x q[1];
rz(-2.8393306) q[2];
sx q[2];
rz(-2.5455591) q[2];
sx q[2];
rz(-0.1421393) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14256515) q[1];
sx q[1];
rz(-1.5836645) q[1];
sx q[1];
rz(-0.60728118) q[1];
rz(2.2315352) q[3];
sx q[3];
rz(-2.4168192) q[3];
sx q[3];
rz(-1.7006601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(0.78655085) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.440381) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(2.4386141) q[1];
rz(-pi/2) q[2];
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
rz(1.2146644) q[2];
sx q[2];
rz(-2.0424358) q[2];
sx q[2];
rz(-0.65442649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8689879) q[1];
sx q[1];
rz(-1.6801949) q[1];
sx q[1];
rz(-1.617856) q[1];
x q[2];
rz(-2.1941575) q[3];
sx q[3];
rz(-1.7121592) q[3];
sx q[3];
rz(-0.75368222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(0.91397816) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(-2.3655868) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(-2.5783096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.027772) q[0];
sx q[0];
rz(-1.8437244) q[0];
sx q[0];
rz(2.2233637) q[0];
x q[1];
rz(-1.9070542) q[2];
sx q[2];
rz(-0.90869892) q[2];
sx q[2];
rz(1.3627571) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56983419) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(0.024072577) q[1];
rz(-pi) q[2];
rz(1.8654278) q[3];
sx q[3];
rz(-2.2664824) q[3];
sx q[3];
rz(0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(-2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1253818) q[0];
sx q[0];
rz(-1.0291161) q[0];
sx q[0];
rz(1.1075695) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0931307) q[2];
sx q[2];
rz(-1.9198717) q[2];
sx q[2];
rz(-2.366684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1197549) q[1];
sx q[1];
rz(-0.86958414) q[1];
sx q[1];
rz(-2.5874058) q[1];
rz(-1.5557628) q[3];
sx q[3];
rz(-0.84926499) q[3];
sx q[3];
rz(2.815747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(-3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86876774) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.437498) q[0];
sx q[0];
rz(-1.9537582) q[0];
sx q[0];
rz(-1.0136481) q[0];
rz(-pi) q[1];
rz(-1.7734217) q[2];
sx q[2];
rz(-2.0261923) q[2];
sx q[2];
rz(2.1855598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.4218876) q[3];
sx q[3];
rz(-0.9871261) q[3];
sx q[3];
rz(2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(-2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-2.0948998) q[0];
rz(1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(2.7244862) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0102651) q[0];
sx q[0];
rz(-1.350704) q[0];
sx q[0];
rz(-3.1165645) q[0];
rz(-pi) q[1];
rz(1.8646556) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(-0.35640946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6228094) q[1];
sx q[1];
rz(-0.86474027) q[1];
sx q[1];
rz(0.6219567) q[1];
rz(-pi) q[2];
rz(0.38325558) q[3];
sx q[3];
rz(-0.68985046) q[3];
sx q[3];
rz(-2.3243429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(-1.7741514) q[2];
rz(0.55316365) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(-0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32507867) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(-1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(-0.11238012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988797) q[0];
sx q[0];
rz(-0.50857022) q[0];
sx q[0];
rz(1.2099427) q[0];
x q[1];
rz(-0.24765315) q[2];
sx q[2];
rz(-2.7981749) q[2];
sx q[2];
rz(-2.0695956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6205794) q[1];
sx q[1];
rz(-2.3304906) q[1];
sx q[1];
rz(-1.0902507) q[1];
rz(-1.1376082) q[3];
sx q[3];
rz(-0.28372753) q[3];
sx q[3];
rz(1.9445436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(1.2049234) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444721) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-2.2299178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.864894) q[0];
sx q[0];
rz(-1.8407341) q[0];
sx q[0];
rz(2.9093268) q[0];
x q[1];
rz(1.4439092) q[2];
sx q[2];
rz(-2.1107026) q[2];
sx q[2];
rz(1.1586231) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94757838) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(2.6737763) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73260143) q[3];
sx q[3];
rz(-0.67389518) q[3];
sx q[3];
rz(0.037308824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-0.5138548) q[2];
sx q[2];
rz(0.24924499) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
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
rz(2.7695079) q[0];
rz(-0.58139873) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6045195) q[0];
sx q[0];
rz(-1.3086638) q[0];
sx q[0];
rz(-1.3634691) q[0];
rz(-pi) q[1];
rz(-1.0572817) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(-1.8884115) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9006173) q[1];
sx q[1];
rz(-1.4957349) q[1];
sx q[1];
rz(1.6402628) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98202242) q[3];
sx q[3];
rz(-2.5128799) q[3];
sx q[3];
rz(-2.2002937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32594484) q[2];
sx q[2];
rz(-0.14020136) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(1.5851371) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(-1.8150868) q[2];
sx q[2];
rz(-2.6752224) q[2];
sx q[2];
rz(0.35717076) q[2];
rz(1.4428044) q[3];
sx q[3];
rz(-2.5656869) q[3];
sx q[3];
rz(2.2681469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
