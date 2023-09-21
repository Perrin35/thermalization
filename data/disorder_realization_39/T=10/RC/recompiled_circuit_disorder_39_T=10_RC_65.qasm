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
rz(-2.0733641) q[0];
sx q[0];
rz(0.46407035) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(-0.067117604) q[1];
sx q[1];
rz(1.2844515) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0811631) q[0];
sx q[0];
rz(-1.8627177) q[0];
sx q[0];
rz(1.97776) q[0];
x q[1];
rz(0.061878248) q[2];
sx q[2];
rz(-2.6510694) q[2];
sx q[2];
rz(-2.2762736) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1661108) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(-2.6544177) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0331217) q[3];
sx q[3];
rz(-2.9900108) q[3];
sx q[3];
rz(2.3173995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(2.822067) q[2];
rz(-0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-0.67392504) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(2.3449576) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731199) q[0];
sx q[0];
rz(-0.94808775) q[0];
sx q[0];
rz(2.7362583) q[0];
rz(2.566922) q[2];
sx q[2];
rz(-1.7386912) q[2];
sx q[2];
rz(1.4603953) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14256515) q[1];
sx q[1];
rz(-1.5836645) q[1];
sx q[1];
rz(-2.5343115) q[1];
rz(-pi) q[2];
rz(0.91005743) q[3];
sx q[3];
rz(-0.72477341) q[3];
sx q[3];
rz(1.4409325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18156302) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(1.440381) q[0];
rz(0.72021833) q[1];
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
rz(-2.7754606) q[0];
sx q[0];
rz(-1.6876939) q[0];
sx q[0];
rz(-1.2890105) q[0];
x q[1];
rz(2.643232) q[2];
sx q[2];
rz(-1.2549855) q[2];
sx q[2];
rz(-2.3926546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2762201) q[1];
sx q[1];
rz(-3.0225388) q[1];
sx q[1];
rz(2.7369376) q[1];
rz(-pi) q[2];
rz(1.8099144) q[3];
sx q[3];
rz(-2.5044887) q[3];
sx q[3];
rz(-1.0106196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(2.3655868) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(-0.56328303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3455428) q[0];
sx q[0];
rz(-2.4420218) q[0];
sx q[0];
rz(2.0027341) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4513676) q[2];
sx q[2];
rz(-1.3075271) q[2];
sx q[2];
rz(-0.0036247591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5717585) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(-3.1175201) q[1];
rz(-pi) q[2];
rz(2.8068845) q[3];
sx q[3];
rz(-2.3957806) q[3];
sx q[3];
rz(-2.0832182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48878601) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-2.9768067) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8673458) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(0.85246032) q[0];
rz(-0.35119855) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1253818) q[0];
sx q[0];
rz(-2.1124766) q[0];
sx q[0];
rz(1.1075695) q[0];
rz(-pi) q[1];
rz(-0.94050546) q[2];
sx q[2];
rz(-0.61912196) q[2];
sx q[2];
rz(1.3319912) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.1197549) q[1];
sx q[1];
rz(-2.2720085) q[1];
sx q[1];
rz(2.5874058) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4200053) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(1.23502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44624415) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.86876774) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(0.75063467) q[0];
rz(2.0320832) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.437498) q[0];
sx q[0];
rz(-1.9537582) q[0];
sx q[0];
rz(-1.0136481) q[0];
rz(2.7517031) q[2];
sx q[2];
rz(-2.6460558) q[2];
sx q[2];
rz(-1.748566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8514511) q[1];
sx q[1];
rz(-0.72039225) q[1];
sx q[1];
rz(-0.011522567) q[1];
rz(-1.4218876) q[3];
sx q[3];
rz(-0.9871261) q[3];
sx q[3];
rz(-0.35051171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(-0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951185) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(-1.0466928) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(0.41710645) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55506599) q[0];
sx q[0];
rz(-1.546372) q[0];
sx q[0];
rz(-1.7909554) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8646556) q[2];
sx q[2];
rz(-1.2760578) q[2];
sx q[2];
rz(0.35640946) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
rz(0.38325558) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(-0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(-2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32507867) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(-0.11238012) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6067137) q[0];
sx q[0];
rz(-2.0438072) q[0];
sx q[0];
rz(2.9472449) q[0];
x q[1];
rz(-2.8939395) q[2];
sx q[2];
rz(-2.7981749) q[2];
sx q[2];
rz(-1.0719971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8466859) q[1];
sx q[1];
rz(-1.9125738) q[1];
sx q[1];
rz(0.81975598) q[1];
rz(-pi) q[2];
rz(-1.3120679) q[3];
sx q[3];
rz(-1.4530164) q[3];
sx q[3];
rz(3.0974914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0743951) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-3.1066185) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-2.2299178) q[1];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1107026) q[2];
sx q[2];
rz(1.1586231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72494353) q[1];
sx q[1];
rz(-2.0282201) q[1];
sx q[1];
rz(-1.7979421) q[1];
rz(-pi) q[2];
rz(-0.73260143) q[3];
sx q[3];
rz(-0.67389518) q[3];
sx q[3];
rz(-0.037308824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
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
rz(-0.77722469) q[1];
sx q[1];
rz(1.7262329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088179528) q[0];
sx q[0];
rz(-1.3706494) q[0];
sx q[0];
rz(-2.8739909) q[0];
rz(0.66843372) q[2];
sx q[2];
rz(-0.738315) q[2];
sx q[2];
rz(2.7065606) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2409754) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(-1.6402628) q[1];
rz(2.7578027) q[3];
sx q[3];
rz(-2.0818315) q[3];
sx q[3];
rz(2.8904861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32594484) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-0.20467219) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50080147) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(1.1164222) q[2];
sx q[2];
rz(-1.6797671) q[2];
sx q[2];
rz(-1.4327008) q[2];
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