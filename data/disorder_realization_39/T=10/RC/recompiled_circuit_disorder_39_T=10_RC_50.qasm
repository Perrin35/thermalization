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
rz(1.9595454) q[1];
sx q[1];
rz(6.2160677) q[1];
sx q[1];
rz(10.709229) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5285437) q[0];
sx q[0];
rz(-1.9595946) q[0];
sx q[0];
rz(-0.31624985) q[0];
rz(-pi) q[1];
rz(-0.061878248) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(0.86531901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7317036) q[1];
sx q[1];
rz(-2.0239502) q[1];
sx q[1];
rz(1.1660006) q[1];
rz(-0.15070446) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(-2.5022262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(-2.3449576) q[1];
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
rz(0.57467069) q[2];
sx q[2];
rz(-1.4029014) q[2];
sx q[2];
rz(-1.6811973) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9990275) q[1];
sx q[1];
rz(-1.5579281) q[1];
sx q[1];
rz(-2.5343115) q[1];
x q[2];
rz(0.96062406) q[3];
sx q[3];
rz(-1.989813) q[3];
sx q[3];
rz(-0.39715365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18156302) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(0.49318796) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(-1.440381) q[0];
rz(2.4213743) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709135) q[0];
sx q[0];
rz(-1.2909856) q[0];
sx q[0];
rz(-3.0199416) q[0];
rz(-pi) q[1];
rz(-1.2146644) q[2];
sx q[2];
rz(-2.0424358) q[2];
sx q[2];
rz(-2.4871662) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2762201) q[1];
sx q[1];
rz(-3.0225388) q[1];
sx q[1];
rz(0.40465506) q[1];
rz(-pi) q[2];
rz(1.3316783) q[3];
sx q[3];
rz(-0.63710391) q[3];
sx q[3];
rz(2.1309731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(0.91153574) q[2];
rz(2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-0.56328303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11382064) q[0];
sx q[0];
rz(-1.2978683) q[0];
sx q[0];
rz(2.2233637) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9070542) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(1.7788356) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.024151) q[1];
sx q[1];
rz(-1.5643331) q[1];
sx q[1];
rz(-0.27177377) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71743439) q[3];
sx q[3];
rz(-1.3460025) q[3];
sx q[3];
rz(-2.3790529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48878601) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(-2.9131043) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(2.1951108) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(0.85246032) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441372) q[0];
sx q[0];
rz(-1.9636969) q[0];
sx q[0];
rz(-2.5494954) q[0];
x q[1];
rz(2.0931307) q[2];
sx q[2];
rz(-1.9198717) q[2];
sx q[2];
rz(0.77490865) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0218378) q[1];
sx q[1];
rz(-2.2720085) q[1];
sx q[1];
rz(-0.55418684) q[1];
x q[2];
rz(1.5858298) q[3];
sx q[3];
rz(-0.84926499) q[3];
sx q[3];
rz(-0.3258457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.8943141) q[2];
rz(-3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-0.22918992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40697843) q[0];
sx q[0];
rz(-2.4771871) q[0];
sx q[0];
rz(0.91974308) q[0];
rz(-2.6779804) q[2];
sx q[2];
rz(-1.7525275) q[2];
sx q[2];
rz(0.52464991) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8667824) q[1];
sx q[1];
rz(-2.2911303) q[1];
sx q[1];
rz(-1.5606828) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2209729) q[3];
sx q[3];
rz(-0.60022012) q[3];
sx q[3];
rz(-2.5252987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(-1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-2.7244862) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0171623) q[0];
sx q[0];
rz(-2.9201047) q[0];
sx q[0];
rz(1.682196) q[0];
rz(1.8646556) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(-0.35640946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4579791) q[1];
sx q[1];
rz(-2.2375467) q[1];
sx q[1];
rz(0.97138202) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38325558) q[3];
sx q[3];
rz(-0.68985046) q[3];
sx q[3];
rz(0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(-2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
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
rz(-2.863527) q[1];
sx q[1];
rz(-0.11238012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9427129) q[0];
sx q[0];
rz(-0.50857022) q[0];
sx q[0];
rz(-1.2099427) q[0];
rz(-2.8078812) q[2];
sx q[2];
rz(-1.4881655) q[2];
sx q[2];
rz(2.4090648) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29490678) q[1];
sx q[1];
rz(-1.2290188) q[1];
sx q[1];
rz(-0.81975598) q[1];
rz(-pi) q[2];
rz(-0.12179575) q[3];
sx q[3];
rz(-1.8276916) q[3];
sx q[3];
rz(-1.4956054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0743951) q[2];
sx q[2];
rz(-1.9125568) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444721) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(-3.1066185) q[0];
rz(0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-2.2299178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2766987) q[0];
sx q[0];
rz(-1.8407341) q[0];
sx q[0];
rz(-2.9093268) q[0];
rz(-pi) q[1];
rz(0.54347221) q[2];
sx q[2];
rz(-1.4620355) q[2];
sx q[2];
rz(0.47765884) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4166491) q[1];
sx q[1];
rz(-1.1133725) q[1];
sx q[1];
rz(1.7979421) q[1];
rz(2.4089912) q[3];
sx q[3];
rz(-0.67389518) q[3];
sx q[3];
rz(-0.037308824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-2.8923477) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(2.7695079) q[0];
rz(-2.5601939) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(1.4153597) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5370731) q[0];
sx q[0];
rz(-1.3086638) q[0];
sx q[0];
rz(-1.3634691) q[0];
rz(-1.0572817) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(1.2531812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6485939) q[1];
sx q[1];
rz(-0.10222888) q[1];
sx q[1];
rz(-0.74536721) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7578027) q[3];
sx q[3];
rz(-2.0818315) q[3];
sx q[3];
rz(2.8904861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-0.20467219) q[2];
rz(-1.7278016) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(1.5564556) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(2.0251705) q[2];
sx q[2];
rz(-1.4618256) q[2];
sx q[2];
rz(1.7088919) q[2];
rz(-1.6987883) q[3];
sx q[3];
rz(-2.5656869) q[3];
sx q[3];
rz(2.2681469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
