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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-2.6455542) q[0];
sx q[0];
rz(0.92143671) q[0];
x q[1];
rz(3.0797144) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(-2.2762736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97548188) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(0.48717498) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9908882) q[3];
sx q[3];
rz(-1.5544484) q[3];
sx q[3];
rz(-0.63936641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39711943) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(0.31952566) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(0.67392504) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.5698313) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(-0.79663509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731199) q[0];
sx q[0];
rz(-0.94808775) q[0];
sx q[0];
rz(-0.40533439) q[0];
rz(-2.566922) q[2];
sx q[2];
rz(-1.7386912) q[2];
sx q[2];
rz(-1.4603953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7318774) q[1];
sx q[1];
rz(-0.60740031) q[1];
sx q[1];
rz(0.022547988) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49780952) q[3];
sx q[3];
rz(-1.0199162) q[3];
sx q[3];
rz(2.2450972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(-0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86984533) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7754606) q[0];
sx q[0];
rz(-1.6876939) q[0];
sx q[0];
rz(1.2890105) q[0];
rz(-pi) q[1];
rz(-1.2146644) q[2];
sx q[2];
rz(-2.0424358) q[2];
sx q[2];
rz(0.65442649) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8382593) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(0.10951885) q[1];
x q[2];
rz(-2.1941575) q[3];
sx q[3];
rz(-1.7121592) q[3];
sx q[3];
rz(-0.75368222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(0.91397816) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-2.3655868) q[0];
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
rz(3.027772) q[0];
sx q[0];
rz(-1.8437244) q[0];
sx q[0];
rz(0.91822894) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40043719) q[2];
sx q[2];
rz(-0.73095989) q[2];
sx q[2];
rz(1.2618582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54484425) q[1];
sx q[1];
rz(-1.8425643) q[1];
sx q[1];
rz(1.5640869) q[1];
rz(-1.8654278) q[3];
sx q[3];
rz(-2.2664824) q[3];
sx q[3];
rz(-0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(-0.35119855) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974385) q[0];
sx q[0];
rz(-0.69735202) q[0];
sx q[0];
rz(-2.5028412) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94050546) q[2];
sx q[2];
rz(-2.5224707) q[2];
sx q[2];
rz(1.8096015) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1197549) q[1];
sx q[1];
rz(-2.2720085) q[1];
sx q[1];
rz(-0.55418684) q[1];
rz(-pi) q[2];
rz(0.72158738) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(-1.9065726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86876774) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(1.3060588) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7795777) q[0];
sx q[0];
rz(-1.0581731) q[0];
sx q[0];
rz(-0.4431475) q[0];
rz(-pi) q[1];
rz(-0.38988955) q[2];
sx q[2];
rz(-0.49553686) q[2];
sx q[2];
rz(1.748566) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27481025) q[1];
sx q[1];
rz(-2.2911303) q[1];
sx q[1];
rz(-1.5606828) q[1];
rz(-pi) q[2];
rz(-2.9206198) q[3];
sx q[3];
rz(-0.60022012) q[3];
sx q[3];
rz(0.61629399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7541472) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(-0.35282648) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-1.0466928) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-2.7244862) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0171623) q[0];
sx q[0];
rz(-0.22148795) q[0];
sx q[0];
rz(1.4593967) q[0];
rz(2.8344526) q[2];
sx q[2];
rz(-1.8516314) q[2];
sx q[2];
rz(1.3020696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6542146) q[1];
sx q[1];
rz(-2.0301135) q[1];
sx q[1];
rz(-2.3801801) q[1];
rz(-pi) q[2];
rz(1.8700637) q[3];
sx q[3];
rz(-2.2021658) q[3];
sx q[3];
rz(-2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(-2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0880786) q[0];
sx q[0];
rz(-1.7435762) q[0];
sx q[0];
rz(-2.0515576) q[0];
x q[1];
rz(-0.33371146) q[2];
sx q[2];
rz(-1.6534272) q[2];
sx q[2];
rz(-0.73252788) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5210133) q[1];
sx q[1];
rz(-0.81110209) q[1];
sx q[1];
rz(-1.0902507) q[1];
rz(-pi) q[2];
rz(-2.0039844) q[3];
sx q[3];
rz(-0.28372753) q[3];
sx q[3];
rz(-1.9445436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0743951) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(-1.2049234) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(0.91167489) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.864894) q[0];
sx q[0];
rz(-1.8407341) q[0];
sx q[0];
rz(-0.2322659) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54347221) q[2];
sx q[2];
rz(-1.4620355) q[2];
sx q[2];
rz(-2.6639338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4166491) q[1];
sx q[1];
rz(-1.1133725) q[1];
sx q[1];
rz(1.3436505) q[1];
rz(2.4089912) q[3];
sx q[3];
rz(-0.67389518) q[3];
sx q[3];
rz(-0.037308824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(0.24924499) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0534131) q[0];
sx q[0];
rz(-1.7709433) q[0];
sx q[0];
rz(-0.26760177) q[0];
rz(-1.0572817) q[2];
sx q[2];
rz(-2.1272749) q[2];
sx q[2];
rz(-1.2531812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2409754) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(1.6402628) q[1];
rz(-pi) q[2];
rz(-2.1146718) q[3];
sx q[3];
rz(-1.2380935) q[3];
sx q[3];
rz(1.1247016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32594484) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(-1.7278016) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6407912) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(1.5851371) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(1.3265058) q[2];
sx q[2];
rz(-2.6752224) q[2];
sx q[2];
rz(0.35717076) q[2];
rz(0.08269357) q[3];
sx q[3];
rz(-2.1413998) q[3];
sx q[3];
rz(-1.0257046) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];