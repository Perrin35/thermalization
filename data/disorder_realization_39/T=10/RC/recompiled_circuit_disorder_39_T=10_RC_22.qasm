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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0425134) q[0];
sx q[0];
rz(-2.6455542) q[0];
sx q[0];
rz(2.2201559) q[0];
x q[1];
rz(3.0797144) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(-2.2762736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9570934) q[1];
sx q[1];
rz(-0.59809369) q[1];
sx q[1];
rz(-2.4615272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9908882) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(-0.63936641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39711943) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-0.47839034) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4085061) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(2.615036) q[0];
rz(2.5780442) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(2.3449576) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8684727) q[0];
sx q[0];
rz(-0.94808775) q[0];
sx q[0];
rz(2.7362583) q[0];
rz(2.8393306) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(2.9994534) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4097152) q[1];
sx q[1];
rz(-2.5341923) q[1];
sx q[1];
rz(3.1190447) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2315352) q[3];
sx q[3];
rz(-2.4168192) q[3];
sx q[3];
rz(1.7006601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(0.78655085) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.440381) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(2.4386141) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876578) q[0];
sx q[0];
rz(-2.8371187) q[0];
sx q[0];
rz(1.1712043) q[0];
rz(-pi) q[1];
rz(2.5419652) q[2];
sx q[2];
rz(-0.58279524) q[2];
sx q[2];
rz(-1.3404913) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2762201) q[1];
sx q[1];
rz(-3.0225388) q[1];
sx q[1];
rz(-2.7369376) q[1];
x q[2];
rz(1.8099144) q[3];
sx q[3];
rz(-2.5044887) q[3];
sx q[3];
rz(2.1309731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(0.91153574) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-0.77600586) q[0];
rz(-1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8877836) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(2.8028691) q[0];
rz(-pi) q[1];
rz(2.7411555) q[2];
sx q[2];
rz(-2.4106328) q[2];
sx q[2];
rz(-1.2618582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5967484) q[1];
sx q[1];
rz(-1.8425643) q[1];
sx q[1];
rz(1.5775058) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71743439) q[3];
sx q[3];
rz(-1.3460025) q[3];
sx q[3];
rz(2.3790529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48878601) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(2.9768067) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(0.85246032) q[0];
rz(-2.7903941) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(-1.16211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441372) q[0];
sx q[0];
rz(-1.1778957) q[0];
sx q[0];
rz(-2.5494954) q[0];
x q[1];
rz(0.94050546) q[2];
sx q[2];
rz(-0.61912196) q[2];
sx q[2];
rz(1.8096015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8309161) q[1];
sx q[1];
rz(-1.984593) q[1];
sx q[1];
rz(0.78891854) q[1];
rz(-pi) q[2];
x q[2];
rz(0.01708548) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(-0.34860308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.8943141) q[2];
rz(-3.128483) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(-2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86876774) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-2.390958) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7795777) q[0];
sx q[0];
rz(-1.0581731) q[0];
sx q[0];
rz(0.4431475) q[0];
rz(-pi) q[1];
rz(-2.7517031) q[2];
sx q[2];
rz(-0.49553686) q[2];
sx q[2];
rz(-1.748566) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8667824) q[1];
sx q[1];
rz(-2.2911303) q[1];
sx q[1];
rz(-1.5809098) q[1];
x q[2];
rz(-1.7197051) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-0.35051171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(2.5409017) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(-2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(-1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(0.41710645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(2.8344526) q[2];
sx q[2];
rz(-1.8516314) q[2];
sx q[2];
rz(-1.839523) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4579791) q[1];
sx q[1];
rz(-2.2375467) q[1];
sx q[1];
rz(-0.97138202) q[1];
x q[2];
rz(-2.7583371) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(-0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816514) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(-1.1897855) q[0];
rz(1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(0.11238012) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.534879) q[0];
sx q[0];
rz(-1.0977854) q[0];
sx q[0];
rz(-2.9472449) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24765315) q[2];
sx q[2];
rz(-0.34341771) q[2];
sx q[2];
rz(-1.0719971) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97265128) q[1];
sx q[1];
rz(-0.87247889) q[1];
sx q[1];
rz(0.45291839) q[1];
rz(-3.0197969) q[3];
sx q[3];
rz(-1.8276916) q[3];
sx q[3];
rz(-1.6459873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444721) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-2.2299178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23110403) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(1.2937806) q[0];
x q[1];
rz(0.54347221) q[2];
sx q[2];
rz(-1.6795571) q[2];
sx q[2];
rz(2.6639338) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1940143) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(2.6737763) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6058042) q[3];
sx q[3];
rz(-1.1402604) q[3];
sx q[3];
rz(0.99539869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(-0.58139873) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.4153597) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0534131) q[0];
sx q[0];
rz(-1.7709433) q[0];
sx q[0];
rz(0.26760177) q[0];
rz(-2.5214228) q[2];
sx q[2];
rz(-2.0010741) q[2];
sx q[2];
rz(-0.60713965) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9006173) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(1.6402628) q[1];
rz(0.98202242) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(-0.94129896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32594484) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(2.9369205) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(-2.0457101) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-2.0251705) q[2];
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
