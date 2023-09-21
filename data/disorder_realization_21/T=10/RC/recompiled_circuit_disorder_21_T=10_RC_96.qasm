OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0857467) q[0];
sx q[0];
rz(-0.081781713) q[0];
sx q[0];
rz(-2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0599521) q[0];
sx q[0];
rz(-1.8917221) q[0];
sx q[0];
rz(0.089846213) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88959496) q[2];
sx q[2];
rz(-1.0729562) q[2];
sx q[2];
rz(1.3537458) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1134125) q[1];
sx q[1];
rz(-0.98787381) q[1];
sx q[1];
rz(2.7229573) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93482165) q[3];
sx q[3];
rz(-1.1879731) q[3];
sx q[3];
rz(0.32360199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6364608) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(-2.5855529) q[2];
rz(-2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(2.1957943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(0.15727501) q[0];
rz(-2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(3.0325586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63576525) q[0];
sx q[0];
rz(-0.286239) q[0];
sx q[0];
rz(0.92606996) q[0];
x q[1];
rz(1.7686339) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(2.3193662) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82169689) q[1];
sx q[1];
rz(-2.3488455) q[1];
sx q[1];
rz(-2.4577623) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3974959) q[3];
sx q[3];
rz(-2.0850075) q[3];
sx q[3];
rz(-1.7931995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(-2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(-1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064421244) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(1.3431312) q[0];
rz(-0.24761565) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(2.6599191) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0945064) q[0];
sx q[0];
rz(-1.1990093) q[0];
sx q[0];
rz(1.0466172) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8086497) q[2];
sx q[2];
rz(-0.98368401) q[2];
sx q[2];
rz(-1.3981896) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18322769) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(-2.4135713) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62044575) q[3];
sx q[3];
rz(-0.9471604) q[3];
sx q[3];
rz(2.7321531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3383011) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-2.5894077) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(-1.2447371) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23403215) q[0];
sx q[0];
rz(-1.910277) q[0];
sx q[0];
rz(1.5690804) q[0];
rz(-pi) q[1];
rz(-2.142698) q[2];
sx q[2];
rz(-2.9036387) q[2];
sx q[2];
rz(-1.9902802) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2326395) q[1];
sx q[1];
rz(-0.5852355) q[1];
sx q[1];
rz(-1.1802243) q[1];
x q[2];
rz(0.57085412) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(1.966147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(1.4771279) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(0.081469014) q[0];
rz(3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.5030456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016184729) q[0];
sx q[0];
rz(-2.5446919) q[0];
sx q[0];
rz(1.7993268) q[0];
rz(0.083085255) q[2];
sx q[2];
rz(-1.6787046) q[2];
sx q[2];
rz(1.2667058) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45080966) q[1];
sx q[1];
rz(-0.329031) q[1];
sx q[1];
rz(-1.8915218) q[1];
rz(0.53346975) q[3];
sx q[3];
rz(-2.5662078) q[3];
sx q[3];
rz(1.7076147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(-1.2379237) q[2];
rz(1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5313107) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(0.26671985) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(2.3430603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0621588) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(0.14607231) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59136765) q[2];
sx q[2];
rz(-2.5632576) q[2];
sx q[2];
rz(2.8665198) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9763899) q[1];
sx q[1];
rz(-0.573728) q[1];
sx q[1];
rz(1.6673253) q[1];
x q[2];
rz(-0.800662) q[3];
sx q[3];
rz(-2.0293529) q[3];
sx q[3];
rz(0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(-1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(-0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325608) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(-0.60638705) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(0.46404776) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0021792) q[0];
sx q[0];
rz(-0.78939775) q[0];
sx q[0];
rz(-2.1887357) q[0];
x q[1];
rz(-1.8115933) q[2];
sx q[2];
rz(-1.5511998) q[2];
sx q[2];
rz(-2.7407676) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0276427) q[1];
sx q[1];
rz(-1.0877891) q[1];
sx q[1];
rz(-2.7999858) q[1];
rz(-pi) q[2];
rz(-2.8125698) q[3];
sx q[3];
rz(-0.24917069) q[3];
sx q[3];
rz(2.3014625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2074034) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(-2.8835473) q[2];
rz(1.1856273) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(-2.7602957) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.7000748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9863319) q[0];
sx q[0];
rz(-3.0763456) q[0];
sx q[0];
rz(0.23373993) q[0];
rz(-pi) q[1];
rz(2.5289815) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(2.999246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4691094) q[1];
sx q[1];
rz(-1.6961349) q[1];
sx q[1];
rz(0.13087665) q[1];
rz(-3.0184047) q[3];
sx q[3];
rz(-2.8052969) q[3];
sx q[3];
rz(-1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1429446) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(0.22658919) q[2];
rz(2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.3990078) q[0];
rz(2.8245068) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(2.1549966) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164455) q[0];
sx q[0];
rz(-1.4129606) q[0];
sx q[0];
rz(-2.6437003) q[0];
rz(0.87395845) q[2];
sx q[2];
rz(-1.9947589) q[2];
sx q[2];
rz(1.6800113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36531891) q[1];
sx q[1];
rz(-0.96818189) q[1];
sx q[1];
rz(-1.3990632) q[1];
rz(-pi) q[2];
rz(1.1986198) q[3];
sx q[3];
rz(-2.6937006) q[3];
sx q[3];
rz(0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2150779) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(0.40965664) q[2];
rz(2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(-2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7423994) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(-1.4051399) q[0];
rz(-2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(1.7260889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50354276) q[0];
sx q[0];
rz(-1.8068131) q[0];
sx q[0];
rz(-0.34061265) q[0];
x q[1];
rz(-2.892898) q[2];
sx q[2];
rz(-1.2173614) q[2];
sx q[2];
rz(2.8949708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43511697) q[1];
sx q[1];
rz(-0.26285989) q[1];
sx q[1];
rz(-0.95107066) q[1];
x q[2];
rz(2.9726082) q[3];
sx q[3];
rz(-1.0645234) q[3];
sx q[3];
rz(-0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71904174) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(0.96578807) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(-0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.60349764) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(-2.8339236) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(1.7514501) q[2];
sx q[2];
rz(-1.8271108) q[2];
sx q[2];
rz(1.9432632) q[2];
rz(1.7210759) q[3];
sx q[3];
rz(-1.0192623) q[3];
sx q[3];
rz(-2.7912959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];