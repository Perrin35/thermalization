OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(1.8811037) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0215065) q[0];
sx q[0];
rz(-1.9061631) q[0];
sx q[0];
rz(-2.04727) q[0];
x q[1];
rz(-1.107723) q[2];
sx q[2];
rz(-0.31422868) q[2];
sx q[2];
rz(-2.3067834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94581374) q[1];
sx q[1];
rz(-1.073277) q[1];
sx q[1];
rz(-1.0832018) q[1];
rz(1.9507017) q[3];
sx q[3];
rz(-1.1067179) q[3];
sx q[3];
rz(2.3604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0682003) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(1.1967999) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(1.4555567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14917063) q[0];
sx q[0];
rz(-0.99827168) q[0];
sx q[0];
rz(0.19897977) q[0];
rz(0.20632867) q[2];
sx q[2];
rz(-0.38197877) q[2];
sx q[2];
rz(-2.1260335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93745366) q[1];
sx q[1];
rz(-2.0935645) q[1];
sx q[1];
rz(0.015603113) q[1];
rz(-pi) q[2];
rz(-1.5442113) q[3];
sx q[3];
rz(-1.9823325) q[3];
sx q[3];
rz(-2.6708024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.3519752) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(-0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(-0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(-1.172539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55484178) q[0];
sx q[0];
rz(-1.8790073) q[0];
sx q[0];
rz(0.59535938) q[0];
rz(-1.2433979) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(-0.74795216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.4634358) q[1];
sx q[1];
rz(-2.0883745) q[1];
sx q[1];
rz(0.28097681) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3200687) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(2.1123321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(2.6180843) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2142221) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(-2.5582696) q[0];
rz(1.4857616) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(-1.0852244) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88672968) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(0.32943326) q[1];
rz(-pi) q[2];
x q[2];
rz(1.773049) q[3];
sx q[3];
rz(-0.60086717) q[3];
sx q[3];
rz(-2.3893389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(0.564044) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(-0.55571663) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6600835) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-0.99194828) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6742453) q[0];
sx q[0];
rz(-1.464198) q[0];
sx q[0];
rz(-2.9664413) q[0];
rz(1.8736585) q[2];
sx q[2];
rz(-1.5793243) q[2];
sx q[2];
rz(2.0451343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76584133) q[1];
sx q[1];
rz(-1.6110794) q[1];
sx q[1];
rz(2.8326616) q[1];
rz(3.0117412) q[3];
sx q[3];
rz(-2.3240945) q[3];
sx q[3];
rz(2.0714456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.7927992) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(-1.7165002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7871008) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(1.3160734) q[0];
rz(-pi) q[1];
rz(-1.3155977) q[2];
sx q[2];
rz(-1.4824502) q[2];
sx q[2];
rz(1.7905854) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6188366) q[1];
sx q[1];
rz(-0.4546051) q[1];
sx q[1];
rz(1.418581) q[1];
rz(-pi) q[2];
rz(-0.53374966) q[3];
sx q[3];
rz(-1.0373877) q[3];
sx q[3];
rz(0.48983869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(0.41931835) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-2.3197876) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9200631) q[0];
sx q[0];
rz(-2.1863345) q[0];
sx q[0];
rz(-0.91381844) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0242545) q[2];
sx q[2];
rz(-0.80596906) q[2];
sx q[2];
rz(-2.1794127) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8824749) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(-2.4005753) q[1];
rz(0.58737289) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(1.3164933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(-1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(2.3186671) q[0];
rz(0.30934632) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(1.8364505) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6179498) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(2.6108517) q[0];
x q[1];
rz(0.70456409) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(-1.2003843) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9040363) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(-2.0380286) q[1];
rz(0.023530258) q[3];
sx q[3];
rz(-1.9100034) q[3];
sx q[3];
rz(-0.48285218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(-2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.3867144) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(0.3219147) q[0];
rz(1.5362668) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(-0.70294356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212595) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(2.3938177) q[0];
rz(1.0614971) q[2];
sx q[2];
rz(-1.6054389) q[2];
sx q[2];
rz(0.73355567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1785537) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(-1.8815243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1144936) q[3];
sx q[3];
rz(-2.4745686) q[3];
sx q[3];
rz(-1.4240571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(1.3396324) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43384957) q[0];
sx q[0];
rz(-1.187547) q[0];
sx q[0];
rz(-1.3462523) q[0];
rz(-pi) q[1];
rz(-1.956316) q[2];
sx q[2];
rz(-1.2670994) q[2];
sx q[2];
rz(-2.7591443) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2965282) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(-1.8370085) q[1];
rz(-pi) q[2];
rz(2.3234899) q[3];
sx q[3];
rz(-2.5460498) q[3];
sx q[3];
rz(0.18225741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.1432077) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951915) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(-0.68998228) q[2];
sx q[2];
rz(-2.1524515) q[2];
sx q[2];
rz(-0.099302789) q[2];
rz(1.4680396) q[3];
sx q[3];
rz(-2.4874874) q[3];
sx q[3];
rz(-2.7174674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];