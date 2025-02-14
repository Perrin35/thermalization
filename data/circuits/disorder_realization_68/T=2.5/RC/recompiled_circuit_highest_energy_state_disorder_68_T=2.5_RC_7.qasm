OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6990868) q[0];
sx q[0];
rz(1.2966172) q[0];
sx q[0];
rz(7.9820493) q[0];
rz(-0.35248414) q[1];
sx q[1];
rz(-0.52630693) q[1];
sx q[1];
rz(1.7736645) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5718846) q[0];
sx q[0];
rz(-1.8284599) q[0];
sx q[0];
rz(2.7710157) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5330731) q[2];
sx q[2];
rz(-1.9262656) q[2];
sx q[2];
rz(1.0743574) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.96208159) q[1];
sx q[1];
rz(-1.700811) q[1];
sx q[1];
rz(-0.52501734) q[1];
rz(-0.75482224) q[3];
sx q[3];
rz(-1.1140545) q[3];
sx q[3];
rz(0.43908248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0736531) q[2];
sx q[2];
rz(-1.0509793) q[2];
sx q[2];
rz(1.2037207) q[2];
rz(-0.97483557) q[3];
sx q[3];
rz(-2.1741512) q[3];
sx q[3];
rz(1.2459285) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75342733) q[0];
sx q[0];
rz(-2.0797256) q[0];
sx q[0];
rz(0.88428307) q[0];
rz(2.4321804) q[1];
sx q[1];
rz(-1.8953036) q[1];
sx q[1];
rz(0.90219227) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3002533) q[0];
sx q[0];
rz(-2.640199) q[0];
sx q[0];
rz(-2.7031913) q[0];
x q[1];
rz(-1.4231176) q[2];
sx q[2];
rz(-2.5294249) q[2];
sx q[2];
rz(-1.3308805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1556405) q[1];
sx q[1];
rz(-1.8690171) q[1];
sx q[1];
rz(1.8267426) q[1];
x q[2];
rz(1.7092) q[3];
sx q[3];
rz(-1.6507727) q[3];
sx q[3];
rz(-0.13025912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0238637) q[2];
sx q[2];
rz(-0.29568299) q[2];
sx q[2];
rz(-1.6531061) q[2];
rz(-1.9272517) q[3];
sx q[3];
rz(-1.4597273) q[3];
sx q[3];
rz(1.2795718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0073256) q[0];
sx q[0];
rz(-1.7420344) q[0];
sx q[0];
rz(-2.9752327) q[0];
rz(2.4436489) q[1];
sx q[1];
rz(-2.3614387) q[1];
sx q[1];
rz(1.4770329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1159984) q[0];
sx q[0];
rz(-0.4420155) q[0];
sx q[0];
rz(1.5777977) q[0];
rz(2.8700902) q[2];
sx q[2];
rz(-1.2997441) q[2];
sx q[2];
rz(-2.0663911) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0863667) q[1];
sx q[1];
rz(-0.7387195) q[1];
sx q[1];
rz(2.2768873) q[1];
rz(-pi) q[2];
rz(0.03607492) q[3];
sx q[3];
rz(-2.6429308) q[3];
sx q[3];
rz(-2.7447678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6537689) q[2];
sx q[2];
rz(-1.1710125) q[2];
sx q[2];
rz(2.5269395) q[2];
rz(1.4608308) q[3];
sx q[3];
rz(-1.5644282) q[3];
sx q[3];
rz(3.0998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6093269) q[0];
sx q[0];
rz(-1.8121239) q[0];
sx q[0];
rz(1.5186658) q[0];
rz(-2.3329349) q[1];
sx q[1];
rz(-1.8452019) q[1];
sx q[1];
rz(-3.0928968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2653462) q[0];
sx q[0];
rz(-2.2314203) q[0];
sx q[0];
rz(-1.9991849) q[0];
x q[1];
rz(-0.38191958) q[2];
sx q[2];
rz(-1.6798229) q[2];
sx q[2];
rz(0.38709059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21203707) q[1];
sx q[1];
rz(-2.0494645) q[1];
sx q[1];
rz(-2.5655866) q[1];
rz(-2.7475098) q[3];
sx q[3];
rz(-1.0021804) q[3];
sx q[3];
rz(2.4134825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0398756) q[2];
sx q[2];
rz(-2.3287435) q[2];
sx q[2];
rz(1.2561049) q[2];
rz(0.060955437) q[3];
sx q[3];
rz(-2.239949) q[3];
sx q[3];
rz(-2.2140293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561279) q[0];
sx q[0];
rz(-2.0138854) q[0];
sx q[0];
rz(-0.10966478) q[0];
rz(-1.0288641) q[1];
sx q[1];
rz(-1.2867462) q[1];
sx q[1];
rz(-1.5493772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3105811) q[0];
sx q[0];
rz(-0.54052522) q[0];
sx q[0];
rz(-2.5700918) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36269168) q[2];
sx q[2];
rz(-1.0981993) q[2];
sx q[2];
rz(-2.9031799) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8126545) q[1];
sx q[1];
rz(-0.65954406) q[1];
sx q[1];
rz(0.75700827) q[1];
rz(0.83170931) q[3];
sx q[3];
rz(-0.050938923) q[3];
sx q[3];
rz(2.4866631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5900383) q[2];
sx q[2];
rz(-2.8503214) q[2];
sx q[2];
rz(0.48198286) q[2];
rz(2.7216952) q[3];
sx q[3];
rz(-1.0651383) q[3];
sx q[3];
rz(-0.70986748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14966203) q[0];
sx q[0];
rz(-0.88352942) q[0];
sx q[0];
rz(2.4290207) q[0];
rz(-2.271671) q[1];
sx q[1];
rz(-0.37418071) q[1];
sx q[1];
rz(0.023699997) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837304) q[0];
sx q[0];
rz(-1.0548158) q[0];
sx q[0];
rz(-1.3829154) q[0];
rz(-pi) q[1];
rz(-2.4871102) q[2];
sx q[2];
rz(-2.3238782) q[2];
sx q[2];
rz(0.47374642) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5795676) q[1];
sx q[1];
rz(-1.4050086) q[1];
sx q[1];
rz(3.1400984) q[1];
x q[2];
rz(0.48051254) q[3];
sx q[3];
rz(-1.7578205) q[3];
sx q[3];
rz(-2.9399237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80453834) q[2];
sx q[2];
rz(-0.64547515) q[2];
sx q[2];
rz(-1.1402593) q[2];
rz(1.987847) q[3];
sx q[3];
rz(-1.2146344) q[3];
sx q[3];
rz(1.8203189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7158647) q[0];
sx q[0];
rz(-1.7540997) q[0];
sx q[0];
rz(2.7737889) q[0];
rz(1.6416719) q[1];
sx q[1];
rz(-2.2190614) q[1];
sx q[1];
rz(-2.0466764) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2832868) q[0];
sx q[0];
rz(-1.0890766) q[0];
sx q[0];
rz(2.879309) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7000654) q[2];
sx q[2];
rz(-1.0241) q[2];
sx q[2];
rz(1.136029) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81568372) q[1];
sx q[1];
rz(-1.5885873) q[1];
sx q[1];
rz(-0.23577006) q[1];
rz(0.021281042) q[3];
sx q[3];
rz(-1.9742734) q[3];
sx q[3];
rz(-1.9627987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0992004) q[2];
sx q[2];
rz(-1.7902057) q[2];
sx q[2];
rz(-3.0858827) q[2];
rz(2.4646711) q[3];
sx q[3];
rz(-2.5261295) q[3];
sx q[3];
rz(0.87853471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95213503) q[0];
sx q[0];
rz(-0.5624693) q[0];
sx q[0];
rz(2.1761555) q[0];
rz(1.2646487) q[1];
sx q[1];
rz(-2.3955884) q[1];
sx q[1];
rz(2.8146578) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.801689) q[0];
sx q[0];
rz(-1.2210232) q[0];
sx q[0];
rz(-3.0236493) q[0];
x q[1];
rz(-1.5816273) q[2];
sx q[2];
rz(-1.1248661) q[2];
sx q[2];
rz(1.2496834) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2911243) q[1];
sx q[1];
rz(-1.3936491) q[1];
sx q[1];
rz(1.7825104) q[1];
rz(-pi) q[2];
rz(3.0013004) q[3];
sx q[3];
rz(-2.6888118) q[3];
sx q[3];
rz(1.0990395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7586729) q[2];
sx q[2];
rz(-0.70432538) q[2];
sx q[2];
rz(-0.79603377) q[2];
rz(-0.087513611) q[3];
sx q[3];
rz(-2.5329866) q[3];
sx q[3];
rz(-2.2039738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7079119) q[0];
sx q[0];
rz(-1.7681363) q[0];
sx q[0];
rz(-0.32980907) q[0];
rz(3.0682796) q[1];
sx q[1];
rz(-0.70711702) q[1];
sx q[1];
rz(1.0940394) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1739376) q[0];
sx q[0];
rz(-2.6693845) q[0];
sx q[0];
rz(-3.0595991) q[0];
x q[1];
rz(1.4118926) q[2];
sx q[2];
rz(-2.7393118) q[2];
sx q[2];
rz(-0.91644588) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7787012) q[1];
sx q[1];
rz(-0.56462949) q[1];
sx q[1];
rz(2.8933011) q[1];
rz(-pi) q[2];
x q[2];
rz(0.039765914) q[3];
sx q[3];
rz(-1.13969) q[3];
sx q[3];
rz(-2.414961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0988934) q[2];
sx q[2];
rz(-2.3697479) q[2];
sx q[2];
rz(1.137255) q[2];
rz(2.5134261) q[3];
sx q[3];
rz(-0.38574949) q[3];
sx q[3];
rz(1.0073193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6657418) q[0];
sx q[0];
rz(-2.6306212) q[0];
sx q[0];
rz(-2.0487336) q[0];
rz(1.8765556) q[1];
sx q[1];
rz(-1.8362412) q[1];
sx q[1];
rz(-2.1078033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7395679) q[0];
sx q[0];
rz(-1.0453859) q[0];
sx q[0];
rz(-0.55243405) q[0];
rz(1.3273507) q[2];
sx q[2];
rz(-0.82591354) q[2];
sx q[2];
rz(-1.5856397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7180192) q[1];
sx q[1];
rz(-1.9818273) q[1];
sx q[1];
rz(0.6842821) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.035461144) q[3];
sx q[3];
rz(-1.4322114) q[3];
sx q[3];
rz(1.1690804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8826302) q[2];
sx q[2];
rz(-2.4179103) q[2];
sx q[2];
rz(-2.9222729) q[2];
rz(-0.80260459) q[3];
sx q[3];
rz(-1.8290627) q[3];
sx q[3];
rz(2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6312859) q[0];
sx q[0];
rz(-1.318537) q[0];
sx q[0];
rz(1.0986811) q[0];
rz(0.17768271) q[1];
sx q[1];
rz(-2.1251353) q[1];
sx q[1];
rz(2.9716117) q[1];
rz(2.2732757) q[2];
sx q[2];
rz(-1.5219896) q[2];
sx q[2];
rz(0.513214) q[2];
rz(-1.7614638) q[3];
sx q[3];
rz(-1.4545961) q[3];
sx q[3];
rz(-2.0300099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
