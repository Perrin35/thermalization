OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8283451) q[0];
sx q[0];
rz(-0.77594835) q[0];
sx q[0];
rz(-1.1021855) q[0];
rz(-1.518353) q[1];
sx q[1];
rz(-2.0145388) q[1];
sx q[1];
rz(2.9161646) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0505053) q[0];
sx q[0];
rz(-0.69982547) q[0];
sx q[0];
rz(-0.096629337) q[0];
rz(-pi) q[1];
rz(1.7182452) q[2];
sx q[2];
rz(-2.1814697) q[2];
sx q[2];
rz(-1.5491279) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0590093) q[1];
sx q[1];
rz(-1.2066926) q[1];
sx q[1];
rz(-0.28739448) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5479911) q[3];
sx q[3];
rz(-2.5656325) q[3];
sx q[3];
rz(-2.0570448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48308358) q[2];
sx q[2];
rz(-2.4675214) q[2];
sx q[2];
rz(0.020641208) q[2];
rz(-0.189273) q[3];
sx q[3];
rz(-0.34957379) q[3];
sx q[3];
rz(1.6674204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78057688) q[0];
sx q[0];
rz(-2.3605232) q[0];
sx q[0];
rz(2.1625157) q[0];
rz(3.0304404) q[1];
sx q[1];
rz(-2.4039098) q[1];
sx q[1];
rz(2.0806064) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1887061) q[0];
sx q[0];
rz(-2.0659819) q[0];
sx q[0];
rz(-2.1006891) q[0];
x q[1];
rz(0.13745549) q[2];
sx q[2];
rz(-1.571621) q[2];
sx q[2];
rz(-1.0975791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0506405) q[1];
sx q[1];
rz(-2.0555858) q[1];
sx q[1];
rz(-1.8043955) q[1];
rz(-pi) q[2];
rz(1.4270328) q[3];
sx q[3];
rz(-1.9562436) q[3];
sx q[3];
rz(0.072162554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4336808) q[2];
sx q[2];
rz(-1.3206864) q[2];
sx q[2];
rz(-0.39805463) q[2];
rz(1.4237283) q[3];
sx q[3];
rz(-2.8545696) q[3];
sx q[3];
rz(-2.6557693) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61755359) q[0];
sx q[0];
rz(-2.6300639) q[0];
sx q[0];
rz(1.0523354) q[0];
rz(0.44318336) q[1];
sx q[1];
rz(-2.1801703) q[1];
sx q[1];
rz(0.65394941) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4674462) q[0];
sx q[0];
rz(-0.83195126) q[0];
sx q[0];
rz(1.3988926) q[0];
rz(-pi) q[1];
rz(2.1939799) q[2];
sx q[2];
rz(-1.6296367) q[2];
sx q[2];
rz(0.72966444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2005077) q[1];
sx q[1];
rz(-1.9882747) q[1];
sx q[1];
rz(-0.79065506) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4846333) q[3];
sx q[3];
rz(-0.28626075) q[3];
sx q[3];
rz(-0.036542758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37086481) q[2];
sx q[2];
rz(-2.3720522) q[2];
sx q[2];
rz(2.9097843) q[2];
rz(1.2109463) q[3];
sx q[3];
rz(-1.959266) q[3];
sx q[3];
rz(0.32538357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66832191) q[0];
sx q[0];
rz(-2.7372161) q[0];
sx q[0];
rz(-2.7647198) q[0];
rz(-0.51141524) q[1];
sx q[1];
rz(-1.2579974) q[1];
sx q[1];
rz(2.2976141) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653042) q[0];
sx q[0];
rz(-1.6805959) q[0];
sx q[0];
rz(3.1240433) q[0];
x q[1];
rz(-1.2449712) q[2];
sx q[2];
rz(-2.5375536) q[2];
sx q[2];
rz(0.86261311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7755345) q[1];
sx q[1];
rz(-1.0909799) q[1];
sx q[1];
rz(-2.5262031) q[1];
x q[2];
rz(1.549349) q[3];
sx q[3];
rz(-2.0923695) q[3];
sx q[3];
rz(2.4873268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9722998) q[2];
sx q[2];
rz(-2.1268714) q[2];
sx q[2];
rz(-2.443552) q[2];
rz(-1.1997148) q[3];
sx q[3];
rz(-0.90568393) q[3];
sx q[3];
rz(-2.5419295) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9508764) q[0];
sx q[0];
rz(-0.27442014) q[0];
sx q[0];
rz(-0.12572591) q[0];
rz(-1.0267286) q[1];
sx q[1];
rz(-1.7544361) q[1];
sx q[1];
rz(0.52621192) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95740283) q[0];
sx q[0];
rz(-2.0483575) q[0];
sx q[0];
rz(1.6543426) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7650928) q[2];
sx q[2];
rz(-0.70617968) q[2];
sx q[2];
rz(-2.6841109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6263585) q[1];
sx q[1];
rz(-0.39455104) q[1];
sx q[1];
rz(2.872353) q[1];
x q[2];
rz(0.69198841) q[3];
sx q[3];
rz(-0.92433911) q[3];
sx q[3];
rz(2.5916168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0478263) q[2];
sx q[2];
rz(-0.83790773) q[2];
sx q[2];
rz(-2.8313336) q[2];
rz(3.1066762) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(-0.78970277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0652311) q[0];
sx q[0];
rz(-1.4820453) q[0];
sx q[0];
rz(2.7705627) q[0];
rz(-3.0767483) q[1];
sx q[1];
rz(-0.82227451) q[1];
sx q[1];
rz(-1.8745905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93343914) q[0];
sx q[0];
rz(-1.9133718) q[0];
sx q[0];
rz(-2.2517363) q[0];
rz(0.94159884) q[2];
sx q[2];
rz(-1.8591188) q[2];
sx q[2];
rz(1.1278111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1888926) q[1];
sx q[1];
rz(-1.5363534) q[1];
sx q[1];
rz(0.39826213) q[1];
rz(2.606845) q[3];
sx q[3];
rz(-2.6125997) q[3];
sx q[3];
rz(2.3803866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0176257) q[2];
sx q[2];
rz(-1.4767246) q[2];
sx q[2];
rz(1.5383447) q[2];
rz(-2.7247143) q[3];
sx q[3];
rz(-0.50691253) q[3];
sx q[3];
rz(0.1703593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012731) q[0];
sx q[0];
rz(-0.23389255) q[0];
sx q[0];
rz(-0.42863578) q[0];
rz(1.1242695) q[1];
sx q[1];
rz(-1.602403) q[1];
sx q[1];
rz(0.7695778) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3807629) q[0];
sx q[0];
rz(-1.5676982) q[0];
sx q[0];
rz(-3.1397545) q[0];
rz(-pi) q[1];
rz(2.6181302) q[2];
sx q[2];
rz(-2.0085196) q[2];
sx q[2];
rz(1.7132393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0391072) q[1];
sx q[1];
rz(-1.277248) q[1];
sx q[1];
rz(2.5204646) q[1];
x q[2];
rz(2.8002987) q[3];
sx q[3];
rz(-1.0847632) q[3];
sx q[3];
rz(-2.7529972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4838685) q[2];
sx q[2];
rz(-1.8572073) q[2];
sx q[2];
rz(3.0460975) q[2];
rz(-0.5419845) q[3];
sx q[3];
rz(-0.46613765) q[3];
sx q[3];
rz(-0.12921648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1456881) q[0];
sx q[0];
rz(-0.93038428) q[0];
sx q[0];
rz(-3.0297739) q[0];
rz(-1.8226786) q[1];
sx q[1];
rz(-1.2448064) q[1];
sx q[1];
rz(1.3056614) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5163286) q[0];
sx q[0];
rz(-2.5808307) q[0];
sx q[0];
rz(-1.8050342) q[0];
x q[1];
rz(0.28152604) q[2];
sx q[2];
rz(-1.7067995) q[2];
sx q[2];
rz(-2.1815262) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9255786) q[1];
sx q[1];
rz(-2.6655201) q[1];
sx q[1];
rz(-2.4904597) q[1];
rz(-pi) q[2];
rz(1.1505125) q[3];
sx q[3];
rz(-2.0788361) q[3];
sx q[3];
rz(-0.86364323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32885113) q[2];
sx q[2];
rz(-1.988827) q[2];
sx q[2];
rz(-1.1197155) q[2];
rz(2.8730734) q[3];
sx q[3];
rz(-1.7374141) q[3];
sx q[3];
rz(1.1855116) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28065228) q[0];
sx q[0];
rz(-2.379731) q[0];
sx q[0];
rz(-0.60894388) q[0];
rz(0.87743419) q[1];
sx q[1];
rz(-2.139822) q[1];
sx q[1];
rz(-2.529349) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37892482) q[0];
sx q[0];
rz(-0.84419709) q[0];
sx q[0];
rz(-0.55171497) q[0];
rz(-pi) q[1];
rz(0.52895542) q[2];
sx q[2];
rz(-0.63807708) q[2];
sx q[2];
rz(0.70032373) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22837199) q[1];
sx q[1];
rz(-2.2971584) q[1];
sx q[1];
rz(-2.5102763) q[1];
x q[2];
rz(2.6004535) q[3];
sx q[3];
rz(-1.5039467) q[3];
sx q[3];
rz(-0.56136405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0683384) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(-2.9541328) q[2];
rz(-1.1120262) q[3];
sx q[3];
rz(-2.2234962) q[3];
sx q[3];
rz(2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40792313) q[0];
sx q[0];
rz(-2.3275571) q[0];
sx q[0];
rz(-2.8975876) q[0];
rz(-1.6581274) q[1];
sx q[1];
rz(-1.6226059) q[1];
sx q[1];
rz(0.77267486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059390776) q[0];
sx q[0];
rz(-2.0652949) q[0];
sx q[0];
rz(-0.54532651) q[0];
rz(-pi) q[1];
rz(2.8171854) q[2];
sx q[2];
rz(-1.7008616) q[2];
sx q[2];
rz(2.7489894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6405787) q[1];
sx q[1];
rz(-1.2440885) q[1];
sx q[1];
rz(-1.4540637) q[1];
rz(2.3546702) q[3];
sx q[3];
rz(-0.95810181) q[3];
sx q[3];
rz(2.4502692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8030168) q[2];
sx q[2];
rz(-0.96119857) q[2];
sx q[2];
rz(-2.4147066) q[2];
rz(-1.1072985) q[3];
sx q[3];
rz(-2.6328583) q[3];
sx q[3];
rz(-2.1342962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0016639391) q[0];
sx q[0];
rz(-1.7353084) q[0];
sx q[0];
rz(-1.1575862) q[0];
rz(2.6955556) q[1];
sx q[1];
rz(-1.0855433) q[1];
sx q[1];
rz(-0.6378508) q[1];
rz(-2.73861) q[2];
sx q[2];
rz(-1.4558515) q[2];
sx q[2];
rz(0.2414138) q[2];
rz(0.31728716) q[3];
sx q[3];
rz(-2.696456) q[3];
sx q[3];
rz(1.0831931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
