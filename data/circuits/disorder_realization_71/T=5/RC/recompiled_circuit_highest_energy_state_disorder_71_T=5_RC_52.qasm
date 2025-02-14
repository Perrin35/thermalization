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
rz(1.1310391) q[0];
sx q[0];
rz(-0.57589632) q[0];
sx q[0];
rz(-2.0146712) q[0];
rz(-0.95070401) q[1];
sx q[1];
rz(-0.60125142) q[1];
sx q[1];
rz(2.8079005) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5883049) q[0];
sx q[0];
rz(-1.9899564) q[0];
sx q[0];
rz(2.9774211) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1154217) q[2];
sx q[2];
rz(-1.8118333) q[2];
sx q[2];
rz(-3.0731887) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1378977) q[1];
sx q[1];
rz(-2.7337791) q[1];
sx q[1];
rz(0.56038709) q[1];
x q[2];
rz(-1.2980311) q[3];
sx q[3];
rz(-2.354151) q[3];
sx q[3];
rz(-1.8026601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8814964) q[2];
sx q[2];
rz(-1.0755971) q[2];
sx q[2];
rz(-1.0502226) q[2];
rz(-0.29551926) q[3];
sx q[3];
rz(-2.3727356) q[3];
sx q[3];
rz(-1.3692921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722026) q[0];
sx q[0];
rz(-2.1222332) q[0];
sx q[0];
rz(-0.66725677) q[0];
rz(2.9908906) q[1];
sx q[1];
rz(-2.0867911) q[1];
sx q[1];
rz(-0.31164718) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7512322) q[0];
sx q[0];
rz(-1.5003107) q[0];
sx q[0];
rz(-2.0930392) q[0];
x q[1];
rz(-1.1318593) q[2];
sx q[2];
rz(-1.8510185) q[2];
sx q[2];
rz(2.9228022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63285349) q[1];
sx q[1];
rz(-1.753974) q[1];
sx q[1];
rz(-1.4289209) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17201356) q[3];
sx q[3];
rz(-0.40296754) q[3];
sx q[3];
rz(-0.95517077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47505891) q[2];
sx q[2];
rz(-2.2481613) q[2];
sx q[2];
rz(0.7676355) q[2];
rz(2.660699) q[3];
sx q[3];
rz(-2.3700263) q[3];
sx q[3];
rz(-1.5546999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84592205) q[0];
sx q[0];
rz(-1.6039811) q[0];
sx q[0];
rz(2.3620102) q[0];
rz(-2.7698611) q[1];
sx q[1];
rz(-2.1340243) q[1];
sx q[1];
rz(0.76098162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8801988) q[0];
sx q[0];
rz(-2.7028599) q[0];
sx q[0];
rz(-0.83476557) q[0];
x q[1];
rz(2.2566608) q[2];
sx q[2];
rz(-1.2239309) q[2];
sx q[2];
rz(-0.82277966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5058712) q[1];
sx q[1];
rz(-2.2835915) q[1];
sx q[1];
rz(2.9221623) q[1];
rz(-pi) q[2];
rz(-2.2888921) q[3];
sx q[3];
rz(-1.5670683) q[3];
sx q[3];
rz(2.0223597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80428213) q[2];
sx q[2];
rz(-2.2411942) q[2];
sx q[2];
rz(-2.6864181) q[2];
rz(2.4496487) q[3];
sx q[3];
rz(-2.4271836) q[3];
sx q[3];
rz(2.3327995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1728425) q[0];
sx q[0];
rz(-1.3469561) q[0];
sx q[0];
rz(-2.8853048) q[0];
rz(-1.434727) q[1];
sx q[1];
rz(-1.7211569) q[1];
sx q[1];
rz(0.11875471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6390566) q[0];
sx q[0];
rz(-2.0169741) q[0];
sx q[0];
rz(1.7429551) q[0];
rz(-1.9740749) q[2];
sx q[2];
rz(-2.8244947) q[2];
sx q[2];
rz(2.4479579) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.465047) q[1];
sx q[1];
rz(-0.80261723) q[1];
sx q[1];
rz(1.8957183) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3975735) q[3];
sx q[3];
rz(-2.3459917) q[3];
sx q[3];
rz(-0.52594409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28748301) q[2];
sx q[2];
rz(-2.867925) q[2];
sx q[2];
rz(-2.0859094) q[2];
rz(-1.6915197) q[3];
sx q[3];
rz(-1.4472716) q[3];
sx q[3];
rz(-2.292574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5797822) q[0];
sx q[0];
rz(-0.17794839) q[0];
sx q[0];
rz(-1.5280888) q[0];
rz(-1.9644507) q[1];
sx q[1];
rz(-1.2700932) q[1];
sx q[1];
rz(-1.5783763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7684473) q[0];
sx q[0];
rz(-2.5119436) q[0];
sx q[0];
rz(-0.11778732) q[0];
x q[1];
rz(1.8991532) q[2];
sx q[2];
rz(-1.17982) q[2];
sx q[2];
rz(-2.7308488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.745143) q[1];
sx q[1];
rz(-2.1985594) q[1];
sx q[1];
rz(1.3954074) q[1];
rz(-pi) q[2];
rz(-1.568119) q[3];
sx q[3];
rz(-1.9802046) q[3];
sx q[3];
rz(-0.97312991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13581181) q[2];
sx q[2];
rz(-1.9848738) q[2];
sx q[2];
rz(1.8394252) q[2];
rz(0.33995134) q[3];
sx q[3];
rz(-2.3348742) q[3];
sx q[3];
rz(2.4792041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.133701) q[0];
sx q[0];
rz(-1.5971203) q[0];
sx q[0];
rz(1.1997461) q[0];
rz(1.3613191) q[1];
sx q[1];
rz(-0.97313762) q[1];
sx q[1];
rz(0.57788411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91980714) q[0];
sx q[0];
rz(-3.1263906) q[0];
sx q[0];
rz(2.0073246) q[0];
rz(-pi) q[1];
rz(0.33249929) q[2];
sx q[2];
rz(-2.6869171) q[2];
sx q[2];
rz(-0.60078675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1513704) q[1];
sx q[1];
rz(-1.7066907) q[1];
sx q[1];
rz(-2.6841867) q[1];
x q[2];
rz(1.6541566) q[3];
sx q[3];
rz(-2.8417086) q[3];
sx q[3];
rz(0.56902992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.63152385) q[2];
sx q[2];
rz(-0.57454595) q[2];
sx q[2];
rz(-0.4734545) q[2];
rz(-1.8691285) q[3];
sx q[3];
rz(-1.3926287) q[3];
sx q[3];
rz(-0.22245358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1943787) q[0];
sx q[0];
rz(-2.4838303) q[0];
sx q[0];
rz(-2.8084602) q[0];
rz(-2.9130452) q[1];
sx q[1];
rz(-1.8014149) q[1];
sx q[1];
rz(0.57565912) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822214) q[0];
sx q[0];
rz(-2.3239808) q[0];
sx q[0];
rz(-0.14540093) q[0];
rz(-pi) q[1];
rz(-2.6645053) q[2];
sx q[2];
rz(-1.2735575) q[2];
sx q[2];
rz(1.5971668) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.043283431) q[1];
sx q[1];
rz(-0.38615882) q[1];
sx q[1];
rz(-2.1257867) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6603465) q[3];
sx q[3];
rz(-1.3433045) q[3];
sx q[3];
rz(2.9862822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.26611844) q[2];
sx q[2];
rz(-2.5901651) q[2];
sx q[2];
rz(-1.4761338) q[2];
rz(1.1009781) q[3];
sx q[3];
rz(-1.063238) q[3];
sx q[3];
rz(-2.6235918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5440893) q[0];
sx q[0];
rz(-0.88633716) q[0];
sx q[0];
rz(-0.30712095) q[0];
rz(-0.33044526) q[1];
sx q[1];
rz(-1.5978866) q[1];
sx q[1];
rz(-1.7764567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70676409) q[0];
sx q[0];
rz(-1.4894823) q[0];
sx q[0];
rz(1.5887194) q[0];
x q[1];
rz(-0.34806378) q[2];
sx q[2];
rz(-1.3937147) q[2];
sx q[2];
rz(-0.40906104) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.04330895) q[1];
sx q[1];
rz(-1.829147) q[1];
sx q[1];
rz(-1.1275379) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7264113) q[3];
sx q[3];
rz(-0.74374108) q[3];
sx q[3];
rz(-1.2228325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8660628) q[2];
sx q[2];
rz(-2.5538462) q[2];
sx q[2];
rz(-2.8928939) q[2];
rz(-1.9281049) q[3];
sx q[3];
rz(-1.4444193) q[3];
sx q[3];
rz(0.061633751) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1520749) q[0];
sx q[0];
rz(-2.2353421) q[0];
sx q[0];
rz(-2.4549947) q[0];
rz(1.1129414) q[1];
sx q[1];
rz(-2.3390892) q[1];
sx q[1];
rz(-0.31563219) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9292727) q[0];
sx q[0];
rz(-1.2480191) q[0];
sx q[0];
rz(-0.048521532) q[0];
rz(-2.8807345) q[2];
sx q[2];
rz(-2.3645698) q[2];
sx q[2];
rz(2.379247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0084969) q[1];
sx q[1];
rz(-1.6649455) q[1];
sx q[1];
rz(-0.90742438) q[1];
x q[2];
rz(-0.078048869) q[3];
sx q[3];
rz(-2.4025318) q[3];
sx q[3];
rz(-2.259425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.687872) q[2];
sx q[2];
rz(-2.5747955) q[2];
sx q[2];
rz(2.4570214) q[2];
rz(-1.941393) q[3];
sx q[3];
rz(-2.0122416) q[3];
sx q[3];
rz(2.9620192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0803364) q[0];
sx q[0];
rz(-2.6604524) q[0];
sx q[0];
rz(3.1367593) q[0];
rz(-1.5176516) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(2.8878816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4986202) q[0];
sx q[0];
rz(-1.7812087) q[0];
sx q[0];
rz(3.0543113) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0982622) q[2];
sx q[2];
rz(-0.8923549) q[2];
sx q[2];
rz(2.0022165) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1349134) q[1];
sx q[1];
rz(-2.2354534) q[1];
sx q[1];
rz(-2.8921739) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7697479) q[3];
sx q[3];
rz(-2.1124438) q[3];
sx q[3];
rz(-2.1218079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1634875) q[2];
sx q[2];
rz(-1.9065964) q[2];
sx q[2];
rz(1.9592436) q[2];
rz(-2.4649418) q[3];
sx q[3];
rz(-0.38656056) q[3];
sx q[3];
rz(2.1277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61426281) q[0];
sx q[0];
rz(-2.3557721) q[0];
sx q[0];
rz(-2.8489805) q[0];
rz(2.09265) q[1];
sx q[1];
rz(-1.7221778) q[1];
sx q[1];
rz(-2.4377951) q[1];
rz(0.60381966) q[2];
sx q[2];
rz(-0.55470822) q[2];
sx q[2];
rz(1.8798238) q[2];
rz(-0.098930704) q[3];
sx q[3];
rz(-1.9361648) q[3];
sx q[3];
rz(1.786644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
