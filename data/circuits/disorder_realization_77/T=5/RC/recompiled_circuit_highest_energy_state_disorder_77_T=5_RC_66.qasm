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
rz(2.0692985) q[0];
sx q[0];
rz(-2.9579853) q[0];
sx q[0];
rz(-2.973383) q[0];
rz(-1.3104562) q[1];
sx q[1];
rz(1.1163196) q[1];
sx q[1];
rz(7.4012227) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2471825) q[0];
sx q[0];
rz(-1.6798816) q[0];
sx q[0];
rz(-2.8139958) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9498212) q[2];
sx q[2];
rz(-2.2281335) q[2];
sx q[2];
rz(2.71703) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1217482) q[1];
sx q[1];
rz(-2.3308874) q[1];
sx q[1];
rz(-1.5365063) q[1];
rz(-pi) q[2];
rz(-3.0039677) q[3];
sx q[3];
rz(-1.3018908) q[3];
sx q[3];
rz(2.5816504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2844124) q[2];
sx q[2];
rz(-0.99049157) q[2];
sx q[2];
rz(-1.8641137) q[2];
rz(-0.60753456) q[3];
sx q[3];
rz(-1.7888864) q[3];
sx q[3];
rz(1.7551306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69433895) q[0];
sx q[0];
rz(-2.4632833) q[0];
sx q[0];
rz(0.69764486) q[0];
rz(0.33277008) q[1];
sx q[1];
rz(-1.1079301) q[1];
sx q[1];
rz(0.3322126) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8150114) q[0];
sx q[0];
rz(-1.4403116) q[0];
sx q[0];
rz(-1.4490118) q[0];
x q[1];
rz(0.9941446) q[2];
sx q[2];
rz(-1.3704164) q[2];
sx q[2];
rz(2.8746614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7113641) q[1];
sx q[1];
rz(-1.1540742) q[1];
sx q[1];
rz(2.8695089) q[1];
rz(2.4139514) q[3];
sx q[3];
rz(-1.8944358) q[3];
sx q[3];
rz(-1.0689842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.701) q[2];
sx q[2];
rz(-1.168246) q[2];
sx q[2];
rz(-2.8458703) q[2];
rz(3.131007) q[3];
sx q[3];
rz(-0.9897832) q[3];
sx q[3];
rz(-2.4081965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6538786) q[0];
sx q[0];
rz(-0.45844498) q[0];
sx q[0];
rz(2.606875) q[0];
rz(-1.1010822) q[1];
sx q[1];
rz(-1.4357166) q[1];
sx q[1];
rz(-0.91317493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3951626) q[0];
sx q[0];
rz(-2.5341153) q[0];
sx q[0];
rz(2.1588401) q[0];
rz(-1.1712672) q[2];
sx q[2];
rz(-0.59430423) q[2];
sx q[2];
rz(2.8565796) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7353834) q[1];
sx q[1];
rz(-2.2158631) q[1];
sx q[1];
rz(2.2570133) q[1];
rz(-pi) q[2];
rz(-1.6524773) q[3];
sx q[3];
rz(-0.73058999) q[3];
sx q[3];
rz(-0.047103492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59213006) q[2];
sx q[2];
rz(-2.2495146) q[2];
sx q[2];
rz(-2.7590052) q[2];
rz(1.654918) q[3];
sx q[3];
rz(-0.29522172) q[3];
sx q[3];
rz(2.2216643) q[3];
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
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53052467) q[0];
sx q[0];
rz(-1.684573) q[0];
sx q[0];
rz(2.0953505) q[0];
rz(0.56398448) q[1];
sx q[1];
rz(-1.0020741) q[1];
sx q[1];
rz(-1.8246338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038669) q[0];
sx q[0];
rz(-2.236294) q[0];
sx q[0];
rz(-0.72030495) q[0];
x q[1];
rz(0.91805075) q[2];
sx q[2];
rz(-1.3155283) q[2];
sx q[2];
rz(1.7328615) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7822588) q[1];
sx q[1];
rz(-0.31681199) q[1];
sx q[1];
rz(-0.13765361) q[1];
rz(1.9924367) q[3];
sx q[3];
rz(-2.3083335) q[3];
sx q[3];
rz(1.6374388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4468533) q[2];
sx q[2];
rz(-1.4444192) q[2];
sx q[2];
rz(0.99855885) q[2];
rz(0.57991943) q[3];
sx q[3];
rz(-1.6603575) q[3];
sx q[3];
rz(1.8370321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9272083) q[0];
sx q[0];
rz(-1.913232) q[0];
sx q[0];
rz(-2.4217915) q[0];
rz(0.7431227) q[1];
sx q[1];
rz(-1.9235976) q[1];
sx q[1];
rz(-0.057083759) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9225153) q[0];
sx q[0];
rz(-1.2033899) q[0];
sx q[0];
rz(1.4266582) q[0];
x q[1];
rz(-2.7774335) q[2];
sx q[2];
rz(-2.3500867) q[2];
sx q[2];
rz(3.0730545) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6634906) q[1];
sx q[1];
rz(-1.2896207) q[1];
sx q[1];
rz(1.737514) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43566849) q[3];
sx q[3];
rz(-1.0188802) q[3];
sx q[3];
rz(-2.4727522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7032949) q[2];
sx q[2];
rz(-1.1763828) q[2];
sx q[2];
rz(1.8717742) q[2];
rz(0.54369175) q[3];
sx q[3];
rz(-2.4589977) q[3];
sx q[3];
rz(1.3683176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0755587) q[0];
sx q[0];
rz(-0.37519255) q[0];
sx q[0];
rz(-0.22848465) q[0];
rz(0.98555073) q[1];
sx q[1];
rz(-1.5864213) q[1];
sx q[1];
rz(1.2062581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.120283) q[0];
sx q[0];
rz(-0.37719676) q[0];
sx q[0];
rz(-0.46964236) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94088631) q[2];
sx q[2];
rz(-0.36782757) q[2];
sx q[2];
rz(-1.6339677) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9282189) q[1];
sx q[1];
rz(-2.5210533) q[1];
sx q[1];
rz(-1.7804407) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5980232) q[3];
sx q[3];
rz(-1.479621) q[3];
sx q[3];
rz(-1.4379355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0050521684) q[2];
sx q[2];
rz(-1.2753992) q[2];
sx q[2];
rz(1.8955815) q[2];
rz(-0.15240845) q[3];
sx q[3];
rz(-0.14901769) q[3];
sx q[3];
rz(-0.77308956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8361255) q[0];
sx q[0];
rz(-1.2141328) q[0];
sx q[0];
rz(0.25492302) q[0];
rz(0.98006025) q[1];
sx q[1];
rz(-1.3420811) q[1];
sx q[1];
rz(-2.9068388) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7388044) q[0];
sx q[0];
rz(-1.0287893) q[0];
sx q[0];
rz(0.2861342) q[0];
x q[1];
rz(-0.94257109) q[2];
sx q[2];
rz(-2.8543575) q[2];
sx q[2];
rz(-0.38610215) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3153441) q[1];
sx q[1];
rz(-0.49737206) q[1];
sx q[1];
rz(-2.5647463) q[1];
x q[2];
rz(1.4244377) q[3];
sx q[3];
rz(-0.68181935) q[3];
sx q[3];
rz(0.54659971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.047082575) q[2];
sx q[2];
rz(-0.42028704) q[2];
sx q[2];
rz(-1.1959929) q[2];
rz(-2.1434873) q[3];
sx q[3];
rz(-1.866303) q[3];
sx q[3];
rz(-2.3534145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3274479) q[0];
sx q[0];
rz(-2.271551) q[0];
sx q[0];
rz(0.63234627) q[0];
rz(1.5517722) q[1];
sx q[1];
rz(-1.3464144) q[1];
sx q[1];
rz(-1.2987035) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1237268) q[0];
sx q[0];
rz(-1.3204823) q[0];
sx q[0];
rz(-1.1160442) q[0];
rz(-1.5517637) q[2];
sx q[2];
rz(-0.90523273) q[2];
sx q[2];
rz(-1.9766146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4711705) q[1];
sx q[1];
rz(-1.6123796) q[1];
sx q[1];
rz(-0.50703185) q[1];
x q[2];
rz(1.5428144) q[3];
sx q[3];
rz(-1.4194427) q[3];
sx q[3];
rz(1.3642814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87045264) q[2];
sx q[2];
rz(-0.91338921) q[2];
sx q[2];
rz(0.75922981) q[2];
rz(-2.7325654) q[3];
sx q[3];
rz(-1.4141915) q[3];
sx q[3];
rz(-3.0294688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.2893534) q[0];
sx q[0];
rz(-1.8426581) q[0];
sx q[0];
rz(-2.7759743) q[0];
rz(-3.1386555) q[1];
sx q[1];
rz(-1.6117088) q[1];
sx q[1];
rz(-2.073435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6474175) q[0];
sx q[0];
rz(-1.023698) q[0];
sx q[0];
rz(1.8965333) q[0];
rz(2.3084776) q[2];
sx q[2];
rz(-2.1016869) q[2];
sx q[2];
rz(-2.7754727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4959855) q[1];
sx q[1];
rz(-1.9580972) q[1];
sx q[1];
rz(3.1003968) q[1];
rz(0.13309388) q[3];
sx q[3];
rz(-0.20340098) q[3];
sx q[3];
rz(1.4730367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3785765) q[2];
sx q[2];
rz(-2.2282659) q[2];
sx q[2];
rz(2.1545289) q[2];
rz(-2.3297564) q[3];
sx q[3];
rz(-2.2631009) q[3];
sx q[3];
rz(-1.2442376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(3.0768452) q[0];
sx q[0];
rz(-0.61926121) q[0];
sx q[0];
rz(-1.4641807) q[0];
rz(0.96425104) q[1];
sx q[1];
rz(-1.827927) q[1];
sx q[1];
rz(-1.3826694) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.156517) q[0];
sx q[0];
rz(-1.1958953) q[0];
sx q[0];
rz(0.84701726) q[0];
x q[1];
rz(-1.8455454) q[2];
sx q[2];
rz(-1.2128069) q[2];
sx q[2];
rz(1.4520979) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5359491) q[1];
sx q[1];
rz(-1.1472196) q[1];
sx q[1];
rz(-0.98137318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85282214) q[3];
sx q[3];
rz(-2.3896165) q[3];
sx q[3];
rz(-2.9319921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9619523) q[2];
sx q[2];
rz(-1.8018758) q[2];
sx q[2];
rz(1.4947653) q[2];
rz(-1.1163813) q[3];
sx q[3];
rz(-0.79018441) q[3];
sx q[3];
rz(-0.60379973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8981001) q[0];
sx q[0];
rz(-1.4673163) q[0];
sx q[0];
rz(-1.2661288) q[0];
rz(2.948214) q[1];
sx q[1];
rz(-1.0693751) q[1];
sx q[1];
rz(1.3042915) q[1];
rz(0.73430772) q[2];
sx q[2];
rz(-1.3398017) q[2];
sx q[2];
rz(-0.2504713) q[2];
rz(0.51283097) q[3];
sx q[3];
rz(-0.75779946) q[3];
sx q[3];
rz(2.1661314) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
