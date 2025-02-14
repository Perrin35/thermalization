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
rz(0.69819063) q[0];
sx q[0];
rz(-0.31261045) q[0];
sx q[0];
rz(-0.99553776) q[0];
rz(1.1065464) q[1];
sx q[1];
rz(-0.76835978) q[1];
sx q[1];
rz(-0.33222693) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8790008) q[0];
sx q[0];
rz(-0.79918062) q[0];
sx q[0];
rz(-2.5627717) q[0];
x q[1];
rz(-2.034445) q[2];
sx q[2];
rz(-1.6939179) q[2];
sx q[2];
rz(2.5277918) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39980934) q[1];
sx q[1];
rz(-2.4917549) q[1];
sx q[1];
rz(0.8292966) q[1];
rz(-pi) q[2];
rz(2.1463582) q[3];
sx q[3];
rz(-1.1249295) q[3];
sx q[3];
rz(0.80328548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.60653162) q[2];
sx q[2];
rz(-2.1273095) q[2];
sx q[2];
rz(-0.052058546) q[2];
rz(-2.0590797) q[3];
sx q[3];
rz(-2.1942997) q[3];
sx q[3];
rz(-1.989971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.624619) q[0];
sx q[0];
rz(-3.1095412) q[0];
sx q[0];
rz(0.33729851) q[0];
rz(-0.15268046) q[1];
sx q[1];
rz(-2.3744507) q[1];
sx q[1];
rz(-1.618636) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82459506) q[0];
sx q[0];
rz(-0.64205161) q[0];
sx q[0];
rz(-2.8606728) q[0];
rz(-pi) q[1];
rz(2.8711393) q[2];
sx q[2];
rz(-2.4861397) q[2];
sx q[2];
rz(-1.7846817) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3831269) q[1];
sx q[1];
rz(-0.43769203) q[1];
sx q[1];
rz(1.0785224) q[1];
x q[2];
rz(-1.7256152) q[3];
sx q[3];
rz(-0.972675) q[3];
sx q[3];
rz(-3.104676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0308257) q[2];
sx q[2];
rz(-0.36588565) q[2];
sx q[2];
rz(-2.7552674) q[2];
rz(-1.857916) q[3];
sx q[3];
rz(-1.176703) q[3];
sx q[3];
rz(1.9234689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.074742643) q[0];
sx q[0];
rz(-2.061494) q[0];
sx q[0];
rz(2.1225488) q[0];
rz(-2.3661803) q[1];
sx q[1];
rz(-1.2762028) q[1];
sx q[1];
rz(-2.6554241) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2371444) q[0];
sx q[0];
rz(-1.5394566) q[0];
sx q[0];
rz(-1.5737359) q[0];
rz(-pi) q[1];
rz(-1.08467) q[2];
sx q[2];
rz(-2.297431) q[2];
sx q[2];
rz(-1.655533) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8067161) q[1];
sx q[1];
rz(-1.7893605) q[1];
sx q[1];
rz(0.49592917) q[1];
rz(-pi) q[2];
rz(2.4047732) q[3];
sx q[3];
rz(-0.42356682) q[3];
sx q[3];
rz(-1.119348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7458618) q[2];
sx q[2];
rz(-1.9755325) q[2];
sx q[2];
rz(2.3365848) q[2];
rz(0.65195596) q[3];
sx q[3];
rz(-1.1019573) q[3];
sx q[3];
rz(0.93418795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4790799) q[0];
sx q[0];
rz(-1.3211687) q[0];
sx q[0];
rz(-1.3385734) q[0];
rz(-1.5486807) q[1];
sx q[1];
rz(-2.4376696) q[1];
sx q[1];
rz(-2.7159363) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5558472) q[0];
sx q[0];
rz(-0.092215538) q[0];
sx q[0];
rz(-3.0526428) q[0];
rz(-1.6792308) q[2];
sx q[2];
rz(-0.30955704) q[2];
sx q[2];
rz(-2.5029687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7158) q[1];
sx q[1];
rz(-1.1669901) q[1];
sx q[1];
rz(2.3753662) q[1];
rz(-pi) q[2];
rz(-2.8208977) q[3];
sx q[3];
rz(-2.3983722) q[3];
sx q[3];
rz(2.7666465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8694596) q[2];
sx q[2];
rz(-1.3674066) q[2];
sx q[2];
rz(-0.41932219) q[2];
rz(-1.9648633) q[3];
sx q[3];
rz(-2.4265225) q[3];
sx q[3];
rz(2.5756605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45593539) q[0];
sx q[0];
rz(-2.3929907) q[0];
sx q[0];
rz(-1.5022044) q[0];
rz(1.4337076) q[1];
sx q[1];
rz(-1.2805484) q[1];
sx q[1];
rz(0.40275231) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086632816) q[0];
sx q[0];
rz(-1.827359) q[0];
sx q[0];
rz(-2.199928) q[0];
rz(-pi) q[1];
rz(1.5115159) q[2];
sx q[2];
rz(-1.0067954) q[2];
sx q[2];
rz(-0.21280542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6357543) q[1];
sx q[1];
rz(-0.16316667) q[1];
sx q[1];
rz(-0.6357622) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66688315) q[3];
sx q[3];
rz(-1.6630238) q[3];
sx q[3];
rz(0.86574829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5002284) q[2];
sx q[2];
rz(-2.0945175) q[2];
sx q[2];
rz(2.8889528) q[2];
rz(0.80444515) q[3];
sx q[3];
rz(-0.52144709) q[3];
sx q[3];
rz(-0.49853244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.8135524) q[0];
sx q[0];
rz(-0.10843065) q[0];
sx q[0];
rz(2.257708) q[0];
rz(2.6416595) q[1];
sx q[1];
rz(-1.6729313) q[1];
sx q[1];
rz(-0.51441851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5926625) q[0];
sx q[0];
rz(-0.99931006) q[0];
sx q[0];
rz(0.40056374) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53192682) q[2];
sx q[2];
rz(-1.3585603) q[2];
sx q[2];
rz(-1.797582) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41082327) q[1];
sx q[1];
rz(-1.1699153) q[1];
sx q[1];
rz(-2.904179) q[1];
rz(3.1153999) q[3];
sx q[3];
rz(-2.3520654) q[3];
sx q[3];
rz(1.4007555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9875235) q[2];
sx q[2];
rz(-1.46393) q[2];
sx q[2];
rz(-0.12518159) q[2];
rz(2.3212738) q[3];
sx q[3];
rz(-1.9671974) q[3];
sx q[3];
rz(-2.7952747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5527375) q[0];
sx q[0];
rz(-0.6468361) q[0];
sx q[0];
rz(0.11189017) q[0];
rz(-2.0018068) q[1];
sx q[1];
rz(-0.21509376) q[1];
sx q[1];
rz(1.8591759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70045602) q[0];
sx q[0];
rz(-0.64986357) q[0];
sx q[0];
rz(2.031424) q[0];
x q[1];
rz(0.23156802) q[2];
sx q[2];
rz(-1.7121669) q[2];
sx q[2];
rz(-1.8166775) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3610745) q[1];
sx q[1];
rz(-0.86336923) q[1];
sx q[1];
rz(2.646628) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68422785) q[3];
sx q[3];
rz(-1.6155532) q[3];
sx q[3];
rz(1.0563204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.68975) q[2];
sx q[2];
rz(-0.30561438) q[2];
sx q[2];
rz(-1.8243194) q[2];
rz(-3.1210323) q[3];
sx q[3];
rz(-0.33187425) q[3];
sx q[3];
rz(2.1338972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2059712) q[0];
sx q[0];
rz(-2.1501849) q[0];
sx q[0];
rz(0.29888612) q[0];
rz(-2.0805953) q[1];
sx q[1];
rz(-1.7540878) q[1];
sx q[1];
rz(-1.908173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80432928) q[0];
sx q[0];
rz(-1.6676039) q[0];
sx q[0];
rz(1.2362572) q[0];
rz(-pi) q[1];
rz(-1.226108) q[2];
sx q[2];
rz(-0.12154254) q[2];
sx q[2];
rz(-2.8555088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2005991) q[1];
sx q[1];
rz(-0.72807136) q[1];
sx q[1];
rz(1.1173825) q[1];
rz(-0.66090195) q[3];
sx q[3];
rz(-1.1422302) q[3];
sx q[3];
rz(2.1235583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79248205) q[2];
sx q[2];
rz(-1.0976378) q[2];
sx q[2];
rz(-0.20723542) q[2];
rz(2.4810897) q[3];
sx q[3];
rz(-2.1410172) q[3];
sx q[3];
rz(1.2787261) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5527363) q[0];
sx q[0];
rz(-1.739946) q[0];
sx q[0];
rz(1.2218342) q[0];
rz(0.51512042) q[1];
sx q[1];
rz(-0.11038596) q[1];
sx q[1];
rz(2.5882904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.57896) q[0];
sx q[0];
rz(-1.4427358) q[0];
sx q[0];
rz(2.1306031) q[0];
rz(-1.1498575) q[2];
sx q[2];
rz(-1.5113665) q[2];
sx q[2];
rz(0.42039117) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33722222) q[1];
sx q[1];
rz(-2.9644199) q[1];
sx q[1];
rz(1.6044264) q[1];
x q[2];
rz(-0.030144088) q[3];
sx q[3];
rz(-0.6693535) q[3];
sx q[3];
rz(0.02756674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33636609) q[2];
sx q[2];
rz(-2.0311425) q[2];
sx q[2];
rz(2.6160348) q[2];
rz(1.6622274) q[3];
sx q[3];
rz(-0.95366228) q[3];
sx q[3];
rz(0.76013887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084759921) q[0];
sx q[0];
rz(-2.3469717) q[0];
sx q[0];
rz(-0.42924616) q[0];
rz(1.6191354) q[1];
sx q[1];
rz(-1.2702962) q[1];
sx q[1];
rz(2.2116908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7741755) q[0];
sx q[0];
rz(-0.40433592) q[0];
sx q[0];
rz(-1.0003759) q[0];
x q[1];
rz(-1.11082) q[2];
sx q[2];
rz(-0.24402741) q[2];
sx q[2];
rz(-2.8619253) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9318951) q[1];
sx q[1];
rz(-1.731809) q[1];
sx q[1];
rz(1.2847394) q[1];
x q[2];
rz(2.5885196) q[3];
sx q[3];
rz(-1.6216842) q[3];
sx q[3];
rz(1.0158051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.736019) q[2];
sx q[2];
rz(-0.47601998) q[2];
sx q[2];
rz(-0.19301566) q[2];
rz(-0.080282601) q[3];
sx q[3];
rz(-0.86997) q[3];
sx q[3];
rz(-1.75753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3748462) q[0];
sx q[0];
rz(-1.9369047) q[0];
sx q[0];
rz(1.9932224) q[0];
rz(2.171352) q[1];
sx q[1];
rz(-1.2087676) q[1];
sx q[1];
rz(-1.2420775) q[1];
rz(2.8497981) q[2];
sx q[2];
rz(-0.84635432) q[2];
sx q[2];
rz(1.9683471) q[2];
rz(-2.8945782) q[3];
sx q[3];
rz(-0.09709662) q[3];
sx q[3];
rz(0.63963565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
