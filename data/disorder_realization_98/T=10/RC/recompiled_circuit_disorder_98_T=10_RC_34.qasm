OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(-0.87061849) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4275442) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(0.9057522) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6354191) q[2];
sx q[2];
rz(-1.0867599) q[2];
sx q[2];
rz(1.5278221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.183061) q[1];
sx q[1];
rz(-1.5325938) q[1];
sx q[1];
rz(-1.3903635) q[1];
x q[2];
rz(0.19102328) q[3];
sx q[3];
rz(-2.1130307) q[3];
sx q[3];
rz(0.99381002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(2.4893563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61440496) q[0];
sx q[0];
rz(-1.5584649) q[0];
sx q[0];
rz(0.028610882) q[0];
x q[1];
rz(-1.2781497) q[2];
sx q[2];
rz(-1.3424982) q[2];
sx q[2];
rz(-1.6739068) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7420885) q[1];
sx q[1];
rz(-1.794145) q[1];
sx q[1];
rz(0.29466596) q[1];
rz(-pi) q[2];
rz(-2.2964301) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(1.0172539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(-1.9799505) q[2];
rz(1.15796) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(-2.2913349) q[0];
rz(2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(-1.3495548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36950668) q[0];
sx q[0];
rz(-1.8096576) q[0];
sx q[0];
rz(-0.68349616) q[0];
x q[1];
rz(0.037319855) q[2];
sx q[2];
rz(-1.63675) q[2];
sx q[2];
rz(-1.0730336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68223665) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(0.45046803) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1317741) q[3];
sx q[3];
rz(-1.1447284) q[3];
sx q[3];
rz(0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(-0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(3.1052123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2934389) q[0];
sx q[0];
rz(-1.358195) q[0];
sx q[0];
rz(2.2178177) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5600299) q[2];
sx q[2];
rz(-1.8459324) q[2];
sx q[2];
rz(-0.046422596) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42602793) q[1];
sx q[1];
rz(-1.702311) q[1];
sx q[1];
rz(-2.2711666) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7500238) q[3];
sx q[3];
rz(-2.8885926) q[3];
sx q[3];
rz(-0.60245017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-2.9072445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7097965) q[0];
sx q[0];
rz(-1.0456107) q[0];
sx q[0];
rz(-1.1774506) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0767897) q[2];
sx q[2];
rz(-1.6550118) q[2];
sx q[2];
rz(-0.46352026) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8824132) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(-1.773136) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8946394) q[3];
sx q[3];
rz(-1.4959335) q[3];
sx q[3];
rz(2.4002241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(0.22053545) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41546145) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(0.56610402) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9691539) q[0];
sx q[0];
rz(-1.5307431) q[0];
sx q[0];
rz(-2.7705631) q[0];
rz(-pi) q[1];
rz(0.67310682) q[2];
sx q[2];
rz(-1.7761201) q[2];
sx q[2];
rz(-0.38773195) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5551344) q[1];
sx q[1];
rz(-2.0124334) q[1];
sx q[1];
rz(0.099979062) q[1];
x q[2];
rz(2.6236344) q[3];
sx q[3];
rz(-1.8149788) q[3];
sx q[3];
rz(1.1740008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(0.4894408) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.8547159) q[0];
rz(0.6634179) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(1.9082665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7466) q[0];
sx q[0];
rz(-1.9548423) q[0];
sx q[0];
rz(-3.1296484) q[0];
rz(-pi) q[1];
rz(0.6255409) q[2];
sx q[2];
rz(-0.95226804) q[2];
sx q[2];
rz(2.2369838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.308097) q[1];
sx q[1];
rz(-1.4596246) q[1];
sx q[1];
rz(-1.584946) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4864278) q[3];
sx q[3];
rz(-1.2573811) q[3];
sx q[3];
rz(0.13669554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(-1.0160149) q[2];
rz(-0.070090381) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163517) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(-0.21487543) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.8833556) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90342605) q[0];
sx q[0];
rz(-1.705372) q[0];
sx q[0];
rz(1.1358791) q[0];
rz(-0.9332946) q[2];
sx q[2];
rz(-1.0273232) q[2];
sx q[2];
rz(-0.18804929) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82842365) q[1];
sx q[1];
rz(-1.0272044) q[1];
sx q[1];
rz(1.6924752) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97790896) q[3];
sx q[3];
rz(-0.80855723) q[3];
sx q[3];
rz(0.18329328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(0.56274596) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(0.02773157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8160307) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-0.82657878) q[0];
x q[1];
rz(2.8743923) q[2];
sx q[2];
rz(-2.2520817) q[2];
sx q[2];
rz(2.6097678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8393644) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(-0.60217963) q[1];
rz(3.1157007) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(-2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(-0.85723248) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(0.25751105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637708) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(-1.9985984) q[0];
rz(0.54458877) q[2];
sx q[2];
rz(-2.0047744) q[2];
sx q[2];
rz(-1.4101654) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7911437) q[1];
sx q[1];
rz(-0.87462438) q[1];
sx q[1];
rz(-1.5894366) q[1];
x q[2];
rz(-2.8994843) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(0.47484067) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941147) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(2.9150302) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(-2.6735641) q[2];
sx q[2];
rz(-1.2266908) q[2];
sx q[2];
rz(-2.5897349) q[2];
rz(-2.0541035) q[3];
sx q[3];
rz(-1.8046422) q[3];
sx q[3];
rz(1.3840152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
