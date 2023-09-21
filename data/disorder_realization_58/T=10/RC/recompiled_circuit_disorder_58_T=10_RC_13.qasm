OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(-2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18618628) q[0];
sx q[0];
rz(-0.49476981) q[0];
sx q[0];
rz(-0.68899378) q[0];
x q[1];
rz(2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(-2.8505461) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4315011) q[1];
sx q[1];
rz(-2.3570865) q[1];
sx q[1];
rz(-2.0246519) q[1];
rz(-pi) q[2];
rz(0.72676267) q[3];
sx q[3];
rz(-0.57075497) q[3];
sx q[3];
rz(-1.4445514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(2.2297915) q[0];
rz(0.78951019) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2355736) q[0];
sx q[0];
rz(-2.2740907) q[0];
sx q[0];
rz(2.2554382) q[0];
x q[1];
rz(1.1268483) q[2];
sx q[2];
rz(-1.7146646) q[2];
sx q[2];
rz(0.680188) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5527199) q[1];
sx q[1];
rz(-0.96541222) q[1];
sx q[1];
rz(2.6436716) q[1];
x q[2];
rz(-0.88300206) q[3];
sx q[3];
rz(-1.6371173) q[3];
sx q[3];
rz(1.6460713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(2.8116995) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(2.0842016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111423) q[0];
sx q[0];
rz(-1.550195) q[0];
sx q[0];
rz(1.9111454) q[0];
rz(-pi) q[1];
rz(-2.8003642) q[2];
sx q[2];
rz(-1.382302) q[2];
sx q[2];
rz(-2.598912) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2705546) q[1];
sx q[1];
rz(-1.9005617) q[1];
sx q[1];
rz(2.2689181) q[1];
rz(-1.7626761) q[3];
sx q[3];
rz(-1.1550316) q[3];
sx q[3];
rz(0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(0.81956285) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574922) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9388158) q[0];
sx q[0];
rz(-1.446615) q[0];
sx q[0];
rz(1.7569655) q[0];
x q[1];
rz(1.7521162) q[2];
sx q[2];
rz(-1.4099979) q[2];
sx q[2];
rz(2.5374075) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0927825) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(-2.584143) q[1];
rz(-pi) q[2];
rz(1.1449279) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(-1.4262517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(-1.127634) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(-1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(0.23194557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1629099) q[0];
sx q[0];
rz(-1.1497578) q[0];
sx q[0];
rz(-2.1091503) q[0];
x q[1];
rz(2.4460692) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(0.34851375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1272808) q[1];
sx q[1];
rz(-2.8297272) q[1];
sx q[1];
rz(0.77906268) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0911345) q[3];
sx q[3];
rz(-2.7726463) q[3];
sx q[3];
rz(1.7794533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(-2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0971138) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(0.99745497) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4076685) q[0];
sx q[0];
rz(-1.7051538) q[0];
sx q[0];
rz(2.2785447) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41001292) q[2];
sx q[2];
rz(-0.7747246) q[2];
sx q[2];
rz(-2.0395525) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84038489) q[1];
sx q[1];
rz(-1.1494698) q[1];
sx q[1];
rz(-2.8788484) q[1];
rz(-pi) q[2];
rz(0.84647471) q[3];
sx q[3];
rz(-1.9123565) q[3];
sx q[3];
rz(-0.70111707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(-0.4513936) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(0.61378941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44871556) q[0];
sx q[0];
rz(-2.2971417) q[0];
sx q[0];
rz(-2.8918299) q[0];
x q[1];
rz(2.8579312) q[2];
sx q[2];
rz(-1.6322871) q[2];
sx q[2];
rz(0.21360699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9871414) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(-2.6277072) q[1];
rz(3.0244163) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(-1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(2.1571295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280633) q[0];
sx q[0];
rz(-0.54740471) q[0];
sx q[0];
rz(-3.0696966) q[0];
rz(-0.41495277) q[2];
sx q[2];
rz(-2.3256362) q[2];
sx q[2];
rz(2.6331537) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8967646) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(2.5823309) q[3];
sx q[3];
rz(-2.5945633) q[3];
sx q[3];
rz(-0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(0.53517503) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(0.92591441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6322964) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(2.1584828) q[0];
rz(-0.70897734) q[2];
sx q[2];
rz(-0.93009863) q[2];
sx q[2];
rz(2.7880653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33854252) q[1];
sx q[1];
rz(-2.0085137) q[1];
sx q[1];
rz(-2.9796757) q[1];
rz(-pi) q[2];
rz(1.6894475) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(0.53393902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(-1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-2.1283456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37975882) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(-0.14302111) q[0];
rz(-pi) q[1];
rz(-2.7878694) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(1.8281787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80876795) q[1];
sx q[1];
rz(-1.1507532) q[1];
sx q[1];
rz(1.4648449) q[1];
rz(-2.2447484) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(0.53900063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-2.3821793) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(-1.4105994) q[2];
sx q[2];
rz(-0.44993958) q[2];
sx q[2];
rz(-1.1100563) q[2];
rz(-0.063354062) q[3];
sx q[3];
rz(-2.2676716) q[3];
sx q[3];
rz(0.36361658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];