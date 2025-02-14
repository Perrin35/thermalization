OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0754492) q[0];
sx q[0];
rz(-1.0439405) q[0];
sx q[0];
rz(-3.1312842) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(4.5865321) q[1];
sx q[1];
rz(10.001339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2808997) q[0];
sx q[0];
rz(-1.6917949) q[0];
sx q[0];
rz(-0.31009407) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85421487) q[2];
sx q[2];
rz(-0.81939261) q[2];
sx q[2];
rz(-1.508282) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.623201) q[1];
sx q[1];
rz(-1.745959) q[1];
sx q[1];
rz(-0.53304146) q[1];
rz(0.83769862) q[3];
sx q[3];
rz(-1.4464966) q[3];
sx q[3];
rz(-1.4364786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52446857) q[2];
sx q[2];
rz(-1.4602129) q[2];
sx q[2];
rz(1.8560393) q[2];
rz(-1.5517976) q[3];
sx q[3];
rz(-2.0701305) q[3];
sx q[3];
rz(1.0514528) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0153506) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(-2.8958564) q[0];
rz(1.0579146) q[1];
sx q[1];
rz(-1.3163047) q[1];
sx q[1];
rz(-0.33448514) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4870105) q[0];
sx q[0];
rz(-2.062487) q[0];
sx q[0];
rz(-2.0915178) q[0];
x q[1];
rz(-2.3804383) q[2];
sx q[2];
rz(-1.8705767) q[2];
sx q[2];
rz(1.3083991) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86395553) q[1];
sx q[1];
rz(-2.5647128) q[1];
sx q[1];
rz(2.3052892) q[1];
rz(-pi) q[2];
rz(1.8556375) q[3];
sx q[3];
rz(-0.99203101) q[3];
sx q[3];
rz(-2.9841686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.9574531) q[2];
sx q[2];
rz(-0.89749557) q[2];
sx q[2];
rz(-1.3857566) q[2];
rz(-0.98006788) q[3];
sx q[3];
rz(-2.6762784) q[3];
sx q[3];
rz(0.050962713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362713) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(-1.3901688) q[0];
rz(-0.35274371) q[1];
sx q[1];
rz(-2.2309062) q[1];
sx q[1];
rz(2.5211451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0095795) q[0];
sx q[0];
rz(-1.4963829) q[0];
sx q[0];
rz(0.090601765) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0376126) q[2];
sx q[2];
rz(-2.1562088) q[2];
sx q[2];
rz(0.81987655) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.57261833) q[1];
sx q[1];
rz(-2.6959722) q[1];
sx q[1];
rz(2.3072412) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87832344) q[3];
sx q[3];
rz(-1.9832423) q[3];
sx q[3];
rz(0.45468047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1227526) q[2];
sx q[2];
rz(-0.41308013) q[2];
sx q[2];
rz(-1.8943465) q[2];
rz(0.16417575) q[3];
sx q[3];
rz(-1.7069867) q[3];
sx q[3];
rz(2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4267047) q[0];
sx q[0];
rz(-2.1737104) q[0];
sx q[0];
rz(1.8545275) q[0];
rz(0.60028589) q[1];
sx q[1];
rz(-1.7900107) q[1];
sx q[1];
rz(0.90369019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0080338399) q[0];
sx q[0];
rz(-1.5045591) q[0];
sx q[0];
rz(-0.107482) q[0];
x q[1];
rz(0.16924567) q[2];
sx q[2];
rz(-2.7566559) q[2];
sx q[2];
rz(-0.16916179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7304222) q[1];
sx q[1];
rz(-2.8287585) q[1];
sx q[1];
rz(2.5223047) q[1];
rz(2.729506) q[3];
sx q[3];
rz(-1.2743605) q[3];
sx q[3];
rz(2.3282098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51957447) q[2];
sx q[2];
rz(-2.3081686) q[2];
sx q[2];
rz(0.77345094) q[2];
rz(-1.2889688) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(2.4066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2499823) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(2.6494359) q[0];
rz(-2.6026169) q[1];
sx q[1];
rz(-1.69918) q[1];
sx q[1];
rz(-2.0416562) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625222) q[0];
sx q[0];
rz(-1.8305792) q[0];
sx q[0];
rz(0.35210877) q[0];
rz(-pi) q[1];
x q[1];
rz(0.01077588) q[2];
sx q[2];
rz(-2.1824565) q[2];
sx q[2];
rz(1.7623368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27515181) q[1];
sx q[1];
rz(-0.76443618) q[1];
sx q[1];
rz(1.4363097) q[1];
x q[2];
rz(0.57657974) q[3];
sx q[3];
rz(-2.5825273) q[3];
sx q[3];
rz(-0.41043974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8010572) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(-2.578793) q[2];
rz(-2.4417012) q[3];
sx q[3];
rz(-1.2908582) q[3];
sx q[3];
rz(-2.8904397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3444779) q[0];
sx q[0];
rz(-2.4176702) q[0];
sx q[0];
rz(-2.7864454) q[0];
rz(-1.2190602) q[1];
sx q[1];
rz(-1.9025758) q[1];
sx q[1];
rz(0.452279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4773524) q[0];
sx q[0];
rz(-1.1731804) q[0];
sx q[0];
rz(-0.94066046) q[0];
x q[1];
rz(2.3269749) q[2];
sx q[2];
rz(-2.6067408) q[2];
sx q[2];
rz(-1.2592821) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7004343) q[1];
sx q[1];
rz(-1.5467073) q[1];
sx q[1];
rz(-1.248293) q[1];
rz(-pi) q[2];
rz(1.7132492) q[3];
sx q[3];
rz(-0.94678426) q[3];
sx q[3];
rz(-0.54572661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63048116) q[2];
sx q[2];
rz(-0.52383542) q[2];
sx q[2];
rz(3.0044978) q[2];
rz(1.5591722) q[3];
sx q[3];
rz(-1.1538006) q[3];
sx q[3];
rz(-2.3649575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.1906076) q[0];
sx q[0];
rz(-0.76222104) q[0];
sx q[0];
rz(-1.5964339) q[0];
rz(0.51721382) q[1];
sx q[1];
rz(-1.1511753) q[1];
sx q[1];
rz(-2.1545765) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1008489) q[0];
sx q[0];
rz(-0.13253875) q[0];
sx q[0];
rz(2.9805471) q[0];
rz(-2.1783243) q[2];
sx q[2];
rz(-1.844256) q[2];
sx q[2];
rz(0.69404049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.034711866) q[1];
sx q[1];
rz(-1.3415404) q[1];
sx q[1];
rz(2.8696612) q[1];
rz(-pi) q[2];
rz(0.084264755) q[3];
sx q[3];
rz(-1.2532506) q[3];
sx q[3];
rz(-1.0061044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3551657) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(-0.008358566) q[2];
rz(0.23056325) q[3];
sx q[3];
rz(-1.9356666) q[3];
sx q[3];
rz(-1.4768538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.01345988) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(1.6101366) q[0];
rz(1.5229335) q[1];
sx q[1];
rz(-1.5459272) q[1];
sx q[1];
rz(1.5628372) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42199907) q[0];
sx q[0];
rz(-0.61924705) q[0];
sx q[0];
rz(-2.8887755) q[0];
rz(-pi) q[1];
rz(-0.35308102) q[2];
sx q[2];
rz(-2.7432979) q[2];
sx q[2];
rz(-0.19206599) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4304428) q[1];
sx q[1];
rz(-2.4949269) q[1];
sx q[1];
rz(2.4025687) q[1];
rz(-pi) q[2];
rz(0.39852972) q[3];
sx q[3];
rz(-0.70793286) q[3];
sx q[3];
rz(1.7205451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83850399) q[2];
sx q[2];
rz(-1.1469301) q[2];
sx q[2];
rz(3.1040891) q[2];
rz(-0.23669067) q[3];
sx q[3];
rz(-2.762837) q[3];
sx q[3];
rz(1.2934575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732052) q[0];
sx q[0];
rz(-0.78998843) q[0];
sx q[0];
rz(-1.0733806) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-0.57917246) q[1];
sx q[1];
rz(-2.5198708) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17162308) q[0];
sx q[0];
rz(-1.8299654) q[0];
sx q[0];
rz(0.39879946) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66017588) q[2];
sx q[2];
rz(-1.4731506) q[2];
sx q[2];
rz(2.4843189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.67111193) q[1];
sx q[1];
rz(-0.59146229) q[1];
sx q[1];
rz(0.5670814) q[1];
x q[2];
rz(2.5863566) q[3];
sx q[3];
rz(-0.16517565) q[3];
sx q[3];
rz(-0.092008807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5198034) q[2];
sx q[2];
rz(-0.18814627) q[2];
sx q[2];
rz(0.56358799) q[2];
rz(1.311519) q[3];
sx q[3];
rz(-1.2389641) q[3];
sx q[3];
rz(-0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9594864) q[0];
sx q[0];
rz(-2.2725548) q[0];
sx q[0];
rz(-2.2667789) q[0];
rz(-1.733571) q[1];
sx q[1];
rz(-2.4751016) q[1];
sx q[1];
rz(0.63546884) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498813) q[0];
sx q[0];
rz(-0.81560613) q[0];
sx q[0];
rz(0.41378234) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0032015) q[2];
sx q[2];
rz(-0.48241189) q[2];
sx q[2];
rz(2.8495827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6136352) q[1];
sx q[1];
rz(-1.0571169) q[1];
sx q[1];
rz(-0.46544816) q[1];
rz(1.0069153) q[3];
sx q[3];
rz(-1.0936519) q[3];
sx q[3];
rz(-1.7686219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1982939) q[2];
sx q[2];
rz(-1.0903) q[2];
sx q[2];
rz(2.3401006) q[2];
rz(1.9258026) q[3];
sx q[3];
rz(-0.91167584) q[3];
sx q[3];
rz(-0.041778684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53032482) q[0];
sx q[0];
rz(-1.4401191) q[0];
sx q[0];
rz(-2.8339207) q[0];
rz(1.8511741) q[1];
sx q[1];
rz(-2.1967874) q[1];
sx q[1];
rz(-2.8312942) q[1];
rz(-0.63772884) q[2];
sx q[2];
rz(-0.22508937) q[2];
sx q[2];
rz(-1.3300016) q[2];
rz(0.17388969) q[3];
sx q[3];
rz(-0.66413838) q[3];
sx q[3];
rz(-0.037353368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
