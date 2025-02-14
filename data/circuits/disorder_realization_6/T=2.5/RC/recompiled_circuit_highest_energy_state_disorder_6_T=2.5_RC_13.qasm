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
rz(-1.8974798) q[0];
sx q[0];
rz(-0.31960684) q[0];
sx q[0];
rz(2.6407916) q[0];
rz(0.29419857) q[1];
sx q[1];
rz(3.8837641) q[1];
sx q[1];
rz(9.972218) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476737) q[0];
sx q[0];
rz(-1.8800056) q[0];
sx q[0];
rz(-1.4178314) q[0];
x q[1];
rz(2.6771678) q[2];
sx q[2];
rz(-1.8193442) q[2];
sx q[2];
rz(2.3119702) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55996229) q[1];
sx q[1];
rz(-2.0764372) q[1];
sx q[1];
rz(1.4675114) q[1];
rz(2.1087113) q[3];
sx q[3];
rz(-1.8807372) q[3];
sx q[3];
rz(0.021837318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1529634) q[2];
sx q[2];
rz(-0.84501481) q[2];
sx q[2];
rz(2.5353234) q[2];
rz(-1.9340949) q[3];
sx q[3];
rz(-1.2946125) q[3];
sx q[3];
rz(-1.5906364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.9015053) q[0];
sx q[0];
rz(-1.6143687) q[0];
sx q[0];
rz(2.4497581) q[0];
rz(-1.2884619) q[1];
sx q[1];
rz(-1.144578) q[1];
sx q[1];
rz(-2.8537234) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3436056) q[0];
sx q[0];
rz(-1.6988861) q[0];
sx q[0];
rz(-1.8592181) q[0];
rz(-2.4601654) q[2];
sx q[2];
rz(-0.81290258) q[2];
sx q[2];
rz(0.63804783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58093666) q[1];
sx q[1];
rz(-1.1252778) q[1];
sx q[1];
rz(1.6220644) q[1];
rz(-pi) q[2];
rz(3.0045142) q[3];
sx q[3];
rz(-1.6774639) q[3];
sx q[3];
rz(2.9003473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3466779) q[2];
sx q[2];
rz(-1.4823464) q[2];
sx q[2];
rz(-2.1686926) q[2];
rz(-1.3341303) q[3];
sx q[3];
rz(-2.3502246) q[3];
sx q[3];
rz(2.5202675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7767104) q[0];
sx q[0];
rz(-0.24672395) q[0];
sx q[0];
rz(-0.2488939) q[0];
rz(1.6913255) q[1];
sx q[1];
rz(-1.1509117) q[1];
sx q[1];
rz(2.2379564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8189175) q[0];
sx q[0];
rz(-1.5455017) q[0];
sx q[0];
rz(2.5324244) q[0];
rz(0.10714178) q[2];
sx q[2];
rz(-1.2276548) q[2];
sx q[2];
rz(-1.8245152) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8775656) q[1];
sx q[1];
rz(-2.2912044) q[1];
sx q[1];
rz(-2.6191465) q[1];
x q[2];
rz(0.60402212) q[3];
sx q[3];
rz(-1.9088749) q[3];
sx q[3];
rz(-3.1228408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9768208) q[2];
sx q[2];
rz(-0.73126078) q[2];
sx q[2];
rz(-0.33052793) q[2];
rz(2.9183689) q[3];
sx q[3];
rz(-1.1130755) q[3];
sx q[3];
rz(-0.37842512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44499236) q[0];
sx q[0];
rz(-1.3560504) q[0];
sx q[0];
rz(-0.14183216) q[0];
rz(0.089135535) q[1];
sx q[1];
rz(-0.24239692) q[1];
sx q[1];
rz(2.1884122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87751276) q[0];
sx q[0];
rz(-2.3126214) q[0];
sx q[0];
rz(2.6116051) q[0];
rz(1.6905027) q[2];
sx q[2];
rz(-1.9485932) q[2];
sx q[2];
rz(-2.1816513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.61342) q[1];
sx q[1];
rz(-1.0549567) q[1];
sx q[1];
rz(-2.1891761) q[1];
x q[2];
rz(-2.6450883) q[3];
sx q[3];
rz(-1.3937573) q[3];
sx q[3];
rz(-2.4868709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7715093) q[2];
sx q[2];
rz(-0.13088317) q[2];
sx q[2];
rz(-2.6123602) q[2];
rz(1.0445163) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(-2.9862826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20045497) q[0];
sx q[0];
rz(-1.8998242) q[0];
sx q[0];
rz(-0.45503765) q[0];
rz(-0.014668839) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(-0.98247772) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331086) q[0];
sx q[0];
rz(-2.8245375) q[0];
sx q[0];
rz(0.32755537) q[0];
rz(2.2682954) q[2];
sx q[2];
rz(-1.1050537) q[2];
sx q[2];
rz(2.1312775) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3396709) q[1];
sx q[1];
rz(-1.675444) q[1];
sx q[1];
rz(-2.8644979) q[1];
rz(-1.3339046) q[3];
sx q[3];
rz(-1.8474527) q[3];
sx q[3];
rz(1.2048126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5386397) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(0.61720103) q[2];
rz(-1.3868825) q[3];
sx q[3];
rz(-0.91104561) q[3];
sx q[3];
rz(2.9673747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0844326) q[0];
sx q[0];
rz(-2.3800157) q[0];
sx q[0];
rz(-0.97849751) q[0];
rz(-1.017336) q[1];
sx q[1];
rz(-2.8637736) q[1];
sx q[1];
rz(2.691958) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89116865) q[0];
sx q[0];
rz(-0.40968597) q[0];
sx q[0];
rz(-2.4158359) q[0];
rz(-pi) q[1];
rz(1.1005747) q[2];
sx q[2];
rz(-0.82056773) q[2];
sx q[2];
rz(0.596753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3510531) q[1];
sx q[1];
rz(-1.0996523) q[1];
sx q[1];
rz(-3.0912249) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95532598) q[3];
sx q[3];
rz(-2.1571473) q[3];
sx q[3];
rz(-2.9406656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52266878) q[2];
sx q[2];
rz(-1.7974682) q[2];
sx q[2];
rz(1.3783003) q[2];
rz(1.1528692) q[3];
sx q[3];
rz(-1.1533573) q[3];
sx q[3];
rz(-2.8809179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5744837) q[0];
sx q[0];
rz(-0.672988) q[0];
sx q[0];
rz(-2.8048977) q[0];
rz(2.4163272) q[1];
sx q[1];
rz(-1.5545574) q[1];
sx q[1];
rz(2.7560962) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3025997) q[0];
sx q[0];
rz(-1.2205048) q[0];
sx q[0];
rz(-2.4953609) q[0];
x q[1];
rz(1.0573559) q[2];
sx q[2];
rz(-0.2845736) q[2];
sx q[2];
rz(-1.3651266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.817877) q[1];
sx q[1];
rz(-0.95600946) q[1];
sx q[1];
rz(2.7842658) q[1];
rz(-2.518744) q[3];
sx q[3];
rz(-1.9241352) q[3];
sx q[3];
rz(1.6652157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54726279) q[2];
sx q[2];
rz(-1.2196093) q[2];
sx q[2];
rz(-2.8426113) q[2];
rz(-0.79214054) q[3];
sx q[3];
rz(-2.7463089) q[3];
sx q[3];
rz(1.7637174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2250273) q[0];
sx q[0];
rz(-1.4341609) q[0];
sx q[0];
rz(2.3634971) q[0];
rz(1.3968702) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(3.0527589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7185811) q[0];
sx q[0];
rz(-1.6751583) q[0];
sx q[0];
rz(2.7947438) q[0];
x q[1];
rz(-0.5702604) q[2];
sx q[2];
rz(-2.2024254) q[2];
sx q[2];
rz(2.1960902) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.962304) q[1];
sx q[1];
rz(-1.3109164) q[1];
sx q[1];
rz(1.13729) q[1];
rz(-pi) q[2];
rz(-1.6501637) q[3];
sx q[3];
rz(-0.36874394) q[3];
sx q[3];
rz(-2.6811932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5953411) q[2];
sx q[2];
rz(-0.83117861) q[2];
sx q[2];
rz(-1.1031995) q[2];
rz(0.57175076) q[3];
sx q[3];
rz(-2.5751556) q[3];
sx q[3];
rz(1.5248689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53200805) q[0];
sx q[0];
rz(-0.12140618) q[0];
sx q[0];
rz(-2.4914361) q[0];
rz(2.0434168) q[1];
sx q[1];
rz(-1.8868586) q[1];
sx q[1];
rz(0.065902725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1797721) q[0];
sx q[0];
rz(-2.0393436) q[0];
sx q[0];
rz(2.382032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6292105) q[2];
sx q[2];
rz(-1.3810535) q[2];
sx q[2];
rz(2.3924215) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53532584) q[1];
sx q[1];
rz(-2.3542778) q[1];
sx q[1];
rz(1.3523065) q[1];
rz(-pi) q[2];
rz(0.66146694) q[3];
sx q[3];
rz(-1.3693878) q[3];
sx q[3];
rz(-0.097518343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99964511) q[2];
sx q[2];
rz(-1.9530752) q[2];
sx q[2];
rz(2.578242) q[2];
rz(-2.3422286) q[3];
sx q[3];
rz(-1.0146217) q[3];
sx q[3];
rz(0.25580251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19485165) q[0];
sx q[0];
rz(-0.020337157) q[0];
sx q[0];
rz(2.643423) q[0];
rz(2.8794471) q[1];
sx q[1];
rz(-1.405895) q[1];
sx q[1];
rz(1.3462876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57980832) q[0];
sx q[0];
rz(-1.1129079) q[0];
sx q[0];
rz(2.3696128) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5813229) q[2];
sx q[2];
rz(-1.5731647) q[2];
sx q[2];
rz(-2.6244768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87835364) q[1];
sx q[1];
rz(-1.2290956) q[1];
sx q[1];
rz(-0.71818761) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78018434) q[3];
sx q[3];
rz(-1.7276075) q[3];
sx q[3];
rz(-3.1136857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6654309) q[2];
sx q[2];
rz(-0.68887812) q[2];
sx q[2];
rz(0.8821744) q[2];
rz(2.4836508) q[3];
sx q[3];
rz(-2.0177757) q[3];
sx q[3];
rz(-0.9995802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5049725) q[0];
sx q[0];
rz(-1.5047147) q[0];
sx q[0];
rz(2.104105) q[0];
rz(2.6344015) q[1];
sx q[1];
rz(-2.3585944) q[1];
sx q[1];
rz(2.0814887) q[1];
rz(1.3846332) q[2];
sx q[2];
rz(-1.2825479) q[2];
sx q[2];
rz(-2.1946236) q[2];
rz(0.67340452) q[3];
sx q[3];
rz(-2.2670204) q[3];
sx q[3];
rz(-1.9011998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
