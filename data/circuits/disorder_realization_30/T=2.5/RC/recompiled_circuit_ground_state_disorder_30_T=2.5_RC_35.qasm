OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.423288) q[0];
sx q[0];
rz(-2.365132) q[0];
sx q[0];
rz(-2.7756696) q[0];
rz(0.34691063) q[1];
sx q[1];
rz(-0.28471714) q[1];
sx q[1];
rz(-1.3887583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4735218) q[0];
sx q[0];
rz(-1.6374303) q[0];
sx q[0];
rz(-0.20247831) q[0];
rz(-pi) q[1];
rz(-0.35438673) q[2];
sx q[2];
rz(-2.4172108) q[2];
sx q[2];
rz(-2.5222833) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5468816) q[1];
sx q[1];
rz(-1.1818019) q[1];
sx q[1];
rz(-1.0884029) q[1];
x q[2];
rz(0.71309375) q[3];
sx q[3];
rz(-0.72848195) q[3];
sx q[3];
rz(-0.99458867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95691386) q[2];
sx q[2];
rz(-2.2214486) q[2];
sx q[2];
rz(-2.206395) q[2];
rz(2.8397371) q[3];
sx q[3];
rz(-1.5422834) q[3];
sx q[3];
rz(-0.39366084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97717706) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(2.6432977) q[0];
rz(1.4277108) q[1];
sx q[1];
rz(-1.9352813) q[1];
sx q[1];
rz(-2.0548342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55288494) q[0];
sx q[0];
rz(-2.4037139) q[0];
sx q[0];
rz(0.35564519) q[0];
rz(-pi) q[1];
rz(0.19379036) q[2];
sx q[2];
rz(-1.8326129) q[2];
sx q[2];
rz(-1.7270391) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6689008) q[1];
sx q[1];
rz(-1.7311061) q[1];
sx q[1];
rz(2.7008788) q[1];
rz(-1.5828176) q[3];
sx q[3];
rz(-1.2526027) q[3];
sx q[3];
rz(2.1592922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9429417) q[2];
sx q[2];
rz(-0.31060878) q[2];
sx q[2];
rz(1.3322213) q[2];
rz(-1.9940935) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(-1.8089627) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9320817) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(2.984356) q[0];
rz(-1.2278185) q[1];
sx q[1];
rz(-2.9324052) q[1];
sx q[1];
rz(0.42207178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2651208) q[0];
sx q[0];
rz(-0.19959627) q[0];
sx q[0];
rz(1.100698) q[0];
rz(-pi) q[1];
rz(-1.7447168) q[2];
sx q[2];
rz(-1.2971767) q[2];
sx q[2];
rz(-1.5762941) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2873309) q[1];
sx q[1];
rz(-2.7497083) q[1];
sx q[1];
rz(-1.803483) q[1];
x q[2];
rz(2.3396569) q[3];
sx q[3];
rz(-1.7953331) q[3];
sx q[3];
rz(0.79659407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.813039) q[2];
sx q[2];
rz(-0.96031323) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(-1.0775393) q[3];
sx q[3];
rz(-2.4534093) q[3];
sx q[3];
rz(0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587191) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(3.1251113) q[0];
rz(-0.074542848) q[1];
sx q[1];
rz(-0.45769474) q[1];
sx q[1];
rz(-0.25513908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82127386) q[0];
sx q[0];
rz(-2.6373568) q[0];
sx q[0];
rz(2.8404499) q[0];
rz(-pi) q[1];
rz(-2.8066977) q[2];
sx q[2];
rz(-2.1112295) q[2];
sx q[2];
rz(1.188736) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.093016) q[1];
sx q[1];
rz(-0.79461351) q[1];
sx q[1];
rz(-0.62127415) q[1];
rz(1.9202407) q[3];
sx q[3];
rz(-1.6779052) q[3];
sx q[3];
rz(-2.8529594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5269346) q[2];
sx q[2];
rz(-2.4444828) q[2];
sx q[2];
rz(0.69980168) q[2];
rz(-0.091863306) q[3];
sx q[3];
rz(-1.9242761) q[3];
sx q[3];
rz(-2.6868668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62234539) q[0];
sx q[0];
rz(-2.3212101) q[0];
sx q[0];
rz(-1.1118332) q[0];
rz(2.9513997) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(1.9817188) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1721508) q[0];
sx q[0];
rz(-1.5363201) q[0];
sx q[0];
rz(1.9146754) q[0];
rz(1.3746475) q[2];
sx q[2];
rz(-2.4939459) q[2];
sx q[2];
rz(0.37774936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76471699) q[1];
sx q[1];
rz(-0.60827601) q[1];
sx q[1];
rz(2.6981955) q[1];
rz(-pi) q[2];
rz(-1.492736) q[3];
sx q[3];
rz(-0.52320899) q[3];
sx q[3];
rz(-1.4735492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5708892) q[2];
sx q[2];
rz(-2.3184226) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(1.0620091) q[3];
sx q[3];
rz(-1.8657203) q[3];
sx q[3];
rz(1.6331858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.688711) q[0];
sx q[0];
rz(-2.0827561) q[0];
sx q[0];
rz(1.3558615) q[0];
rz(-0.52976766) q[1];
sx q[1];
rz(-0.7822839) q[1];
sx q[1];
rz(-1.9897602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6731176) q[0];
sx q[0];
rz(-1.6694607) q[0];
sx q[0];
rz(2.7844546) q[0];
x q[1];
rz(1.4131613) q[2];
sx q[2];
rz(-0.54716483) q[2];
sx q[2];
rz(2.5131651) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63910145) q[1];
sx q[1];
rz(-1.7646953) q[1];
sx q[1];
rz(-2.7706233) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21422503) q[3];
sx q[3];
rz(-1.6856632) q[3];
sx q[3];
rz(-1.8095995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6666224) q[2];
sx q[2];
rz(-2.5002561) q[2];
sx q[2];
rz(2.6744911) q[2];
rz(-1.0111672) q[3];
sx q[3];
rz(-2.7979388) q[3];
sx q[3];
rz(-1.7721734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8491299) q[0];
sx q[0];
rz(-2.9865773) q[0];
sx q[0];
rz(-2.9169061) q[0];
rz(-0.16381964) q[1];
sx q[1];
rz(-2.8024709) q[1];
sx q[1];
rz(2.4952369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3098329) q[0];
sx q[0];
rz(-2.0947959) q[0];
sx q[0];
rz(2.8519467) q[0];
rz(-pi) q[1];
rz(0.21804131) q[2];
sx q[2];
rz(-1.8865117) q[2];
sx q[2];
rz(1.174508) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.028370628) q[1];
sx q[1];
rz(-0.54873768) q[1];
sx q[1];
rz(-0.1446677) q[1];
rz(2.7541646) q[3];
sx q[3];
rz(-1.4447851) q[3];
sx q[3];
rz(-0.84539094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9859163) q[2];
sx q[2];
rz(-1.6935657) q[2];
sx q[2];
rz(0.4099561) q[2];
rz(0.83034596) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.603867) q[0];
sx q[0];
rz(-2.4483838) q[0];
sx q[0];
rz(0.24630462) q[0];
rz(-0.82659563) q[1];
sx q[1];
rz(-1.3984503) q[1];
sx q[1];
rz(2.8776317) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0279585) q[0];
sx q[0];
rz(-1.6662681) q[0];
sx q[0];
rz(-0.0072633538) q[0];
x q[1];
rz(-2.0461778) q[2];
sx q[2];
rz(-2.2959297) q[2];
sx q[2];
rz(2.1755168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3726793) q[1];
sx q[1];
rz(-2.2376334) q[1];
sx q[1];
rz(-0.36075488) q[1];
rz(1.1353536) q[3];
sx q[3];
rz(-2.1135984) q[3];
sx q[3];
rz(2.8538728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0366514) q[2];
sx q[2];
rz(-1.8737326) q[2];
sx q[2];
rz(-2.4264753) q[2];
rz(0.07130833) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(0.74688545) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1320553) q[0];
sx q[0];
rz(-1.4893463) q[0];
sx q[0];
rz(2.706053) q[0];
rz(2.9938193) q[1];
sx q[1];
rz(-1.7099893) q[1];
sx q[1];
rz(1.6730283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4184121) q[0];
sx q[0];
rz(-1.3910798) q[0];
sx q[0];
rz(0.43850684) q[0];
x q[1];
rz(1.4917489) q[2];
sx q[2];
rz(-0.87495308) q[2];
sx q[2];
rz(0.33338132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0316511) q[1];
sx q[1];
rz(-1.1835956) q[1];
sx q[1];
rz(-2.9613858) q[1];
rz(-pi) q[2];
rz(1.3545669) q[3];
sx q[3];
rz(-1.9659967) q[3];
sx q[3];
rz(-0.54602269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74464166) q[2];
sx q[2];
rz(-1.7640742) q[2];
sx q[2];
rz(2.2486539) q[2];
rz(1.2785771) q[3];
sx q[3];
rz(-1.9847001) q[3];
sx q[3];
rz(1.4634092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8466012) q[0];
sx q[0];
rz(-0.33116594) q[0];
sx q[0];
rz(-2.5355329) q[0];
rz(3.1312969) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(-3.0616679) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3256667) q[0];
sx q[0];
rz(-1.0116966) q[0];
sx q[0];
rz(-2.8432106) q[0];
rz(-pi) q[1];
rz(1.2508008) q[2];
sx q[2];
rz(-2.3872132) q[2];
sx q[2];
rz(2.2571074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8680758) q[1];
sx q[1];
rz(-2.7339327) q[1];
sx q[1];
rz(3.0855387) q[1];
rz(-pi) q[2];
rz(0.20584835) q[3];
sx q[3];
rz(-1.7966812) q[3];
sx q[3];
rz(-3.1354648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14110485) q[2];
sx q[2];
rz(-0.92685574) q[2];
sx q[2];
rz(-2.1263988) q[2];
rz(-2.5403533) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(-2.0950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4973608) q[0];
sx q[0];
rz(-1.624122) q[0];
sx q[0];
rz(-2.4540785) q[0];
rz(-0.58912206) q[1];
sx q[1];
rz(-2.4644869) q[1];
sx q[1];
rz(-0.45225515) q[1];
rz(2.6964006) q[2];
sx q[2];
rz(-1.6576851) q[2];
sx q[2];
rz(0.97752006) q[2];
rz(2.3169869) q[3];
sx q[3];
rz(-1.6578703) q[3];
sx q[3];
rz(-0.35005611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
