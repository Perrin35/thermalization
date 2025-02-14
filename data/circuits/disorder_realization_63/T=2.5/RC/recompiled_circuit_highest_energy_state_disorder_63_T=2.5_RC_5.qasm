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
rz(1.0032049) q[0];
sx q[0];
rz(1.9411074) q[0];
sx q[0];
rz(8.3984126) q[0];
rz(0.29844555) q[1];
sx q[1];
rz(-1.6331853) q[1];
sx q[1];
rz(1.0025947) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110785) q[0];
sx q[0];
rz(-1.7745695) q[0];
sx q[0];
rz(3.1183) q[0];
x q[1];
rz(0.070656405) q[2];
sx q[2];
rz(-0.66412726) q[2];
sx q[2];
rz(-0.16985591) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94002515) q[1];
sx q[1];
rz(-2.7301412) q[1];
sx q[1];
rz(0.33951851) q[1];
rz(-2.3839124) q[3];
sx q[3];
rz(-0.79059383) q[3];
sx q[3];
rz(-2.9942715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5528494) q[2];
sx q[2];
rz(-2.1407949) q[2];
sx q[2];
rz(3.106485) q[2];
rz(-1.9378174) q[3];
sx q[3];
rz(-1.8193918) q[3];
sx q[3];
rz(-0.28425851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0725919) q[0];
sx q[0];
rz(-0.44206107) q[0];
sx q[0];
rz(2.1373855) q[0];
rz(-0.12241157) q[1];
sx q[1];
rz(-1.6740084) q[1];
sx q[1];
rz(2.4845128) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5733684) q[0];
sx q[0];
rz(-1.4889297) q[0];
sx q[0];
rz(-2.7533145) q[0];
rz(-pi) q[1];
rz(-2.0900656) q[2];
sx q[2];
rz(-2.8247571) q[2];
sx q[2];
rz(2.6749382) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66533684) q[1];
sx q[1];
rz(-1.0968535) q[1];
sx q[1];
rz(0.83935763) q[1];
rz(1.8155926) q[3];
sx q[3];
rz(-0.80202937) q[3];
sx q[3];
rz(2.0662465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.55679) q[2];
sx q[2];
rz(-0.57791296) q[2];
sx q[2];
rz(-2.54971) q[2];
rz(-1.6651734) q[3];
sx q[3];
rz(-1.6074601) q[3];
sx q[3];
rz(2.2679451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0037435) q[0];
sx q[0];
rz(-2.9326404) q[0];
sx q[0];
rz(-1.9269706) q[0];
rz(2.1851723) q[1];
sx q[1];
rz(-1.950187) q[1];
sx q[1];
rz(-1.1716243) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3879165) q[0];
sx q[0];
rz(-1.9646586) q[0];
sx q[0];
rz(2.0772165) q[0];
x q[1];
rz(2.2713216) q[2];
sx q[2];
rz(-1.4129801) q[2];
sx q[2];
rz(-1.4627139) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4332648) q[1];
sx q[1];
rz(-1.9076132) q[1];
sx q[1];
rz(-2.5418806) q[1];
rz(-pi) q[2];
rz(-0.089853386) q[3];
sx q[3];
rz(-1.615534) q[3];
sx q[3];
rz(2.5344332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9413635) q[2];
sx q[2];
rz(-2.0912632) q[2];
sx q[2];
rz(-1.7639147) q[2];
rz(-2.7348943) q[3];
sx q[3];
rz(-0.64928693) q[3];
sx q[3];
rz(-2.6810834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3520626) q[0];
sx q[0];
rz(-1.8822414) q[0];
sx q[0];
rz(1.9011185) q[0];
rz(-2.4234096) q[1];
sx q[1];
rz(-1.8156464) q[1];
sx q[1];
rz(2.9026418) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093042308) q[0];
sx q[0];
rz(-1.6918381) q[0];
sx q[0];
rz(0.28994513) q[0];
rz(-pi) q[1];
rz(2.3283655) q[2];
sx q[2];
rz(-1.0095749) q[2];
sx q[2];
rz(2.6475069) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0480369) q[1];
sx q[1];
rz(-0.29680064) q[1];
sx q[1];
rz(1.0637299) q[1];
x q[2];
rz(1.7051058) q[3];
sx q[3];
rz(-2.1013936) q[3];
sx q[3];
rz(2.1710896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6255528) q[2];
sx q[2];
rz(-1.7480787) q[2];
sx q[2];
rz(1.3429406) q[2];
rz(0.88137734) q[3];
sx q[3];
rz(-2.9727029) q[3];
sx q[3];
rz(-0.78099293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7278904) q[0];
sx q[0];
rz(-0.24402937) q[0];
sx q[0];
rz(-1.4843041) q[0];
rz(-0.65385747) q[1];
sx q[1];
rz(-1.6203208) q[1];
sx q[1];
rz(-0.13793764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50820275) q[0];
sx q[0];
rz(-1.422056) q[0];
sx q[0];
rz(-1.9277907) q[0];
rz(0.65278585) q[2];
sx q[2];
rz(-1.372027) q[2];
sx q[2];
rz(2.7592289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5522462) q[1];
sx q[1];
rz(-1.4728496) q[1];
sx q[1];
rz(-1.1393273) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3119427) q[3];
sx q[3];
rz(-1.6248684) q[3];
sx q[3];
rz(1.8340221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58978927) q[2];
sx q[2];
rz(-1.1729596) q[2];
sx q[2];
rz(-0.087184437) q[2];
rz(-0.11463556) q[3];
sx q[3];
rz(-0.81511027) q[3];
sx q[3];
rz(0.4172999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7755985) q[0];
sx q[0];
rz(-1.5394779) q[0];
sx q[0];
rz(-0.48496801) q[0];
rz(0.94351774) q[1];
sx q[1];
rz(-0.8129932) q[1];
sx q[1];
rz(1.9224723) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2526357) q[0];
sx q[0];
rz(-0.85272861) q[0];
sx q[0];
rz(-3.0120993) q[0];
rz(-pi) q[1];
rz(1.2487683) q[2];
sx q[2];
rz(-1.5862984) q[2];
sx q[2];
rz(-0.44336927) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18029848) q[1];
sx q[1];
rz(-1.4438259) q[1];
sx q[1];
rz(-3.0650223) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0659473) q[3];
sx q[3];
rz(-2.2643914) q[3];
sx q[3];
rz(-1.1923196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.24011) q[2];
sx q[2];
rz(-1.0604475) q[2];
sx q[2];
rz(-0.94129747) q[2];
rz(2.8955722) q[3];
sx q[3];
rz(-1.4964208) q[3];
sx q[3];
rz(-1.7814319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0435903) q[0];
sx q[0];
rz(-1.101838) q[0];
sx q[0];
rz(-1.164042) q[0];
rz(-1.7835435) q[1];
sx q[1];
rz(-0.86654228) q[1];
sx q[1];
rz(1.39303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355125) q[0];
sx q[0];
rz(-2.7580259) q[0];
sx q[0];
rz(-1.7696898) q[0];
rz(1.8556526) q[2];
sx q[2];
rz(-1.0394888) q[2];
sx q[2];
rz(1.667995) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6026615) q[1];
sx q[1];
rz(-1.4653413) q[1];
sx q[1];
rz(0.47172539) q[1];
rz(-pi) q[2];
rz(0.22402899) q[3];
sx q[3];
rz(-0.90440905) q[3];
sx q[3];
rz(1.9206604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67036575) q[2];
sx q[2];
rz(-1.7134075) q[2];
sx q[2];
rz(2.8426389) q[2];
rz(2.7361338) q[3];
sx q[3];
rz(-0.36646989) q[3];
sx q[3];
rz(2.5772742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7120755) q[0];
sx q[0];
rz(-0.27292621) q[0];
sx q[0];
rz(2.3522229) q[0];
rz(-2.7366267) q[1];
sx q[1];
rz(-1.7908432) q[1];
sx q[1];
rz(0.49496034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4433561) q[0];
sx q[0];
rz(-2.3670417) q[0];
sx q[0];
rz(0.86408918) q[0];
rz(-pi) q[1];
rz(-2.8414824) q[2];
sx q[2];
rz(-1.0213739) q[2];
sx q[2];
rz(0.97893366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72515762) q[1];
sx q[1];
rz(-0.817653) q[1];
sx q[1];
rz(-2.8055951) q[1];
rz(2.1405419) q[3];
sx q[3];
rz(-2.3214139) q[3];
sx q[3];
rz(0.4789744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1627545) q[2];
sx q[2];
rz(-1.2252204) q[2];
sx q[2];
rz(0.2552574) q[2];
rz(3.1066331) q[3];
sx q[3];
rz(-1.5233663) q[3];
sx q[3];
rz(-2.8953654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73664767) q[0];
sx q[0];
rz(-1.0628137) q[0];
sx q[0];
rz(-2.431562) q[0];
rz(2.1404449) q[1];
sx q[1];
rz(-2.2444057) q[1];
sx q[1];
rz(0.62090105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9644755) q[0];
sx q[0];
rz(-1.6311446) q[0];
sx q[0];
rz(1.9326769) q[0];
rz(-pi) q[1];
x q[1];
rz(2.484637) q[2];
sx q[2];
rz(-1.1816506) q[2];
sx q[2];
rz(2.3537113) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94540751) q[1];
sx q[1];
rz(-1.2239478) q[1];
sx q[1];
rz(1.6849243) q[1];
rz(-pi) q[2];
rz(2.8594703) q[3];
sx q[3];
rz(-2.3446313) q[3];
sx q[3];
rz(-2.7180501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8093449) q[2];
sx q[2];
rz(-1.6184018) q[2];
sx q[2];
rz(0.28918949) q[2];
rz(2.0293763) q[3];
sx q[3];
rz(-0.29505348) q[3];
sx q[3];
rz(1.0203863) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6638829) q[0];
sx q[0];
rz(-1.7998671) q[0];
sx q[0];
rz(-2.2148602) q[0];
rz(2.0345188) q[1];
sx q[1];
rz(-1.5137545) q[1];
sx q[1];
rz(-0.56799299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04161373) q[0];
sx q[0];
rz(-1.8868514) q[0];
sx q[0];
rz(0.63412447) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6935518) q[2];
sx q[2];
rz(-2.5841641) q[2];
sx q[2];
rz(-0.008226062) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5458249) q[1];
sx q[1];
rz(-0.95965702) q[1];
sx q[1];
rz(0.45054884) q[1];
rz(2.5015478) q[3];
sx q[3];
rz(-2.0012104) q[3];
sx q[3];
rz(-2.0295582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0893112) q[2];
sx q[2];
rz(-2.3312882) q[2];
sx q[2];
rz(1.7029765) q[2];
rz(2.8449521) q[3];
sx q[3];
rz(-0.34160015) q[3];
sx q[3];
rz(-2.6945485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68375568) q[0];
sx q[0];
rz(-2.044027) q[0];
sx q[0];
rz(-2.5720163) q[0];
rz(2.8029022) q[1];
sx q[1];
rz(-1.6117922) q[1];
sx q[1];
rz(1.502996) q[1];
rz(2.4581134) q[2];
sx q[2];
rz(-0.31961315) q[2];
sx q[2];
rz(1.4994958) q[2];
rz(-1.8453311) q[3];
sx q[3];
rz(-1.408315) q[3];
sx q[3];
rz(-2.3111865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
