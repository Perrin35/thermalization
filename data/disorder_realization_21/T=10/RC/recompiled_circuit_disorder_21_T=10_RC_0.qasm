OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(6.364967) q[0];
sx q[0];
rz(9.9262417) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8034536) q[0];
sx q[0];
rz(-0.33284602) q[0];
sx q[0];
rz(1.8344318) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6104923) q[2];
sx q[2];
rz(-2.1571026) q[2];
sx q[2];
rz(0.58639484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7079561) q[1];
sx q[1];
rz(-2.4383713) q[1];
sx q[1];
rz(-1.0183079) q[1];
rz(-pi) q[2];
rz(-2.6775042) q[3];
sx q[3];
rz(-2.1543192) q[3];
sx q[3];
rz(0.97809631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(0.15727501) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(-0.10903407) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5058274) q[0];
sx q[0];
rz(-0.286239) q[0];
sx q[0];
rz(0.92606996) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26308665) q[2];
sx q[2];
rz(-1.7619942) q[2];
sx q[2];
rz(-2.4441602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1034531) q[1];
sx q[1];
rz(-2.1557169) q[1];
sx q[1];
rz(1.000688) q[1];
rz(-pi) q[2];
rz(1.3974959) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(1.3483931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2382425) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.2634574) q[2];
rz(-2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.3431312) q[0];
rz(0.24761565) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-2.6599191) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7307229) q[0];
sx q[0];
rz(-2.0559089) q[0];
sx q[0];
rz(2.7184125) q[0];
rz(-pi) q[1];
rz(-2.1842314) q[2];
sx q[2];
rz(-1.295225) q[2];
sx q[2];
rz(-0.36188175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3232705) q[1];
sx q[1];
rz(-0.96410492) q[1];
sx q[1];
rz(0.89241772) q[1];
x q[2];
rz(-0.89110156) q[3];
sx q[3];
rz(-0.84935969) q[3];
sx q[3];
rz(0.47618714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3383011) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(-2.6417007) q[2];
rz(-2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-0.552185) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(-1.2447371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(-3.1367338) q[0];
rz(-pi) q[1];
rz(-2.142698) q[2];
sx q[2];
rz(-0.23795393) q[2];
sx q[2];
rz(-1.1513125) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4501805) q[1];
sx q[1];
rz(-1.0346864) q[1];
sx q[1];
rz(-2.8944573) q[1];
x q[2];
rz(2.7987715) q[3];
sx q[3];
rz(-2.5431513) q[3];
sx q[3];
rz(3.0330021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(1.9909031) q[2];
rz(-1.4771279) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(2.6586444) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0634336) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(0.081469014) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.6385471) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1254079) q[0];
sx q[0];
rz(-0.59690079) q[0];
sx q[0];
rz(-1.3422658) q[0];
x q[1];
rz(0.083085255) q[2];
sx q[2];
rz(-1.462888) q[2];
sx q[2];
rz(1.8748869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11322184) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(-0.10722864) q[1];
x q[2];
rz(1.8893858) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(0.82061758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5082671) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(-1.903669) q[2];
rz(2.0189019) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(-2.5807014) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(-2.3430603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0621588) q[0];
sx q[0];
rz(-1.845813) q[0];
sx q[0];
rz(-2.9955203) q[0];
x q[1];
rz(-0.59136765) q[2];
sx q[2];
rz(-0.57833507) q[2];
sx q[2];
rz(-0.27507281) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27998566) q[1];
sx q[1];
rz(-2.1415188) q[1];
sx q[1];
rz(-0.062203783) q[1];
rz(-0.800662) q[3];
sx q[3];
rz(-2.0293529) q[3];
sx q[3];
rz(0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325608) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(0.60638705) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(2.6775449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0021792) q[0];
sx q[0];
rz(-2.3521949) q[0];
sx q[0];
rz(-0.95285691) q[0];
x q[1];
rz(1.3299994) q[2];
sx q[2];
rz(-1.5511998) q[2];
sx q[2];
rz(0.400825) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5210919) q[1];
sx q[1];
rz(-1.2695841) q[1];
sx q[1];
rz(2.0786933) q[1];
x q[2];
rz(0.23630996) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(2.7304756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93418926) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(0.25804538) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946452) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(2.7602957) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.7000748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521116) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(-1.585929) q[0];
x q[1];
rz(0.61261119) q[2];
sx q[2];
rz(-1.6022041) q[2];
sx q[2];
rz(-0.14234662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.085233363) q[1];
sx q[1];
rz(-1.7006405) q[1];
sx q[1];
rz(-1.4443881) q[1];
rz(1.5278682) q[3];
sx q[3];
rz(-1.2371484) q[3];
sx q[3];
rz(-1.7121401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1429446) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(-2.9150035) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-2.3118238) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3806234) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(-2.1549966) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2664468) q[0];
sx q[0];
rz(-2.6212924) q[0];
sx q[0];
rz(-0.32169028) q[0];
x q[1];
rz(0.53194745) q[2];
sx q[2];
rz(-2.1956458) q[2];
sx q[2];
rz(2.9192386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.36531891) q[1];
sx q[1];
rz(-0.96818189) q[1];
sx q[1];
rz(-1.3990632) q[1];
rz(-pi) q[2];
rz(-1.1986198) q[3];
sx q[3];
rz(-0.44789207) q[3];
sx q[3];
rz(0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(-0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(-1.4051399) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(1.4155037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6575359) q[0];
sx q[0];
rz(-0.41175479) q[0];
sx q[0];
rz(-0.62396892) q[0];
rz(-pi) q[1];
rz(2.892898) q[2];
sx q[2];
rz(-1.2173614) q[2];
sx q[2];
rz(-2.8949708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43511697) q[1];
sx q[1];
rz(-0.26285989) q[1];
sx q[1];
rz(-2.190522) q[1];
rz(-pi) q[2];
rz(1.2763001) q[3];
sx q[3];
rz(-0.53139549) q[3];
sx q[3];
rz(-2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(-0.96578807) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(1.3901426) q[2];
sx q[2];
rz(-1.3144819) q[2];
sx q[2];
rz(-1.1983295) q[2];
rz(2.9028805) q[3];
sx q[3];
rz(-0.56959) q[3];
sx q[3];
rz(0.068985229) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];