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
rz(-1.8981847) q[0];
sx q[0];
rz(-1.564448) q[0];
sx q[0];
rz(2.2301883) q[0];
rz(-0.66943327) q[1];
sx q[1];
rz(-0.27639204) q[1];
sx q[1];
rz(-2.3120094) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4938717) q[0];
sx q[0];
rz(-0.52807489) q[0];
sx q[0];
rz(-2.2654669) q[0];
rz(-1.4777415) q[2];
sx q[2];
rz(-1.3246127) q[2];
sx q[2];
rz(-0.024178084) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1948283) q[1];
sx q[1];
rz(-0.84673893) q[1];
sx q[1];
rz(-1.1652693) q[1];
rz(0.43609332) q[3];
sx q[3];
rz(-1.7096409) q[3];
sx q[3];
rz(0.88263765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36898819) q[2];
sx q[2];
rz(-2.172894) q[2];
sx q[2];
rz(1.1645092) q[2];
rz(-2.4146967) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(-0.070778457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(0.29139274) q[0];
rz(2.5413051) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(-2.9765863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5801738) q[0];
sx q[0];
rz(-0.2720755) q[0];
sx q[0];
rz(1.3035357) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8012456) q[2];
sx q[2];
rz(-2.682095) q[2];
sx q[2];
rz(2.8682617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8540604) q[1];
sx q[1];
rz(-2.2496576) q[1];
sx q[1];
rz(-1.6599342) q[1];
x q[2];
rz(0.24834541) q[3];
sx q[3];
rz(-1.6237061) q[3];
sx q[3];
rz(-1.9236881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0642285) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(-2.9020818) q[2];
rz(0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(-0.094836205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25642446) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(1.5727795) q[0];
rz(2.7478711) q[1];
sx q[1];
rz(-0.55501333) q[1];
sx q[1];
rz(0.47600019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25903738) q[0];
sx q[0];
rz(-2.5406484) q[0];
sx q[0];
rz(-1.1968234) q[0];
rz(-pi) q[1];
rz(1.7577241) q[2];
sx q[2];
rz(-1.6003967) q[2];
sx q[2];
rz(1.6254978) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6913298) q[1];
sx q[1];
rz(-0.99585497) q[1];
sx q[1];
rz(-2.6460669) q[1];
x q[2];
rz(-1.2331486) q[3];
sx q[3];
rz(-1.8051355) q[3];
sx q[3];
rz(2.4837787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3289045) q[2];
sx q[2];
rz(-2.4220971) q[2];
sx q[2];
rz(2.1997814) q[2];
rz(1.4224667) q[3];
sx q[3];
rz(-2.7673281) q[3];
sx q[3];
rz(-2.4678738) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.5330676) q[0];
rz(2.7948921) q[1];
sx q[1];
rz(-2.0997212) q[1];
sx q[1];
rz(0.18992058) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3118211) q[0];
sx q[0];
rz(-1.30952) q[0];
sx q[0];
rz(-0.7423865) q[0];
rz(-2.0702576) q[2];
sx q[2];
rz(-0.53218666) q[2];
sx q[2];
rz(0.70390534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20224883) q[1];
sx q[1];
rz(-0.96572438) q[1];
sx q[1];
rz(3.0720965) q[1];
x q[2];
rz(-1.1433825) q[3];
sx q[3];
rz(-1.9460925) q[3];
sx q[3];
rz(-2.0209598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3612264) q[2];
sx q[2];
rz(-0.85005886) q[2];
sx q[2];
rz(-0.38193199) q[2];
rz(3.0319038) q[3];
sx q[3];
rz(-1.0332801) q[3];
sx q[3];
rz(0.1194574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-0.1031951) q[0];
sx q[0];
rz(-1.9586451) q[0];
sx q[0];
rz(-0.84684816) q[0];
rz(3.0021744) q[1];
sx q[1];
rz(-0.81652313) q[1];
sx q[1];
rz(2.3925508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3232305) q[0];
sx q[0];
rz(-1.4233755) q[0];
sx q[0];
rz(-0.69578232) q[0];
rz(-1.4290358) q[2];
sx q[2];
rz(-2.0516011) q[2];
sx q[2];
rz(2.9713809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74816675) q[1];
sx q[1];
rz(-1.6977377) q[1];
sx q[1];
rz(2.9462161) q[1];
rz(-pi) q[2];
rz(1.7827634) q[3];
sx q[3];
rz(-2.6304768) q[3];
sx q[3];
rz(1.2455452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2060812) q[2];
sx q[2];
rz(-1.8847909) q[2];
sx q[2];
rz(1.067266) q[2];
rz(2.4308128) q[3];
sx q[3];
rz(-1.4830517) q[3];
sx q[3];
rz(1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.34273219) q[0];
sx q[0];
rz(-1.4787759) q[0];
sx q[0];
rz(0.96858281) q[0];
rz(-0.69106483) q[1];
sx q[1];
rz(-1.3422809) q[1];
sx q[1];
rz(-1.8341281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027031487) q[0];
sx q[0];
rz(-1.3163211) q[0];
sx q[0];
rz(-0.64395321) q[0];
x q[1];
rz(2.5619823) q[2];
sx q[2];
rz(-2.7977849) q[2];
sx q[2];
rz(2.0909617) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4279514) q[1];
sx q[1];
rz(-2.7507857) q[1];
sx q[1];
rz(1.4321838) q[1];
rz(-2.9581466) q[3];
sx q[3];
rz(-1.8357092) q[3];
sx q[3];
rz(1.9062756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4751733) q[2];
sx q[2];
rz(-2.9066777) q[2];
sx q[2];
rz(1.8443745) q[2];
rz(1.746486) q[3];
sx q[3];
rz(-1.3336983) q[3];
sx q[3];
rz(1.6102128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.7995826) q[0];
sx q[0];
rz(-2.3260249) q[0];
sx q[0];
rz(-2.2014501) q[0];
rz(-0.6483342) q[1];
sx q[1];
rz(-1.9377361) q[1];
sx q[1];
rz(-2.8727093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20737442) q[0];
sx q[0];
rz(-1.3257799) q[0];
sx q[0];
rz(2.5859358) q[0];
x q[1];
rz(-0.70272311) q[2];
sx q[2];
rz(-1.3303192) q[2];
sx q[2];
rz(2.3287077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2317048) q[1];
sx q[1];
rz(-0.99867601) q[1];
sx q[1];
rz(2.8673299) q[1];
x q[2];
rz(-0.42952764) q[3];
sx q[3];
rz(-1.8801573) q[3];
sx q[3];
rz(-1.3573028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0073283422) q[2];
sx q[2];
rz(-1.2790044) q[2];
sx q[2];
rz(-2.4556665) q[2];
rz(-0.736233) q[3];
sx q[3];
rz(-2.1015621) q[3];
sx q[3];
rz(-2.9161016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.4418884) q[0];
sx q[0];
rz(-1.9206973) q[0];
sx q[0];
rz(0.7487444) q[0];
rz(-0.092983149) q[1];
sx q[1];
rz(-2.3817606) q[1];
sx q[1];
rz(1.071113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1919928) q[0];
sx q[0];
rz(-1.9849536) q[0];
sx q[0];
rz(-0.55121514) q[0];
rz(-1.1125715) q[2];
sx q[2];
rz(-2.1734218) q[2];
sx q[2];
rz(-0.17652179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6735355) q[1];
sx q[1];
rz(-2.2186465) q[1];
sx q[1];
rz(-2.7537529) q[1];
x q[2];
rz(-1.981826) q[3];
sx q[3];
rz(-1.5229406) q[3];
sx q[3];
rz(-1.4519297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6894655) q[2];
sx q[2];
rz(-0.63673821) q[2];
sx q[2];
rz(0.20991906) q[2];
rz(0.12957761) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(1.5642222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.121948) q[0];
sx q[0];
rz(-0.26109281) q[0];
sx q[0];
rz(-0.013539465) q[0];
rz(-2.6059222) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(0.013669107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9122152) q[0];
sx q[0];
rz(-1.3446864) q[0];
sx q[0];
rz(2.3588728) q[0];
rz(-2.2315908) q[2];
sx q[2];
rz(-2.5507413) q[2];
sx q[2];
rz(0.94205085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0422349) q[1];
sx q[1];
rz(-2.3546948) q[1];
sx q[1];
rz(2.6762647) q[1];
rz(0.39710703) q[3];
sx q[3];
rz(-0.79063168) q[3];
sx q[3];
rz(-1.1943898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51836625) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(-1.3631442) q[2];
rz(2.3220883) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(-2.4578186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0226456) q[0];
sx q[0];
rz(-2.2830257) q[0];
sx q[0];
rz(-1.488142) q[0];
rz(-0.34291521) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(2.3905579) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35151154) q[0];
sx q[0];
rz(-1.292391) q[0];
sx q[0];
rz(-0.06601575) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1079587) q[2];
sx q[2];
rz(-1.3909391) q[2];
sx q[2];
rz(3.0015869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0768834) q[1];
sx q[1];
rz(-2.5297477) q[1];
sx q[1];
rz(-0.20670273) q[1];
rz(-pi) q[2];
rz(1.4806101) q[3];
sx q[3];
rz(-1.6038365) q[3];
sx q[3];
rz(-2.9905013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47120961) q[2];
sx q[2];
rz(-2.1977916) q[2];
sx q[2];
rz(2.75441) q[2];
rz(0.47532982) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(0.79132426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8925856) q[0];
sx q[0];
rz(-0.96207608) q[0];
sx q[0];
rz(0.67208653) q[0];
rz(-0.36920209) q[1];
sx q[1];
rz(-2.3875356) q[1];
sx q[1];
rz(1.748132) q[1];
rz(-2.9291991) q[2];
sx q[2];
rz(-2.4133189) q[2];
sx q[2];
rz(2.6475785) q[2];
rz(1.9561097) q[3];
sx q[3];
rz(-0.29373071) q[3];
sx q[3];
rz(3.1391524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
