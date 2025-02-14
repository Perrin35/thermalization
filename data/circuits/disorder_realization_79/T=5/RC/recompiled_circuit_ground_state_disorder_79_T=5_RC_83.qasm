OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.76685846) q[0];
sx q[0];
rz(2.2280966) q[0];
sx q[0];
rz(10.812617) q[0];
rz(0.18985441) q[1];
sx q[1];
rz(-1.8228276) q[1];
sx q[1];
rz(0.28611046) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090534276) q[0];
sx q[0];
rz(-1.5468452) q[0];
sx q[0];
rz(-2.6811203) q[0];
rz(-pi) q[1];
rz(-1.453581) q[2];
sx q[2];
rz(-2.5273537) q[2];
sx q[2];
rz(0.87488824) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.58264153) q[1];
sx q[1];
rz(-2.1182334) q[1];
sx q[1];
rz(1.2657341) q[1];
x q[2];
rz(1.1299575) q[3];
sx q[3];
rz(-1.5209271) q[3];
sx q[3];
rz(0.7383371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.046595786) q[2];
sx q[2];
rz(-1.4904212) q[2];
sx q[2];
rz(2.4563834) q[2];
rz(-1.6738711) q[3];
sx q[3];
rz(-2.0592212) q[3];
sx q[3];
rz(-2.8732324) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14107038) q[0];
sx q[0];
rz(-1.594161) q[0];
sx q[0];
rz(0.94456124) q[0];
rz(-3.0682849) q[1];
sx q[1];
rz(-2.3088375) q[1];
sx q[1];
rz(1.2578957) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2680141) q[0];
sx q[0];
rz(-1.1167545) q[0];
sx q[0];
rz(0.94768967) q[0];
rz(2.804685) q[2];
sx q[2];
rz(-2.5818129) q[2];
sx q[2];
rz(2.155456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98698178) q[1];
sx q[1];
rz(-1.4263337) q[1];
sx q[1];
rz(-2.7491436) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6173094) q[3];
sx q[3];
rz(-1.5336972) q[3];
sx q[3];
rz(1.0853037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32626095) q[2];
sx q[2];
rz(-1.3911284) q[2];
sx q[2];
rz(2.1168671) q[2];
rz(-2.4552086) q[3];
sx q[3];
rz(-2.4095583) q[3];
sx q[3];
rz(-0.91275233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3150385) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(0.81845534) q[0];
rz(1.2089027) q[1];
sx q[1];
rz(-1.5070288) q[1];
sx q[1];
rz(2.2817629) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6085109) q[0];
sx q[0];
rz(-1.9010592) q[0];
sx q[0];
rz(0.94655605) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6099714) q[2];
sx q[2];
rz(-0.79382703) q[2];
sx q[2];
rz(-0.17228157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2261994) q[1];
sx q[1];
rz(-2.2186154) q[1];
sx q[1];
rz(-0.97572216) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1016261) q[3];
sx q[3];
rz(-0.2519603) q[3];
sx q[3];
rz(0.89706883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39765057) q[2];
sx q[2];
rz(-2.8068779) q[2];
sx q[2];
rz(-3.0261377) q[2];
rz(0.3913106) q[3];
sx q[3];
rz(-1.0262998) q[3];
sx q[3];
rz(-1.9976043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35599577) q[0];
sx q[0];
rz(-1.9556671) q[0];
sx q[0];
rz(-2.9010229) q[0];
rz(1.7754414) q[1];
sx q[1];
rz(-1.6484478) q[1];
sx q[1];
rz(1.8919401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8842953) q[0];
sx q[0];
rz(-2.4438071) q[0];
sx q[0];
rz(1.1691514) q[0];
rz(-1.0860822) q[2];
sx q[2];
rz(-0.67469413) q[2];
sx q[2];
rz(-2.7759068) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0084997) q[1];
sx q[1];
rz(-1.3357107) q[1];
sx q[1];
rz(2.4753768) q[1];
rz(-pi) q[2];
rz(-0.68393882) q[3];
sx q[3];
rz(-1.9762282) q[3];
sx q[3];
rz(-0.33567521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25632855) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(-1.7735927) q[2];
rz(-0.3041501) q[3];
sx q[3];
rz(-1.5658028) q[3];
sx q[3];
rz(-2.4450541) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91365415) q[0];
sx q[0];
rz(-0.21622394) q[0];
sx q[0];
rz(-2.6044593) q[0];
rz(2.3917603) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(-0.73712635) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110187) q[0];
sx q[0];
rz(-1.6665742) q[0];
sx q[0];
rz(-1.6349313) q[0];
rz(-3.0546435) q[2];
sx q[2];
rz(-1.7450404) q[2];
sx q[2];
rz(-0.50424313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29038844) q[1];
sx q[1];
rz(-1.3121288) q[1];
sx q[1];
rz(1.7028536) q[1];
x q[2];
rz(0.26925663) q[3];
sx q[3];
rz(-2.807121) q[3];
sx q[3];
rz(2.0967029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4445112) q[2];
sx q[2];
rz(-0.35327521) q[2];
sx q[2];
rz(0.14661655) q[2];
rz(1.692449) q[3];
sx q[3];
rz(-1.2796947) q[3];
sx q[3];
rz(2.6176738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61831063) q[0];
sx q[0];
rz(-1.016541) q[0];
sx q[0];
rz(-1.7200394) q[0];
rz(-1.8824185) q[1];
sx q[1];
rz(-0.6764532) q[1];
sx q[1];
rz(-0.4037942) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97966563) q[0];
sx q[0];
rz(-2.296519) q[0];
sx q[0];
rz(-0.37958522) q[0];
x q[1];
rz(0.25713276) q[2];
sx q[2];
rz(-1.3916573) q[2];
sx q[2];
rz(2.9828097) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67298573) q[1];
sx q[1];
rz(-2.1588209) q[1];
sx q[1];
rz(0.67317669) q[1];
rz(-2.4533692) q[3];
sx q[3];
rz(-1.7587708) q[3];
sx q[3];
rz(1.926762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8518565) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(0.16656052) q[2];
rz(1.6377431) q[3];
sx q[3];
rz(-2.6681191) q[3];
sx q[3];
rz(-3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.903776) q[0];
sx q[0];
rz(-0.38919583) q[0];
sx q[0];
rz(-2.26407) q[0];
rz(0.96218836) q[1];
sx q[1];
rz(-2.018237) q[1];
sx q[1];
rz(0.35433623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0577045) q[0];
sx q[0];
rz(-1.5530927) q[0];
sx q[0];
rz(0.23772765) q[0];
rz(-1.7155859) q[2];
sx q[2];
rz(-1.5939404) q[2];
sx q[2];
rz(1.0095846) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6207032) q[1];
sx q[1];
rz(-2.1701309) q[1];
sx q[1];
rz(0.66461892) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6316492) q[3];
sx q[3];
rz(-1.6367607) q[3];
sx q[3];
rz(0.13291453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3366036) q[2];
sx q[2];
rz(-2.4038834) q[2];
sx q[2];
rz(2.0965516) q[2];
rz(-0.40237829) q[3];
sx q[3];
rz(-1.1674403) q[3];
sx q[3];
rz(-1.6654525) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104915) q[0];
sx q[0];
rz(-3.0248108) q[0];
sx q[0];
rz(-2.2905599) q[0];
rz(-2.7878413) q[1];
sx q[1];
rz(-1.7314792) q[1];
sx q[1];
rz(2.2869349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8847722) q[0];
sx q[0];
rz(-1.8516292) q[0];
sx q[0];
rz(-1.3815085) q[0];
rz(-pi) q[1];
rz(2.607367) q[2];
sx q[2];
rz(-2.8534128) q[2];
sx q[2];
rz(-2.2792918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9766337) q[1];
sx q[1];
rz(-1.7384496) q[1];
sx q[1];
rz(0.58604764) q[1];
rz(-pi) q[2];
rz(0.1885957) q[3];
sx q[3];
rz(-0.40939399) q[3];
sx q[3];
rz(0.90367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9510368) q[2];
sx q[2];
rz(-2.4943116) q[2];
sx q[2];
rz(2.5212042) q[2];
rz(0.22647151) q[3];
sx q[3];
rz(-1.7445931) q[3];
sx q[3];
rz(2.5605719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.9613551) q[0];
sx q[0];
rz(-1.9397475) q[0];
sx q[0];
rz(3.1264547) q[0];
rz(-1.0221647) q[1];
sx q[1];
rz(-2.731555) q[1];
sx q[1];
rz(-1.2000363) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2101636) q[0];
sx q[0];
rz(-2.4151122) q[0];
sx q[0];
rz(1.3036215) q[0];
rz(-pi) q[1];
rz(1.6300417) q[2];
sx q[2];
rz(-2.1631835) q[2];
sx q[2];
rz(2.5424344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85644279) q[1];
sx q[1];
rz(-1.0593482) q[1];
sx q[1];
rz(0.034354547) q[1];
x q[2];
rz(-2.3288378) q[3];
sx q[3];
rz(-0.65633196) q[3];
sx q[3];
rz(-1.1876653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6156561) q[2];
sx q[2];
rz(-2.2286712) q[2];
sx q[2];
rz(-2.9021662) q[2];
rz(-3.0827403) q[3];
sx q[3];
rz(-2.3299496) q[3];
sx q[3];
rz(2.8656901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1926159) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(-0.97491997) q[0];
rz(-2.4661567) q[1];
sx q[1];
rz(-2.2223739) q[1];
sx q[1];
rz(-2.1598037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73327195) q[0];
sx q[0];
rz(-1.3106723) q[0];
sx q[0];
rz(-0.97113804) q[0];
rz(-0.15530972) q[2];
sx q[2];
rz(-2.7255645) q[2];
sx q[2];
rz(1.1890026) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2416545) q[1];
sx q[1];
rz(-1.3119427) q[1];
sx q[1];
rz(-2.3090825) q[1];
rz(-pi) q[2];
rz(1.6493787) q[3];
sx q[3];
rz(-0.44789568) q[3];
sx q[3];
rz(3.0495897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5741817) q[2];
sx q[2];
rz(-2.1705706) q[2];
sx q[2];
rz(1.2249472) q[2];
rz(2.775906) q[3];
sx q[3];
rz(-0.94958011) q[3];
sx q[3];
rz(-1.1398116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097261978) q[0];
sx q[0];
rz(-1.2428357) q[0];
sx q[0];
rz(-1.6999929) q[0];
rz(-3.0615831) q[1];
sx q[1];
rz(-1.3792104) q[1];
sx q[1];
rz(3.1389799) q[1];
rz(0.44543191) q[2];
sx q[2];
rz(-1.6966698) q[2];
sx q[2];
rz(2.0296283) q[2];
rz(-3.0946685) q[3];
sx q[3];
rz(-0.45294807) q[3];
sx q[3];
rz(3.1035085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
