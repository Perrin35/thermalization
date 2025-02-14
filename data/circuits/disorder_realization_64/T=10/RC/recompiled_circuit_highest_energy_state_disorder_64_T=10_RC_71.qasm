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
rz(2.0848701) q[0];
sx q[0];
rz(4.4216006) q[0];
sx q[0];
rz(13.292424) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(2.1093624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5846378) q[0];
sx q[0];
rz(-2.224486) q[0];
sx q[0];
rz(-1.4435221) q[0];
x q[1];
rz(-0.039041877) q[2];
sx q[2];
rz(-1.0532951) q[2];
sx q[2];
rz(1.7797178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91019708) q[1];
sx q[1];
rz(-0.33581844) q[1];
sx q[1];
rz(0.37792716) q[1];
x q[2];
rz(-0.49555669) q[3];
sx q[3];
rz(-0.86641387) q[3];
sx q[3];
rz(0.92951194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7008179) q[2];
sx q[2];
rz(-2.5371607) q[2];
sx q[2];
rz(-0.087892858) q[2];
rz(1.6739738) q[3];
sx q[3];
rz(-1.2940977) q[3];
sx q[3];
rz(0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217594) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(-1.649296) q[0];
rz(0.78798931) q[1];
sx q[1];
rz(-1.9580656) q[1];
sx q[1];
rz(-2.1764887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3639587) q[0];
sx q[0];
rz(-2.2978362) q[0];
sx q[0];
rz(1.4476089) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4282706) q[2];
sx q[2];
rz(-0.83907849) q[2];
sx q[2];
rz(0.97396641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31307855) q[1];
sx q[1];
rz(-2.0277436) q[1];
sx q[1];
rz(-2.9345153) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4697462) q[3];
sx q[3];
rz(-1.489991) q[3];
sx q[3];
rz(-1.684379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1284156) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(1.817912) q[2];
rz(-2.2549818) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(-0.19865856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-2.1710945) q[0];
sx q[0];
rz(0.6955198) q[0];
rz(0.8356525) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(2.4998891) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4280689) q[0];
sx q[0];
rz(-1.4260573) q[0];
sx q[0];
rz(-2.3109396) q[0];
x q[1];
rz(-1.4029866) q[2];
sx q[2];
rz(-2.2750759) q[2];
sx q[2];
rz(1.86098) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0271488) q[1];
sx q[1];
rz(-0.58747298) q[1];
sx q[1];
rz(2.1617166) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99567497) q[3];
sx q[3];
rz(-2.174587) q[3];
sx q[3];
rz(0.010494516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(1.3345435) q[2];
rz(1.7415907) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(2.002772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982933) q[0];
sx q[0];
rz(-2.2704953) q[0];
sx q[0];
rz(-2.368108) q[0];
rz(3.0553014) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(-0.48826826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7542758) q[0];
sx q[0];
rz(-1.5026706) q[0];
sx q[0];
rz(-0.46275567) q[0];
x q[1];
rz(-0.62520844) q[2];
sx q[2];
rz(-2.3397988) q[2];
sx q[2];
rz(-1.2259353) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99566702) q[1];
sx q[1];
rz(-2.2856376) q[1];
sx q[1];
rz(-2.8413474) q[1];
x q[2];
rz(-2.7375324) q[3];
sx q[3];
rz(-1.9480138) q[3];
sx q[3];
rz(2.4665143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46828541) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(2.9735273) q[2];
rz(0.42260653) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(1.9510829) q[0];
rz(-0.52781421) q[1];
sx q[1];
rz(-2.0595136) q[1];
sx q[1];
rz(0.0091008069) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4531698) q[0];
sx q[0];
rz(-0.97081447) q[0];
sx q[0];
rz(-2.1196516) q[0];
rz(-pi) q[1];
rz(2.3528655) q[2];
sx q[2];
rz(-0.085914748) q[2];
sx q[2];
rz(-0.14363657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.121351) q[1];
sx q[1];
rz(-2.1727736) q[1];
sx q[1];
rz(-0.51687981) q[1];
rz(-pi) q[2];
rz(2.1808257) q[3];
sx q[3];
rz(-0.72566635) q[3];
sx q[3];
rz(2.8270367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8268383) q[2];
sx q[2];
rz(-0.49586168) q[2];
sx q[2];
rz(-2.4763988) q[2];
rz(-1.0264171) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219532) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(0.39805472) q[0];
rz(-1.6710501) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(0.7801396) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0288643) q[0];
sx q[0];
rz(-2.5327589) q[0];
sx q[0];
rz(2.2236149) q[0];
rz(-2.4498135) q[2];
sx q[2];
rz(-0.31111141) q[2];
sx q[2];
rz(1.9115314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8079905) q[1];
sx q[1];
rz(-1.2871075) q[1];
sx q[1];
rz(2.4391443) q[1];
rz(-pi) q[2];
rz(1.1885181) q[3];
sx q[3];
rz(-2.4197289) q[3];
sx q[3];
rz(0.75306276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.857343) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(-0.16172376) q[2];
rz(-0.11180793) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(0.73981729) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(0.69851843) q[0];
rz(-2.0493719) q[1];
sx q[1];
rz(-2.3869546) q[1];
sx q[1];
rz(0.05365595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221401) q[0];
sx q[0];
rz(-0.86556731) q[0];
sx q[0];
rz(0.34366519) q[0];
x q[1];
rz(1.1321819) q[2];
sx q[2];
rz(-1.2234067) q[2];
sx q[2];
rz(1.5127104) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0567017) q[1];
sx q[1];
rz(-0.47353256) q[1];
sx q[1];
rz(3.0242306) q[1];
rz(-2.9916006) q[3];
sx q[3];
rz(-0.92954554) q[3];
sx q[3];
rz(0.36442001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1412389) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(-2.9290283) q[2];
rz(2.8138748) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(-1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89606365) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(2.8118706) q[0];
rz(0.7715191) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(0.82569295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77149174) q[0];
sx q[0];
rz(-2.2326062) q[0];
sx q[0];
rz(-1.6812117) q[0];
x q[1];
rz(-2.1671381) q[2];
sx q[2];
rz(-1.505841) q[2];
sx q[2];
rz(-2.7686053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18620423) q[1];
sx q[1];
rz(-2.3729747) q[1];
sx q[1];
rz(0.15665084) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65123043) q[3];
sx q[3];
rz(-0.31933258) q[3];
sx q[3];
rz(-2.8791219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19994911) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(-1.6597718) q[2];
rz(3.0485349) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2411497) q[0];
sx q[0];
rz(-0.41541442) q[0];
sx q[0];
rz(0.39733091) q[0];
rz(0.98995248) q[1];
sx q[1];
rz(-1.7918469) q[1];
sx q[1];
rz(2.1053402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2275804) q[0];
sx q[0];
rz(-1.5591994) q[0];
sx q[0];
rz(0.40475233) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35276995) q[2];
sx q[2];
rz(-1.4401541) q[2];
sx q[2];
rz(-0.58707159) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9246722) q[1];
sx q[1];
rz(-2.9045871) q[1];
sx q[1];
rz(2.7525178) q[1];
rz(-pi) q[2];
rz(-1.051722) q[3];
sx q[3];
rz(-1.241893) q[3];
sx q[3];
rz(-0.84612209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1194666) q[2];
sx q[2];
rz(-0.72212044) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(-0.25257603) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13755688) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(0.23605119) q[0];
rz(0.099418489) q[1];
sx q[1];
rz(-2.3853018) q[1];
sx q[1];
rz(1.8831467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041708) q[0];
sx q[0];
rz(-1.4930471) q[0];
sx q[0];
rz(2.581152) q[0];
rz(-1.0706606) q[2];
sx q[2];
rz(-0.96958435) q[2];
sx q[2];
rz(-0.76155969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.395605) q[1];
sx q[1];
rz(-1.6853991) q[1];
sx q[1];
rz(-2.458861) q[1];
rz(-2.7332337) q[3];
sx q[3];
rz(-1.1757188) q[3];
sx q[3];
rz(1.4407106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9749757) q[2];
sx q[2];
rz(-2.6585343) q[2];
sx q[2];
rz(1.0413337) q[2];
rz(2.5552022) q[3];
sx q[3];
rz(-2.6643463) q[3];
sx q[3];
rz(0.81048107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36622421) q[0];
sx q[0];
rz(-2.0851705) q[0];
sx q[0];
rz(2.0023517) q[0];
rz(2.1495023) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(0.61983776) q[2];
sx q[2];
rz(-0.38417338) q[2];
sx q[2];
rz(-1.2585121) q[2];
rz(0.6057779) q[3];
sx q[3];
rz(-0.68690261) q[3];
sx q[3];
rz(1.5920832) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
