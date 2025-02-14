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
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(-2.9297096) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(2.6091726) q[1];
sx q[1];
rz(6.6611023) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.657562) q[0];
sx q[0];
rz(-0.0095417984) q[0];
sx q[0];
rz(1.5242759) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67899668) q[2];
sx q[2];
rz(-1.8583991) q[2];
sx q[2];
rz(-1.2305413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18641414) q[1];
sx q[1];
rz(-0.70611533) q[1];
sx q[1];
rz(-0.43118711) q[1];
x q[2];
rz(-0.8182977) q[3];
sx q[3];
rz(-2.1967271) q[3];
sx q[3];
rz(-1.675954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2892896) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(0.079785384) q[2];
rz(2.1763109) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(0.7307581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47736436) q[0];
sx q[0];
rz(-2.1381162) q[0];
sx q[0];
rz(-0.088951237) q[0];
rz(-1.3720007) q[1];
sx q[1];
rz(-1.7141432) q[1];
sx q[1];
rz(-1.1044097) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0029582214) q[0];
sx q[0];
rz(-1.3413652) q[0];
sx q[0];
rz(2.5087439) q[0];
rz(0.78271336) q[2];
sx q[2];
rz(-1.209895) q[2];
sx q[2];
rz(0.95907839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.96530487) q[1];
sx q[1];
rz(-1.7200279) q[1];
sx q[1];
rz(2.01841) q[1];
rz(-0.16544754) q[3];
sx q[3];
rz(-2.5962127) q[3];
sx q[3];
rz(2.489413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.264512) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(-1.6748927) q[2];
rz(0.029021164) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(-0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324209) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(1.4311283) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(-2.7412282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6872266) q[0];
sx q[0];
rz(-0.80601826) q[0];
sx q[0];
rz(2.5044652) q[0];
rz(-2.4552271) q[2];
sx q[2];
rz(-1.371742) q[2];
sx q[2];
rz(-0.65716568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0770994) q[1];
sx q[1];
rz(-1.1005741) q[1];
sx q[1];
rz(-1.749586) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8767897) q[3];
sx q[3];
rz(-1.8974432) q[3];
sx q[3];
rz(-1.6040925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.532734) q[2];
sx q[2];
rz(0.45480967) q[2];
rz(-1.2416035) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(-0.72030592) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860745) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(-3.1410826) q[0];
rz(-0.60091248) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(-0.14437637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7847608) q[0];
sx q[0];
rz(-0.5934754) q[0];
sx q[0];
rz(-1.9463825) q[0];
rz(-pi) q[1];
rz(0.25344376) q[2];
sx q[2];
rz(-1.4254693) q[2];
sx q[2];
rz(0.0090473024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3444654) q[1];
sx q[1];
rz(-0.68736156) q[1];
sx q[1];
rz(-0.3631773) q[1];
x q[2];
rz(-0.5471121) q[3];
sx q[3];
rz(-1.4792076) q[3];
sx q[3];
rz(-1.6325337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.942975) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(-0.96251881) q[2];
rz(-1.5396384) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(2.8008154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2365504) q[0];
sx q[0];
rz(-1.6860697) q[0];
sx q[0];
rz(1.1849674) q[0];
rz(-0.22625893) q[1];
sx q[1];
rz(-2.2626651) q[1];
sx q[1];
rz(-2.102898) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13030355) q[0];
sx q[0];
rz(-1.652392) q[0];
sx q[0];
rz(-0.58873436) q[0];
rz(-pi) q[1];
rz(2.2009497) q[2];
sx q[2];
rz(-2.2046208) q[2];
sx q[2];
rz(1.188082) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24327899) q[1];
sx q[1];
rz(-0.89363511) q[1];
sx q[1];
rz(-2.2449298) q[1];
rz(-pi) q[2];
rz(-3.076055) q[3];
sx q[3];
rz(-2.0100694) q[3];
sx q[3];
rz(0.36079839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.431939) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(2.8254438) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-2.0114055) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6370711) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(2.0380518) q[0];
rz(2.6612813) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(2.5441817) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3854073) q[0];
sx q[0];
rz(-0.28235897) q[0];
sx q[0];
rz(0.83448164) q[0];
rz(-pi) q[1];
rz(2.3961847) q[2];
sx q[2];
rz(-1.3920203) q[2];
sx q[2];
rz(-2.488236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.589041) q[1];
sx q[1];
rz(-1.3197761) q[1];
sx q[1];
rz(1.7415206) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3378998) q[3];
sx q[3];
rz(-2.7613556) q[3];
sx q[3];
rz(-1.5986795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(-1.9988029) q[3];
sx q[3];
rz(-2.3484774) q[3];
sx q[3];
rz(-2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(3.1031188) q[0];
rz(3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(0.23385349) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9486987) q[0];
sx q[0];
rz(-1.7561551) q[0];
sx q[0];
rz(-2.9167487) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31189708) q[2];
sx q[2];
rz(-1.2205178) q[2];
sx q[2];
rz(0.52679449) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.13713) q[1];
sx q[1];
rz(-2.4871768) q[1];
sx q[1];
rz(-0.14738247) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7247808) q[3];
sx q[3];
rz(-2.9997065) q[3];
sx q[3];
rz(2.8519423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48876277) q[2];
sx q[2];
rz(-1.5831999) q[2];
sx q[2];
rz(-0.59824198) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(-2.2940476) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4050196) q[0];
sx q[0];
rz(-2.2990655) q[0];
sx q[0];
rz(-2.8952428) q[0];
rz(1.2975289) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(-2.6447703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778567) q[0];
sx q[0];
rz(-1.9851159) q[0];
sx q[0];
rz(0.68092771) q[0];
x q[1];
rz(2.970817) q[2];
sx q[2];
rz(-2.429212) q[2];
sx q[2];
rz(0.26584372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5649181) q[1];
sx q[1];
rz(-2.0510625) q[1];
sx q[1];
rz(2.8757921) q[1];
x q[2];
rz(-1.2874576) q[3];
sx q[3];
rz(-0.50931286) q[3];
sx q[3];
rz(0.71374245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42314998) q[2];
sx q[2];
rz(-0.53360525) q[2];
sx q[2];
rz(-1.144484) q[2];
rz(1.9874969) q[3];
sx q[3];
rz(-1.7172979) q[3];
sx q[3];
rz(0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474739) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(-3.0754454) q[0];
rz(1.1478395) q[1];
sx q[1];
rz(-1.7639953) q[1];
sx q[1];
rz(0.60417169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19949958) q[0];
sx q[0];
rz(-1.8460994) q[0];
sx q[0];
rz(1.3574187) q[0];
rz(-2.6148817) q[2];
sx q[2];
rz(-1.4350495) q[2];
sx q[2];
rz(1.2413687) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79163247) q[1];
sx q[1];
rz(-0.59795982) q[1];
sx q[1];
rz(-2.7954742) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3202032) q[3];
sx q[3];
rz(-2.2940972) q[3];
sx q[3];
rz(-2.5625474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26312795) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(0.43295941) q[2];
rz(1.2290907) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(-1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94806725) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(1.8684335) q[0];
rz(-2.676447) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(2.926631) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8037655) q[0];
sx q[0];
rz(-1.318232) q[0];
sx q[0];
rz(2.4830677) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6647756) q[2];
sx q[2];
rz(-0.6419581) q[2];
sx q[2];
rz(-2.4017815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82483208) q[1];
sx q[1];
rz(-0.90710282) q[1];
sx q[1];
rz(0.65726991) q[1];
x q[2];
rz(-1.3137903) q[3];
sx q[3];
rz(-2.1177835) q[3];
sx q[3];
rz(1.6099324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2994069) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(-0.95883933) q[2];
rz(-1.1622608) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(-1.4452665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5398298) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(-0.70855793) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(0.56253994) q[2];
sx q[2];
rz(-1.8481135) q[2];
sx q[2];
rz(-2.6279966) q[2];
rz(-1.3832573) q[3];
sx q[3];
rz(-2.2839727) q[3];
sx q[3];
rz(-1.4061389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
