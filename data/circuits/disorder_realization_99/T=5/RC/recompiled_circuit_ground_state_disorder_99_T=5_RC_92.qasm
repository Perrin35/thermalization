OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(4.0532589) q[0];
sx q[0];
rz(3.4687923) q[0];
sx q[0];
rz(11.240525) q[0];
rz(1.3658547) q[1];
sx q[1];
rz(-0.39443016) q[1];
sx q[1];
rz(-2.3529513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74371007) q[0];
sx q[0];
rz(-1.2628444) q[0];
sx q[0];
rz(0.90713199) q[0];
rz(-0.6782669) q[2];
sx q[2];
rz(-0.94647289) q[2];
sx q[2];
rz(-2.2140142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0732735) q[1];
sx q[1];
rz(-1.8536659) q[1];
sx q[1];
rz(1.3702931) q[1];
x q[2];
rz(1.7472505) q[3];
sx q[3];
rz(-0.6354699) q[3];
sx q[3];
rz(0.42144767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6476562) q[2];
sx q[2];
rz(-1.6948573) q[2];
sx q[2];
rz(2.9392865) q[2];
rz(-0.10937396) q[3];
sx q[3];
rz(-0.60905639) q[3];
sx q[3];
rz(-3.0561395) q[3];
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
rz(-1.0626471) q[0];
sx q[0];
rz(-0.92342347) q[0];
sx q[0];
rz(-1.1751291) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-2.1390095) q[1];
sx q[1];
rz(-2.7820803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26665154) q[0];
sx q[0];
rz(-1.734899) q[0];
sx q[0];
rz(0.070673857) q[0];
rz(-0.55745947) q[2];
sx q[2];
rz(-1.2458846) q[2];
sx q[2];
rz(0.2211472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4103315) q[1];
sx q[1];
rz(-0.71440694) q[1];
sx q[1];
rz(2.6055715) q[1];
x q[2];
rz(2.3399598) q[3];
sx q[3];
rz(-2.3838861) q[3];
sx q[3];
rz(3.0996292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8466865) q[2];
sx q[2];
rz(-1.0470231) q[2];
sx q[2];
rz(-0.27493757) q[2];
rz(1.8168195) q[3];
sx q[3];
rz(-1.6144269) q[3];
sx q[3];
rz(1.8435318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0637958) q[0];
sx q[0];
rz(-0.15811385) q[0];
sx q[0];
rz(2.7445444) q[0];
rz(2.5322757) q[1];
sx q[1];
rz(-1.807223) q[1];
sx q[1];
rz(1.5416001) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7924431) q[0];
sx q[0];
rz(-0.31013784) q[0];
sx q[0];
rz(1.22195) q[0];
rz(1.6360388) q[2];
sx q[2];
rz(-1.5582623) q[2];
sx q[2];
rz(-1.5312851) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0761281) q[1];
sx q[1];
rz(-2.4148126) q[1];
sx q[1];
rz(1.5962334) q[1];
rz(-2.479781) q[3];
sx q[3];
rz(-2.1490009) q[3];
sx q[3];
rz(0.44826642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2043173) q[2];
sx q[2];
rz(-1.5466362) q[2];
sx q[2];
rz(-2.9079962) q[2];
rz(-1.2967845) q[3];
sx q[3];
rz(-1.1917043) q[3];
sx q[3];
rz(-0.72062033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5066756) q[0];
sx q[0];
rz(-1.1577865) q[0];
sx q[0];
rz(-1.7764212) q[0];
rz(1.8695976) q[1];
sx q[1];
rz(-1.1993661) q[1];
sx q[1];
rz(-1.2619527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.711447) q[0];
sx q[0];
rz(-1.0017348) q[0];
sx q[0];
rz(1.7993159) q[0];
rz(-pi) q[1];
rz(1.8795525) q[2];
sx q[2];
rz(-1.595904) q[2];
sx q[2];
rz(2.0659735) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8183444) q[1];
sx q[1];
rz(-1.7533315) q[1];
sx q[1];
rz(-1.6365461) q[1];
rz(-2.9449894) q[3];
sx q[3];
rz(-0.8340618) q[3];
sx q[3];
rz(1.3722668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9103553) q[2];
sx q[2];
rz(-2.1114025) q[2];
sx q[2];
rz(0.78872284) q[2];
rz(1.8606868) q[3];
sx q[3];
rz(-1.7773881) q[3];
sx q[3];
rz(1.0219215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3874409) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(-0.66106853) q[0];
rz(-1.1152274) q[1];
sx q[1];
rz(-1.3122357) q[1];
sx q[1];
rz(-2.1818395) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5710053) q[0];
sx q[0];
rz(-0.59413213) q[0];
sx q[0];
rz(-2.6460365) q[0];
x q[1];
rz(1.8371253) q[2];
sx q[2];
rz(-2.5779471) q[2];
sx q[2];
rz(0.061390419) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58127357) q[1];
sx q[1];
rz(-1.4380102) q[1];
sx q[1];
rz(-1.2788692) q[1];
x q[2];
rz(0.19790217) q[3];
sx q[3];
rz(-2.2060135) q[3];
sx q[3];
rz(-2.9445348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5614718) q[2];
sx q[2];
rz(-0.46130195) q[2];
sx q[2];
rz(-1.0599773) q[2];
rz(-0.75016841) q[3];
sx q[3];
rz(-1.7629905) q[3];
sx q[3];
rz(1.9157971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2724514) q[0];
sx q[0];
rz(-1.5831818) q[0];
sx q[0];
rz(-2.9172752) q[0];
rz(-1.51651) q[1];
sx q[1];
rz(-2.5254011) q[1];
sx q[1];
rz(0.075627653) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4556548) q[0];
sx q[0];
rz(-1.1854396) q[0];
sx q[0];
rz(2.1598201) q[0];
rz(-pi) q[1];
rz(-1.2709898) q[2];
sx q[2];
rz(-0.3170911) q[2];
sx q[2];
rz(0.89978774) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.87958) q[1];
sx q[1];
rz(-2.0063489) q[1];
sx q[1];
rz(1.0578937) q[1];
x q[2];
rz(-0.84556112) q[3];
sx q[3];
rz(-1.9048094) q[3];
sx q[3];
rz(-0.090868252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.533796) q[2];
sx q[2];
rz(-2.3888612) q[2];
sx q[2];
rz(-3.094574) q[2];
rz(-0.67695824) q[3];
sx q[3];
rz(-1.2983863) q[3];
sx q[3];
rz(-3.1162139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15636477) q[0];
sx q[0];
rz(-1.6522836) q[0];
sx q[0];
rz(1.4439247) q[0];
rz(-2.2185183) q[1];
sx q[1];
rz(-1.0052899) q[1];
sx q[1];
rz(0.18327555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8240487) q[0];
sx q[0];
rz(-2.5185761) q[0];
sx q[0];
rz(0.43806324) q[0];
x q[1];
rz(-0.44281339) q[2];
sx q[2];
rz(-1.9573136) q[2];
sx q[2];
rz(-0.3511951) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35872719) q[1];
sx q[1];
rz(-2.065383) q[1];
sx q[1];
rz(-2.7703157) q[1];
x q[2];
rz(0.99162001) q[3];
sx q[3];
rz(-1.3282816) q[3];
sx q[3];
rz(0.33580599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.034059374) q[2];
sx q[2];
rz(-1.4564161) q[2];
sx q[2];
rz(-3.122186) q[2];
rz(0.16128811) q[3];
sx q[3];
rz(-2.1262157) q[3];
sx q[3];
rz(-1.2279145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0374544) q[0];
sx q[0];
rz(-2.7327974) q[0];
sx q[0];
rz(3.0492875) q[0];
rz(-1.9654988) q[1];
sx q[1];
rz(-0.25830019) q[1];
sx q[1];
rz(1.4091122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9285959) q[0];
sx q[0];
rz(-0.11830506) q[0];
sx q[0];
rz(-2.0483093) q[0];
rz(-pi) q[1];
rz(0.95887948) q[2];
sx q[2];
rz(-1.1772708) q[2];
sx q[2];
rz(-2.3589691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8586052) q[1];
sx q[1];
rz(-1.3027362) q[1];
sx q[1];
rz(0.40045935) q[1];
x q[2];
rz(0.73958379) q[3];
sx q[3];
rz(-0.7639262) q[3];
sx q[3];
rz(-0.52810625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1000503) q[2];
sx q[2];
rz(-1.5644194) q[2];
sx q[2];
rz(1.9680295) q[2];
rz(-2.7581577) q[3];
sx q[3];
rz(-1.8337367) q[3];
sx q[3];
rz(1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0122871) q[0];
sx q[0];
rz(-1.6467935) q[0];
sx q[0];
rz(-2.9360085) q[0];
rz(2.3221817) q[1];
sx q[1];
rz(-1.9267547) q[1];
sx q[1];
rz(-0.49088556) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6554479) q[0];
sx q[0];
rz(-1.1974338) q[0];
sx q[0];
rz(2.4400737) q[0];
rz(-0.41342469) q[2];
sx q[2];
rz(-0.67081645) q[2];
sx q[2];
rz(-0.49395032) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5583239) q[1];
sx q[1];
rz(-1.0029536) q[1];
sx q[1];
rz(0.056599157) q[1];
x q[2];
rz(1.396046) q[3];
sx q[3];
rz(-0.51769231) q[3];
sx q[3];
rz(2.2153436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72863355) q[2];
sx q[2];
rz(-2.1042991) q[2];
sx q[2];
rz(-1.884985) q[2];
rz(2.7058153) q[3];
sx q[3];
rz(-1.1083009) q[3];
sx q[3];
rz(-0.88424879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4202145) q[0];
sx q[0];
rz(-0.90737897) q[0];
sx q[0];
rz(-2.9086928) q[0];
rz(0.13889343) q[1];
sx q[1];
rz(-1.1537617) q[1];
sx q[1];
rz(-2.6331666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8247037) q[0];
sx q[0];
rz(-1.034786) q[0];
sx q[0];
rz(2.8404923) q[0];
rz(-pi) q[1];
rz(2.0783689) q[2];
sx q[2];
rz(-1.3295577) q[2];
sx q[2];
rz(-2.9335748) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.057699) q[1];
sx q[1];
rz(-0.13421276) q[1];
sx q[1];
rz(1.6300409) q[1];
rz(0.87457652) q[3];
sx q[3];
rz(-1.1406058) q[3];
sx q[3];
rz(-0.065635292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7306708) q[2];
sx q[2];
rz(-1.6795009) q[2];
sx q[2];
rz(-0.81653583) q[2];
rz(-2.7517892) q[3];
sx q[3];
rz(-1.0023508) q[3];
sx q[3];
rz(-1.7635112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.447406) q[0];
sx q[0];
rz(-0.92300713) q[0];
sx q[0];
rz(-0.22966455) q[0];
rz(-1.6387088) q[1];
sx q[1];
rz(-1.8062183) q[1];
sx q[1];
rz(0.91581215) q[1];
rz(1.7433132) q[2];
sx q[2];
rz(-1.1926706) q[2];
sx q[2];
rz(-0.43422912) q[2];
rz(1.7671723) q[3];
sx q[3];
rz(-0.71153258) q[3];
sx q[3];
rz(1.1179433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
