OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(0.84258643) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9561477) q[0];
sx q[0];
rz(-1.7831793) q[0];
sx q[0];
rz(0.74622112) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(-0.27054271) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7588501) q[1];
sx q[1];
rz(-1.9470125) q[1];
sx q[1];
rz(1.43169) q[1];
x q[2];
rz(1.4487212) q[3];
sx q[3];
rz(-0.97994084) q[3];
sx q[3];
rz(1.5199682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(-0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-0.32546145) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017529537) q[0];
sx q[0];
rz(-1.5520099) q[0];
sx q[0];
rz(-1.5879052) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94475586) q[2];
sx q[2];
rz(-1.2465887) q[2];
sx q[2];
rz(-0.39694436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5408389) q[1];
sx q[1];
rz(-2.3327017) q[1];
sx q[1];
rz(1.6879338) q[1];
x q[2];
rz(-0.18493821) q[3];
sx q[3];
rz(-1.3122845) q[3];
sx q[3];
rz(-1.0227433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(-2.3538891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31323355) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.5154243) q[0];
rz(0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-1.0916969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9804304) q[0];
sx q[0];
rz(-0.83320252) q[0];
sx q[0];
rz(-2.3479793) q[0];
rz(-0.91471471) q[2];
sx q[2];
rz(-1.2351742) q[2];
sx q[2];
rz(2.6732973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4746689) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(-2.073642) q[1];
rz(3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(-1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(0.88095218) q[2];
rz(-1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(-0.4367035) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-2.8312347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.045517) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(-1.7128574) q[0];
rz(-0.68508673) q[2];
sx q[2];
rz(-1.6685467) q[2];
sx q[2];
rz(-0.098066559) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0879285) q[1];
sx q[1];
rz(-1.2128608) q[1];
sx q[1];
rz(3.1217561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.035590812) q[3];
sx q[3];
rz(-1.4471874) q[3];
sx q[3];
rz(-1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(-3.0854026) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.189165) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(0.24965832) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0349883) q[0];
sx q[0];
rz(-0.37811324) q[0];
sx q[0];
rz(-2.5614221) q[0];
x q[1];
rz(-0.22360794) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(-0.86415926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89228499) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(2.409163) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3051885) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(0.30361787) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34981397) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0598037) q[0];
sx q[0];
rz(-0.44133082) q[0];
sx q[0];
rz(2.9102737) q[0];
rz(-pi) q[1];
rz(-0.50478023) q[2];
sx q[2];
rz(-2.3961888) q[2];
sx q[2];
rz(2.8842852) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41960934) q[1];
sx q[1];
rz(-1.0227385) q[1];
sx q[1];
rz(2.4649058) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3098573) q[3];
sx q[3];
rz(-2.76537) q[3];
sx q[3];
rz(-0.46686831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-0.84364676) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(-0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(0.86589083) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(-0.033989865) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2798529) q[0];
sx q[0];
rz(-1.4089157) q[0];
sx q[0];
rz(0.78761657) q[0];
x q[1];
rz(-2.6007973) q[2];
sx q[2];
rz(-1.0058837) q[2];
sx q[2];
rz(2.5164547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6778292) q[1];
sx q[1];
rz(-0.67968183) q[1];
sx q[1];
rz(1.6512647) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5266685) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(-1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6035446) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(-0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.4454909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60364265) q[0];
sx q[0];
rz(-1.1374439) q[0];
sx q[0];
rz(-3.135878) q[0];
x q[1];
rz(0.58480279) q[2];
sx q[2];
rz(-1.0324761) q[2];
sx q[2];
rz(-2.1540097) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6519424) q[1];
sx q[1];
rz(-1.4133269) q[1];
sx q[1];
rz(0.55649346) q[1];
rz(-pi) q[2];
rz(-2.3950855) q[3];
sx q[3];
rz(-2.3254447) q[3];
sx q[3];
rz(1.5619123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-2.7344446) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(2.8318185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6471841) q[0];
sx q[0];
rz(-1.4443047) q[0];
sx q[0];
rz(-2.6148318) q[0];
rz(-pi) q[1];
rz(-1.7477112) q[2];
sx q[2];
rz(-1.5801016) q[2];
sx q[2];
rz(-2.6975346) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5619547) q[1];
sx q[1];
rz(-1.5127752) q[1];
sx q[1];
rz(-1.7768363) q[1];
rz(-pi) q[2];
rz(0.36848948) q[3];
sx q[3];
rz(-0.62928761) q[3];
sx q[3];
rz(3.1143509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.2072198) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(0.65565482) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(-1.9357095) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.6419798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033894) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(3.035726) q[0];
rz(3.0707804) q[2];
sx q[2];
rz(-0.76528463) q[2];
sx q[2];
rz(2.8634957) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8780898) q[1];
sx q[1];
rz(-0.51551688) q[1];
sx q[1];
rz(1.131119) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5564735) q[3];
sx q[3];
rz(-0.17281547) q[3];
sx q[3];
rz(2.1567791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88400921) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-2.1949027) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832227) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(-0.46438607) q[2];
sx q[2];
rz(-2.0225564) q[2];
sx q[2];
rz(0.49973942) q[2];
rz(0.55340135) q[3];
sx q[3];
rz(-2.3407866) q[3];
sx q[3];
rz(1.8368807) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
