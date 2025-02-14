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
rz(-1.0492078) q[0];
sx q[0];
rz(-2.3823491) q[0];
sx q[0];
rz(-2.7651751) q[0];
rz(1.5578101) q[1];
sx q[1];
rz(-1.0695142) q[1];
sx q[1];
rz(-2.1853316) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.555684) q[0];
sx q[0];
rz(-2.1901185) q[0];
sx q[0];
rz(0.39062087) q[0];
rz(-1.7643433) q[2];
sx q[2];
rz(-1.9864533) q[2];
sx q[2];
rz(-0.035482835) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6105683) q[1];
sx q[1];
rz(-2.3914365) q[1];
sx q[1];
rz(1.7498739) q[1];
x q[2];
rz(0.27122916) q[3];
sx q[3];
rz(-1.3406957) q[3];
sx q[3];
rz(2.7903737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.817953) q[2];
sx q[2];
rz(-0.91922593) q[2];
sx q[2];
rz(-1.9383355) q[2];
rz(2.1008927) q[3];
sx q[3];
rz(-0.87472707) q[3];
sx q[3];
rz(3.021595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0284718) q[0];
sx q[0];
rz(-2.5682243) q[0];
sx q[0];
rz(2.0290802) q[0];
rz(-1.6587229) q[1];
sx q[1];
rz(-2.5506546) q[1];
sx q[1];
rz(-2.9842751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4388889) q[0];
sx q[0];
rz(-1.3957141) q[0];
sx q[0];
rz(0.69467993) q[0];
rz(-pi) q[1];
rz(2.9221465) q[2];
sx q[2];
rz(-0.38139653) q[2];
sx q[2];
rz(1.4666605) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29780218) q[1];
sx q[1];
rz(-2.4239967) q[1];
sx q[1];
rz(1.1949431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4955161) q[3];
sx q[3];
rz(-1.3485678) q[3];
sx q[3];
rz(0.27929515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6436254) q[2];
sx q[2];
rz(-0.47933856) q[2];
sx q[2];
rz(2.7304926) q[2];
rz(0.57605612) q[3];
sx q[3];
rz(-1.0141625) q[3];
sx q[3];
rz(-2.1435553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75329798) q[0];
sx q[0];
rz(-1.0030712) q[0];
sx q[0];
rz(2.4004747) q[0];
rz(1.6916212) q[1];
sx q[1];
rz(-1.8284109) q[1];
sx q[1];
rz(0.57201874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2021671) q[0];
sx q[0];
rz(-1.8644445) q[0];
sx q[0];
rz(-2.9333326) q[0];
rz(-pi) q[1];
rz(0.26511221) q[2];
sx q[2];
rz(-1.7705317) q[2];
sx q[2];
rz(2.0743362) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95840824) q[1];
sx q[1];
rz(-0.60380492) q[1];
sx q[1];
rz(0.81012695) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23046646) q[3];
sx q[3];
rz(-1.8982045) q[3];
sx q[3];
rz(-2.5203506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.584562) q[2];
sx q[2];
rz(-1.9494373) q[2];
sx q[2];
rz(-0.7977879) q[2];
rz(-1.8320463) q[3];
sx q[3];
rz(-1.8833501) q[3];
sx q[3];
rz(1.7358739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0107467) q[0];
sx q[0];
rz(-2.2444785) q[0];
sx q[0];
rz(-2.77453) q[0];
rz(1.2182073) q[1];
sx q[1];
rz(-1.7299078) q[1];
sx q[1];
rz(0.22148111) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9629134) q[0];
sx q[0];
rz(-1.7811462) q[0];
sx q[0];
rz(-1.5014929) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34299739) q[2];
sx q[2];
rz(-0.38831899) q[2];
sx q[2];
rz(1.8754387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1180956) q[1];
sx q[1];
rz(-2.6129236) q[1];
sx q[1];
rz(0.65331991) q[1];
x q[2];
rz(-2.6148197) q[3];
sx q[3];
rz(-1.793141) q[3];
sx q[3];
rz(0.46565817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9941142) q[2];
sx q[2];
rz(-1.9693815) q[2];
sx q[2];
rz(1.3383024) q[2];
rz(-0.5021247) q[3];
sx q[3];
rz(-3.0410671) q[3];
sx q[3];
rz(-2.8978735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.5153656) q[0];
sx q[0];
rz(-2.5318662) q[0];
sx q[0];
rz(1.6060265) q[0];
rz(-3.0961127) q[1];
sx q[1];
rz(-1.7605503) q[1];
sx q[1];
rz(2.6754726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9392747) q[0];
sx q[0];
rz(-0.28616787) q[0];
sx q[0];
rz(-1.0658468) q[0];
x q[1];
rz(0.32987288) q[2];
sx q[2];
rz(-0.95392841) q[2];
sx q[2];
rz(2.5170779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.921195) q[1];
sx q[1];
rz(-0.092893727) q[1];
sx q[1];
rz(-2.9542406) q[1];
rz(0.65069549) q[3];
sx q[3];
rz(-2.0632207) q[3];
sx q[3];
rz(0.7975815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8972299) q[2];
sx q[2];
rz(-1.676061) q[2];
sx q[2];
rz(-0.7828632) q[2];
rz(-0.016294567) q[3];
sx q[3];
rz(-1.8330845) q[3];
sx q[3];
rz(1.0619987) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29430729) q[0];
sx q[0];
rz(-2.6562302) q[0];
sx q[0];
rz(-2.7622727) q[0];
rz(0.30917057) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(-2.9177623) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9681518) q[0];
sx q[0];
rz(-1.8814527) q[0];
sx q[0];
rz(2.2947427) q[0];
x q[1];
rz(1.3794704) q[2];
sx q[2];
rz(-1.680964) q[2];
sx q[2];
rz(-0.7502816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31948446) q[1];
sx q[1];
rz(-2.9908239) q[1];
sx q[1];
rz(-2.3586078) q[1];
rz(1.0204487) q[3];
sx q[3];
rz(-0.942505) q[3];
sx q[3];
rz(-0.8250784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0763187) q[2];
sx q[2];
rz(-1.167401) q[2];
sx q[2];
rz(-0.25263146) q[2];
rz(-0.97918716) q[3];
sx q[3];
rz(-2.2604209) q[3];
sx q[3];
rz(-0.34669909) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52794367) q[0];
sx q[0];
rz(-0.76512965) q[0];
sx q[0];
rz(1.2921523) q[0];
rz(-2.6076803) q[1];
sx q[1];
rz(-1.4521867) q[1];
sx q[1];
rz(0.36010489) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2947073) q[0];
sx q[0];
rz(-0.36723235) q[0];
sx q[0];
rz(-1.7647554) q[0];
rz(-pi) q[1];
rz(1.8890284) q[2];
sx q[2];
rz(-1.1063853) q[2];
sx q[2];
rz(1.1030674) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26038995) q[1];
sx q[1];
rz(-2.5298928) q[1];
sx q[1];
rz(-2.7124497) q[1];
rz(-pi) q[2];
rz(1.259616) q[3];
sx q[3];
rz(-3.0336685) q[3];
sx q[3];
rz(0.51801658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1503633) q[2];
sx q[2];
rz(-1.5183307) q[2];
sx q[2];
rz(2.4156127) q[2];
rz(0.43068543) q[3];
sx q[3];
rz(-0.78023282) q[3];
sx q[3];
rz(-0.45346692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78366572) q[0];
sx q[0];
rz(-1.6538606) q[0];
sx q[0];
rz(-0.7861535) q[0];
rz(0.17732009) q[1];
sx q[1];
rz(-2.0568078) q[1];
sx q[1];
rz(-2.1254983) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35619104) q[0];
sx q[0];
rz(-1.2895661) q[0];
sx q[0];
rz(-0.10939557) q[0];
rz(-pi) q[1];
rz(-1.4403604) q[2];
sx q[2];
rz(-1.36048) q[2];
sx q[2];
rz(0.73784251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7238771) q[1];
sx q[1];
rz(-0.9068562) q[1];
sx q[1];
rz(-1.5491897) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75288624) q[3];
sx q[3];
rz(-1.2883923) q[3];
sx q[3];
rz(1.2644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59755406) q[2];
sx q[2];
rz(-0.95557135) q[2];
sx q[2];
rz(2.3717144) q[2];
rz(0.17635135) q[3];
sx q[3];
rz(-1.4498815) q[3];
sx q[3];
rz(2.8281853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4835994) q[0];
sx q[0];
rz(-0.083881065) q[0];
sx q[0];
rz(2.0954848) q[0];
rz(-1.5931891) q[1];
sx q[1];
rz(-1.4175339) q[1];
sx q[1];
rz(1.2215325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6283327) q[0];
sx q[0];
rz(-1.9258537) q[0];
sx q[0];
rz(3.0299788) q[0];
rz(-0.1510001) q[2];
sx q[2];
rz(-1.2556453) q[2];
sx q[2];
rz(1.1282819) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9323349) q[1];
sx q[1];
rz(-0.29587165) q[1];
sx q[1];
rz(2.3274755) q[1];
rz(-pi) q[2];
rz(3.0980403) q[3];
sx q[3];
rz(-1.5659025) q[3];
sx q[3];
rz(2.2088449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.017642411) q[2];
sx q[2];
rz(-1.7240179) q[2];
sx q[2];
rz(-1.1851912) q[2];
rz(-2.8017398) q[3];
sx q[3];
rz(-0.9797107) q[3];
sx q[3];
rz(3.0237696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3621984) q[0];
sx q[0];
rz(-1.8147991) q[0];
sx q[0];
rz(0.69778824) q[0];
rz(-0.70480529) q[1];
sx q[1];
rz(-1.1230527) q[1];
sx q[1];
rz(2.8885081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20257178) q[0];
sx q[0];
rz(-0.24469412) q[0];
sx q[0];
rz(1.2608725) q[0];
x q[1];
rz(1.0519876) q[2];
sx q[2];
rz(-2.1588783) q[2];
sx q[2];
rz(-1.1747509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29927897) q[1];
sx q[1];
rz(-0.83228534) q[1];
sx q[1];
rz(-0.93361093) q[1];
x q[2];
rz(0.730597) q[3];
sx q[3];
rz(-1.3775148) q[3];
sx q[3];
rz(2.766942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7207328) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(-2.4534658) q[2];
rz(-2.1963035) q[3];
sx q[3];
rz(-1.9663845) q[3];
sx q[3];
rz(-2.2139886) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59104334) q[0];
sx q[0];
rz(-1.3591546) q[0];
sx q[0];
rz(-1.2426283) q[0];
rz(-0.91724829) q[1];
sx q[1];
rz(-1.6684253) q[1];
sx q[1];
rz(-2.565276) q[1];
rz(-1.9333712) q[2];
sx q[2];
rz(-2.1167663) q[2];
sx q[2];
rz(2.703985) q[2];
rz(3.0455899) q[3];
sx q[3];
rz(-1.1461555) q[3];
sx q[3];
rz(-2.8050527) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
