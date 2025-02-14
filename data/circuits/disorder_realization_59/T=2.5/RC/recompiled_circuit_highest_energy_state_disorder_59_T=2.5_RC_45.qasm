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
rz(-0.72379273) q[0];
sx q[0];
rz(-2.0279217) q[0];
sx q[0];
rz(0.96473515) q[0];
rz(-0.073702987) q[1];
sx q[1];
rz(-0.6414203) q[1];
sx q[1];
rz(1.6471656) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9091563) q[0];
sx q[0];
rz(-0.80003827) q[0];
sx q[0];
rz(2.590488) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2531075) q[2];
sx q[2];
rz(-1.471638) q[2];
sx q[2];
rz(0.74867189) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.6801939) q[1];
sx q[1];
rz(-2.0909266) q[1];
sx q[1];
rz(1.0375711) q[1];
x q[2];
rz(2.1001746) q[3];
sx q[3];
rz(-1.6106859) q[3];
sx q[3];
rz(-2.1969469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.73704314) q[2];
sx q[2];
rz(-1.5207542) q[2];
sx q[2];
rz(-0.21657319) q[2];
rz(0.51566044) q[3];
sx q[3];
rz(-2.0629081) q[3];
sx q[3];
rz(-0.083958538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0292922) q[0];
sx q[0];
rz(-3.0700505) q[0];
sx q[0];
rz(-1.4805502) q[0];
rz(-2.5484565) q[1];
sx q[1];
rz(-0.99855223) q[1];
sx q[1];
rz(-2.8592529) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4078699) q[0];
sx q[0];
rz(-1.5417011) q[0];
sx q[0];
rz(-3.0139774) q[0];
x q[1];
rz(2.1611358) q[2];
sx q[2];
rz(-0.23721248) q[2];
sx q[2];
rz(2.1179167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.39277252) q[1];
sx q[1];
rz(-1.7085679) q[1];
sx q[1];
rz(2.5974136) q[1];
x q[2];
rz(-2.2766791) q[3];
sx q[3];
rz(-0.54482116) q[3];
sx q[3];
rz(3.0044905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4234408) q[2];
sx q[2];
rz(-2.580692) q[2];
sx q[2];
rz(-2.1924428) q[2];
rz(3.0632422) q[3];
sx q[3];
rz(-1.6716985) q[3];
sx q[3];
rz(1.9115492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345728) q[0];
sx q[0];
rz(-0.0581352) q[0];
sx q[0];
rz(0.19244254) q[0];
rz(-2.1678534) q[1];
sx q[1];
rz(-2.1111919) q[1];
sx q[1];
rz(1.2724426) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5010492) q[0];
sx q[0];
rz(-1.328381) q[0];
sx q[0];
rz(0.69037504) q[0];
rz(-pi) q[1];
rz(-1.8143623) q[2];
sx q[2];
rz(-2.3911016) q[2];
sx q[2];
rz(1.0433973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7682926) q[1];
sx q[1];
rz(-1.321861) q[1];
sx q[1];
rz(1.8087923) q[1];
rz(1.4951823) q[3];
sx q[3];
rz(-2.6845884) q[3];
sx q[3];
rz(-1.6058894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5729313) q[2];
sx q[2];
rz(-0.91721407) q[2];
sx q[2];
rz(0.60079637) q[2];
rz(-2.3066547) q[3];
sx q[3];
rz(-2.2227414) q[3];
sx q[3];
rz(2.6981573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610953) q[0];
sx q[0];
rz(-1.1325855) q[0];
sx q[0];
rz(1.4196716) q[0];
rz(-0.28600606) q[1];
sx q[1];
rz(-1.1347457) q[1];
sx q[1];
rz(2.6223415) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88548286) q[0];
sx q[0];
rz(-1.7392495) q[0];
sx q[0];
rz(1.662111) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6250766) q[2];
sx q[2];
rz(-2.9331547) q[2];
sx q[2];
rz(-2.1303653) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1468526) q[1];
sx q[1];
rz(-2.6255872) q[1];
sx q[1];
rz(0.13928646) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.052922225) q[3];
sx q[3];
rz(-0.54659802) q[3];
sx q[3];
rz(1.6569759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1661735) q[2];
sx q[2];
rz(-0.68048802) q[2];
sx q[2];
rz(-0.55230459) q[2];
rz(0.78260261) q[3];
sx q[3];
rz(-1.8389333) q[3];
sx q[3];
rz(2.9569614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2833589) q[0];
sx q[0];
rz(-0.62223804) q[0];
sx q[0];
rz(1.2828113) q[0];
rz(-1.917631) q[1];
sx q[1];
rz(-1.6061648) q[1];
sx q[1];
rz(2.431869) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26742753) q[0];
sx q[0];
rz(-1.5418959) q[0];
sx q[0];
rz(-1.548882) q[0];
rz(1.5915839) q[2];
sx q[2];
rz(-1.1026303) q[2];
sx q[2];
rz(0.028108953) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2096055) q[1];
sx q[1];
rz(-2.8199204) q[1];
sx q[1];
rz(-1.7965921) q[1];
x q[2];
rz(-0.17931314) q[3];
sx q[3];
rz(-1.9783522) q[3];
sx q[3];
rz(-0.61791622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2964786) q[2];
sx q[2];
rz(-0.60135403) q[2];
sx q[2];
rz(2.5884886) q[2];
rz(-0.84135711) q[3];
sx q[3];
rz(-1.0480806) q[3];
sx q[3];
rz(-1.1355737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76310054) q[0];
sx q[0];
rz(-1.0569514) q[0];
sx q[0];
rz(-0.77504778) q[0];
rz(1.4729602) q[1];
sx q[1];
rz(-2.5170363) q[1];
sx q[1];
rz(-0.83736173) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178029) q[0];
sx q[0];
rz(-0.63321165) q[0];
sx q[0];
rz(1.4073611) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7532127) q[2];
sx q[2];
rz(-1.1216394) q[2];
sx q[2];
rz(2.7936008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9186975) q[1];
sx q[1];
rz(-0.59894005) q[1];
sx q[1];
rz(-0.75798613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0927518) q[3];
sx q[3];
rz(-0.47419729) q[3];
sx q[3];
rz(-1.3124581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2396635) q[2];
sx q[2];
rz(-1.1666802) q[2];
sx q[2];
rz(-2.3789294) q[2];
rz(2.934382) q[3];
sx q[3];
rz(-1.8800294) q[3];
sx q[3];
rz(-2.5337849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36829456) q[0];
sx q[0];
rz(-0.93450707) q[0];
sx q[0];
rz(-1.4542798) q[0];
rz(1.8621209) q[1];
sx q[1];
rz(-2.2984633) q[1];
sx q[1];
rz(-1.7284733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2103677) q[0];
sx q[0];
rz(-2.702091) q[0];
sx q[0];
rz(2.6511433) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94788289) q[2];
sx q[2];
rz(-2.3051734) q[2];
sx q[2];
rz(-1.9854058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3154732) q[1];
sx q[1];
rz(-2.5971088) q[1];
sx q[1];
rz(-2.524091) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2388494) q[3];
sx q[3];
rz(-2.6029426) q[3];
sx q[3];
rz(-2.1421632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1450682) q[2];
sx q[2];
rz(-2.2729496) q[2];
sx q[2];
rz(2.9800912) q[2];
rz(0.7343556) q[3];
sx q[3];
rz(-0.86638325) q[3];
sx q[3];
rz(-1.8374779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.9228009) q[0];
sx q[0];
rz(-2.6113593) q[0];
sx q[0];
rz(-0.23736048) q[0];
rz(-2.3573719) q[1];
sx q[1];
rz(-1.3619245) q[1];
sx q[1];
rz(-1.721419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7331384) q[0];
sx q[0];
rz(-1.8693921) q[0];
sx q[0];
rz(-0.51513735) q[0];
x q[1];
rz(2.2591822) q[2];
sx q[2];
rz(-1.7441445) q[2];
sx q[2];
rz(-1.2899931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8964149) q[1];
sx q[1];
rz(-0.62508622) q[1];
sx q[1];
rz(-2.6253789) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5077614) q[3];
sx q[3];
rz(-2.0326893) q[3];
sx q[3];
rz(0.92538242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64804849) q[2];
sx q[2];
rz(-1.250114) q[2];
sx q[2];
rz(1.9165215) q[2];
rz(-1.8152292) q[3];
sx q[3];
rz(-2.1882961) q[3];
sx q[3];
rz(1.1939322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61106435) q[0];
sx q[0];
rz(-0.079212991) q[0];
sx q[0];
rz(-2.5031669) q[0];
rz(-3.0910659) q[1];
sx q[1];
rz(-1.2082929) q[1];
sx q[1];
rz(-0.60727492) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716403) q[0];
sx q[0];
rz(-1.1881614) q[0];
sx q[0];
rz(0.34374551) q[0];
rz(-pi) q[1];
rz(2.2556317) q[2];
sx q[2];
rz(-1.6185624) q[2];
sx q[2];
rz(-1.1807962) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.133195) q[1];
sx q[1];
rz(-1.6411157) q[1];
sx q[1];
rz(-1.2446872) q[1];
rz(-pi) q[2];
rz(2.9565892) q[3];
sx q[3];
rz(-1.7599938) q[3];
sx q[3];
rz(-0.69533492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66024572) q[2];
sx q[2];
rz(-1.9504184) q[2];
sx q[2];
rz(-0.78835362) q[2];
rz(-2.3039019) q[3];
sx q[3];
rz(-0.82101429) q[3];
sx q[3];
rz(-1.8590417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7343219) q[0];
sx q[0];
rz(-0.73198524) q[0];
sx q[0];
rz(0.5087854) q[0];
rz(1.58163) q[1];
sx q[1];
rz(-1.0204693) q[1];
sx q[1];
rz(-0.4090974) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74241766) q[0];
sx q[0];
rz(-1.8054726) q[0];
sx q[0];
rz(-0.15813903) q[0];
x q[1];
rz(1.0739543) q[2];
sx q[2];
rz(-1.1449305) q[2];
sx q[2];
rz(-1.9857621) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4933192) q[1];
sx q[1];
rz(-1.5602605) q[1];
sx q[1];
rz(0.22064836) q[1];
rz(-pi) q[2];
rz(-1.9473829) q[3];
sx q[3];
rz(-1.9344182) q[3];
sx q[3];
rz(-2.0998433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6009377) q[2];
sx q[2];
rz(-1.6161852) q[2];
sx q[2];
rz(1.8589004) q[2];
rz(-3.1278074) q[3];
sx q[3];
rz(-1.1252334) q[3];
sx q[3];
rz(-1.8085326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29671058) q[0];
sx q[0];
rz(-1.409335) q[0];
sx q[0];
rz(0.62406337) q[0];
rz(0.86112549) q[1];
sx q[1];
rz(-0.54665165) q[1];
sx q[1];
rz(0.13269592) q[1];
rz(0.3202318) q[2];
sx q[2];
rz(-1.474517) q[2];
sx q[2];
rz(1.3040645) q[2];
rz(2.5100711) q[3];
sx q[3];
rz(-1.0167663) q[3];
sx q[3];
rz(2.9192703) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
