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
rz(-1.9788111) q[0];
sx q[0];
rz(-2.2012308) q[0];
sx q[0];
rz(2.9475687) q[0];
rz(-0.54430517) q[1];
sx q[1];
rz(-1.5679918) q[1];
sx q[1];
rz(-0.1967217) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17795632) q[0];
sx q[0];
rz(-2.5546409) q[0];
sx q[0];
rz(1.3474083) q[0];
rz(-pi) q[1];
x q[1];
rz(0.098539515) q[2];
sx q[2];
rz(-0.218501) q[2];
sx q[2];
rz(-0.69536415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3595092) q[1];
sx q[1];
rz(-1.8386786) q[1];
sx q[1];
rz(-0.33133502) q[1];
x q[2];
rz(2.9858304) q[3];
sx q[3];
rz(-2.1418396) q[3];
sx q[3];
rz(0.076249853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0611614) q[2];
sx q[2];
rz(-1.4208527) q[2];
sx q[2];
rz(-0.013193456) q[2];
rz(2.0773928) q[3];
sx q[3];
rz(-2.1049757) q[3];
sx q[3];
rz(-0.17087759) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5203633) q[0];
sx q[0];
rz(-2.5907488) q[0];
sx q[0];
rz(-2.6589822) q[0];
rz(-2.6699325) q[1];
sx q[1];
rz(-0.99716798) q[1];
sx q[1];
rz(0.68798033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88408282) q[0];
sx q[0];
rz(-0.36195746) q[0];
sx q[0];
rz(0.56364001) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2763763) q[2];
sx q[2];
rz(-2.3465119) q[2];
sx q[2];
rz(-2.9816422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5045484) q[1];
sx q[1];
rz(-0.47330644) q[1];
sx q[1];
rz(-0.096122336) q[1];
rz(-pi) q[2];
rz(-1.7741995) q[3];
sx q[3];
rz(-1.6546094) q[3];
sx q[3];
rz(1.0871241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.017435) q[2];
sx q[2];
rz(-1.0729125) q[2];
sx q[2];
rz(-2.9166481) q[2];
rz(-2.0624835) q[3];
sx q[3];
rz(-2.2033117) q[3];
sx q[3];
rz(2.6911531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17063046) q[0];
sx q[0];
rz(-2.5522975) q[0];
sx q[0];
rz(1.2901837) q[0];
rz(-1.8633441) q[1];
sx q[1];
rz(-1.9854913) q[1];
sx q[1];
rz(1.2826756) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0874003) q[0];
sx q[0];
rz(-2.4303959) q[0];
sx q[0];
rz(-1.0048701) q[0];
rz(-2.569246) q[2];
sx q[2];
rz(-0.1466912) q[2];
sx q[2];
rz(1.0425422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1206818) q[1];
sx q[1];
rz(-1.5889165) q[1];
sx q[1];
rz(-1.5906035) q[1];
rz(-0.77529737) q[3];
sx q[3];
rz(-0.92758815) q[3];
sx q[3];
rz(-2.624529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67905417) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(2.6666717) q[2];
rz(1.9654407) q[3];
sx q[3];
rz(-0.85634309) q[3];
sx q[3];
rz(-0.28158751) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7684286) q[0];
sx q[0];
rz(-1.7730862) q[0];
sx q[0];
rz(2.9737293) q[0];
rz(-0.31790512) q[1];
sx q[1];
rz(-1.8524421) q[1];
sx q[1];
rz(-2.1071404) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2892923) q[0];
sx q[0];
rz(-1.5358791) q[0];
sx q[0];
rz(-1.8330397) q[0];
rz(-pi) q[1];
rz(0.6851721) q[2];
sx q[2];
rz(-1.2735441) q[2];
sx q[2];
rz(-2.9135421) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2047836) q[1];
sx q[1];
rz(-2.2310871) q[1];
sx q[1];
rz(-2.5461063) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4698704) q[3];
sx q[3];
rz(-0.96138326) q[3];
sx q[3];
rz(0.17796365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0340524) q[2];
sx q[2];
rz(-2.5120021) q[2];
sx q[2];
rz(-0.26629392) q[2];
rz(0.1693503) q[3];
sx q[3];
rz(-1.063187) q[3];
sx q[3];
rz(0.17899409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2696445) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(1.2931152) q[0];
rz(3.0710908) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(0.806113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9490813) q[0];
sx q[0];
rz(-1.7777788) q[0];
sx q[0];
rz(-3.0144318) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.057439645) q[2];
sx q[2];
rz(-1.4220793) q[2];
sx q[2];
rz(1.1401389) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.986113) q[1];
sx q[1];
rz(-0.9782741) q[1];
sx q[1];
rz(1.3377124) q[1];
rz(-pi) q[2];
rz(-0.53907303) q[3];
sx q[3];
rz(-1.381449) q[3];
sx q[3];
rz(2.2089778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9301794) q[2];
sx q[2];
rz(-1.6665919) q[2];
sx q[2];
rz(-0.30964568) q[2];
rz(-2.6304701) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(-2.3206594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1593793) q[0];
sx q[0];
rz(-2.5143304) q[0];
sx q[0];
rz(2.7045767) q[0];
rz(-2.6868195) q[1];
sx q[1];
rz(-2.3450856) q[1];
sx q[1];
rz(-2.2589267) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1293212) q[0];
sx q[0];
rz(-2.8573749) q[0];
sx q[0];
rz(1.2491262) q[0];
rz(-pi) q[1];
rz(1.4910554) q[2];
sx q[2];
rz(-0.3796595) q[2];
sx q[2];
rz(-0.14282957) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5875579) q[1];
sx q[1];
rz(-1.9022501) q[1];
sx q[1];
rz(1.0629088) q[1];
rz(-2.168247) q[3];
sx q[3];
rz(-2.4472467) q[3];
sx q[3];
rz(2.2996817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5060045) q[2];
sx q[2];
rz(-1.7134066) q[2];
sx q[2];
rz(-0.44090718) q[2];
rz(-2.0697557) q[3];
sx q[3];
rz(-0.25522885) q[3];
sx q[3];
rz(2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014024409) q[0];
sx q[0];
rz(-1.9381645) q[0];
sx q[0];
rz(1.2822275) q[0];
rz(1.0857238) q[1];
sx q[1];
rz(-1.5241357) q[1];
sx q[1];
rz(-1.6283183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075628932) q[0];
sx q[0];
rz(-2.3018078) q[0];
sx q[0];
rz(-2.9828067) q[0];
x q[1];
rz(-0.70019987) q[2];
sx q[2];
rz(-1.321903) q[2];
sx q[2];
rz(0.85485103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3320184) q[1];
sx q[1];
rz(-2.2032713) q[1];
sx q[1];
rz(-2.6676086) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40005584) q[3];
sx q[3];
rz(-0.7981185) q[3];
sx q[3];
rz(1.9824073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6699803) q[2];
sx q[2];
rz(-2.4358304) q[2];
sx q[2];
rz(-0.82028779) q[2];
rz(-0.39508501) q[3];
sx q[3];
rz(-0.79212752) q[3];
sx q[3];
rz(2.5274966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87913269) q[0];
sx q[0];
rz(-1.6430055) q[0];
sx q[0];
rz(-2.5347128) q[0];
rz(-2.795769) q[1];
sx q[1];
rz(-1.9551829) q[1];
sx q[1];
rz(0.80723673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8223752) q[0];
sx q[0];
rz(-1.2273047) q[0];
sx q[0];
rz(0.10256711) q[0];
rz(-pi) q[1];
rz(0.99844653) q[2];
sx q[2];
rz(-0.46914161) q[2];
sx q[2];
rz(2.1690705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2507627) q[1];
sx q[1];
rz(-1.3034262) q[1];
sx q[1];
rz(1.3469338) q[1];
x q[2];
rz(-0.96514197) q[3];
sx q[3];
rz(-2.0803841) q[3];
sx q[3];
rz(-1.4450362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3240933) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(-0.49638805) q[2];
rz(0.084970623) q[3];
sx q[3];
rz(-0.43046633) q[3];
sx q[3];
rz(-0.5808723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519415) q[0];
sx q[0];
rz(-1.4284416) q[0];
sx q[0];
rz(0.38452837) q[0];
rz(-2.8893068) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(-1.5834454) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2551832) q[0];
sx q[0];
rz(-2.6758868) q[0];
sx q[0];
rz(0.021115594) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27416147) q[2];
sx q[2];
rz(-1.4574852) q[2];
sx q[2];
rz(2.1459412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0817683) q[1];
sx q[1];
rz(-1.5835973) q[1];
sx q[1];
rz(-2.1600385) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6739613) q[3];
sx q[3];
rz(-2.1836851) q[3];
sx q[3];
rz(-0.4919258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(1.1400878) q[2];
rz(1.7831066) q[3];
sx q[3];
rz(-1.8890817) q[3];
sx q[3];
rz(2.8235249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110905) q[0];
sx q[0];
rz(-1.4416114) q[0];
sx q[0];
rz(-2.7464113) q[0];
rz(1.2218062) q[1];
sx q[1];
rz(-2.2876574) q[1];
sx q[1];
rz(-0.88776678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319755) q[0];
sx q[0];
rz(-1.8450301) q[0];
sx q[0];
rz(0.72276177) q[0];
rz(1.5776921) q[2];
sx q[2];
rz(-0.47499945) q[2];
sx q[2];
rz(1.4003225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.58037591) q[1];
sx q[1];
rz(-1.2151011) q[1];
sx q[1];
rz(3.1144193) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10668175) q[3];
sx q[3];
rz(-1.8802207) q[3];
sx q[3];
rz(2.9984561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6937574) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(-0.57175222) q[2];
rz(-1.6311496) q[3];
sx q[3];
rz(-1.7084728) q[3];
sx q[3];
rz(-1.800764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0157264) q[0];
sx q[0];
rz(-1.8342352) q[0];
sx q[0];
rz(-2.9974708) q[0];
rz(1.9563328) q[1];
sx q[1];
rz(-0.85694255) q[1];
sx q[1];
rz(-1.8640929) q[1];
rz(0.55276362) q[2];
sx q[2];
rz(-2.6949099) q[2];
sx q[2];
rz(0.75134122) q[2];
rz(-0.64438728) q[3];
sx q[3];
rz(-1.3105262) q[3];
sx q[3];
rz(-2.5986828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
