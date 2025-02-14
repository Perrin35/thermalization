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
rz(1.1627816) q[0];
sx q[0];
rz(5.3428234) q[0];
sx q[0];
rz(12.760395) q[0];
rz(-0.54430517) q[1];
sx q[1];
rz(-1.5679918) q[1];
sx q[1];
rz(-0.1967217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5618043) q[0];
sx q[0];
rz(-1.4477948) q[0];
sx q[0];
rz(0.99539006) q[0];
rz(-pi) q[1];
rz(-2.9241184) q[2];
sx q[2];
rz(-1.5921235) q[2];
sx q[2];
rz(2.3623717) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87952033) q[1];
sx q[1];
rz(-1.889887) q[1];
sx q[1];
rz(-1.2882922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.807735) q[3];
sx q[3];
rz(-0.58962017) q[3];
sx q[3];
rz(-2.9350964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0611614) q[2];
sx q[2];
rz(-1.4208527) q[2];
sx q[2];
rz(-0.013193456) q[2];
rz(-2.0773928) q[3];
sx q[3];
rz(-1.0366169) q[3];
sx q[3];
rz(-0.17087759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6212293) q[0];
sx q[0];
rz(-0.55084387) q[0];
sx q[0];
rz(-2.6589822) q[0];
rz(-2.6699325) q[1];
sx q[1];
rz(-2.1444247) q[1];
sx q[1];
rz(2.4536123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4783966) q[0];
sx q[0];
rz(-1.8747878) q[0];
sx q[0];
rz(-1.770397) q[0];
rz(-pi) q[1];
rz(1.2763763) q[2];
sx q[2];
rz(-0.79508077) q[2];
sx q[2];
rz(0.15995041) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.019362018) q[1];
sx q[1];
rz(-1.6145585) q[1];
sx q[1];
rz(-2.6701609) q[1];
rz(-pi) q[2];
rz(1.7741995) q[3];
sx q[3];
rz(-1.6546094) q[3];
sx q[3];
rz(2.0544685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.017435) q[2];
sx q[2];
rz(-2.0686801) q[2];
sx q[2];
rz(0.2249445) q[2];
rz(2.0624835) q[3];
sx q[3];
rz(-0.93828097) q[3];
sx q[3];
rz(2.6911531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17063046) q[0];
sx q[0];
rz(-0.58929515) q[0];
sx q[0];
rz(1.851409) q[0];
rz(-1.2782485) q[1];
sx q[1];
rz(-1.1561013) q[1];
sx q[1];
rz(1.2826756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0874003) q[0];
sx q[0];
rz(-2.4303959) q[0];
sx q[0];
rz(1.0048701) q[0];
rz(-pi) q[1];
rz(-3.01802) q[2];
sx q[2];
rz(-1.6500435) q[2];
sx q[2];
rz(1.0956956) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8508529) q[1];
sx q[1];
rz(-3.1147482) q[1];
sx q[1];
rz(-0.82976262) q[1];
rz(-pi) q[2];
rz(-0.81961378) q[3];
sx q[3];
rz(-0.96246877) q[3];
sx q[3];
rz(1.5386594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67905417) q[2];
sx q[2];
rz(-1.6109698) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3731641) q[0];
sx q[0];
rz(-1.3685065) q[0];
sx q[0];
rz(2.9737293) q[0];
rz(-0.31790512) q[1];
sx q[1];
rz(-1.2891506) q[1];
sx q[1];
rz(-1.0344523) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8694591) q[0];
sx q[0];
rz(-1.3087166) q[0];
sx q[0];
rz(-0.036152187) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1940766) q[2];
sx q[2];
rz(-2.2206306) q[2];
sx q[2];
rz(1.5638994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2047836) q[1];
sx q[1];
rz(-2.2310871) q[1];
sx q[1];
rz(0.59548638) q[1];
rz(2.529783) q[3];
sx q[3];
rz(-1.4880848) q[3];
sx q[3];
rz(-1.8066607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1075403) q[2];
sx q[2];
rz(-0.62959051) q[2];
sx q[2];
rz(0.26629392) q[2];
rz(2.9722424) q[3];
sx q[3];
rz(-2.0784056) q[3];
sx q[3];
rz(0.17899409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87194815) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(-1.2931152) q[0];
rz(0.070501892) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(2.3354796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370395) q[0];
sx q[0];
rz(-1.4463639) q[0];
sx q[0];
rz(1.7794153) q[0];
rz(-pi) q[1];
x q[1];
rz(0.057439645) q[2];
sx q[2];
rz(-1.4220793) q[2];
sx q[2];
rz(-1.1401389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5841585) q[1];
sx q[1];
rz(-0.63156742) q[1];
sx q[1];
rz(-2.8110792) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6025196) q[3];
sx q[3];
rz(-1.381449) q[3];
sx q[3];
rz(0.93261485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21141323) q[2];
sx q[2];
rz(-1.6665919) q[2];
sx q[2];
rz(-2.831947) q[2];
rz(-0.51112255) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(2.3206594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822134) q[0];
sx q[0];
rz(-0.62726227) q[0];
sx q[0];
rz(0.43701592) q[0];
rz(0.45477319) q[1];
sx q[1];
rz(-0.79650703) q[1];
sx q[1];
rz(2.2589267) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1293212) q[0];
sx q[0];
rz(-0.28421775) q[0];
sx q[0];
rz(1.8924665) q[0];
x q[1];
rz(-1.1922311) q[2];
sx q[2];
rz(-1.541271) q[2];
sx q[2];
rz(1.6395417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5875579) q[1];
sx q[1];
rz(-1.2393426) q[1];
sx q[1];
rz(-2.0786839) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96787937) q[3];
sx q[3];
rz(-1.9390188) q[3];
sx q[3];
rz(0.24711025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5060045) q[2];
sx q[2];
rz(-1.4281861) q[2];
sx q[2];
rz(0.44090718) q[2];
rz(1.0718369) q[3];
sx q[3];
rz(-0.25522885) q[3];
sx q[3];
rz(2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014024409) q[0];
sx q[0];
rz(-1.2034282) q[0];
sx q[0];
rz(1.2822275) q[0];
rz(-1.0857238) q[1];
sx q[1];
rz(-1.5241357) q[1];
sx q[1];
rz(-1.5132743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8305539) q[0];
sx q[0];
rz(-2.3966602) q[0];
sx q[0];
rz(5*pi/9) q[0];
x q[1];
rz(2.4413928) q[2];
sx q[2];
rz(-1.321903) q[2];
sx q[2];
rz(0.85485103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.46678695) q[1];
sx q[1];
rz(-1.9478096) q[1];
sx q[1];
rz(2.2598355) q[1];
rz(2.3845838) q[3];
sx q[3];
rz(-1.2881713) q[3];
sx q[3];
rz(0.69863897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6699803) q[2];
sx q[2];
rz(-0.70576224) q[2];
sx q[2];
rz(-0.82028779) q[2];
rz(2.7465076) q[3];
sx q[3];
rz(-0.79212752) q[3];
sx q[3];
rz(2.5274966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87913269) q[0];
sx q[0];
rz(-1.4985871) q[0];
sx q[0];
rz(0.60687989) q[0];
rz(0.34582368) q[1];
sx q[1];
rz(-1.1864097) q[1];
sx q[1];
rz(-0.80723673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31921747) q[0];
sx q[0];
rz(-1.2273047) q[0];
sx q[0];
rz(-0.10256711) q[0];
rz(2.1431461) q[2];
sx q[2];
rz(-0.46914161) q[2];
sx q[2];
rz(-2.1690705) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2507627) q[1];
sx q[1];
rz(-1.8381665) q[1];
sx q[1];
rz(1.7946589) q[1];
x q[2];
rz(-0.96514197) q[3];
sx q[3];
rz(-2.0803841) q[3];
sx q[3];
rz(1.6965564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81749934) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(0.49638805) q[2];
rz(-3.056622) q[3];
sx q[3];
rz(-0.43046633) q[3];
sx q[3];
rz(2.5607204) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896511) q[0];
sx q[0];
rz(-1.4284416) q[0];
sx q[0];
rz(-0.38452837) q[0];
rz(0.25228581) q[1];
sx q[1];
rz(-1.7583881) q[1];
sx q[1];
rz(1.5834454) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66551946) q[0];
sx q[0];
rz(-1.5613149) q[0];
sx q[0];
rz(-0.46561636) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6884638) q[2];
sx q[2];
rz(-1.8431547) q[2];
sx q[2];
rz(0.54335574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0817683) q[1];
sx q[1];
rz(-1.5835973) q[1];
sx q[1];
rz(-0.98155419) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83797601) q[3];
sx q[3];
rz(-2.1065578) q[3];
sx q[3];
rz(2.4934078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3462476) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(2.0015049) q[2];
rz(1.358486) q[3];
sx q[3];
rz(-1.2525109) q[3];
sx q[3];
rz(-0.31806773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43050218) q[0];
sx q[0];
rz(-1.6999812) q[0];
sx q[0];
rz(2.7464113) q[0];
rz(-1.9197865) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(0.88776678) q[1];
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
rz(-1.5776921) q[2];
sx q[2];
rz(-0.47499945) q[2];
sx q[2];
rz(-1.4003225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.58037591) q[1];
sx q[1];
rz(-1.2151011) q[1];
sx q[1];
rz(-0.027173398) q[1];
x q[2];
rz(1.8818782) q[3];
sx q[3];
rz(-1.4691989) q[3];
sx q[3];
rz(1.7465308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6937574) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(2.5698404) q[2];
rz(-1.5104431) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(-1.800764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1258662) q[0];
sx q[0];
rz(-1.8342352) q[0];
sx q[0];
rz(-2.9974708) q[0];
rz(-1.9563328) q[1];
sx q[1];
rz(-2.2846501) q[1];
sx q[1];
rz(1.2774998) q[1];
rz(-1.324426) q[2];
sx q[2];
rz(-1.194321) q[2];
sx q[2];
rz(1.3512667) q[2];
rz(2.4972054) q[3];
sx q[3];
rz(-1.3105262) q[3];
sx q[3];
rz(-2.5986828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
