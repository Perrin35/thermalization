OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(-0.88859963) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(-1.5421901) q[1];
sx q[1];
rz(1.3808274) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61533538) q[0];
sx q[0];
rz(-1.3852296) q[0];
sx q[0];
rz(-0.19677563) q[0];
x q[1];
rz(1.0384473) q[2];
sx q[2];
rz(-2.2962036) q[2];
sx q[2];
rz(-0.40276819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4145248) q[1];
sx q[1];
rz(-1.5611851) q[1];
sx q[1];
rz(-0.33551402) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9350151) q[3];
sx q[3];
rz(-1.7813588) q[3];
sx q[3];
rz(-0.59277804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8218653) q[2];
sx q[2];
rz(-1.8165908) q[2];
sx q[2];
rz(-0.74903178) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(-0.099893959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(1.9248167) q[0];
rz(2.079839) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(-2.7144576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7593133) q[0];
sx q[0];
rz(-2.1370892) q[0];
sx q[0];
rz(1.9840389) q[0];
x q[1];
rz(-2.8982694) q[2];
sx q[2];
rz(-2.2017225) q[2];
sx q[2];
rz(0.48829809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.037405326) q[1];
sx q[1];
rz(-1.4302642) q[1];
sx q[1];
rz(-1.8611845) q[1];
x q[2];
rz(-1.8187567) q[3];
sx q[3];
rz(-2.9313847) q[3];
sx q[3];
rz(1.5402286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(-1.742935) q[2];
rz(0.22377293) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021521213) q[0];
sx q[0];
rz(-1.7912309) q[0];
sx q[0];
rz(-3.0960826) q[0];
rz(-1.9728164) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(-0.57560903) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110014) q[0];
sx q[0];
rz(-1.5596034) q[0];
sx q[0];
rz(2.5024947) q[0];
rz(0.91012886) q[2];
sx q[2];
rz(-1.5843387) q[2];
sx q[2];
rz(-2.8000268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3176885) q[1];
sx q[1];
rz(-1.0270734) q[1];
sx q[1];
rz(0.25441092) q[1];
x q[2];
rz(-3.1166385) q[3];
sx q[3];
rz(-2.2234699) q[3];
sx q[3];
rz(3.0252473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(-2.1843145) q[2];
rz(0.57473985) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(-1.3055698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15605536) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(1.8804469) q[0];
rz(3.0592697) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(1.7877158) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1151617) q[0];
sx q[0];
rz(-1.8055834) q[0];
sx q[0];
rz(2.874766) q[0];
x q[1];
rz(-0.33892314) q[2];
sx q[2];
rz(-0.51033516) q[2];
sx q[2];
rz(1.5938544) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8435858) q[1];
sx q[1];
rz(-1.516696) q[1];
sx q[1];
rz(0.49523103) q[1];
rz(-pi) q[2];
rz(2.6526101) q[3];
sx q[3];
rz(-2.1900898) q[3];
sx q[3];
rz(-1.796738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7202619) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(1.8850231) q[2];
rz(2.4300857) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(0.28212696) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.310815) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(0.34314439) q[0];
rz(-3.0798196) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(-1.7247346) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2767216) q[0];
sx q[0];
rz(-1.4007995) q[0];
sx q[0];
rz(-1.8525339) q[0];
rz(0.70059641) q[2];
sx q[2];
rz(-1.8237517) q[2];
sx q[2];
rz(-0.41159901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8843536) q[1];
sx q[1];
rz(-1.4352918) q[1];
sx q[1];
rz(1.0865023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9067578) q[3];
sx q[3];
rz(-1.9855301) q[3];
sx q[3];
rz(1.0339586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5938277) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(-2.055638) q[2];
rz(-0.91935277) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(-0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74349657) q[0];
sx q[0];
rz(-1.9323823) q[0];
sx q[0];
rz(-2.3486775) q[0];
rz(2.4010557) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(-0.83121306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75126264) q[0];
sx q[0];
rz(-2.7075276) q[0];
sx q[0];
rz(-2.084226) q[0];
rz(-pi) q[1];
rz(-1.5969876) q[2];
sx q[2];
rz(-0.27715836) q[2];
sx q[2];
rz(1.3326534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97396353) q[1];
sx q[1];
rz(-0.87461738) q[1];
sx q[1];
rz(0.40145282) q[1];
rz(1.5300203) q[3];
sx q[3];
rz(-2.6705461) q[3];
sx q[3];
rz(0.99179964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.071216019) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(-1.1191204) q[2];
rz(-1.8708771) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374461) q[0];
sx q[0];
rz(-0.16682145) q[0];
sx q[0];
rz(1.5995837) q[0];
rz(2.1265325) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(0.17280811) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240066) q[0];
sx q[0];
rz(-2.3476971) q[0];
sx q[0];
rz(2.2909597) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44266959) q[2];
sx q[2];
rz(-0.73390642) q[2];
sx q[2];
rz(1.1709605) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1041278) q[1];
sx q[1];
rz(-1.3224673) q[1];
sx q[1];
rz(-1.1315956) q[1];
x q[2];
rz(-1.0181581) q[3];
sx q[3];
rz(-0.88200906) q[3];
sx q[3];
rz(-1.8398726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66199866) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(0.97770989) q[2];
rz(-2.5665723) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(-1.2303111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32790312) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(-3.1258702) q[0];
rz(1.188259) q[1];
sx q[1];
rz(-2.5091722) q[1];
sx q[1];
rz(1.1994919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9415814) q[0];
sx q[0];
rz(-1.1304378) q[0];
sx q[0];
rz(2.8811127) q[0];
rz(-pi) q[1];
rz(-1.6686317) q[2];
sx q[2];
rz(-0.49287686) q[2];
sx q[2];
rz(-2.5081483) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4839061) q[1];
sx q[1];
rz(-0.71345854) q[1];
sx q[1];
rz(2.3785527) q[1];
x q[2];
rz(2.0301129) q[3];
sx q[3];
rz(-2.5190341) q[3];
sx q[3];
rz(-2.7186269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99379313) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(1.7564868) q[2];
rz(-2.5177453) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(-1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.1856336) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(-0.69325915) q[0];
rz(0.26501003) q[1];
sx q[1];
rz(-2.3131504) q[1];
sx q[1];
rz(1.5230491) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5361621) q[0];
sx q[0];
rz(-2.3443065) q[0];
sx q[0];
rz(1.5862341) q[0];
rz(-pi) q[1];
rz(-0.62317836) q[2];
sx q[2];
rz(-1.0747391) q[2];
sx q[2];
rz(-0.11962275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5793188) q[1];
sx q[1];
rz(-0.34046158) q[1];
sx q[1];
rz(-0.071694362) q[1];
x q[2];
rz(1.4632439) q[3];
sx q[3];
rz(-2.5578024) q[3];
sx q[3];
rz(-1.8369758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30274621) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(-1.6839074) q[2];
rz(2.1863106) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(-0.83520755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38462287) q[0];
sx q[0];
rz(-0.56461016) q[0];
sx q[0];
rz(2.6589174) q[0];
rz(-2.2118498) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(2.7453056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38448725) q[0];
sx q[0];
rz(-1.1570017) q[0];
sx q[0];
rz(1.9648457) q[0];
x q[1];
rz(2.9626289) q[2];
sx q[2];
rz(-1.71993) q[2];
sx q[2];
rz(2.3863132) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5037646) q[1];
sx q[1];
rz(-2.2866268) q[1];
sx q[1];
rz(-0.24055918) q[1];
rz(-1.9307053) q[3];
sx q[3];
rz(-1.9511838) q[3];
sx q[3];
rz(-1.4402657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3489909) q[2];
sx q[2];
rz(-1.4395809) q[2];
sx q[2];
rz(2.8144042) q[2];
rz(0.78091019) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(1.6646632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0384211) q[0];
sx q[0];
rz(-0.99089834) q[0];
sx q[0];
rz(0.25767576) q[0];
rz(-2.7720263) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(1.5255899) q[2];
sx q[2];
rz(-2.2635985) q[2];
sx q[2];
rz(0.23512693) q[2];
rz(1.0450324) q[3];
sx q[3];
rz(-0.67787328) q[3];
sx q[3];
rz(-1.9280435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
