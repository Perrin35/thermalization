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
rz(0.64463717) q[0];
sx q[0];
rz(3.7312464) q[0];
sx q[0];
rz(10.51242) q[0];
rz(-0.015406869) q[1];
sx q[1];
rz(-2.7518123) q[1];
sx q[1];
rz(1.9773693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.919214) q[0];
sx q[0];
rz(-1.5755249) q[0];
sx q[0];
rz(1.7008002) q[0];
rz(-pi) q[1];
rz(2.4504775) q[2];
sx q[2];
rz(-1.8370312) q[2];
sx q[2];
rz(-2.63111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1401745) q[1];
sx q[1];
rz(-1.3571897) q[1];
sx q[1];
rz(1.0402388) q[1];
rz(-1.6887929) q[3];
sx q[3];
rz(-1.495365) q[3];
sx q[3];
rz(2.3153967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9001793) q[2];
sx q[2];
rz(-0.56168491) q[2];
sx q[2];
rz(2.4955595) q[2];
rz(0.25520405) q[3];
sx q[3];
rz(-2.6845158) q[3];
sx q[3];
rz(1.7469762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141041) q[0];
sx q[0];
rz(-0.33559594) q[0];
sx q[0];
rz(-2.3345729) q[0];
rz(-1.457816) q[1];
sx q[1];
rz(-2.5648263) q[1];
sx q[1];
rz(-1.3438276) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9066853) q[0];
sx q[0];
rz(-0.32706383) q[0];
sx q[0];
rz(0.48717888) q[0];
x q[1];
rz(1.0033002) q[2];
sx q[2];
rz(-1.9229182) q[2];
sx q[2];
rz(-2.9825236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7585425) q[1];
sx q[1];
rz(-0.82393007) q[1];
sx q[1];
rz(-2.2609018) q[1];
rz(-pi) q[2];
x q[2];
rz(2.299661) q[3];
sx q[3];
rz(-1.6042097) q[3];
sx q[3];
rz(0.15671003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1809711) q[2];
sx q[2];
rz(-2.0714859) q[2];
sx q[2];
rz(0.19372678) q[2];
rz(-1.5242029) q[3];
sx q[3];
rz(-0.75810713) q[3];
sx q[3];
rz(-2.962842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6206361) q[0];
sx q[0];
rz(-0.23796029) q[0];
sx q[0];
rz(-0.97610193) q[0];
rz(0.8887662) q[1];
sx q[1];
rz(-0.74405324) q[1];
sx q[1];
rz(-0.1730473) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3510378) q[0];
sx q[0];
rz(-1.5332744) q[0];
sx q[0];
rz(-3.13878) q[0];
rz(-pi) q[1];
rz(-2.9251418) q[2];
sx q[2];
rz(-2.6107288) q[2];
sx q[2];
rz(2.5420839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4394697) q[1];
sx q[1];
rz(-2.6071848) q[1];
sx q[1];
rz(-0.32976697) q[1];
x q[2];
rz(-0.31555303) q[3];
sx q[3];
rz(-2.8036661) q[3];
sx q[3];
rz(-0.38697836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.139107) q[2];
sx q[2];
rz(-1.8128914) q[2];
sx q[2];
rz(-2.4277182) q[2];
rz(2.930323) q[3];
sx q[3];
rz(-2.1043089) q[3];
sx q[3];
rz(-0.56536388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7738889) q[0];
sx q[0];
rz(-1.1834894) q[0];
sx q[0];
rz(-1.1327889) q[0];
rz(2.6702113) q[1];
sx q[1];
rz(-1.2939204) q[1];
sx q[1];
rz(0.43101355) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9760808) q[0];
sx q[0];
rz(-1.3595194) q[0];
sx q[0];
rz(-1.6863281) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9507031) q[2];
sx q[2];
rz(-1.2326733) q[2];
sx q[2];
rz(0.1491216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1872594) q[1];
sx q[1];
rz(-2.1012839) q[1];
sx q[1];
rz(-1.2198636) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1924344) q[3];
sx q[3];
rz(-1.7931058) q[3];
sx q[3];
rz(0.75321001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6212578) q[2];
sx q[2];
rz(-0.73953491) q[2];
sx q[2];
rz(0.4062824) q[2];
rz(1.2064365) q[3];
sx q[3];
rz(-2.0483569) q[3];
sx q[3];
rz(2.5549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5093812) q[0];
sx q[0];
rz(-2.2659232) q[0];
sx q[0];
rz(-0.32633728) q[0];
rz(-2.1874766) q[1];
sx q[1];
rz(-2.4931144) q[1];
sx q[1];
rz(-3.1180678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3259567) q[0];
sx q[0];
rz(-1.1679497) q[0];
sx q[0];
rz(-2.7340552) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65803501) q[2];
sx q[2];
rz(-1.3172704) q[2];
sx q[2];
rz(2.4352698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70290804) q[1];
sx q[1];
rz(-2.6881235) q[1];
sx q[1];
rz(-2.2681342) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.708528) q[3];
sx q[3];
rz(-0.32363656) q[3];
sx q[3];
rz(1.8214846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24713369) q[2];
sx q[2];
rz(-0.45866141) q[2];
sx q[2];
rz(2.4776754) q[2];
rz(-0.89020056) q[3];
sx q[3];
rz(-1.592344) q[3];
sx q[3];
rz(-3.1288872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3506055) q[0];
sx q[0];
rz(-2.7101639) q[0];
sx q[0];
rz(-2.8107693) q[0];
rz(-0.36258969) q[1];
sx q[1];
rz(-1.9622842) q[1];
sx q[1];
rz(-0.51574743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88721758) q[0];
sx q[0];
rz(-1.8083113) q[0];
sx q[0];
rz(-0.22254469) q[0];
rz(-2.8644805) q[2];
sx q[2];
rz(-2.5293969) q[2];
sx q[2];
rz(-3.1316568) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6192316) q[1];
sx q[1];
rz(-0.69612487) q[1];
sx q[1];
rz(6*pi/13) q[1];
x q[2];
rz(2.3456818) q[3];
sx q[3];
rz(-1.6636968) q[3];
sx q[3];
rz(-1.2248298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7911239) q[2];
sx q[2];
rz(-1.4501269) q[2];
sx q[2];
rz(-0.80751944) q[2];
rz(3.0430072) q[3];
sx q[3];
rz(-2.2435296) q[3];
sx q[3];
rz(1.0928104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7035141) q[0];
sx q[0];
rz(-2.5683537) q[0];
sx q[0];
rz(2.692063) q[0];
rz(2.4564157) q[1];
sx q[1];
rz(-0.40760577) q[1];
sx q[1];
rz(2.5221241) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44815608) q[0];
sx q[0];
rz(-1.5991028) q[0];
sx q[0];
rz(-2.8319915) q[0];
rz(-0.52823587) q[2];
sx q[2];
rz(-1.5291457) q[2];
sx q[2];
rz(1.3170751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5743235) q[1];
sx q[1];
rz(-0.78227121) q[1];
sx q[1];
rz(-2.6560654) q[1];
x q[2];
rz(-1.3656627) q[3];
sx q[3];
rz(-1.7516802) q[3];
sx q[3];
rz(0.10553979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56457907) q[2];
sx q[2];
rz(-1.4489633) q[2];
sx q[2];
rz(-2.9435834) q[2];
rz(-1.7791344) q[3];
sx q[3];
rz(-0.30006108) q[3];
sx q[3];
rz(0.70039606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9971767) q[0];
sx q[0];
rz(-0.62189019) q[0];
sx q[0];
rz(-1.2855592) q[0];
rz(0.27664912) q[1];
sx q[1];
rz(-1.3686907) q[1];
sx q[1];
rz(-0.10861529) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1254107) q[0];
sx q[0];
rz(-1.6526395) q[0];
sx q[0];
rz(-1.0684408) q[0];
rz(3.0622903) q[2];
sx q[2];
rz(-0.9482884) q[2];
sx q[2];
rz(2.9626737) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7759571) q[1];
sx q[1];
rz(-1.7851019) q[1];
sx q[1];
rz(2.2909095) q[1];
rz(-1.2694938) q[3];
sx q[3];
rz(-2.7387894) q[3];
sx q[3];
rz(2.3294152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30839977) q[2];
sx q[2];
rz(-1.1966285) q[2];
sx q[2];
rz(1.2742554) q[2];
rz(-1.7878923) q[3];
sx q[3];
rz(-2.7022868) q[3];
sx q[3];
rz(-0.86658365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15247791) q[0];
sx q[0];
rz(-0.72565961) q[0];
sx q[0];
rz(3.0210378) q[0];
rz(2.4001135) q[1];
sx q[1];
rz(-1.2040141) q[1];
sx q[1];
rz(-0.075686879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78223373) q[0];
sx q[0];
rz(-1.4067211) q[0];
sx q[0];
rz(2.8316281) q[0];
x q[1];
rz(1.7255177) q[2];
sx q[2];
rz(-2.675229) q[2];
sx q[2];
rz(1.6365964) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8119252) q[1];
sx q[1];
rz(-0.65663785) q[1];
sx q[1];
rz(3.1001904) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5406557) q[3];
sx q[3];
rz(-1.7225725) q[3];
sx q[3];
rz(2.5545718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0691444) q[2];
sx q[2];
rz(-2.2587903) q[2];
sx q[2];
rz(2.7952588) q[2];
rz(0.94541466) q[3];
sx q[3];
rz(-1.3225222) q[3];
sx q[3];
rz(2.1944428) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866975) q[0];
sx q[0];
rz(-1.9434384) q[0];
sx q[0];
rz(0.44788885) q[0];
rz(2.2121494) q[1];
sx q[1];
rz(-1.7421236) q[1];
sx q[1];
rz(0.67451745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69792992) q[0];
sx q[0];
rz(-2.3259458) q[0];
sx q[0];
rz(-1.8471414) q[0];
rz(-2.0777763) q[2];
sx q[2];
rz(-0.21179767) q[2];
sx q[2];
rz(1.3174881) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50215534) q[1];
sx q[1];
rz(-1.5696671) q[1];
sx q[1];
rz(-1.3974031) q[1];
rz(-pi) q[2];
rz(1.4025979) q[3];
sx q[3];
rz(-0.4240331) q[3];
sx q[3];
rz(-0.096120983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21005361) q[2];
sx q[2];
rz(-1.209582) q[2];
sx q[2];
rz(-2.8741969) q[2];
rz(-1.1072655) q[3];
sx q[3];
rz(-0.78247726) q[3];
sx q[3];
rz(-2.2926008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4103107) q[0];
sx q[0];
rz(-1.5189497) q[0];
sx q[0];
rz(2.0547163) q[0];
rz(0.80401737) q[1];
sx q[1];
rz(-1.6537279) q[1];
sx q[1];
rz(-2.0860685) q[1];
rz(-1.7179391) q[2];
sx q[2];
rz(-1.6416807) q[2];
sx q[2];
rz(0.15105187) q[2];
rz(-2.238672) q[3];
sx q[3];
rz(-1.5048448) q[3];
sx q[3];
rz(2.0497143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
