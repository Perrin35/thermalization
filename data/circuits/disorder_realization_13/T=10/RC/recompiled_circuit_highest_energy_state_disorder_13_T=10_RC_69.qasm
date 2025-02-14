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
rz(-2.4969555) q[0];
sx q[0];
rz(-0.58965373) q[0];
sx q[0];
rz(-1.0876422) q[0];
rz(3.1261858) q[1];
sx q[1];
rz(-0.3897804) q[1];
sx q[1];
rz(1.1642233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7925568) q[0];
sx q[0];
rz(-1.4407939) q[0];
sx q[0];
rz(-3.1368238) q[0];
rz(-0.40428031) q[2];
sx q[2];
rz(-2.4089185) q[2];
sx q[2];
rz(2.3892185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.30689209) q[1];
sx q[1];
rz(-2.0880719) q[1];
sx q[1];
rz(0.24637988) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4527997) q[3];
sx q[3];
rz(-1.495365) q[3];
sx q[3];
rz(-0.82619595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9001793) q[2];
sx q[2];
rz(-0.56168491) q[2];
sx q[2];
rz(-2.4955595) q[2];
rz(-2.8863886) q[3];
sx q[3];
rz(-0.45707688) q[3];
sx q[3];
rz(-1.7469762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141041) q[0];
sx q[0];
rz(-0.33559594) q[0];
sx q[0];
rz(-2.3345729) q[0];
rz(-1.6837766) q[1];
sx q[1];
rz(-0.57676637) q[1];
sx q[1];
rz(1.797765) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80088491) q[0];
sx q[0];
rz(-1.721764) q[0];
sx q[0];
rz(-0.29125352) q[0];
rz(-pi) q[1];
rz(2.1704111) q[2];
sx q[2];
rz(-2.4840725) q[2];
sx q[2];
rz(-2.2258141) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8763153) q[1];
sx q[1];
rz(-0.96926276) q[1];
sx q[1];
rz(-2.5391891) q[1];
rz(-pi) q[2];
rz(1.5206512) q[3];
sx q[3];
rz(-0.72948958) q[3];
sx q[3];
rz(1.6901085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1809711) q[2];
sx q[2];
rz(-1.0701067) q[2];
sx q[2];
rz(-0.19372678) q[2];
rz(1.6173897) q[3];
sx q[3];
rz(-0.75810713) q[3];
sx q[3];
rz(-2.962842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6206361) q[0];
sx q[0];
rz(-0.23796029) q[0];
sx q[0];
rz(-0.97610193) q[0];
rz(-2.2528265) q[1];
sx q[1];
rz(-0.74405324) q[1];
sx q[1];
rz(2.9685453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.219653) q[0];
sx q[0];
rz(-1.573607) q[0];
sx q[0];
rz(-1.5332743) q[0];
x q[1];
rz(-1.6962181) q[2];
sx q[2];
rz(-2.0880155) q[2];
sx q[2];
rz(2.7917535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.58223078) q[1];
sx q[1];
rz(-1.4051062) q[1];
sx q[1];
rz(-0.5104602) q[1];
rz(-pi) q[2];
rz(0.32239989) q[3];
sx q[3];
rz(-1.4677257) q[3];
sx q[3];
rz(0.88503982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0024857) q[2];
sx q[2];
rz(-1.8128914) q[2];
sx q[2];
rz(2.4277182) q[2];
rz(-0.21126963) q[3];
sx q[3];
rz(-1.0372838) q[3];
sx q[3];
rz(-2.5762288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7738889) q[0];
sx q[0];
rz(-1.1834894) q[0];
sx q[0];
rz(2.0088038) q[0];
rz(-0.47138131) q[1];
sx q[1];
rz(-1.2939204) q[1];
sx q[1];
rz(0.43101355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16551183) q[0];
sx q[0];
rz(-1.3595194) q[0];
sx q[0];
rz(-1.6863281) q[0];
rz(-2.3296146) q[2];
sx q[2];
rz(-0.50305191) q[2];
sx q[2];
rz(0.72847086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20038988) q[1];
sx q[1];
rz(-1.2697436) q[1];
sx q[1];
rz(-2.5832504) q[1];
x q[2];
rz(-2.1924344) q[3];
sx q[3];
rz(-1.3484869) q[3];
sx q[3];
rz(2.3883826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5203349) q[2];
sx q[2];
rz(-0.73953491) q[2];
sx q[2];
rz(0.4062824) q[2];
rz(1.9351561) q[3];
sx q[3];
rz(-2.0483569) q[3];
sx q[3];
rz(-2.5549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63221145) q[0];
sx q[0];
rz(-0.87566942) q[0];
sx q[0];
rz(2.8152554) q[0];
rz(-0.95411602) q[1];
sx q[1];
rz(-0.64847821) q[1];
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
rz(0.40753741) q[0];
x q[1];
rz(-0.65803501) q[2];
sx q[2];
rz(-1.3172704) q[2];
sx q[2];
rz(0.70632284) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.919173) q[1];
sx q[1];
rz(-1.2856163) q[1];
sx q[1];
rz(1.213277) q[1];
x q[2];
rz(2.8460449) q[3];
sx q[3];
rz(-1.4369399) q[3];
sx q[3];
rz(-2.9791711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.894459) q[2];
sx q[2];
rz(-0.45866141) q[2];
sx q[2];
rz(2.4776754) q[2];
rz(-2.2513921) q[3];
sx q[3];
rz(-1.5492487) q[3];
sx q[3];
rz(0.012705407) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3506055) q[0];
sx q[0];
rz(-2.7101639) q[0];
sx q[0];
rz(0.33082333) q[0];
rz(2.779003) q[1];
sx q[1];
rz(-1.9622842) q[1];
sx q[1];
rz(-0.51574743) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12130485) q[0];
sx q[0];
rz(-2.8175531) q[0];
sx q[0];
rz(2.3100467) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3810034) q[2];
sx q[2];
rz(-2.1564473) q[2];
sx q[2];
rz(0.3444258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52236103) q[1];
sx q[1];
rz(-0.69612487) q[1];
sx q[1];
rz(-6*pi/13) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79591085) q[3];
sx q[3];
rz(-1.6636968) q[3];
sx q[3];
rz(1.9167629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35046878) q[2];
sx q[2];
rz(-1.4501269) q[2];
sx q[2];
rz(2.3340732) q[2];
rz(-0.098585419) q[3];
sx q[3];
rz(-0.8980631) q[3];
sx q[3];
rz(2.0487823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43807855) q[0];
sx q[0];
rz(-0.57323891) q[0];
sx q[0];
rz(2.692063) q[0];
rz(-0.685177) q[1];
sx q[1];
rz(-0.40760577) q[1];
sx q[1];
rz(2.5221241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9307053) q[0];
sx q[0];
rz(-2.8307417) q[0];
sx q[0];
rz(0.092664552) q[0];
rz(-1.5225822) q[2];
sx q[2];
rz(-1.0430665) q[2];
sx q[2];
rz(2.9121648) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56726914) q[1];
sx q[1];
rz(-0.78227121) q[1];
sx q[1];
rz(0.48552728) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3024019) q[3];
sx q[3];
rz(-0.27264914) q[3];
sx q[3];
rz(-2.1781937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56457907) q[2];
sx q[2];
rz(-1.4489633) q[2];
sx q[2];
rz(-2.9435834) q[2];
rz(-1.7791344) q[3];
sx q[3];
rz(-2.8415316) q[3];
sx q[3];
rz(2.4411966) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.144416) q[0];
sx q[0];
rz(-2.5197025) q[0];
sx q[0];
rz(1.8560334) q[0];
rz(-0.27664912) q[1];
sx q[1];
rz(-1.7729019) q[1];
sx q[1];
rz(-0.10861529) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40688694) q[0];
sx q[0];
rz(-2.6331775) q[0];
sx q[0];
rz(-1.4020573) q[0];
rz(-3.0622903) q[2];
sx q[2];
rz(-0.9482884) q[2];
sx q[2];
rz(0.17891893) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.120879) q[1];
sx q[1];
rz(-0.87057897) q[1];
sx q[1];
rz(2.8597699) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1843508) q[3];
sx q[3];
rz(-1.4542011) q[3];
sx q[3];
rz(-2.1045096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8331929) q[2];
sx q[2];
rz(-1.9449642) q[2];
sx q[2];
rz(-1.2742554) q[2];
rz(1.3537004) q[3];
sx q[3];
rz(-2.7022868) q[3];
sx q[3];
rz(2.275009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9891147) q[0];
sx q[0];
rz(-2.415933) q[0];
sx q[0];
rz(3.0210378) q[0];
rz(-0.74147916) q[1];
sx q[1];
rz(-1.9375786) q[1];
sx q[1];
rz(0.075686879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73629433) q[0];
sx q[0];
rz(-1.8764623) q[0];
sx q[0];
rz(1.7429245) q[0];
rz(-0.077421256) q[2];
sx q[2];
rz(-1.110437) q[2];
sx q[2];
rz(1.3321239) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32966742) q[1];
sx q[1];
rz(-2.4849548) q[1];
sx q[1];
rz(-0.041402264) q[1];
rz(-pi) q[2];
rz(-2.5406557) q[3];
sx q[3];
rz(-1.4190201) q[3];
sx q[3];
rz(-2.5545718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0691444) q[2];
sx q[2];
rz(-2.2587903) q[2];
sx q[2];
rz(-2.7952588) q[2];
rz(-0.94541466) q[3];
sx q[3];
rz(-1.3225222) q[3];
sx q[3];
rz(0.94714981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866975) q[0];
sx q[0];
rz(-1.1981542) q[0];
sx q[0];
rz(2.6937038) q[0];
rz(2.2121494) q[1];
sx q[1];
rz(-1.7421236) q[1];
sx q[1];
rz(-2.4670752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69792992) q[0];
sx q[0];
rz(-2.3259458) q[0];
sx q[0];
rz(-1.2944512) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3849887) q[2];
sx q[2];
rz(-1.6730434) q[2];
sx q[2];
rz(-0.24413689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0688388) q[1];
sx q[1];
rz(-1.7441894) q[1];
sx q[1];
rz(0.0011464107) q[1];
rz(-1.7389948) q[3];
sx q[3];
rz(-2.7175596) q[3];
sx q[3];
rz(0.096120983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21005361) q[2];
sx q[2];
rz(-1.9320107) q[2];
sx q[2];
rz(-0.26739576) q[2];
rz(-2.0343272) q[3];
sx q[3];
rz(-2.3591154) q[3];
sx q[3];
rz(0.84899181) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4103107) q[0];
sx q[0];
rz(-1.622643) q[0];
sx q[0];
rz(-1.0868764) q[0];
rz(0.80401737) q[1];
sx q[1];
rz(-1.6537279) q[1];
sx q[1];
rz(-2.0860685) q[1];
rz(3.0699365) q[2];
sx q[2];
rz(-1.4240257) q[2];
sx q[2];
rz(1.7113513) q[2];
rz(1.4645529) q[3];
sx q[3];
rz(-0.67062702) q[3];
sx q[3];
rz(0.39556243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
