OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(1.3786432) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8060018) q[0];
sx q[0];
rz(-1.5854892) q[0];
sx q[0];
rz(-0.088784119) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11159201) q[2];
sx q[2];
rz(-2.0473695) q[2];
sx q[2];
rz(-0.0026207844) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8494107) q[1];
sx q[1];
rz(-1.0725478) q[1];
sx q[1];
rz(-1.8241747) q[1];
x q[2];
rz(-2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(2.297304) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(2.8785008) q[0];
rz(-0.94353765) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.1862322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037558) q[0];
sx q[0];
rz(-0.058996011) q[0];
sx q[0];
rz(-2.8179413) q[0];
x q[1];
rz(2.8430014) q[2];
sx q[2];
rz(-0.208374) q[2];
sx q[2];
rz(1.5339472) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9895049) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(1.3701887) q[1];
rz(1.7632742) q[3];
sx q[3];
rz(-2.2566183) q[3];
sx q[3];
rz(-1.457085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0120323) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(1.9821232) q[2];
rz(-0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(-2.3348715) q[0];
rz(2.9280248) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(0.82021964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41528156) q[0];
sx q[0];
rz(-1.4584686) q[0];
sx q[0];
rz(0.88322722) q[0];
x q[1];
rz(1.1782896) q[2];
sx q[2];
rz(-0.5272534) q[2];
sx q[2];
rz(-0.96780992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.080938235) q[1];
sx q[1];
rz(-1.9229691) q[1];
sx q[1];
rz(-1.1206131) q[1];
rz(-pi) q[2];
rz(-3.0235602) q[3];
sx q[3];
rz(-1.1088088) q[3];
sx q[3];
rz(-3.1363917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(-2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(-0.91127515) q[0];
rz(2.7032734) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(-1.320425) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60162773) q[0];
sx q[0];
rz(-2.7318582) q[0];
sx q[0];
rz(-1.5448624) q[0];
x q[1];
rz(-2.7337012) q[2];
sx q[2];
rz(-2.6558999) q[2];
sx q[2];
rz(-2.8913468) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1530694) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(-0.70446976) q[1];
x q[2];
rz(2.3539691) q[3];
sx q[3];
rz(-2.0260603) q[3];
sx q[3];
rz(2.5256707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300366) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(-1.226549) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.3006166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9601701) q[0];
sx q[0];
rz(-1.2167861) q[0];
sx q[0];
rz(-1.5439073) q[0];
rz(-pi) q[1];
rz(-1.6426716) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(2.9938811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.865766) q[1];
sx q[1];
rz(-1.8991125) q[1];
sx q[1];
rz(2.5715716) q[1];
rz(-pi) q[2];
rz(1.424765) q[3];
sx q[3];
rz(-1.8304123) q[3];
sx q[3];
rz(2.8433593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(2.9495083) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(0.011750301) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(1.5884429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9827305) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(1.1984675) q[0];
rz(-pi) q[1];
rz(-1.9366829) q[2];
sx q[2];
rz(-1.5503746) q[2];
sx q[2];
rz(-0.4991971) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3192056) q[1];
sx q[1];
rz(-0.86383312) q[1];
sx q[1];
rz(1.6824526) q[1];
rz(-3.0670777) q[3];
sx q[3];
rz(-2.2722368) q[3];
sx q[3];
rz(2.6889192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(-1.7717308) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5605374) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(0.67725956) q[0];
rz(2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(-2.1645434) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517294) q[0];
sx q[0];
rz(-1.4954733) q[0];
sx q[0];
rz(0.59318869) q[0];
rz(-pi) q[1];
rz(-3.0095519) q[2];
sx q[2];
rz(-1.5706976) q[2];
sx q[2];
rz(1.6824818) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.44927412) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(1.3432137) q[1];
rz(0.3427152) q[3];
sx q[3];
rz(-1.2667155) q[3];
sx q[3];
rz(-2.8976687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7523505) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(-0.69865984) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218924) q[0];
sx q[0];
rz(-1.3601174) q[0];
sx q[0];
rz(1.6181437) q[0];
x q[1];
rz(-2.2970389) q[2];
sx q[2];
rz(-0.65659467) q[2];
sx q[2];
rz(-0.086364634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.318995) q[1];
sx q[1];
rz(-1.5929475) q[1];
sx q[1];
rz(-1.5763361) q[1];
rz(0.28835339) q[3];
sx q[3];
rz(-0.95803146) q[3];
sx q[3];
rz(1.310865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(1.684749) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(-2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.8485803) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(0.59757772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1406527) q[0];
sx q[0];
rz(-1.6093328) q[0];
sx q[0];
rz(1.4356104) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.045268) q[2];
sx q[2];
rz(-1.0768441) q[2];
sx q[2];
rz(-1.8906821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24039195) q[1];
sx q[1];
rz(-1.8933834) q[1];
sx q[1];
rz(-2.8278973) q[1];
rz(-pi) q[2];
rz(-3.0258614) q[3];
sx q[3];
rz(-0.58905187) q[3];
sx q[3];
rz(-0.91050402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22275816) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(-1.2333599) q[2];
rz(-0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(2.4694209) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2490847) q[0];
sx q[0];
rz(-1.5810284) q[0];
sx q[0];
rz(-0.0052878629) q[0];
x q[1];
rz(1.5800843) q[2];
sx q[2];
rz(-2.0822968) q[2];
sx q[2];
rz(1.7363422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2939261) q[1];
sx q[1];
rz(-2.3617509) q[1];
sx q[1];
rz(0.49077175) q[1];
x q[2];
rz(2.7087595) q[3];
sx q[3];
rz(-1.7107309) q[3];
sx q[3];
rz(2.7452552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621915) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-0.72369408) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(1.1719218) q[2];
sx q[2];
rz(-1.4234067) q[2];
sx q[2];
rz(1.1248551) q[2];
rz(2.1914992) q[3];
sx q[3];
rz(-2.6549669) q[3];
sx q[3];
rz(-2.0012729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
