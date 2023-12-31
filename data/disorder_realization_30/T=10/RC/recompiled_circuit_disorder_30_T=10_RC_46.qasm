OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(-0.39443031) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(-1.2000097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.278468) q[0];
sx q[0];
rz(-1.5651363) q[0];
sx q[0];
rz(-1.5229043) q[0];
x q[1];
rz(0.33485246) q[2];
sx q[2];
rz(-1.1455673) q[2];
sx q[2];
rz(1.2644757) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0389509) q[1];
sx q[1];
rz(-2.0248374) q[1];
sx q[1];
rz(2.7137043) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0311534) q[3];
sx q[3];
rz(-0.38023284) q[3];
sx q[3];
rz(-1.5649753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3068984) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-2.7976024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(1.3551691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25027572) q[0];
sx q[0];
rz(-0.96069562) q[0];
sx q[0];
rz(3.0623869) q[0];
rz(-2.9789574) q[2];
sx q[2];
rz(-1.4504823) q[2];
sx q[2];
rz(-1.6441117) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.558074) q[1];
sx q[1];
rz(-1.1221702) q[1];
sx q[1];
rz(-3.0706057) q[1];
x q[2];
rz(-1.1739028) q[3];
sx q[3];
rz(-2.9429884) q[3];
sx q[3];
rz(-1.7484776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051620313) q[0];
sx q[0];
rz(-1.211477) q[0];
sx q[0];
rz(3.1348455) q[0];
rz(2.2674019) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(3.0490321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6667337) q[1];
sx q[1];
rz(-1.6068646) q[1];
sx q[1];
rz(-2.1875847) q[1];
x q[2];
rz(0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.76434) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(-0.68850368) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(-0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(3.1304741) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5866833) q[0];
sx q[0];
rz(-2.0218098) q[0];
sx q[0];
rz(0.45278544) q[0];
x q[1];
rz(0.83799329) q[2];
sx q[2];
rz(-2.0760771) q[2];
sx q[2];
rz(-2.879564) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0032469) q[1];
sx q[1];
rz(-1.0838638) q[1];
sx q[1];
rz(0.37133118) q[1];
rz(-0.53400455) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(-0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6461688) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(-2.4781573) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.3269075) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715695) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(-0.75385401) q[0];
rz(-pi) q[1];
rz(2.3438498) q[2];
sx q[2];
rz(-1.1290871) q[2];
sx q[2];
rz(-1.9517348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8618968) q[1];
sx q[1];
rz(-1.3394636) q[1];
sx q[1];
rz(0.33613236) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1760686) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(0.6814878) q[0];
rz(2.9340414) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(-1.1157657) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9003446) q[0];
sx q[0];
rz(-1.9755409) q[0];
sx q[0];
rz(2.3099398) q[0];
rz(-0.43004604) q[2];
sx q[2];
rz(-2.1950245) q[2];
sx q[2];
rz(2.4732694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0259077) q[1];
sx q[1];
rz(-2.5528862) q[1];
sx q[1];
rz(3.1097417) q[1];
rz(-2.9168105) q[3];
sx q[3];
rz(-0.68192476) q[3];
sx q[3];
rz(-3.0608321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21268022) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(0.3195233) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(-2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(-2.0293503) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(2.5792714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599265) q[0];
sx q[0];
rz(-0.66837464) q[0];
sx q[0];
rz(-1.5947123) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.855905) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(-3.0722741) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6866236) q[1];
sx q[1];
rz(-1.7640055) q[1];
sx q[1];
rz(-1.7186233) q[1];
rz(-pi) q[2];
rz(3.0543442) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(-1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(0.34004655) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44889221) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(0.39644077) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.5213535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0884468) q[0];
sx q[0];
rz(-1.7674812) q[0];
sx q[0];
rz(1.6403857) q[0];
rz(0.40201681) q[2];
sx q[2];
rz(-1.1321628) q[2];
sx q[2];
rz(1.5356262) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8926881) q[1];
sx q[1];
rz(-0.78033328) q[1];
sx q[1];
rz(-0.08463879) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.429854) q[3];
sx q[3];
rz(-2.4574276) q[3];
sx q[3];
rz(-1.9977026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5788995) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(1.1307905) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-2.14595) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(-2.8709581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7662738) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(2.4256698) q[0];
rz(-pi) q[1];
rz(0.80910271) q[2];
sx q[2];
rz(-1.9168233) q[2];
sx q[2];
rz(1.5835294) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2808025) q[1];
sx q[1];
rz(-1.5250912) q[1];
sx q[1];
rz(1.526282) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7089013) q[3];
sx q[3];
rz(-2.0376251) q[3];
sx q[3];
rz(-0.69463581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-0.45544004) q[2];
rz(2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(0.49490067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1586944) q[0];
sx q[0];
rz(-1.6086846) q[0];
sx q[0];
rz(2.057103) q[0];
rz(-pi) q[1];
rz(-2.8137389) q[2];
sx q[2];
rz(-1.5640537) q[2];
sx q[2];
rz(-1.9890131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7086664) q[1];
sx q[1];
rz(-0.64879829) q[1];
sx q[1];
rz(-1.7556346) q[1];
rz(-0.095330843) q[3];
sx q[3];
rz(-0.87626002) q[3];
sx q[3];
rz(1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(-2.1081934) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0192169) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(2.3090251) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-0.20559786) q[2];
sx q[2];
rz(-2.1046706) q[2];
sx q[2];
rz(-0.54866366) q[2];
rz(2.8752747) q[3];
sx q[3];
rz(-1.3238293) q[3];
sx q[3];
rz(-2.7662591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
