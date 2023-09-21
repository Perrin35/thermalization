OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0177512) q[0];
sx q[0];
rz(-3.0525644) q[0];
sx q[0];
rz(2.428399) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1630152) q[2];
sx q[2];
rz(-0.41523146) q[2];
sx q[2];
rz(-0.65939553) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7072304) q[1];
sx q[1];
rz(-2.011236) q[1];
sx q[1];
rz(0.67727725) q[1];
rz(-0.37791032) q[3];
sx q[3];
rz(-2.1034735) q[3];
sx q[3];
rz(2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(2.0430298) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-2.0400955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(1.4651441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629048) q[0];
sx q[0];
rz(-1.573274) q[0];
sx q[0];
rz(-2.4173827) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0287839) q[2];
sx q[2];
rz(-2.6070242) q[2];
sx q[2];
rz(2.2528258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.022545594) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(2.4382298) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74554262) q[3];
sx q[3];
rz(-1.8288444) q[3];
sx q[3];
rz(-2.6170956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96317545) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(0.80054545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313989) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(-1.9171159) q[0];
rz(-pi) q[1];
rz(0.82528798) q[2];
sx q[2];
rz(-2.2824259) q[2];
sx q[2];
rz(-2.3724144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70558207) q[1];
sx q[1];
rz(-1.5758621) q[1];
sx q[1];
rz(2.7286121) q[1];
rz(-2.27182) q[3];
sx q[3];
rz(-1.3372034) q[3];
sx q[3];
rz(1.1886532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(2.676679) q[0];
rz(2.7930296) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(-2.0565313) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59180075) q[0];
sx q[0];
rz(-1.340197) q[0];
sx q[0];
rz(2.6016298) q[0];
x q[1];
rz(0.98767878) q[2];
sx q[2];
rz(-0.43118011) q[2];
sx q[2];
rz(2.0476066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9492053) q[1];
sx q[1];
rz(-0.17377033) q[1];
sx q[1];
rz(0.6178356) q[1];
x q[2];
rz(1.9784847) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(2.98416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(2.4604649) q[2];
rz(0.51182169) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78445804) q[0];
sx q[0];
rz(-1.1103837) q[0];
sx q[0];
rz(2.5255192) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.328674) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(2.9850933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86451929) q[1];
sx q[1];
rz(-0.48950567) q[1];
sx q[1];
rz(-0.31788748) q[1];
x q[2];
rz(-0.97335191) q[3];
sx q[3];
rz(-1.5095599) q[3];
sx q[3];
rz(-0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-0.16528027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0342456) q[0];
sx q[0];
rz(-3.0421071) q[0];
sx q[0];
rz(-0.29883595) q[0];
x q[1];
rz(0.084947305) q[2];
sx q[2];
rz(-2.4521378) q[2];
sx q[2];
rz(-1.8817608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2416934) q[1];
sx q[1];
rz(-2.4447933) q[1];
sx q[1];
rz(2.4844869) q[1];
rz(-1.1214439) q[3];
sx q[3];
rz(-1.7835788) q[3];
sx q[3];
rz(1.7920115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91339397) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(-0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3570324) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1897141) q[0];
sx q[0];
rz(-1.6197546) q[0];
sx q[0];
rz(0.017107054) q[0];
rz(-1.9095124) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(2.2155857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10662096) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(-2.3751395) q[1];
x q[2];
rz(-0.69865366) q[3];
sx q[3];
rz(-0.62056382) q[3];
sx q[3];
rz(-2.8750949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-0.53156701) q[2];
rz(2.452204) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.2954856) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078995973) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-2.2198026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4066276) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(-2.448003) q[0];
rz(-pi) q[1];
rz(2.8554101) q[2];
sx q[2];
rz(-2.2879062) q[2];
sx q[2];
rz(1.66695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8340048) q[1];
sx q[1];
rz(-1.7760135) q[1];
sx q[1];
rz(-2.8819041) q[1];
x q[2];
rz(-0.86254085) q[3];
sx q[3];
rz(-1.7077912) q[3];
sx q[3];
rz(-2.9555637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(1.9699338) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(1.6660447) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(2.4618861) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70521351) q[0];
sx q[0];
rz(-2.0853015) q[0];
sx q[0];
rz(-1.3374469) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3557415) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(-3.0422473) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6913773) q[1];
sx q[1];
rz(-0.9141578) q[1];
sx q[1];
rz(-1.1378098) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71984843) q[3];
sx q[3];
rz(-0.82092972) q[3];
sx q[3];
rz(-1.5035226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1264964) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.05242059) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(1.1766599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.048621) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(-0.049157814) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87877019) q[2];
sx q[2];
rz(-2.3792017) q[2];
sx q[2];
rz(2.9825485) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1365876) q[1];
sx q[1];
rz(-1.0141546) q[1];
sx q[1];
rz(0.088076061) q[1];
x q[2];
rz(-1.482974) q[3];
sx q[3];
rz(-1.0613958) q[3];
sx q[3];
rz(-1.0305962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(0.39623109) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(-2.836543) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(0.35076326) q[3];
sx q[3];
rz(-1.5630994) q[3];
sx q[3];
rz(0.85170436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];