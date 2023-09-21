OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(0.18145951) q[0];
rz(2.060086) q[1];
sx q[1];
rz(-0.67343155) q[1];
sx q[1];
rz(-1.0531309) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3436514) q[0];
sx q[0];
rz(-0.5125674) q[0];
sx q[0];
rz(-2.0128065) q[0];
rz(0.40197576) q[2];
sx q[2];
rz(-0.66265124) q[2];
sx q[2];
rz(2.0768009) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.83208132) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(1.8413715) q[1];
rz(-pi) q[2];
rz(-2.9795322) q[3];
sx q[3];
rz(-1.2352304) q[3];
sx q[3];
rz(-0.1048564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9782605) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(-2.1286428) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098009) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(-2.5464771) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.5140623) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94905108) q[0];
sx q[0];
rz(-0.66883123) q[0];
sx q[0];
rz(-0.87829725) q[0];
x q[1];
rz(1.5397649) q[2];
sx q[2];
rz(-0.73756644) q[2];
sx q[2];
rz(2.8013128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0552057) q[1];
sx q[1];
rz(-1.3761763) q[1];
sx q[1];
rz(-0.17533949) q[1];
rz(1.1176923) q[3];
sx q[3];
rz(-3.0278904) q[3];
sx q[3];
rz(-2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(-0.97186175) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(-2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650836) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(-0.064095108) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(-1.6832738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653939) q[0];
sx q[0];
rz(-1.6178693) q[0];
sx q[0];
rz(-1.7182299) q[0];
rz(-pi) q[1];
rz(-0.96785737) q[2];
sx q[2];
rz(-1.0087011) q[2];
sx q[2];
rz(-2.2130655) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.872936) q[1];
sx q[1];
rz(-1.6246512) q[1];
sx q[1];
rz(-0.40952803) q[1];
x q[2];
rz(-3.0890373) q[3];
sx q[3];
rz(-2.8263546) q[3];
sx q[3];
rz(1.0759575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1594499) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(-1.7689765) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963592) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(-2.0171719) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(-1.4845928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2161908) q[0];
sx q[0];
rz(-1.6214217) q[0];
sx q[0];
rz(3.1111858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.053327) q[2];
sx q[2];
rz(-1.4033068) q[2];
sx q[2];
rz(-2.6094764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2034519) q[1];
sx q[1];
rz(-1.005299) q[1];
sx q[1];
rz(-2.6184665) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7709511) q[3];
sx q[3];
rz(-2.3834043) q[3];
sx q[3];
rz(-1.7565808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1241887) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(2.1253288) q[2];
rz(1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50773412) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(0.39988363) q[0];
rz(-1.9790861) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-2.9679325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4325718) q[0];
sx q[0];
rz(-3.0637494) q[0];
sx q[0];
rz(2.7709333) q[0];
rz(-pi) q[1];
rz(-1.7816254) q[2];
sx q[2];
rz(-1.1846917) q[2];
sx q[2];
rz(-2.2112276) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.793321) q[1];
sx q[1];
rz(-1.6760567) q[1];
sx q[1];
rz(-3.1131016) q[1];
x q[2];
rz(-2.8097314) q[3];
sx q[3];
rz(-1.3734986) q[3];
sx q[3];
rz(-0.33533898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(-0.85401946) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(1.4703898) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(1.1434198) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66577673) q[0];
sx q[0];
rz(-1.4205298) q[0];
sx q[0];
rz(-1.4861602) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5383699) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(3.0657257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1868362) q[1];
sx q[1];
rz(-1.7446767) q[1];
sx q[1];
rz(-1.9435005) q[1];
rz(2.0538123) q[3];
sx q[3];
rz(-2.5452168) q[3];
sx q[3];
rz(3.0697825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4986971) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(1.6112304) q[2];
rz(1.4536084) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(3.0221508) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(2.4123736) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(2.008332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92111174) q[0];
sx q[0];
rz(-0.32520121) q[0];
sx q[0];
rz(0.13334206) q[0];
x q[1];
rz(-1.5845756) q[2];
sx q[2];
rz(-2.2829977) q[2];
sx q[2];
rz(-2.6028002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1843695) q[1];
sx q[1];
rz(-1.8436699) q[1];
sx q[1];
rz(0.89783122) q[1];
rz(-pi) q[2];
rz(-0.53448581) q[3];
sx q[3];
rz(-0.30189461) q[3];
sx q[3];
rz(-1.9369949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(2.6531632) q[2];
rz(-1.3119665) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(-0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(0.07117614) q[0];
rz(-0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.9326928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8457444) q[0];
sx q[0];
rz(-2.309572) q[0];
sx q[0];
rz(-0.31006281) q[0];
rz(2.2279943) q[2];
sx q[2];
rz(-2.845394) q[2];
sx q[2];
rz(1.2072472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3370918) q[1];
sx q[1];
rz(-0.13131222) q[1];
sx q[1];
rz(-1.4485703) q[1];
x q[2];
rz(0.62187059) q[3];
sx q[3];
rz(-1.2688046) q[3];
sx q[3];
rz(-0.57884502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9341087) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(2.7895555) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1677925) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-0.31162509) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.5100381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85743839) q[0];
sx q[0];
rz(-2.7630002) q[0];
sx q[0];
rz(-2.5749102) q[0];
rz(-pi) q[1];
rz(0.19599234) q[2];
sx q[2];
rz(-1.9078983) q[2];
sx q[2];
rz(-2.5981284) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.65731243) q[1];
sx q[1];
rz(-1.3687951) q[1];
sx q[1];
rz(2.5217418) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66442499) q[3];
sx q[3];
rz(-0.29237745) q[3];
sx q[3];
rz(-1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2733549) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.2379439) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8907392) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-0.2302641) q[0];
rz(0.62581217) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(2.4831916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31357665) q[0];
sx q[0];
rz(-0.99452924) q[0];
sx q[0];
rz(-0.97998951) q[0];
rz(-pi) q[1];
rz(1.9071104) q[2];
sx q[2];
rz(-2.2173777) q[2];
sx q[2];
rz(-2.1474349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2645996) q[1];
sx q[1];
rz(-2.0473285) q[1];
sx q[1];
rz(-0.18696733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55024054) q[3];
sx q[3];
rz(-0.60866683) q[3];
sx q[3];
rz(0.88702162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(-2.1181469) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(-1.9256928) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6170549) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(1.8160309) q[2];
sx q[2];
rz(-2.500848) q[2];
sx q[2];
rz(1.4304888) q[2];
rz(-0.099048793) q[3];
sx q[3];
rz(-2.5669813) q[3];
sx q[3];
rz(3.0909227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];