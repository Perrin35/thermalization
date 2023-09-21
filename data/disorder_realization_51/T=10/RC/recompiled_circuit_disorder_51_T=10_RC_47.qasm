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
rz(5.6097538) q[1];
sx q[1];
rz(14.654832) q[1];
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
rz(-2.7396169) q[2];
sx q[2];
rz(-2.4789414) q[2];
sx q[2];
rz(-2.0768009) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44838201) q[1];
sx q[1];
rz(-0.795185) q[1];
sx q[1];
rz(2.8661212) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16206046) q[3];
sx q[3];
rz(-1.2352304) q[3];
sx q[3];
rz(-0.1048564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16333214) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(1.7738316) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098009) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(0.59511551) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(-1.6275303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1925416) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(2.2632954) q[0];
x q[1];
rz(-0.02818429) q[2];
sx q[2];
rz(-2.3079254) q[2];
sx q[2];
rz(-0.29836269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0552057) q[1];
sx q[1];
rz(-1.7654164) q[1];
sx q[1];
rz(-0.17533949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.049948143) q[3];
sx q[3];
rz(-1.6729828) q[3];
sx q[3];
rz(0.12168836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4864768) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(2.1697309) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(-0.064095108) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(-1.4583189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4401535) q[0];
sx q[0];
rz(-2.9868786) q[0];
sx q[0];
rz(-1.2604777) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96785737) q[2];
sx q[2];
rz(-2.1328916) q[2];
sx q[2];
rz(0.92852718) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4255193) q[1];
sx q[1];
rz(-0.41285535) q[1];
sx q[1];
rz(-3.0070261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5879257) q[3];
sx q[3];
rz(-1.2560085) q[3];
sx q[3];
rz(1.1312315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98214275) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(-0.94846559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.8452334) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(-2.0171719) q[0];
rz(-1.9212978) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.6569998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67498818) q[0];
sx q[0];
rz(-0.059048422) q[0];
sx q[0];
rz(1.0303636) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2414615) q[2];
sx q[2];
rz(-0.54154684) q[2];
sx q[2];
rz(-2.3878218) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1160108) q[1];
sx q[1];
rz(-2.3911747) q[1];
sx q[1];
rz(-2.2376899) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7709511) q[3];
sx q[3];
rz(-0.7581884) q[3];
sx q[3];
rz(1.3850118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-2.1253288) q[2];
rz(-1.4034363) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(2.0531634) q[3];
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
rz(0.50773412) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(-1.9790861) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(-0.17366017) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50786103) q[0];
sx q[0];
rz(-1.5989688) q[0];
sx q[0];
rz(0.072576056) q[0];
rz(-pi) q[1];
rz(-0.47541754) q[2];
sx q[2];
rz(-2.7042275) q[2];
sx q[2];
rz(1.4471444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.34827161) q[1];
sx q[1];
rz(-1.6760567) q[1];
sx q[1];
rz(-0.028491032) q[1];
x q[2];
rz(0.55032702) q[3];
sx q[3];
rz(-0.384207) q[3];
sx q[3];
rz(1.7526383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48012039) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(2.373467) q[2];
rz(-0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1059234) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(1.4703898) q[0];
rz(-2.629783) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.1434198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4758159) q[0];
sx q[0];
rz(-1.7210628) q[0];
sx q[0];
rz(-1.4861602) q[0];
rz(-pi) q[1];
rz(-0.60322275) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(3.0657257) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.032635078) q[1];
sx q[1];
rz(-0.40954486) q[1];
sx q[1];
rz(-1.1213379) q[1];
x q[2];
rz(2.8361736) q[3];
sx q[3];
rz(-2.0913887) q[3];
sx q[3];
rz(0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4986971) q[2];
sx q[2];
rz(-0.84474793) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0221508) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(1.4402333) q[0];
rz(-2.4123736) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-2.008332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92111174) q[0];
sx q[0];
rz(-0.32520121) q[0];
sx q[0];
rz(0.13334206) q[0];
rz(-pi) q[1];
rz(-1.5845756) q[2];
sx q[2];
rz(-2.2829977) q[2];
sx q[2];
rz(0.53879246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82517351) q[1];
sx q[1];
rz(-0.9269886) q[1];
sx q[1];
rz(-0.34367798) q[1];
rz(-2.8797638) q[3];
sx q[3];
rz(-1.7228408) q[3];
sx q[3];
rz(0.1482299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40522727) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-2.6531632) q[2];
rz(1.8296261) q[3];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0768123) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(-3.0704165) q[0];
rz(0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.2088998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14825386) q[0];
sx q[0];
rz(-2.3518666) q[0];
sx q[0];
rz(-1.8940311) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2279943) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(1.9343455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4603685) q[1];
sx q[1];
rz(-1.4404693) q[1];
sx q[1];
rz(-3.1254915) q[1];
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
rz(-0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677925) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(1.6315546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841543) q[0];
sx q[0];
rz(-0.3785924) q[0];
sx q[0];
rz(2.5749102) q[0];
x q[1];
rz(0.19599234) q[2];
sx q[2];
rz(-1.9078983) q[2];
sx q[2];
rz(0.54346426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4842802) q[1];
sx q[1];
rz(-1.3687951) q[1];
sx q[1];
rz(-2.5217418) q[1];
x q[2];
rz(1.7543091) q[3];
sx q[3];
rz(-1.7997051) q[3];
sx q[3];
rz(-2.6218417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2733549) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(-2.6319035) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8907392) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(2.9113286) q[0];
rz(-2.5157805) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(2.4831916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9068245) q[0];
sx q[0];
rz(-1.0848197) q[0];
sx q[0];
rz(2.4776239) q[0];
rz(2.7294331) q[2];
sx q[2];
rz(-0.71752749) q[2];
sx q[2];
rz(1.621643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.534348) q[1];
sx q[1];
rz(-1.736728) q[1];
sx q[1];
rz(1.0870618) q[1];
x q[2];
rz(-1.9202616) q[3];
sx q[3];
rz(-1.06171) q[3];
sx q[3];
rz(-1.5290608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-2.6170549) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(-1.2383923) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
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