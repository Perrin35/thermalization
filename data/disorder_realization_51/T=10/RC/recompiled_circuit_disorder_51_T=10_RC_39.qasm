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
rz(-2.9601331) q[0];
rz(2.060086) q[1];
sx q[1];
rz(5.6097538) q[1];
sx q[1];
rz(14.654832) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3436514) q[0];
sx q[0];
rz(-2.6290253) q[0];
sx q[0];
rz(-1.1287862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5187831) q[2];
sx q[2];
rz(-1.8138759) q[2];
sx q[2];
rz(-0.18261766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3095113) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(-1.3002212) q[1];
x q[2];
rz(1.231108) q[3];
sx q[3];
rz(-1.4178483) q[3];
sx q[3];
rz(-1.7294401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16333214) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(-1.7738316) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098009) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(0.59511551) q[0];
rz(-1.0455421) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.6275303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1925416) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(2.2632954) q[0];
x q[1];
rz(-3.1134084) q[2];
sx q[2];
rz(-0.83366725) q[2];
sx q[2];
rz(-0.29836269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6914312) q[1];
sx q[1];
rz(-1.7427923) q[1];
sx q[1];
rz(1.7683692) q[1];
rz(-pi) q[2];
rz(1.6731095) q[3];
sx q[3];
rz(-1.6204837) q[3];
sx q[3];
rz(-1.697584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65511584) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(0.97186175) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(2.9698353) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650836) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(1.4583189) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4401535) q[0];
sx q[0];
rz(-0.15471409) q[0];
sx q[0];
rz(-1.8811149) q[0];
x q[1];
rz(-2.1737353) q[2];
sx q[2];
rz(-1.0087011) q[2];
sx q[2];
rz(-0.92852718) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26865667) q[1];
sx q[1];
rz(-1.6246512) q[1];
sx q[1];
rz(2.7320646) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8267616) q[3];
sx q[3];
rz(-1.5545087) q[3];
sx q[3];
rz(-0.44486886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(0.67908755) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963592) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-1.1244208) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.6569998) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.785448) q[0];
sx q[0];
rz(-1.5404285) q[0];
sx q[0];
rz(1.6214451) q[0];
x q[1];
rz(1.053327) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(0.53211624) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.209219) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(-0.93851628) q[1];
rz(-pi) q[2];
rz(0.72317601) q[3];
sx q[3];
rz(-1.8225065) q[3];
sx q[3];
rz(0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(2.1253288) q[2];
rz(1.4034363) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338585) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(-1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(2.9679325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7090209) q[0];
sx q[0];
rz(-0.077843277) q[0];
sx q[0];
rz(2.7709333) q[0];
x q[1];
rz(2.6661751) q[2];
sx q[2];
rz(-0.43736514) q[2];
sx q[2];
rz(1.6944483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.083399051) q[1];
sx q[1];
rz(-3.0325583) q[1];
sx q[1];
rz(-1.3074247) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55032702) q[3];
sx q[3];
rz(-2.7573857) q[3];
sx q[3];
rz(-1.7526383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48012039) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(-0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.0356692) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(-1.6712028) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.1434198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15018806) q[0];
sx q[0];
rz(-2.9692869) q[0];
sx q[0];
rz(-0.50921391) q[0];
x q[1];
rz(0.60322275) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(0.075866931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1868362) q[1];
sx q[1];
rz(-1.3969159) q[1];
sx q[1];
rz(-1.1980921) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8361736) q[3];
sx q[3];
rz(-1.050204) q[3];
sx q[3];
rz(0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4986971) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.5303622) q[2];
rz(1.4536084) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11944184) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(2.4123736) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(-2.008332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92111174) q[0];
sx q[0];
rz(-0.32520121) q[0];
sx q[0];
rz(3.0082506) q[0];
rz(0.015958162) q[2];
sx q[2];
rz(-0.71231132) q[2];
sx q[2];
rz(-0.51770891) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1843695) q[1];
sx q[1];
rz(-1.8436699) q[1];
sx q[1];
rz(2.2437614) q[1];
rz(0.53448581) q[3];
sx q[3];
rz(-0.30189461) q[3];
sx q[3];
rz(-1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40522727) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-2.6531632) q[2];
rz(-1.8296261) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(0.64129889) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(0.07117614) q[0];
rz(-3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.2088998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9933388) q[0];
sx q[0];
rz(-0.78972602) q[0];
sx q[0];
rz(-1.8940311) q[0];
rz(-pi) q[1];
rz(-2.2279943) q[2];
sx q[2];
rz(-2.845394) q[2];
sx q[2];
rz(-1.2072472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2541131) q[1];
sx q[1];
rz(-1.5548318) q[1];
sx q[1];
rz(-1.70114) q[1];
x q[2];
rz(-1.204793) q[3];
sx q[3];
rz(-0.98099785) q[3];
sx q[3];
rz(2.3596711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9341087) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(2.2682244) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677925) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(-0.82178003) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.6315546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841543) q[0];
sx q[0];
rz(-0.3785924) q[0];
sx q[0];
rz(-0.5666825) q[0];
rz(-pi) q[1];
rz(-2.9456003) q[2];
sx q[2];
rz(-1.2336944) q[2];
sx q[2];
rz(2.5981284) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65731243) q[1];
sx q[1];
rz(-1.3687951) q[1];
sx q[1];
rz(-0.61985086) q[1];
x q[2];
rz(-2.4771677) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(-1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(2.6319035) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508535) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(-0.62581217) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(-0.65840107) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31357665) q[0];
sx q[0];
rz(-2.1470634) q[0];
sx q[0];
rz(2.1616031) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7294331) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(-1.5199496) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.534348) q[1];
sx q[1];
rz(-1.4048647) q[1];
sx q[1];
rz(-2.0545309) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53604605) q[3];
sx q[3];
rz(-1.8744161) q[3];
sx q[3];
rz(0.21751465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(-1.0234458) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(-1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(-2.5692206) q[3];
sx q[3];
rz(-1.5170245) q[3];
sx q[3];
rz(-1.7046884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
