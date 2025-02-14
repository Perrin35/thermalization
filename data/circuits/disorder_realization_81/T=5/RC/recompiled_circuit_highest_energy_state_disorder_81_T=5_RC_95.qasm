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
rz(-2.6851299) q[0];
sx q[0];
rz(-0.87378341) q[0];
sx q[0];
rz(-1.1377347) q[0];
rz(0.046009215) q[1];
sx q[1];
rz(-0.53541056) q[1];
sx q[1];
rz(1.6028264) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7932324) q[0];
sx q[0];
rz(-1.6038415) q[0];
sx q[0];
rz(1.6152302) q[0];
rz(-pi) q[1];
rz(-2.2884772) q[2];
sx q[2];
rz(-1.6632281) q[2];
sx q[2];
rz(-1.5866367) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12373172) q[1];
sx q[1];
rz(-1.3657454) q[1];
sx q[1];
rz(-2.6791632) q[1];
rz(-pi) q[2];
rz(-1.9494667) q[3];
sx q[3];
rz(-0.14992564) q[3];
sx q[3];
rz(-0.42335864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5662745) q[2];
sx q[2];
rz(-2.7988269) q[2];
sx q[2];
rz(2.9052367) q[2];
rz(0.44860873) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(-1.0033222) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3818632) q[0];
sx q[0];
rz(-0.48668447) q[0];
sx q[0];
rz(-0.78616649) q[0];
rz(0.78474125) q[1];
sx q[1];
rz(-1.4737543) q[1];
sx q[1];
rz(-3.0395708) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090793153) q[0];
sx q[0];
rz(-0.57928071) q[0];
sx q[0];
rz(-0.37426484) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3574102) q[2];
sx q[2];
rz(-1.3578556) q[2];
sx q[2];
rz(2.7226457) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1338623) q[1];
sx q[1];
rz(-1.8419767) q[1];
sx q[1];
rz(2.8274963) q[1];
x q[2];
rz(-0.24436538) q[3];
sx q[3];
rz(-0.8042585) q[3];
sx q[3];
rz(1.3028631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1138136) q[2];
sx q[2];
rz(-1.4872097) q[2];
sx q[2];
rz(2.8786744) q[2];
rz(-1.8391838) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(0.7836248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1394434) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(-1.8078467) q[0];
rz(1.7568781) q[1];
sx q[1];
rz(-2.2127071) q[1];
sx q[1];
rz(-0.76470107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36604213) q[0];
sx q[0];
rz(-1.8100272) q[0];
sx q[0];
rz(-1.7704829) q[0];
rz(-pi) q[1];
rz(1.809757) q[2];
sx q[2];
rz(-1.2028678) q[2];
sx q[2];
rz(-2.9604173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6669608) q[1];
sx q[1];
rz(-1.8256011) q[1];
sx q[1];
rz(-0.17837015) q[1];
rz(1.2250879) q[3];
sx q[3];
rz(-0.74784333) q[3];
sx q[3];
rz(-1.2338828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1769522) q[2];
sx q[2];
rz(-1.2531589) q[2];
sx q[2];
rz(-0.18356744) q[2];
rz(2.99672) q[3];
sx q[3];
rz(-1.262007) q[3];
sx q[3];
rz(2.9174793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1862828) q[0];
sx q[0];
rz(-2.0109542) q[0];
sx q[0];
rz(0.282222) q[0];
rz(1.1497633) q[1];
sx q[1];
rz(-1.7825922) q[1];
sx q[1];
rz(-1.5001635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44631413) q[0];
sx q[0];
rz(-1.2124966) q[0];
sx q[0];
rz(-1.3867239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.167539) q[2];
sx q[2];
rz(-2.5563498) q[2];
sx q[2];
rz(-2.5257021) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0974429) q[1];
sx q[1];
rz(-1.9056093) q[1];
sx q[1];
rz(-2.8756932) q[1];
rz(-2.8372466) q[3];
sx q[3];
rz(-1.4831717) q[3];
sx q[3];
rz(0.1803785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1993316) q[2];
sx q[2];
rz(-1.4456238) q[2];
sx q[2];
rz(1.4873827) q[2];
rz(0.88996327) q[3];
sx q[3];
rz(-1.3993989) q[3];
sx q[3];
rz(0.14931211) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1509961) q[0];
sx q[0];
rz(-0.52495933) q[0];
sx q[0];
rz(-1.6816444) q[0];
rz(-1.8563942) q[1];
sx q[1];
rz(-1.7743013) q[1];
sx q[1];
rz(-0.45305124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1465587) q[0];
sx q[0];
rz(-1.6251329) q[0];
sx q[0];
rz(-1.2666235) q[0];
rz(-3.0886623) q[2];
sx q[2];
rz(-0.9890511) q[2];
sx q[2];
rz(-2.5404921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.510524) q[1];
sx q[1];
rz(-1.1981401) q[1];
sx q[1];
rz(-1.371553) q[1];
x q[2];
rz(-2.3052081) q[3];
sx q[3];
rz(-2.5646113) q[3];
sx q[3];
rz(-1.0671237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4840661) q[2];
sx q[2];
rz(-0.47407293) q[2];
sx q[2];
rz(-1.5849812) q[2];
rz(-1.5629684) q[3];
sx q[3];
rz(-1.2789187) q[3];
sx q[3];
rz(1.0274308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22353657) q[0];
sx q[0];
rz(-2.1503088) q[0];
sx q[0];
rz(-1.1728485) q[0];
rz(-2.6808443) q[1];
sx q[1];
rz(-1.8811503) q[1];
sx q[1];
rz(1.6430829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39678412) q[0];
sx q[0];
rz(-0.61273324) q[0];
sx q[0];
rz(1.5694797) q[0];
x q[1];
rz(1.837311) q[2];
sx q[2];
rz(-1.892748) q[2];
sx q[2];
rz(2.6139174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.071049404) q[1];
sx q[1];
rz(-0.26916801) q[1];
sx q[1];
rz(1.7885923) q[1];
x q[2];
rz(0.86885683) q[3];
sx q[3];
rz(-1.464499) q[3];
sx q[3];
rz(-0.93206638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1341165) q[2];
sx q[2];
rz(-1.7191929) q[2];
sx q[2];
rz(-2.9134992) q[2];
rz(-2.5988233) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-2.2293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0246564) q[0];
sx q[0];
rz(-1.0448562) q[0];
sx q[0];
rz(0.60633099) q[0];
rz(-1.2888651) q[1];
sx q[1];
rz(-1.6090798) q[1];
sx q[1];
rz(-2.2859763) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037362075) q[0];
sx q[0];
rz(-2.3920569) q[0];
sx q[0];
rz(-2.5226592) q[0];
x q[1];
rz(1.4855723) q[2];
sx q[2];
rz(-2.3422675) q[2];
sx q[2];
rz(2.8003789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2825477) q[1];
sx q[1];
rz(-0.63524109) q[1];
sx q[1];
rz(2.7905491) q[1];
rz(-pi) q[2];
rz(-2.0315657) q[3];
sx q[3];
rz(-0.52692014) q[3];
sx q[3];
rz(-2.6941018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3768846) q[2];
sx q[2];
rz(-2.7089684) q[2];
sx q[2];
rz(-0.077733668) q[2];
rz(-0.12668315) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(2.3104987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7938101) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(2.1600294) q[0];
rz(1.3579824) q[1];
sx q[1];
rz(-0.49109083) q[1];
sx q[1];
rz(-0.29409274) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3273701) q[0];
sx q[0];
rz(-0.31457385) q[0];
sx q[0];
rz(2.9474763) q[0];
rz(1.2889936) q[2];
sx q[2];
rz(-1.929612) q[2];
sx q[2];
rz(0.11419088) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31042415) q[1];
sx q[1];
rz(-2.4267395) q[1];
sx q[1];
rz(1.1352886) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0915501) q[3];
sx q[3];
rz(-1.4980157) q[3];
sx q[3];
rz(-0.67228729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10245094) q[2];
sx q[2];
rz(-2.4852018) q[2];
sx q[2];
rz(-1.0605109) q[2];
rz(0.55824009) q[3];
sx q[3];
rz(-2.5679936) q[3];
sx q[3];
rz(-0.019088117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7790262) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(-0.78417626) q[0];
rz(-1.342429) q[1];
sx q[1];
rz(-1.289184) q[1];
sx q[1];
rz(0.69650355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36464115) q[0];
sx q[0];
rz(-0.31420194) q[0];
sx q[0];
rz(-0.81885851) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59665859) q[2];
sx q[2];
rz(-0.76603973) q[2];
sx q[2];
rz(0.5853563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1494157) q[1];
sx q[1];
rz(-2.0609629) q[1];
sx q[1];
rz(2.0489245) q[1];
rz(-pi) q[2];
rz(0.83902766) q[3];
sx q[3];
rz(-2.9425042) q[3];
sx q[3];
rz(1.8323048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.33783087) q[2];
sx q[2];
rz(-0.75118128) q[2];
sx q[2];
rz(1.2639698) q[2];
rz(1.7654644) q[3];
sx q[3];
rz(-0.91457808) q[3];
sx q[3];
rz(-1.5862563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2176168) q[0];
sx q[0];
rz(-0.042984977) q[0];
sx q[0];
rz(1.6643583) q[0];
rz(2.2337275) q[1];
sx q[1];
rz(-1.5217179) q[1];
sx q[1];
rz(-0.97000617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15031397) q[0];
sx q[0];
rz(-1.3320001) q[0];
sx q[0];
rz(2.1697609) q[0];
rz(-pi) q[1];
rz(2.8447609) q[2];
sx q[2];
rz(-0.99952159) q[2];
sx q[2];
rz(0.67650696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86783389) q[1];
sx q[1];
rz(-0.9504488) q[1];
sx q[1];
rz(1.378744) q[1];
rz(-2.41955) q[3];
sx q[3];
rz(-2.0835428) q[3];
sx q[3];
rz(1.3861314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19405356) q[2];
sx q[2];
rz(-0.87475646) q[2];
sx q[2];
rz(0.27210316) q[2];
rz(-0.84723204) q[3];
sx q[3];
rz(-0.50744414) q[3];
sx q[3];
rz(-2.0525172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67615164) q[0];
sx q[0];
rz(-2.0068824) q[0];
sx q[0];
rz(-1.6750499) q[0];
rz(-0.33879406) q[1];
sx q[1];
rz(-1.6962961) q[1];
sx q[1];
rz(1.2512339) q[1];
rz(0.39739824) q[2];
sx q[2];
rz(-2.3811679) q[2];
sx q[2];
rz(0.74447348) q[2];
rz(1.7539514) q[3];
sx q[3];
rz(-2.7121501) q[3];
sx q[3];
rz(0.6482758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
