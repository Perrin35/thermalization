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
rz(0.45646271) q[0];
sx q[0];
rz(4.0153761) q[0];
sx q[0];
rz(10.562513) q[0];
rz(-3.0955834) q[1];
sx q[1];
rz(-2.6061821) q[1];
sx q[1];
rz(1.5387662) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3483602) q[0];
sx q[0];
rz(-1.6038415) q[0];
sx q[0];
rz(1.6152302) q[0];
x q[1];
rz(-0.12243226) q[2];
sx q[2];
rz(-2.2847567) q[2];
sx q[2];
rz(-3.0453504) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12373172) q[1];
sx q[1];
rz(-1.3657454) q[1];
sx q[1];
rz(-0.46242949) q[1];
rz(-pi) q[2];
x q[2];
rz(1.192126) q[3];
sx q[3];
rz(-0.14992564) q[3];
sx q[3];
rz(2.718234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5753182) q[2];
sx q[2];
rz(-0.34276572) q[2];
sx q[2];
rz(0.23635593) q[2];
rz(-2.6929839) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(-1.0033222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7597294) q[0];
sx q[0];
rz(-0.48668447) q[0];
sx q[0];
rz(-0.78616649) q[0];
rz(0.78474125) q[1];
sx q[1];
rz(-1.4737543) q[1];
sx q[1];
rz(-3.0395708) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.611972) q[0];
sx q[0];
rz(-1.0361605) q[0];
sx q[0];
rz(1.8055339) q[0];
rz(-1.2743852) q[2];
sx q[2];
rz(-0.80889946) q[2];
sx q[2];
rz(2.1975225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0077304) q[1];
sx q[1];
rz(-1.8419767) q[1];
sx q[1];
rz(2.8274963) q[1];
x q[2];
rz(-2.8972273) q[3];
sx q[3];
rz(-0.8042585) q[3];
sx q[3];
rz(1.8387295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1138136) q[2];
sx q[2];
rz(-1.6543829) q[2];
sx q[2];
rz(0.2629183) q[2];
rz(1.8391838) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(-0.7836248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0021492783) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(-1.333746) q[0];
rz(-1.7568781) q[1];
sx q[1];
rz(-0.92888558) q[1];
sx q[1];
rz(2.3768916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686808) q[0];
sx q[0];
rz(-0.31038767) q[0];
sx q[0];
rz(-2.458802) q[0];
x q[1];
rz(-0.3777067) q[2];
sx q[2];
rz(-1.3481082) q[2];
sx q[2];
rz(1.6645704) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0955268) q[1];
sx q[1];
rz(-2.8316759) q[1];
sx q[1];
rz(0.97280963) q[1];
rz(-0.30454347) q[3];
sx q[3];
rz(-2.2650654) q[3];
sx q[3];
rz(2.3643431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1769522) q[2];
sx q[2];
rz(-1.8884337) q[2];
sx q[2];
rz(-2.9580252) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95530987) q[0];
sx q[0];
rz(-2.0109542) q[0];
sx q[0];
rz(-0.282222) q[0];
rz(1.9918293) q[1];
sx q[1];
rz(-1.7825922) q[1];
sx q[1];
rz(-1.6414292) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.09984) q[0];
sx q[0];
rz(-2.7405993) q[0];
sx q[0];
rz(0.45464154) q[0];
rz(-pi) q[1];
rz(2.1181951) q[2];
sx q[2];
rz(-1.3522902) q[2];
sx q[2];
rz(-2.5283739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0974429) q[1];
sx q[1];
rz(-1.9056093) q[1];
sx q[1];
rz(-0.26589946) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4789739) q[3];
sx q[3];
rz(-1.8739376) q[3];
sx q[3];
rz(1.417899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1993316) q[2];
sx q[2];
rz(-1.4456238) q[2];
sx q[2];
rz(-1.4873827) q[2];
rz(2.2516294) q[3];
sx q[3];
rz(-1.7421937) q[3];
sx q[3];
rz(0.14931211) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99059659) q[0];
sx q[0];
rz(-2.6166333) q[0];
sx q[0];
rz(1.4599482) q[0];
rz(1.8563942) q[1];
sx q[1];
rz(-1.7743013) q[1];
sx q[1];
rz(0.45305124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.734402) q[0];
sx q[0];
rz(-1.2670867) q[0];
sx q[0];
rz(-0.056945368) q[0];
rz(0.052930367) q[2];
sx q[2];
rz(-0.9890511) q[2];
sx q[2];
rz(-2.5404921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1374319) q[1];
sx q[1];
rz(-0.42037005) q[1];
sx q[1];
rz(-2.67291) q[1];
rz(-pi) q[2];
rz(-0.83638453) q[3];
sx q[3];
rz(-0.57698133) q[3];
sx q[3];
rz(2.0744689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4840661) q[2];
sx q[2];
rz(-2.6675197) q[2];
sx q[2];
rz(-1.5849812) q[2];
rz(1.5786242) q[3];
sx q[3];
rz(-1.862674) q[3];
sx q[3];
rz(-1.0274308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22353657) q[0];
sx q[0];
rz(-0.99128381) q[0];
sx q[0];
rz(1.9687442) q[0];
rz(2.6808443) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(-1.4985098) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9665034) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(2.1835292) q[0];
rz(-2.8086964) q[2];
sx q[2];
rz(-1.82331) q[2];
sx q[2];
rz(-1.1292924) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9869796) q[1];
sx q[1];
rz(-1.3081395) q[1];
sx q[1];
rz(-0.059537454) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13883484) q[3];
sx q[3];
rz(-0.87362008) q[3];
sx q[3];
rz(2.4133854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0074761) q[2];
sx q[2];
rz(-1.7191929) q[2];
sx q[2];
rz(0.22809347) q[2];
rz(2.5988233) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-0.91226474) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11693624) q[0];
sx q[0];
rz(-2.0967364) q[0];
sx q[0];
rz(-2.5352617) q[0];
rz(1.2888651) q[1];
sx q[1];
rz(-1.6090798) q[1];
sx q[1];
rz(2.2859763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3324749) q[0];
sx q[0];
rz(-0.98253378) q[0];
sx q[0];
rz(1.0756798) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4855723) q[2];
sx q[2];
rz(-0.79932511) q[2];
sx q[2];
rz(-2.8003789) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42509584) q[1];
sx q[1];
rz(-1.3653058) q[1];
sx q[1];
rz(2.53611) q[1];
rz(-pi) q[2];
rz(2.8884573) q[3];
sx q[3];
rz(-2.0380424) q[3];
sx q[3];
rz(0.96880355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.764708) q[2];
sx q[2];
rz(-2.7089684) q[2];
sx q[2];
rz(0.077733668) q[2];
rz(3.0149095) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(2.3104987) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34778255) q[0];
sx q[0];
rz(-2.7816483) q[0];
sx q[0];
rz(-0.98156324) q[0];
rz(-1.7836102) q[1];
sx q[1];
rz(-2.6505018) q[1];
sx q[1];
rz(-2.8474999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9413853) q[0];
sx q[0];
rz(-1.5110755) q[0];
sx q[0];
rz(0.30903791) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.852599) q[2];
sx q[2];
rz(-1.929612) q[2];
sx q[2];
rz(-3.0274018) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5433054) q[1];
sx q[1];
rz(-1.2906055) q[1];
sx q[1];
rz(2.2375475) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96947396) q[3];
sx q[3];
rz(-3.0532928) q[3];
sx q[3];
rz(1.2754319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0391417) q[2];
sx q[2];
rz(-0.65639085) q[2];
sx q[2];
rz(2.0810818) q[2];
rz(2.5833526) q[3];
sx q[3];
rz(-2.5679936) q[3];
sx q[3];
rz(-3.1225045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3625665) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(-2.3574164) q[0];
rz(-1.7991637) q[1];
sx q[1];
rz(-1.8524086) q[1];
sx q[1];
rz(-2.4450891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7292228) q[0];
sx q[0];
rz(-1.7984838) q[0];
sx q[0];
rz(-0.21842893) q[0];
rz(-2.4694211) q[2];
sx q[2];
rz(-1.1706588) q[2];
sx q[2];
rz(-2.6113685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99217691) q[1];
sx q[1];
rz(-2.0609629) q[1];
sx q[1];
rz(2.0489245) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4217997) q[3];
sx q[3];
rz(-1.4382575) q[3];
sx q[3];
rz(0.46035351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8037618) q[2];
sx q[2];
rz(-0.75118128) q[2];
sx q[2];
rz(-1.2639698) q[2];
rz(-1.7654644) q[3];
sx q[3];
rz(-2.2270146) q[3];
sx q[3];
rz(1.5553364) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92397583) q[0];
sx q[0];
rz(-3.0986077) q[0];
sx q[0];
rz(-1.6643583) q[0];
rz(-2.2337275) q[1];
sx q[1];
rz(-1.5217179) q[1];
sx q[1];
rz(0.97000617) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5805626) q[0];
sx q[0];
rz(-0.99107689) q[0];
sx q[0];
rz(-0.28663488) q[0];
x q[1];
rz(2.1625948) q[2];
sx q[2];
rz(-1.3221957) q[2];
sx q[2];
rz(2.4112041) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1906311) q[1];
sx q[1];
rz(-2.4959557) q[1];
sx q[1];
rz(-0.26107045) q[1];
x q[2];
rz(2.2144058) q[3];
sx q[3];
rz(-0.95706104) q[3];
sx q[3];
rz(-2.9184379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.19405356) q[2];
sx q[2];
rz(-2.2668362) q[2];
sx q[2];
rz(0.27210316) q[2];
rz(2.2943606) q[3];
sx q[3];
rz(-0.50744414) q[3];
sx q[3];
rz(-2.0525172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67615164) q[0];
sx q[0];
rz(-1.1347102) q[0];
sx q[0];
rz(1.4665428) q[0];
rz(-2.8027986) q[1];
sx q[1];
rz(-1.4452965) q[1];
sx q[1];
rz(-1.8903587) q[1];
rz(-1.9235545) q[2];
sx q[2];
rz(-2.259476) q[2];
sx q[2];
rz(-1.8720575) q[2];
rz(-1.7539514) q[3];
sx q[3];
rz(-0.42944254) q[3];
sx q[3];
rz(-2.4933168) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
