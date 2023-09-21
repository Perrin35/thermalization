OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(-0.98288012) q[0];
sx q[0];
rz(-2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(-1.2892105) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50531189) q[0];
sx q[0];
rz(-1.2828151) q[0];
sx q[0];
rz(1.0907537) q[0];
rz(-pi) q[1];
rz(-2.870954) q[2];
sx q[2];
rz(-0.025279609) q[2];
sx q[2];
rz(-1.4852922) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7154555) q[1];
sx q[1];
rz(-1.8278367) q[1];
sx q[1];
rz(0.21547683) q[1];
x q[2];
rz(3.0715277) q[3];
sx q[3];
rz(-1.3032039) q[3];
sx q[3];
rz(-2.23578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-3.0060449) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41483375) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(-2.360789) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(-1.82812) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75269964) q[0];
sx q[0];
rz(-1.7912731) q[0];
sx q[0];
rz(-1.3541143) q[0];
rz(-pi) q[1];
rz(-2.8333089) q[2];
sx q[2];
rz(-1.1051902) q[2];
sx q[2];
rz(-2.3870433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7105512) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(2.980742) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1900858) q[3];
sx q[3];
rz(-0.12976876) q[3];
sx q[3];
rz(-1.6956003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69497067) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(-2.1874645) q[2];
rz(-1.4387087) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(2.4296956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-0.69586786) q[0];
rz(2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(2.9755039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9244708) q[0];
sx q[0];
rz(-0.32453254) q[0];
sx q[0];
rz(-1.0990418) q[0];
x q[1];
rz(2.0735998) q[2];
sx q[2];
rz(-1.0223801) q[2];
sx q[2];
rz(-1.8333679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3661763) q[1];
sx q[1];
rz(-1.1300547) q[1];
sx q[1];
rz(-2.612545) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8576287) q[3];
sx q[3];
rz(-1.3751251) q[3];
sx q[3];
rz(-1.5989725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35963905) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(2.847805) q[0];
rz(-0.44149533) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-0.77484432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.144865) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(0.073547151) q[0];
rz(-pi) q[1];
rz(-2.7735633) q[2];
sx q[2];
rz(-1.5291011) q[2];
sx q[2];
rz(2.7421943) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6619819) q[1];
sx q[1];
rz(-1.4927215) q[1];
sx q[1];
rz(-1.1675203) q[1];
rz(-pi) q[2];
rz(2.1500548) q[3];
sx q[3];
rz(-0.30617985) q[3];
sx q[3];
rz(-1.3907719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0094770771) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-0.57007989) q[2];
rz(-1.305497) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(1.6397887) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(-0.18138012) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(-0.10770527) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5126702) q[0];
sx q[0];
rz(-1.9229917) q[0];
sx q[0];
rz(-3.0432426) q[0];
rz(-pi) q[1];
rz(-0.55107112) q[2];
sx q[2];
rz(-1.0017918) q[2];
sx q[2];
rz(1.0702733) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4957461) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(-0.87162019) q[1];
rz(2.7961736) q[3];
sx q[3];
rz(-1.0343699) q[3];
sx q[3];
rz(0.1766583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(-2.8520612) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(1.6311197) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(1.0169792) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6769584) q[0];
sx q[0];
rz(-0.58566739) q[0];
sx q[0];
rz(0.93924378) q[0];
rz(-pi) q[1];
rz(-0.11634155) q[2];
sx q[2];
rz(-2.939072) q[2];
sx q[2];
rz(-0.6845189) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79343866) q[1];
sx q[1];
rz(-0.70421709) q[1];
sx q[1];
rz(2.3401295) q[1];
rz(-pi) q[2];
rz(2.0820886) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(-2.2389776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-2.5069359) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(-2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.9218146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4635173) q[0];
sx q[0];
rz(-3.0419555) q[0];
sx q[0];
rz(0.3936605) q[0];
x q[1];
rz(-1.9856521) q[2];
sx q[2];
rz(-1.2659237) q[2];
sx q[2];
rz(2.4883303) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4525758) q[1];
sx q[1];
rz(-0.92454443) q[1];
sx q[1];
rz(2.8313682) q[1];
x q[2];
rz(1.9795498) q[3];
sx q[3];
rz(-1.6533274) q[3];
sx q[3];
rz(2.2857034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(-2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379631) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(-1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(-2.7323639) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0497826) q[0];
sx q[0];
rz(-0.70735332) q[0];
sx q[0];
rz(-2.011767) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6312469) q[2];
sx q[2];
rz(-2.5556457) q[2];
sx q[2];
rz(-0.79615359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1815391) q[1];
sx q[1];
rz(-2.1337193) q[1];
sx q[1];
rz(-1.3905418) q[1];
x q[2];
rz(-1.8141361) q[3];
sx q[3];
rz(-1.0662603) q[3];
sx q[3];
rz(-0.53989053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(-3.1316485) q[0];
rz(2.0857281) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(1.6329637) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2752339) q[0];
sx q[0];
rz(-0.70952053) q[0];
sx q[0];
rz(-1.3225681) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6998859) q[2];
sx q[2];
rz(-0.61868762) q[2];
sx q[2];
rz(2.8512851) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3907527) q[1];
sx q[1];
rz(-2.5166593) q[1];
sx q[1];
rz(-2.1363439) q[1];
rz(1.6375332) q[3];
sx q[3];
rz(-2.4373694) q[3];
sx q[3];
rz(3.0610386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(2.1127624) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(2.57634) q[0];
rz(-0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(-1.7601097) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5436343) q[0];
sx q[0];
rz(-1.6126452) q[0];
sx q[0];
rz(2.6393487) q[0];
rz(-2.6673615) q[2];
sx q[2];
rz(-1.6696764) q[2];
sx q[2];
rz(-0.41254166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.46422058) q[1];
sx q[1];
rz(-0.83253011) q[1];
sx q[1];
rz(0.99454753) q[1];
x q[2];
rz(2.284427) q[3];
sx q[3];
rz(-1.790159) q[3];
sx q[3];
rz(0.21881783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(-2.6560442) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.272841) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(1.3960719) q[2];
sx q[2];
rz(-0.25824418) q[2];
sx q[2];
rz(1.2373409) q[2];
rz(-2.6315401) q[3];
sx q[3];
rz(-1.4997713) q[3];
sx q[3];
rz(-0.98239519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];