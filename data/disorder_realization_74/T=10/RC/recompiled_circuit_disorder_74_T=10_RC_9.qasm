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
rz(1.8523822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6362808) q[0];
sx q[0];
rz(-1.8587776) q[0];
sx q[0];
rz(-1.0907537) q[0];
x q[1];
rz(0.024359811) q[2];
sx q[2];
rz(-1.5640386) q[2];
sx q[2];
rz(2.9565405) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7154555) q[1];
sx q[1];
rz(-1.313756) q[1];
sx q[1];
rz(-2.9261158) q[1];
rz(1.8208002) q[3];
sx q[3];
rz(-0.27640009) q[3];
sx q[3];
rz(1.1652511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-0.13554779) q[2];
rz(0.1401976) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(2.8285817) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(-1.82812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1056054) q[0];
sx q[0];
rz(-0.30788883) q[0];
sx q[0];
rz(0.76460989) q[0];
rz(-pi) q[1];
rz(1.08554) q[2];
sx q[2];
rz(-1.2962356) q[2];
sx q[2];
rz(-0.95825125) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5116509) q[1];
sx q[1];
rz(-2.959046) q[1];
sx q[1];
rz(0.49717848) q[1];
x q[2];
rz(-0.075606451) q[3];
sx q[3];
rz(-1.6763655) q[3];
sx q[3];
rz(2.3188863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.69497067) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(-2.4296956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057782877) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(2.7867735) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(-2.9755039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040708) q[0];
sx q[0];
rz(-1.4253758) q[0];
sx q[0];
rz(1.8619596) q[0];
rz(-pi) q[1];
rz(-0.60885749) q[2];
sx q[2];
rz(-1.9945842) q[2];
sx q[2];
rz(-2.599803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7158311) q[1];
sx q[1];
rz(-2.466723) q[1];
sx q[1];
rz(-2.3900044) q[1];
rz(-pi) q[2];
rz(0.8576287) q[3];
sx q[3];
rz(-1.3751251) q[3];
sx q[3];
rz(-1.5989725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-0.98177838) q[2];
rz(1.123547) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(1.3004998) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0825901) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(2.847805) q[0];
rz(-0.44149533) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(0.77484432) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7168717) q[0];
sx q[0];
rz(-1.6443335) q[0];
sx q[0];
rz(-1.5543235) q[0];
x q[1];
rz(-0.36802937) q[2];
sx q[2];
rz(-1.5291011) q[2];
sx q[2];
rz(0.3993984) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0171417) q[1];
sx q[1];
rz(-1.972773) q[1];
sx q[1];
rz(-0.084852858) q[1];
rz(-2.9702441) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(-1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(-0.57007989) q[2];
rz(-1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.6397887) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(-0.10770527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62892249) q[0];
sx q[0];
rz(-1.218601) q[0];
sx q[0];
rz(3.0432426) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88476752) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(-1.2203072) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50385034) q[1];
sx q[1];
rz(-2.3554194) q[1];
sx q[1];
rz(-0.99650683) q[1];
rz(2.0884573) q[3];
sx q[3];
rz(-0.62873757) q[3];
sx q[3];
rz(0.79013463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(0.28953141) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(2.0026813) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(-2.1246134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6769584) q[0];
sx q[0];
rz(-0.58566739) q[0];
sx q[0];
rz(-0.93924378) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20118841) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(0.77229283) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79343866) q[1];
sx q[1];
rz(-0.70421709) q[1];
sx q[1];
rz(-2.3401295) q[1];
x q[2];
rz(-2.6050657) q[3];
sx q[3];
rz(-1.121472) q[3];
sx q[3];
rz(-0.41538737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.25062659) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780753) q[0];
sx q[0];
rz(-3.0419555) q[0];
sx q[0];
rz(-0.3936605) q[0];
rz(-pi) q[1];
rz(0.33118601) q[2];
sx q[2];
rz(-1.965431) q[2];
sx q[2];
rz(2.0926203) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68901686) q[1];
sx q[1];
rz(-0.92454443) q[1];
sx q[1];
rz(-2.8313682) q[1];
x q[2];
rz(-1.1620429) q[3];
sx q[3];
rz(-1.4882653) q[3];
sx q[3];
rz(-2.2857034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62961489) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(2.8249595) q[2];
rz(0.037242446) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(2.8102002) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(0.40922871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054571) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(-0.34988846) q[0];
rz(-0.52493237) q[2];
sx q[2];
rz(-1.8443174) q[2];
sx q[2];
rz(0.33821019) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85247773) q[1];
sx q[1];
rz(-0.58809885) q[1];
sx q[1];
rz(2.8647793) q[1];
rz(-pi) q[2];
rz(-2.6243022) q[3];
sx q[3];
rz(-1.358277) q[3];
sx q[3];
rz(-1.9912491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(-0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(-0.0099442033) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(-1.508629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976534) q[0];
sx q[0];
rz(-0.88730747) q[0];
sx q[0];
rz(-2.9336714) q[0];
rz(-pi) q[1];
rz(-3.0502003) q[2];
sx q[2];
rz(-0.95801991) q[2];
sx q[2];
rz(0.44832715) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29533169) q[1];
sx q[1];
rz(-1.8896855) q[1];
sx q[1];
rz(1.0237414) q[1];
x q[2];
rz(-0.86767254) q[3];
sx q[3];
rz(-1.613986) q[3];
sx q[3];
rz(-1.6004576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(2.7049086) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10903877) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(-1.3814829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917282) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(1.5230595) q[0];
rz(2.9276766) q[2];
sx q[2];
rz(-2.6579318) q[2];
sx q[2];
rz(2.1733401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5188462) q[1];
sx q[1];
rz(-1.985605) q[1];
sx q[1];
rz(-0.82621375) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2425209) q[3];
sx q[3];
rz(-0.74088135) q[3];
sx q[3];
rz(-1.5433943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4505724) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(2.6560442) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.8687517) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(-1.8252774) q[2];
sx q[2];
rz(-1.5263867) q[2];
sx q[2];
rz(2.9771794) q[2];
rz(0.5100526) q[3];
sx q[3];
rz(-1.4997713) q[3];
sx q[3];
rz(-0.98239519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];