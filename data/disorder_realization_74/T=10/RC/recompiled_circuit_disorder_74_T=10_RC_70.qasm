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
rz(1.0097526) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(-0.90254012) q[1];
sx q[1];
rz(-1.8523822) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754817) q[0];
sx q[0];
rz(-0.55395836) q[0];
sx q[0];
rz(-1.0004811) q[0];
rz(-pi) q[1];
rz(1.5775561) q[2];
sx q[2];
rz(-1.5951556) q[2];
sx q[2];
rz(-1.7560132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2002441) q[1];
sx q[1];
rz(-1.779088) q[1];
sx q[1];
rz(-1.3079446) q[1];
x q[2];
rz(-3.0715277) q[3];
sx q[3];
rz(-1.3032039) q[3];
sx q[3];
rz(2.23578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5499128) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(-0.1401976) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(0.31301096) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(-1.82812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.27539) q[0];
sx q[0];
rz(-1.35944) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(-2.8333089) q[2];
sx q[2];
rz(-1.1051902) q[2];
sx q[2];
rz(0.75454933) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0158851) q[1];
sx q[1];
rz(-1.410555) q[1];
sx q[1];
rz(1.6586152) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95150681) q[3];
sx q[3];
rz(-0.12976876) q[3];
sx q[3];
rz(-1.4459923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(0.69586786) q[0];
rz(-0.35481915) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(-2.9755039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21712185) q[0];
sx q[0];
rz(-0.32453254) q[0];
sx q[0];
rz(-1.0990418) q[0];
rz(2.4737349) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(0.49612507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6925466) q[1];
sx q[1];
rz(-1.0967626) q[1];
sx q[1];
rz(-2.0708592) q[1];
rz(-pi) q[2];
rz(-0.2563128) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(-0.13845201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-0.98177838) q[2];
rz(-1.123547) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(2.7000973) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-0.77484432) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7168717) q[0];
sx q[0];
rz(-1.6443335) q[0];
sx q[0];
rz(-1.5543235) q[0];
rz(-pi) q[1];
rz(-1.5261126) q[2];
sx q[2];
rz(-1.2031021) q[2];
sx q[2];
rz(1.9862663) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91037726) q[1];
sx q[1];
rz(-0.41035715) q[1];
sx q[1];
rz(-1.7675722) q[1];
rz(-pi) q[2];
rz(-0.1713486) q[3];
sx q[3];
rz(-1.315794) q[3];
sx q[3];
rz(-1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(2.5715128) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(0.18138012) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(0.10770527) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2337423) q[0];
sx q[0];
rz(-1.478501) q[0];
sx q[0];
rz(-1.2170296) q[0];
x q[1];
rz(0.88476752) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(1.2203072) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6458466) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(-2.2699725) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0072458) q[3];
sx q[3];
rz(-1.8661024) q[3];
sx q[3];
rz(-1.2122648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(0.28953141) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3973629) q[0];
sx q[0];
rz(-2.0332391) q[0];
sx q[0];
rz(-0.37325333) q[0];
rz(-pi) q[1];
rz(-0.20118841) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(-2.3692998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0000227) q[1];
sx q[1];
rz(-1.1035898) q[1];
sx q[1];
rz(-1.0228553) q[1];
rz(-pi) q[2];
rz(-0.75628144) q[3];
sx q[3];
rz(-0.68538266) q[3];
sx q[3];
rz(-1.7862198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54414526) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(2.0641816) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(-0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25062659) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-2.2977258) q[0];
rz(1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4635173) q[0];
sx q[0];
rz(-0.099637195) q[0];
sx q[0];
rz(2.7479322) q[0];
x q[1];
rz(-0.90791038) q[2];
sx q[2];
rz(-2.6320576) q[2];
sx q[2];
rz(0.31928911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9417604) q[1];
sx q[1];
rz(-0.70711771) q[1];
sx q[1];
rz(1.1862399) q[1];
rz(-pi) q[2];
rz(-1.1620429) q[3];
sx q[3];
rz(-1.6533274) q[3];
sx q[3];
rz(-0.85588928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379631) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(-0.40922871) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054571) q[0];
sx q[0];
rz(-2.1989609) q[0];
sx q[0];
rz(0.34988846) q[0];
rz(-pi) q[1];
rz(-2.6312469) q[2];
sx q[2];
rz(-2.5556457) q[2];
sx q[2];
rz(0.79615359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9600535) q[1];
sx q[1];
rz(-1.0078733) q[1];
sx q[1];
rz(1.7510508) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41142923) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(-0.065447741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(-2.5869353) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(-0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(3.1316485) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(1.508629) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976534) q[0];
sx q[0];
rz(-2.2542852) q[0];
sx q[0];
rz(-2.9336714) q[0];
x q[1];
rz(-1.4417068) q[2];
sx q[2];
rz(-0.61868762) q[2];
sx q[2];
rz(-0.29030756) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.75084) q[1];
sx q[1];
rz(-2.5166593) q[1];
sx q[1];
rz(-1.0052488) q[1];
rz(-pi) q[2];
rz(-0.86767254) q[3];
sx q[3];
rz(-1.613986) q[3];
sx q[3];
rz(1.5411351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9987954) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(0.43668401) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(1.0288303) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10903877) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(2.57634) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(-1.7601097) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59795838) q[0];
sx q[0];
rz(-1.6126452) q[0];
sx q[0];
rz(-0.50224395) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21391602) q[2];
sx q[2];
rz(-0.48366085) q[2];
sx q[2];
rz(-0.96825251) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5188462) q[1];
sx q[1];
rz(-1.1559876) q[1];
sx q[1];
rz(0.82621375) q[1];
x q[2];
rz(0.28678203) q[3];
sx q[3];
rz(-0.87773318) q[3];
sx q[3];
rz(-1.1657438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(0.027464494) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.8687517) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(-1.3163153) q[2];
sx q[2];
rz(-1.6152059) q[2];
sx q[2];
rz(-0.16441328) q[2];
rz(0.14470312) q[3];
sx q[3];
rz(-0.51454138) q[3];
sx q[3];
rz(0.71458057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];