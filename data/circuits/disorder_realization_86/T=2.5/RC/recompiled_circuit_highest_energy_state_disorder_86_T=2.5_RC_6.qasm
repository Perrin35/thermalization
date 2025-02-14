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
rz(0.64492172) q[0];
sx q[0];
rz(2.6725197) q[0];
sx q[0];
rz(7.2735431) q[0];
rz(1.9682091) q[1];
sx q[1];
rz(-2.2819509) q[1];
sx q[1];
rz(-0.8492066) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.828091) q[0];
sx q[0];
rz(-1.5679616) q[0];
sx q[0];
rz(1.536973) q[0];
rz(-0.34022944) q[2];
sx q[2];
rz(-0.89982596) q[2];
sx q[2];
rz(0.036341993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.039307825) q[1];
sx q[1];
rz(-1.8370974) q[1];
sx q[1];
rz(-0.92693383) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0931271) q[3];
sx q[3];
rz(-1.4712787) q[3];
sx q[3];
rz(2.8985648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2363756) q[2];
sx q[2];
rz(-1.507501) q[2];
sx q[2];
rz(-0.53984731) q[2];
rz(1.5652462) q[3];
sx q[3];
rz(-2.7326475) q[3];
sx q[3];
rz(1.7319771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09963116) q[0];
sx q[0];
rz(-0.67622447) q[0];
sx q[0];
rz(-1.8145632) q[0];
rz(-0.78500336) q[1];
sx q[1];
rz(-1.7186586) q[1];
sx q[1];
rz(2.6410417) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6645522) q[0];
sx q[0];
rz(-2.121976) q[0];
sx q[0];
rz(1.6001979) q[0];
rz(-pi) q[1];
rz(0.18987101) q[2];
sx q[2];
rz(-1.9315757) q[2];
sx q[2];
rz(-1.9343513) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9119279) q[1];
sx q[1];
rz(-1.2357792) q[1];
sx q[1];
rz(-0.93697164) q[1];
rz(-pi) q[2];
rz(2.0359382) q[3];
sx q[3];
rz(-2.3243679) q[3];
sx q[3];
rz(-1.5508088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0565722) q[2];
sx q[2];
rz(-2.4041921) q[2];
sx q[2];
rz(2.6596587) q[2];
rz(-0.36007544) q[3];
sx q[3];
rz(-1.998338) q[3];
sx q[3];
rz(2.9329407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29927403) q[0];
sx q[0];
rz(-1.1085008) q[0];
sx q[0];
rz(-2.6639248) q[0];
rz(2.281669) q[1];
sx q[1];
rz(-2.0932902) q[1];
sx q[1];
rz(-0.28712505) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0880756) q[0];
sx q[0];
rz(-1.4447803) q[0];
sx q[0];
rz(-0.066670316) q[0];
x q[1];
rz(0.53638228) q[2];
sx q[2];
rz(-0.4066168) q[2];
sx q[2];
rz(-1.0445724) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9707253) q[1];
sx q[1];
rz(-1.5620462) q[1];
sx q[1];
rz(0.88717242) q[1];
rz(-pi) q[2];
rz(1.7805598) q[3];
sx q[3];
rz(-2.752619) q[3];
sx q[3];
rz(2.8613608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3890248) q[2];
sx q[2];
rz(-1.3501046) q[2];
sx q[2];
rz(-0.019850578) q[2];
rz(2.1974473) q[3];
sx q[3];
rz(-0.70725924) q[3];
sx q[3];
rz(-0.39785644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.95553628) q[0];
sx q[0];
rz(-0.63025403) q[0];
sx q[0];
rz(1.8408884) q[0];
rz(0.4862673) q[1];
sx q[1];
rz(-1.174289) q[1];
sx q[1];
rz(-3.1062612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080950532) q[0];
sx q[0];
rz(-1.5968906) q[0];
sx q[0];
rz(-1.9992725) q[0];
rz(-1.8043121) q[2];
sx q[2];
rz(-1.9583251) q[2];
sx q[2];
rz(-3.1336409) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.983528) q[1];
sx q[1];
rz(-3.1175545) q[1];
sx q[1];
rz(-2.9062494) q[1];
x q[2];
rz(1.0159303) q[3];
sx q[3];
rz(-1.0175704) q[3];
sx q[3];
rz(-2.826626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47762927) q[2];
sx q[2];
rz(-2.5783381) q[2];
sx q[2];
rz(-2.8832054) q[2];
rz(-0.7102617) q[3];
sx q[3];
rz(-2.5370772) q[3];
sx q[3];
rz(1.5822423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48359394) q[0];
sx q[0];
rz(-0.76376629) q[0];
sx q[0];
rz(-0.64436954) q[0];
rz(2.1517892) q[1];
sx q[1];
rz(-0.80392307) q[1];
sx q[1];
rz(-1.1092626) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.155494) q[0];
sx q[0];
rz(-1.3751831) q[0];
sx q[0];
rz(1.0367111) q[0];
rz(-pi) q[1];
rz(0.039626683) q[2];
sx q[2];
rz(-2.0840696) q[2];
sx q[2];
rz(-2.4921592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68298816) q[1];
sx q[1];
rz(-0.79958497) q[1];
sx q[1];
rz(-1.7728642) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3886834) q[3];
sx q[3];
rz(-1.782182) q[3];
sx q[3];
rz(-1.1803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41205078) q[2];
sx q[2];
rz(-1.5645626) q[2];
sx q[2];
rz(-2.5315419) q[2];
rz(-0.080502056) q[3];
sx q[3];
rz(-0.1463612) q[3];
sx q[3];
rz(3.1350873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42008156) q[0];
sx q[0];
rz(-2.2048936) q[0];
sx q[0];
rz(1.4790081) q[0];
rz(-1.1889907) q[1];
sx q[1];
rz(-2.7956796) q[1];
sx q[1];
rz(1.4422653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39665297) q[0];
sx q[0];
rz(-1.728894) q[0];
sx q[0];
rz(-0.94435662) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5601588) q[2];
sx q[2];
rz(-1.5155041) q[2];
sx q[2];
rz(-0.96382574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3644818) q[1];
sx q[1];
rz(-1.2419257) q[1];
sx q[1];
rz(2.2996705) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6053794) q[3];
sx q[3];
rz(-1.8701167) q[3];
sx q[3];
rz(-2.0423391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46001616) q[2];
sx q[2];
rz(-0.67395335) q[2];
sx q[2];
rz(-0.67503929) q[2];
rz(-0.33440822) q[3];
sx q[3];
rz(-1.6655917) q[3];
sx q[3];
rz(-0.36791754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6595031) q[0];
sx q[0];
rz(-0.98564321) q[0];
sx q[0];
rz(-0.12568812) q[0];
rz(-0.19371678) q[1];
sx q[1];
rz(-0.88231641) q[1];
sx q[1];
rz(-1.151459) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9917982) q[0];
sx q[0];
rz(-2.3201482) q[0];
sx q[0];
rz(0.69355884) q[0];
rz(1.9415683) q[2];
sx q[2];
rz(-1.8324914) q[2];
sx q[2];
rz(-0.38693025) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1035422) q[1];
sx q[1];
rz(-1.0443496) q[1];
sx q[1];
rz(2.0851233) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5732521) q[3];
sx q[3];
rz(-2.3821444) q[3];
sx q[3];
rz(-2.8041149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8573528) q[2];
sx q[2];
rz(-1.3606631) q[2];
sx q[2];
rz(-0.24687684) q[2];
rz(0.85865584) q[3];
sx q[3];
rz(-2.1571428) q[3];
sx q[3];
rz(0.27913678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3566647) q[0];
sx q[0];
rz(-0.44342884) q[0];
sx q[0];
rz(-2.8163633) q[0];
rz(2.6204956) q[1];
sx q[1];
rz(-0.63322133) q[1];
sx q[1];
rz(-0.31347832) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8756008) q[0];
sx q[0];
rz(-0.82339215) q[0];
sx q[0];
rz(-0.70988795) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.283139) q[2];
sx q[2];
rz(-0.23459841) q[2];
sx q[2];
rz(1.2755659) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4904546) q[1];
sx q[1];
rz(-0.50781194) q[1];
sx q[1];
rz(0.74393028) q[1];
rz(-0.20200396) q[3];
sx q[3];
rz(-2.083255) q[3];
sx q[3];
rz(-0.92512586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7079805) q[2];
sx q[2];
rz(-1.2182451) q[2];
sx q[2];
rz(1.8617967) q[2];
rz(-0.72569877) q[3];
sx q[3];
rz(-1.1596707) q[3];
sx q[3];
rz(-0.80997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63876605) q[0];
sx q[0];
rz(-1.9483197) q[0];
sx q[0];
rz(-2.9916812) q[0];
rz(-1.2431078) q[1];
sx q[1];
rz(-1.5877692) q[1];
sx q[1];
rz(-1.4814203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4325162) q[0];
sx q[0];
rz(-3.0750599) q[0];
sx q[0];
rz(-2.3155022) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3301351) q[2];
sx q[2];
rz(-1.5207371) q[2];
sx q[2];
rz(-0.71739774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.882355) q[1];
sx q[1];
rz(-1.9505525) q[1];
sx q[1];
rz(-1.3391375) q[1];
rz(-pi) q[2];
rz(-2.0444938) q[3];
sx q[3];
rz(-2.8682531) q[3];
sx q[3];
rz(2.7784082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82568613) q[2];
sx q[2];
rz(-1.3797727) q[2];
sx q[2];
rz(-0.97563499) q[2];
rz(2.0129096) q[3];
sx q[3];
rz(-0.55207878) q[3];
sx q[3];
rz(2.8868207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16633701) q[0];
sx q[0];
rz(-1.0003426) q[0];
sx q[0];
rz(-2.9794203) q[0];
rz(1.3316679) q[1];
sx q[1];
rz(-2.5089896) q[1];
sx q[1];
rz(3.0955637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8019288) q[0];
sx q[0];
rz(-1.5886512) q[0];
sx q[0];
rz(-0.0059203832) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.677605) q[2];
sx q[2];
rz(-1.2210238) q[2];
sx q[2];
rz(2.2217563) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.522559) q[1];
sx q[1];
rz(-1.7184034) q[1];
sx q[1];
rz(-1.9760494) q[1];
rz(2.8022553) q[3];
sx q[3];
rz(-1.4095613) q[3];
sx q[3];
rz(1.9766903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8056246) q[2];
sx q[2];
rz(-2.7569572) q[2];
sx q[2];
rz(1.9920721) q[2];
rz(0.51291054) q[3];
sx q[3];
rz(-1.3866813) q[3];
sx q[3];
rz(-0.30470595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49160663) q[0];
sx q[0];
rz(-2.2254324) q[0];
sx q[0];
rz(-1.9944763) q[0];
rz(3.0803549) q[1];
sx q[1];
rz(-1.1920659) q[1];
sx q[1];
rz(1.7658284) q[1];
rz(0.55228615) q[2];
sx q[2];
rz(-0.88197642) q[2];
sx q[2];
rz(-2.3563202) q[2];
rz(2.7313204) q[3];
sx q[3];
rz(-2.3193852) q[3];
sx q[3];
rz(1.074765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
