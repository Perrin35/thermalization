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
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(1.8523822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9292944) q[0];
sx q[0];
rz(-1.1120783) q[0];
sx q[0];
rz(0.32231583) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5640366) q[2];
sx q[2];
rz(-1.5464371) q[2];
sx q[2];
rz(1.3855795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42613712) q[1];
sx q[1];
rz(-1.8278367) q[1];
sx q[1];
rz(2.9261158) q[1];
x q[2];
rz(0.07006499) q[3];
sx q[3];
rz(-1.8383887) q[3];
sx q[3];
rz(-2.23578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5916799) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(3.0060449) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(2.360789) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(1.3134726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.388893) q[0];
sx q[0];
rz(-1.3503195) q[0];
sx q[0];
rz(1.7874784) q[0];
x q[1];
rz(2.0560527) q[2];
sx q[2];
rz(-1.8453571) q[2];
sx q[2];
rz(-0.95825125) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12570757) q[1];
sx q[1];
rz(-1.410555) q[1];
sx q[1];
rz(-1.6586152) q[1];
x q[2];
rz(-0.075606451) q[3];
sx q[3];
rz(-1.6763655) q[3];
sx q[3];
rz(2.3188863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(-0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(2.4457248) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-2.9755039) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9244708) q[0];
sx q[0];
rz(-0.32453254) q[0];
sx q[0];
rz(1.0990418) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60885749) q[2];
sx q[2];
rz(-1.1470084) q[2];
sx q[2];
rz(-2.599803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7754164) q[1];
sx q[1];
rz(-1.1300547) q[1];
sx q[1];
rz(0.52904769) q[1];
x q[2];
rz(0.2563128) q[3];
sx q[3];
rz(-2.267572) q[3];
sx q[3];
rz(-0.13845201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(-1.123547) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(-1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(-0.29378763) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.144865) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(3.0680455) q[0];
rz(-pi) q[1];
rz(-2.7735633) q[2];
sx q[2];
rz(-1.5291011) q[2];
sx q[2];
rz(-0.3993984) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1244509) q[1];
sx q[1];
rz(-1.972773) q[1];
sx q[1];
rz(-0.084852858) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1713486) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(1.1495429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.0094770771) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(-0.57007989) q[2];
rz(-1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(-0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(0.10770527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90785039) q[0];
sx q[0];
rz(-1.6630917) q[0];
sx q[0];
rz(1.924563) q[0];
rz(-2.2147419) q[2];
sx q[2];
rz(-1.1139718) q[2];
sx q[2];
rz(-2.9608375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50385034) q[1];
sx q[1];
rz(-0.78617326) q[1];
sx q[1];
rz(2.1450858) q[1];
rz(-2.7961736) q[3];
sx q[3];
rz(-2.1072227) q[3];
sx q[3];
rz(-2.9649343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058379563) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(-1.6311197) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-2.1246134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3973629) q[0];
sx q[0];
rz(-2.0332391) q[0];
sx q[0];
rz(-2.7683393) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20118841) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(0.77229283) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79343866) q[1];
sx q[1];
rz(-2.4373756) q[1];
sx q[1];
rz(2.3401295) q[1];
x q[2];
rz(2.0820886) q[3];
sx q[3];
rz(-1.0923311) q[3];
sx q[3];
rz(-0.90261501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(0.63465676) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(0.84386688) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-0.73262501) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0680925) q[0];
sx q[0];
rz(-1.4788027) q[0];
sx q[0];
rz(-1.5324701) q[0];
rz(0.33118601) q[2];
sx q[2];
rz(-1.965431) q[2];
sx q[2];
rz(-1.0489724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4525758) q[1];
sx q[1];
rz(-2.2170482) q[1];
sx q[1];
rz(2.8313682) q[1];
x q[2];
rz(1.1620429) q[3];
sx q[3];
rz(-1.4882653) q[3];
sx q[3];
rz(2.2857034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62961489) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(-0.037242446) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(2.7323639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17657875) q[0];
sx q[0];
rz(-1.2897549) q[0];
sx q[0];
rz(-0.9126419) q[0];
x q[1];
rz(-2.6166603) q[2];
sx q[2];
rz(-1.2972752) q[2];
sx q[2];
rz(0.33821019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85247773) q[1];
sx q[1];
rz(-0.58809885) q[1];
sx q[1];
rz(-0.27681338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7301634) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(3.0761449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(-0.0099442033) q[0];
rz(2.0857281) q[1];
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
rz(0.54393923) q[0];
sx q[0];
rz(-2.2542852) q[0];
sx q[0];
rz(-2.9336714) q[0];
x q[1];
rz(0.95605085) q[2];
sx q[2];
rz(-1.6455257) q[2];
sx q[2];
rz(-1.9664619) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29533169) q[1];
sx q[1];
rz(-1.2519072) q[1];
sx q[1];
rz(-2.1178513) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86767254) q[3];
sx q[3];
rz(-1.5276067) q[3];
sx q[3];
rz(-1.5411351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(-2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(2.57634) q[0];
rz(0.25282192) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(-1.7601097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0927267) q[0];
sx q[0];
rz(-2.6377569) q[0];
sx q[0];
rz(3.0548274) q[0];
x q[1];
rz(1.4597458) q[2];
sx q[2];
rz(-1.0990708) q[2];
sx q[2];
rz(-1.2088838) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5188462) q[1];
sx q[1];
rz(-1.985605) q[1];
sx q[1];
rz(-0.82621375) q[1];
rz(-pi) q[2];
rz(1.8990717) q[3];
sx q[3];
rz(-0.74088135) q[3];
sx q[3];
rz(1.5981984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(0.59664574) q[2];
rz(-0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(-0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(-1.272841) q[1];
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
