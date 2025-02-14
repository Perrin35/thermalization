OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.4247704) q[0];
sx q[0];
rz(4.7228887) q[0];
sx q[0];
rz(6.9520998) q[0];
rz(0.34612292) q[1];
sx q[1];
rz(-2.3200413) q[1];
sx q[1];
rz(-0.93924826) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4209983) q[0];
sx q[0];
rz(-1.2828553) q[0];
sx q[0];
rz(1.1714515) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.16667) q[2];
sx q[2];
rz(-1.5277281) q[2];
sx q[2];
rz(-1.7801746) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85744845) q[1];
sx q[1];
rz(-1.3054056) q[1];
sx q[1];
rz(-1.1737203) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40832728) q[3];
sx q[3];
rz(-2.3882867) q[3];
sx q[3];
rz(-1.4391132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5758489) q[2];
sx q[2];
rz(-1.0172903) q[2];
sx q[2];
rz(3.1080833) q[2];
rz(1.9401898) q[3];
sx q[3];
rz(-1.3591432) q[3];
sx q[3];
rz(-2.0638154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(3.1101542) q[0];
sx q[0];
rz(-0.25122508) q[0];
sx q[0];
rz(-2.187619) q[0];
rz(1.6961478) q[1];
sx q[1];
rz(-1.0547538) q[1];
sx q[1];
rz(1.9704069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32098371) q[0];
sx q[0];
rz(-0.96555987) q[0];
sx q[0];
rz(0.73802276) q[0];
x q[1];
rz(1.2759802) q[2];
sx q[2];
rz(-2.1871242) q[2];
sx q[2];
rz(0.64586879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18942434) q[1];
sx q[1];
rz(-1.599424) q[1];
sx q[1];
rz(2.9131123) q[1];
rz(-2.21255) q[3];
sx q[3];
rz(-0.52601846) q[3];
sx q[3];
rz(-1.075901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.235432) q[2];
sx q[2];
rz(-1.4305328) q[2];
sx q[2];
rz(-1.1492427) q[2];
rz(-0.033128459) q[3];
sx q[3];
rz(-1.5002316) q[3];
sx q[3];
rz(-2.5455425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0572492) q[0];
sx q[0];
rz(-0.79279041) q[0];
sx q[0];
rz(-2.1113915) q[0];
rz(1.239981) q[1];
sx q[1];
rz(-1.3571309) q[1];
sx q[1];
rz(-2.1485567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1112615) q[0];
sx q[0];
rz(-1.5227888) q[0];
sx q[0];
rz(-2.6859849) q[0];
rz(-pi) q[1];
rz(-1.4517752) q[2];
sx q[2];
rz(-0.90672158) q[2];
sx q[2];
rz(1.9805976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5883623) q[1];
sx q[1];
rz(-2.2892286) q[1];
sx q[1];
rz(0.78939446) q[1];
rz(-2.7168324) q[3];
sx q[3];
rz(-0.83213193) q[3];
sx q[3];
rz(-0.55603851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21732907) q[2];
sx q[2];
rz(-1.1161048) q[2];
sx q[2];
rz(0.35476157) q[2];
rz(0.99700704) q[3];
sx q[3];
rz(-1.6544147) q[3];
sx q[3];
rz(-0.94064373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16911258) q[0];
sx q[0];
rz(-1.4262154) q[0];
sx q[0];
rz(-2.0347563) q[0];
rz(2.2726982) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(2.4443764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49299846) q[0];
sx q[0];
rz(-1.0169858) q[0];
sx q[0];
rz(1.336344) q[0];
rz(-pi) q[1];
rz(0.72202335) q[2];
sx q[2];
rz(-1.486384) q[2];
sx q[2];
rz(1.9868324) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47095151) q[1];
sx q[1];
rz(-1.1837675) q[1];
sx q[1];
rz(2.7426038) q[1];
rz(-pi) q[2];
rz(-0.79037621) q[3];
sx q[3];
rz(-0.034878313) q[3];
sx q[3];
rz(2.6863475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8640459) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(2.8475658) q[2];
rz(1.0634408) q[3];
sx q[3];
rz(-1.682351) q[3];
sx q[3];
rz(2.3991876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251637) q[0];
sx q[0];
rz(-2.9603781) q[0];
sx q[0];
rz(2.678405) q[0];
rz(-0.40395346) q[1];
sx q[1];
rz(-1.5084167) q[1];
sx q[1];
rz(0.24872669) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8489784) q[0];
sx q[0];
rz(-0.19751829) q[0];
sx q[0];
rz(-1.5621444) q[0];
rz(-pi) q[1];
rz(-1.312448) q[2];
sx q[2];
rz(-2.097192) q[2];
sx q[2];
rz(-2.0451982) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46649964) q[1];
sx q[1];
rz(-2.6394301) q[1];
sx q[1];
rz(-2.2999022) q[1];
x q[2];
rz(2.4234613) q[3];
sx q[3];
rz(-2.2635824) q[3];
sx q[3];
rz(2.3819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2948461) q[2];
sx q[2];
rz(-1.3174026) q[2];
sx q[2];
rz(1.741629) q[2];
rz(1.2830265) q[3];
sx q[3];
rz(-1.7581519) q[3];
sx q[3];
rz(-2.9156901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1840709) q[0];
sx q[0];
rz(-2.3029843) q[0];
sx q[0];
rz(-0.34570178) q[0];
rz(-0.050994571) q[1];
sx q[1];
rz(-1.2268927) q[1];
sx q[1];
rz(-2.3695703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2564023) q[0];
sx q[0];
rz(-1.9213543) q[0];
sx q[0];
rz(1.3023443) q[0];
x q[1];
rz(-1.4441024) q[2];
sx q[2];
rz(-2.4856353) q[2];
sx q[2];
rz(-0.7743338) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.009338) q[1];
sx q[1];
rz(-1.8937832) q[1];
sx q[1];
rz(0.96829523) q[1];
rz(-pi) q[2];
rz(1.0994301) q[3];
sx q[3];
rz(-2.5623218) q[3];
sx q[3];
rz(1.260965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2766075) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(2.9841606) q[2];
rz(-2.7031247) q[3];
sx q[3];
rz(-2.4003568) q[3];
sx q[3];
rz(2.1854775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.081414374) q[0];
sx q[0];
rz(-1.0670476) q[0];
sx q[0];
rz(-9/(11*pi)) q[0];
rz(-0.57834894) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(-2.9885805) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6575847) q[0];
sx q[0];
rz(-1.0562684) q[0];
sx q[0];
rz(2.5637676) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4744841) q[2];
sx q[2];
rz(-1.4862712) q[2];
sx q[2];
rz(0.21749228) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.93219705) q[1];
sx q[1];
rz(-1.160566) q[1];
sx q[1];
rz(-0.61102976) q[1];
x q[2];
rz(1.1122996) q[3];
sx q[3];
rz(-1.8613437) q[3];
sx q[3];
rz(-0.19536881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1795307) q[2];
sx q[2];
rz(-3.0425368) q[2];
sx q[2];
rz(0.2641826) q[2];
rz(-1.3029441) q[3];
sx q[3];
rz(-1.7641726) q[3];
sx q[3];
rz(-2.0343659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1846979) q[0];
sx q[0];
rz(-0.95320025) q[0];
sx q[0];
rz(1.4554998) q[0];
rz(-0.9616583) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(0.89967322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46612211) q[0];
sx q[0];
rz(-0.5044901) q[0];
sx q[0];
rz(-0.44874189) q[0];
rz(-pi) q[1];
rz(-1.3191965) q[2];
sx q[2];
rz(-2.4008022) q[2];
sx q[2];
rz(1.0902001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2921999) q[1];
sx q[1];
rz(-1.6136886) q[1];
sx q[1];
rz(-0.98341839) q[1];
rz(-pi) q[2];
rz(-0.03735383) q[3];
sx q[3];
rz(-1.787623) q[3];
sx q[3];
rz(-1.9268394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4453033) q[2];
sx q[2];
rz(-2.5407007) q[2];
sx q[2];
rz(-0.75330934) q[2];
rz(-0.2937915) q[3];
sx q[3];
rz(-1.1283504) q[3];
sx q[3];
rz(-2.0297348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927354) q[0];
sx q[0];
rz(-0.98965544) q[0];
sx q[0];
rz(-0.85187546) q[0];
rz(1.4082255) q[1];
sx q[1];
rz(-1.4829166) q[1];
sx q[1];
rz(-0.3826938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2792042) q[0];
sx q[0];
rz(-0.48868256) q[0];
sx q[0];
rz(-0.53532289) q[0];
rz(-pi) q[1];
rz(-0.23613249) q[2];
sx q[2];
rz(-2.3575767) q[2];
sx q[2];
rz(1.2003984) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.574461) q[1];
sx q[1];
rz(-1.2700915) q[1];
sx q[1];
rz(-1.6733132) q[1];
rz(-pi) q[2];
rz(2.9858573) q[3];
sx q[3];
rz(-0.49113516) q[3];
sx q[3];
rz(-0.30065003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3771628) q[2];
sx q[2];
rz(-0.96401507) q[2];
sx q[2];
rz(-2.5386179) q[2];
rz(1.0682586) q[3];
sx q[3];
rz(-1.7030741) q[3];
sx q[3];
rz(-0.8409797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0076440796) q[0];
sx q[0];
rz(-1.4643865) q[0];
sx q[0];
rz(2.7879047) q[0];
rz(-2.0091281) q[1];
sx q[1];
rz(-2.1734889) q[1];
sx q[1];
rz(2.7521334) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2089473) q[0];
sx q[0];
rz(-1.07558) q[0];
sx q[0];
rz(1.3071278) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10305291) q[2];
sx q[2];
rz(-1.3259058) q[2];
sx q[2];
rz(-2.6239397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82853854) q[1];
sx q[1];
rz(-1.8450929) q[1];
sx q[1];
rz(-2.6649339) q[1];
rz(-pi) q[2];
rz(-0.52822379) q[3];
sx q[3];
rz(-1.6162795) q[3];
sx q[3];
rz(-0.32873617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13122095) q[2];
sx q[2];
rz(-1.5067357) q[2];
sx q[2];
rz(-1.1466675) q[2];
rz(1.5366588) q[3];
sx q[3];
rz(-2.6585572) q[3];
sx q[3];
rz(-0.35227942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1491886) q[0];
sx q[0];
rz(-0.70269361) q[0];
sx q[0];
rz(0.50444034) q[0];
rz(-3.0814677) q[1];
sx q[1];
rz(-2.6838214) q[1];
sx q[1];
rz(2.6837742) q[1];
rz(-1.4627152) q[2];
sx q[2];
rz(-1.4975794) q[2];
sx q[2];
rz(0.076051408) q[2];
rz(1.6144013) q[3];
sx q[3];
rz(-0.674896) q[3];
sx q[3];
rz(1.9622635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
