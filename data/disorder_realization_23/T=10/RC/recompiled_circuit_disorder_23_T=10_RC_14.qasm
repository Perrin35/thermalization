OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(7.6927778) q[0];
sx q[0];
rz(11.132244) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7941064) q[0];
sx q[0];
rz(-2.5972164) q[0];
sx q[0];
rz(2.0967336) q[0];
x q[1];
rz(1.4346052) q[2];
sx q[2];
rz(-1.681466) q[2];
sx q[2];
rz(-1.1510804) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1474255) q[1];
sx q[1];
rz(-0.52417437) q[1];
sx q[1];
rz(-2.0736573) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22085805) q[3];
sx q[3];
rz(-1.1410032) q[3];
sx q[3];
rz(0.62648279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-2.3036172) q[2];
rz(0.49301246) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(0.078991927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74719602) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.863742) q[0];
rz(2.9648119) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(2.7094254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55789253) q[0];
sx q[0];
rz(-1.0520792) q[0];
sx q[0];
rz(-1.8899263) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.997666) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(-1.2183684) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9464311) q[1];
sx q[1];
rz(-2.6841607) q[1];
sx q[1];
rz(1.8381579) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5315941) q[3];
sx q[3];
rz(-0.56534846) q[3];
sx q[3];
rz(2.431734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5923578) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(1.6529282) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(-0.506385) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8538118) q[0];
sx q[0];
rz(-1.3940485) q[0];
sx q[0];
rz(1.214083) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73451368) q[2];
sx q[2];
rz(-1.7506415) q[2];
sx q[2];
rz(2.432446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.312078) q[1];
sx q[1];
rz(-1.8784338) q[1];
sx q[1];
rz(-2.679146) q[1];
x q[2];
rz(2.0154325) q[3];
sx q[3];
rz(-2.7362842) q[3];
sx q[3];
rz(2.0733548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69904077) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(-0.58602035) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(-1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9716924) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(2.3024094) q[0];
rz(-0.025578586) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(-1.5485839) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73491053) q[0];
sx q[0];
rz(-1.6618177) q[0];
sx q[0];
rz(2.0216366) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3392802) q[2];
sx q[2];
rz(-2.1244086) q[2];
sx q[2];
rz(-2.4137036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5316481) q[1];
sx q[1];
rz(-0.87024401) q[1];
sx q[1];
rz(1.0611666) q[1];
x q[2];
rz(-2.7110093) q[3];
sx q[3];
rz(-1.8292556) q[3];
sx q[3];
rz(-0.57627288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(3.0419066) q[2];
rz(-0.95885197) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-0.25594041) q[0];
rz(2.6804965) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(0.76006132) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1003935) q[0];
sx q[0];
rz(-1.4677591) q[0];
sx q[0];
rz(-3.0242007) q[0];
rz(-1.4270093) q[2];
sx q[2];
rz(-0.41848768) q[2];
sx q[2];
rz(-0.33868044) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1000423) q[1];
sx q[1];
rz(-1.6987545) q[1];
sx q[1];
rz(-0.12771878) q[1];
x q[2];
rz(2.2171668) q[3];
sx q[3];
rz(-1.2741718) q[3];
sx q[3];
rz(2.7588206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5710859) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(-0.6742397) q[2];
rz(0.21480602) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-0.017344346) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6102585) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(-1.1791139) q[0];
rz(-2.9367661) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(2.0746322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63283352) q[0];
sx q[0];
rz(-1.5501223) q[0];
sx q[0];
rz(1.5563006) q[0];
rz(-pi) q[1];
rz(0.42491575) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(-1.0645107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21642099) q[1];
sx q[1];
rz(-0.74837084) q[1];
sx q[1];
rz(1.640635) q[1];
rz(-pi) q[2];
rz(-0.1596998) q[3];
sx q[3];
rz(-2.3387863) q[3];
sx q[3];
rz(-1.6646977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71010464) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(2.5816494) q[2];
rz(-0.7263178) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(-1.7204826) q[0];
rz(0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(-2.2834159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439529) q[0];
sx q[0];
rz(-1.7174935) q[0];
sx q[0];
rz(0.12978817) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0965205) q[2];
sx q[2];
rz(-1.1606693) q[2];
sx q[2];
rz(0.28442597) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8691751) q[1];
sx q[1];
rz(-0.032748001) q[1];
sx q[1];
rz(2.6564024) q[1];
rz(-pi) q[2];
rz(-2.6303597) q[3];
sx q[3];
rz(-2.5066262) q[3];
sx q[3];
rz(2.883203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.3254335) q[2];
rz(-2.251513) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(0.98852283) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637852) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23096202) q[0];
sx q[0];
rz(-2.1719143) q[0];
sx q[0];
rz(1.0966976) q[0];
rz(-pi) q[1];
x q[1];
rz(0.026272341) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(-2.9422613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6264682) q[1];
sx q[1];
rz(-1.5212458) q[1];
sx q[1];
rz(-1.3888166) q[1];
rz(-pi) q[2];
rz(-2.5958063) q[3];
sx q[3];
rz(-1.3457527) q[3];
sx q[3];
rz(-0.033566098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.65138856) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(0.46869579) q[2];
rz(1.9474585) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63012183) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(1.2619031) q[0];
rz(2.966554) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(1.6040241) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14311929) q[0];
sx q[0];
rz(-2.7612918) q[0];
sx q[0];
rz(-2.2024676) q[0];
rz(-0.41861694) q[2];
sx q[2];
rz(-1.6659684) q[2];
sx q[2];
rz(1.0375432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9201587) q[1];
sx q[1];
rz(-2.4281574) q[1];
sx q[1];
rz(0.17381298) q[1];
x q[2];
rz(-2.2213307) q[3];
sx q[3];
rz(-2.5839621) q[3];
sx q[3];
rz(-0.44616163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1016772) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(0.76134479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.3392357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(2.9246869) q[0];
rz(-0.63198173) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(2.1868618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.117884) q[0];
sx q[0];
rz(-0.81239359) q[0];
sx q[0];
rz(-0.448416) q[0];
x q[1];
rz(-2.5920792) q[2];
sx q[2];
rz(-1.9351442) q[2];
sx q[2];
rz(0.074631045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2659457) q[1];
sx q[1];
rz(-1.4946283) q[1];
sx q[1];
rz(-0.45653685) q[1];
rz(-pi) q[2];
rz(-2.2050489) q[3];
sx q[3];
rz(-1.9359971) q[3];
sx q[3];
rz(1.6983502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(-2.1968502) q[2];
rz(0.38481209) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(-2.184536) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5007297) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(0.8846994) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-0.89818556) q[2];
sx q[2];
rz(-0.85815103) q[2];
sx q[2];
rz(1.382538) q[2];
rz(-1.9839722) q[3];
sx q[3];
rz(-2.6631841) q[3];
sx q[3];
rz(-0.14315179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
