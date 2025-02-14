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
rz(1.2394387) q[0];
sx q[0];
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(2.6091726) q[1];
sx q[1];
rz(6.6611023) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53055313) q[0];
sx q[0];
rz(-1.5612649) q[0];
sx q[0];
rz(-0.00044374142) q[0];
rz(0.67899668) q[2];
sx q[2];
rz(-1.2831935) q[2];
sx q[2];
rz(1.2305413) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18641414) q[1];
sx q[1];
rz(-0.70611533) q[1];
sx q[1];
rz(2.7104055) q[1];
rz(-pi) q[2];
rz(2.3610695) q[3];
sx q[3];
rz(-0.98376432) q[3];
sx q[3];
rz(2.7450456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8523031) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(0.079785384) q[2];
rz(-2.1763109) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6642283) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(3.0526414) q[0];
rz(1.7695919) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(1.1044097) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4084642) q[0];
sx q[0];
rz(-0.95703546) q[0];
sx q[0];
rz(1.2888786) q[0];
x q[1];
rz(2.0599819) q[2];
sx q[2];
rz(-2.2913165) q[2];
sx q[2];
rz(2.1921076) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1762878) q[1];
sx q[1];
rz(-1.4215648) q[1];
sx q[1];
rz(-2.01841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6703963) q[3];
sx q[3];
rz(-2.1079113) q[3];
sx q[3];
rz(0.84505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.264512) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(1.6748927) q[2];
rz(3.1125715) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(-2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324209) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(1.4311283) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(0.40036449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63567296) q[0];
sx q[0];
rz(-0.95209661) q[0];
sx q[0];
rz(2.1257504) q[0];
rz(-pi) q[1];
rz(-2.8334191) q[2];
sx q[2];
rz(-2.4314667) q[2];
sx q[2];
rz(-1.1504722) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.4246042) q[1];
sx q[1];
rz(-1.7300055) q[1];
sx q[1];
rz(0.47674322) q[1];
x q[2];
rz(2.4147846) q[3];
sx q[3];
rz(-2.6977959) q[3];
sx q[3];
rz(-0.82647317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4572738) q[2];
sx q[2];
rz(-1.532734) q[2];
sx q[2];
rz(-0.45480967) q[2];
rz(1.8999892) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(-2.4212867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860745) q[0];
sx q[0];
rz(-1.8121413) q[0];
sx q[0];
rz(3.1410826) q[0];
rz(0.60091248) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(-0.14437637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7847608) q[0];
sx q[0];
rz(-2.5481173) q[0];
sx q[0];
rz(-1.9463825) q[0];
x q[1];
rz(-1.7208485) q[2];
sx q[2];
rz(-1.8215107) q[2];
sx q[2];
rz(-1.5423519) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3444654) q[1];
sx q[1];
rz(-0.68736156) q[1];
sx q[1];
rz(-0.3631773) q[1];
rz(-pi) q[2];
rz(2.9668429) q[3];
sx q[3];
rz(-0.55395444) q[3];
sx q[3];
rz(2.9308211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1986177) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(-2.1790738) q[2];
rz(1.6019542) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2365504) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(1.9566253) q[0];
rz(-0.22625893) q[1];
sx q[1];
rz(-2.2626651) q[1];
sx q[1];
rz(-2.102898) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7554692) q[0];
sx q[0];
rz(-2.1573108) q[0];
sx q[0];
rz(1.4727794) q[0];
rz(-pi) q[1];
rz(2.4034385) q[2];
sx q[2];
rz(-2.0655491) q[2];
sx q[2];
rz(-3.11657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24327899) q[1];
sx q[1];
rz(-2.2479575) q[1];
sx q[1];
rz(-0.89666287) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.709278) q[3];
sx q[3];
rz(-2.6977728) q[3];
sx q[3];
rz(0.5139155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.431939) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(0.31614885) q[2];
rz(1.4656674) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(1.3795615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(2.6370711) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(2.0380518) q[0];
rz(2.6612813) q[1];
sx q[1];
rz(-2.4791398) q[1];
sx q[1];
rz(0.59741098) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7561853) q[0];
sx q[0];
rz(-0.28235897) q[0];
sx q[0];
rz(-0.83448164) q[0];
rz(-pi) q[1];
rz(-2.8812203) q[2];
sx q[2];
rz(-0.7625167) q[2];
sx q[2];
rz(1.1077566) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0874071) q[1];
sx q[1];
rz(-0.30255908) q[1];
sx q[1];
rz(2.5564479) q[1];
rz(1.8036929) q[3];
sx q[3];
rz(-0.38023708) q[3];
sx q[3];
rz(1.5429131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20208134) q[2];
sx q[2];
rz(-1.6263522) q[2];
sx q[2];
rz(1.0763947) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(-0.49016652) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550605) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(-0.038473815) q[0];
rz(-0.064727457) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(-2.9077392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8412178) q[0];
sx q[0];
rz(-2.8511957) q[0];
sx q[0];
rz(0.69926326) q[0];
x q[1];
rz(1.2042768) q[2];
sx q[2];
rz(-1.2784174) q[2];
sx q[2];
rz(-2.2077843) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55091399) q[1];
sx q[1];
rz(-1.6603025) q[1];
sx q[1];
rz(-2.4924335) q[1];
rz(-pi) q[2];
rz(1.6285628) q[3];
sx q[3];
rz(-1.7004629) q[3];
sx q[3];
rz(2.4313789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-2.5433507) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(-0.84754506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.7365731) q[0];
sx q[0];
rz(-2.2990655) q[0];
sx q[0];
rz(2.8952428) q[0];
rz(-1.8440638) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(2.6447703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9333736) q[0];
sx q[0];
rz(-2.1849306) q[0];
sx q[0];
rz(1.0557515) q[0];
rz(-pi) q[1];
rz(1.7165347) q[2];
sx q[2];
rz(-2.2706804) q[2];
sx q[2];
rz(0.041797195) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5766746) q[1];
sx q[1];
rz(-1.0905301) q[1];
sx q[1];
rz(0.26580055) q[1];
rz(-pi) q[2];
rz(-1.0786177) q[3];
sx q[3];
rz(-1.7075305) q[3];
sx q[3];
rz(-2.0355952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42314998) q[2];
sx q[2];
rz(-0.53360525) q[2];
sx q[2];
rz(1.9971087) q[2];
rz(-1.9874969) q[3];
sx q[3];
rz(-1.7172979) q[3];
sx q[3];
rz(2.9437039) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474739) q[0];
sx q[0];
rz(-2.9070774) q[0];
sx q[0];
rz(-0.06614729) q[0];
rz(1.1478395) q[1];
sx q[1];
rz(-1.7639953) q[1];
sx q[1];
rz(-2.537421) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19949958) q[0];
sx q[0];
rz(-1.8460994) q[0];
sx q[0];
rz(-1.3574187) q[0];
rz(-pi) q[1];
rz(2.8762875) q[2];
sx q[2];
rz(-2.5992706) q[2];
sx q[2];
rz(-0.10077439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38018885) q[1];
sx q[1];
rz(-1.0127002) q[1];
sx q[1];
rz(-1.7978884) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2630713) q[3];
sx q[3];
rz(-2.1067348) q[3];
sx q[3];
rz(-1.5978447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26312795) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(-0.43295941) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(-1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94806725) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(1.8684335) q[0];
rz(-2.676447) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(2.926631) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8037655) q[0];
sx q[0];
rz(-1.318232) q[0];
sx q[0];
rz(2.4830677) q[0];
x q[1];
rz(-0.6647756) q[2];
sx q[2];
rz(-2.4996346) q[2];
sx q[2];
rz(0.73981111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7219639) q[1];
sx q[1];
rz(-0.89744324) q[1];
sx q[1];
rz(-2.2339905) q[1];
rz(-pi) q[2];
rz(-2.5796579) q[3];
sx q[3];
rz(-1.3519545) q[3];
sx q[3];
rz(-2.9666025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8421858) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(0.95883933) q[2];
rz(1.9793319) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(1.6963262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.5398298) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(0.70855793) q[1];
sx q[1];
rz(-2.9500912) q[1];
sx q[1];
rz(-0.74831829) q[1];
rz(0.49025771) q[2];
sx q[2];
rz(-2.5210862) q[2];
sx q[2];
rz(1.6747337) q[2];
rz(-2.419653) q[3];
sx q[3];
rz(-1.4293213) q[3];
sx q[3];
rz(-2.8534129) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
