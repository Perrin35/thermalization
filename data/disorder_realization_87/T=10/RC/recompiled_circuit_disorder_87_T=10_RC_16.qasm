OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(2.2044866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21298458) q[0];
sx q[0];
rz(-2.1323418) q[0];
sx q[0];
rz(1.0667849) q[0];
rz(2.4891698) q[2];
sx q[2];
rz(-0.7235652) q[2];
sx q[2];
rz(0.10318081) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5697437) q[1];
sx q[1];
rz(-1.7535926) q[1];
sx q[1];
rz(0.079449541) q[1];
rz(1.9107781) q[3];
sx q[3];
rz(-2.8570606) q[3];
sx q[3];
rz(-1.1839379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15930882) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(-0.086159555) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(-1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57698292) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.3264867) q[0];
rz(-1.8857229) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(-0.27145162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46389929) q[0];
sx q[0];
rz(-1.9672183) q[0];
sx q[0];
rz(1.0248313) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2618622) q[2];
sx q[2];
rz(-1.075282) q[2];
sx q[2];
rz(-0.77906424) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27122341) q[1];
sx q[1];
rz(-1.1961812) q[1];
sx q[1];
rz(-2.3323374) q[1];
x q[2];
rz(-1.9430964) q[3];
sx q[3];
rz(-2.3447403) q[3];
sx q[3];
rz(-2.4380395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0469971) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(-0.57717741) q[2];
rz(0.92352891) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24401027) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(0.14257167) q[0];
rz(1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(2.9325063) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3554879) q[0];
sx q[0];
rz(-0.84083637) q[0];
sx q[0];
rz(-2.3441322) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1582583) q[2];
sx q[2];
rz(-2.2670855) q[2];
sx q[2];
rz(-1.0227026) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0774539) q[1];
sx q[1];
rz(-1.9874548) q[1];
sx q[1];
rz(1.0799079) q[1];
rz(1.7861869) q[3];
sx q[3];
rz(-1.8672018) q[3];
sx q[3];
rz(-0.6682369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(1.8557619) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(-2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4362815) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(0.96281111) q[0];
rz(2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-3.1406291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22165933) q[0];
sx q[0];
rz(-2.1533305) q[0];
sx q[0];
rz(-2.4525053) q[0];
rz(-0.0041299303) q[2];
sx q[2];
rz(-0.12400707) q[2];
sx q[2];
rz(-2.4069402) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5597792) q[1];
sx q[1];
rz(-1.0043136) q[1];
sx q[1];
rz(0.15484667) q[1];
x q[2];
rz(-1.4812874) q[3];
sx q[3];
rz(-2.138391) q[3];
sx q[3];
rz(-1.9247418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0659539) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(3.0299419) q[2];
rz(2.3305437) q[3];
sx q[3];
rz(-0.45447293) q[3];
sx q[3];
rz(3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(-0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(2.8809663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97317364) q[0];
sx q[0];
rz(-1.6126313) q[0];
sx q[0];
rz(-2.9996458) q[0];
rz(-1.8243276) q[2];
sx q[2];
rz(-1.8142482) q[2];
sx q[2];
rz(-1.7340811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.534429) q[1];
sx q[1];
rz(-2.5669614) q[1];
sx q[1];
rz(0.032392153) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77633206) q[3];
sx q[3];
rz(-0.20271248) q[3];
sx q[3];
rz(-0.67938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9115209) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(-2.9122706) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(1.649958) q[0];
rz(1.0391327) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(-1.7274436) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6468069) q[0];
sx q[0];
rz(-1.2959058) q[0];
sx q[0];
rz(1.6199714) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3041777) q[2];
sx q[2];
rz(-1.068371) q[2];
sx q[2];
rz(2.7555639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8687739) q[1];
sx q[1];
rz(-1.0827853) q[1];
sx q[1];
rz(0.52656071) q[1];
x q[2];
rz(0.075501637) q[3];
sx q[3];
rz(-1.3550948) q[3];
sx q[3];
rz(-2.2942033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(-0.49267832) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3843) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(-3.0199155) q[0];
rz(-1.1514459) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(2.7391403) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6399022) q[0];
sx q[0];
rz(-0.89986378) q[0];
sx q[0];
rz(-0.34696607) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9629838) q[2];
sx q[2];
rz(-0.046769301) q[2];
sx q[2];
rz(-0.59431078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35343364) q[1];
sx q[1];
rz(-2.806059) q[1];
sx q[1];
rz(-0.56281705) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1137755) q[3];
sx q[3];
rz(-1.2226612) q[3];
sx q[3];
rz(2.5607525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.014331269) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(-2.2793615) q[2];
rz(2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(-0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7858793) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(-1.0345116) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2037172) q[0];
sx q[0];
rz(-2.8563742) q[0];
sx q[0];
rz(2.2144149) q[0];
rz(-1.6413692) q[2];
sx q[2];
rz(-1.6732054) q[2];
sx q[2];
rz(2.9533353) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5518783) q[1];
sx q[1];
rz(-2.1363746) q[1];
sx q[1];
rz(0.043885529) q[1];
rz(-pi) q[2];
rz(2.7233852) q[3];
sx q[3];
rz(-1.6697262) q[3];
sx q[3];
rz(1.6855155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(-0.77587664) q[2];
rz(2.4173229) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(-1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4416606) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(-3.124776) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(0.7787849) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22395615) q[0];
sx q[0];
rz(-0.42376873) q[0];
sx q[0];
rz(-0.64026041) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0428033) q[2];
sx q[2];
rz(-1.325377) q[2];
sx q[2];
rz(0.62308842) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9913745) q[1];
sx q[1];
rz(-0.43979904) q[1];
sx q[1];
rz(0.99779731) q[1];
rz(-pi) q[2];
rz(0.064568297) q[3];
sx q[3];
rz(-1.350292) q[3];
sx q[3];
rz(1.8622423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4497711) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(2.4251535) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.7313892) q[3];
sx q[3];
rz(1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.107782) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(2.9901436) q[0];
rz(-2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(2.418628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3176214) q[0];
sx q[0];
rz(-0.184632) q[0];
sx q[0];
rz(0.95438192) q[0];
x q[1];
rz(-0.16913551) q[2];
sx q[2];
rz(-1.9431056) q[2];
sx q[2];
rz(-1.7793836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1180229) q[1];
sx q[1];
rz(-1.5155063) q[1];
sx q[1];
rz(1.0667332) q[1];
x q[2];
rz(-1.1687752) q[3];
sx q[3];
rz(-2.3178181) q[3];
sx q[3];
rz(-2.0603501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67939776) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(2.6160713) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(-2.7450558) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(-2.6635086) q[2];
sx q[2];
rz(-2.1838084) q[2];
sx q[2];
rz(2.9391391) q[2];
rz(2.5476534) q[3];
sx q[3];
rz(-0.38729061) q[3];
sx q[3];
rz(-0.95872986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
