OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(4.1133147) q[1];
sx q[1];
rz(9.9749554) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101947) q[0];
sx q[0];
rz(-2.9203127) q[0];
sx q[0];
rz(-1.4472603) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43843856) q[2];
sx q[2];
rz(-2.0619832) q[2];
sx q[2];
rz(-0.05471281) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5397415) q[1];
sx q[1];
rz(-1.9508233) q[1];
sx q[1];
rz(0.7976346) q[1];
rz(0.95991858) q[3];
sx q[3];
rz(-1.6392518) q[3];
sx q[3];
rz(-1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(2.7602592) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2663651) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.5126022) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0727901) q[0];
sx q[0];
rz(-1.2501984) q[0];
sx q[0];
rz(-0.94706236) q[0];
rz(-0.28378758) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(-0.46797215) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4334129) q[1];
sx q[1];
rz(-2.4509894) q[1];
sx q[1];
rz(-1.3604926) q[1];
rz(-pi) q[2];
rz(-2.2643331) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(3.0960494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(-0.80336037) q[2];
rz(1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(2.846068) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(0.74584109) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.196516) q[0];
sx q[0];
rz(-1.7917504) q[0];
sx q[0];
rz(-0.093284746) q[0];
x q[1];
rz(0.27772851) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(-0.03253983) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5902139) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(0.53932921) q[1];
x q[2];
rz(2.6882719) q[3];
sx q[3];
rz(-1.3461777) q[3];
sx q[3];
rz(-2.2090705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(-1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(2.8682958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057987) q[0];
sx q[0];
rz(-1.6596284) q[0];
sx q[0];
rz(1.2466794) q[0];
x q[1];
rz(-0.64812135) q[2];
sx q[2];
rz(-1.1037877) q[2];
sx q[2];
rz(-1.8304706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4218688) q[1];
sx q[1];
rz(-2.8641652) q[1];
sx q[1];
rz(-1.006702) q[1];
x q[2];
rz(-0.64951879) q[3];
sx q[3];
rz(-1.7139385) q[3];
sx q[3];
rz(0.79917819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(-3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(0.98714978) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8866681) q[0];
sx q[0];
rz(-2.3912171) q[0];
sx q[0];
rz(2.1819941) q[0];
rz(-pi) q[1];
rz(-0.64381386) q[2];
sx q[2];
rz(-1.5783196) q[2];
sx q[2];
rz(1.3155589) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8122711) q[1];
sx q[1];
rz(-1.2580039) q[1];
sx q[1];
rz(0.66004628) q[1];
rz(-0.042111245) q[3];
sx q[3];
rz(-1.1782421) q[3];
sx q[3];
rz(2.7819355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8905028) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(0.063407272) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(-0.791839) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56023843) q[0];
sx q[0];
rz(-1.2890153) q[0];
sx q[0];
rz(-3.1262585) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2811786) q[2];
sx q[2];
rz(-2.1589303) q[2];
sx q[2];
rz(0.16298018) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8705604) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(1.51103) q[1];
rz(-0.91026129) q[3];
sx q[3];
rz(-2.913774) q[3];
sx q[3];
rz(-1.4753301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(1.6714913) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(2.9260013) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(3.1380222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647033) q[0];
sx q[0];
rz(-1.2029552) q[0];
sx q[0];
rz(-0.44643114) q[0];
x q[1];
rz(1.1793169) q[2];
sx q[2];
rz(-0.85748312) q[2];
sx q[2];
rz(1.5228524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4184985) q[1];
sx q[1];
rz(-1.8685891) q[1];
sx q[1];
rz(-0.94350068) q[1];
x q[2];
rz(1.7690897) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(-2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51628095) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(1.7238808) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-2.4628941) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40020254) q[0];
sx q[0];
rz(-0.42660248) q[0];
sx q[0];
rz(2.461117) q[0];
rz(-pi) q[1];
rz(2.0754201) q[2];
sx q[2];
rz(-2.3373211) q[2];
sx q[2];
rz(1.417516) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8933967) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(3.0688973) q[1];
rz(1.4123165) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(-2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(2.2237681) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(2.6146467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567112) q[0];
sx q[0];
rz(-1.3538133) q[0];
sx q[0];
rz(2.7816539) q[0];
x q[1];
rz(-1.2277463) q[2];
sx q[2];
rz(-0.68708778) q[2];
sx q[2];
rz(-2.3234141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37104169) q[1];
sx q[1];
rz(-2.7275804) q[1];
sx q[1];
rz(-2.4310334) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6960877) q[3];
sx q[3];
rz(-1.8337436) q[3];
sx q[3];
rz(1.6642237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-2.7808166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26786131) q[0];
sx q[0];
rz(-2.972313) q[0];
sx q[0];
rz(1.4219567) q[0];
rz(-pi) q[1];
rz(-3.0214494) q[2];
sx q[2];
rz(-1.7842245) q[2];
sx q[2];
rz(-3.0394768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4611778) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(1.3780891) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1694069) q[3];
sx q[3];
rz(-1.425404) q[3];
sx q[3];
rz(1.3225079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(0.94341755) q[2];
rz(2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671134) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(2.0104682) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(0.75136649) q[3];
sx q[3];
rz(-1.02117) q[3];
sx q[3];
rz(1.1985967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];