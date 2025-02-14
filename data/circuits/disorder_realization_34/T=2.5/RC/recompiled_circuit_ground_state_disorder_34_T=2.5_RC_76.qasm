OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3736149) q[0];
sx q[0];
rz(-0.64426214) q[0];
sx q[0];
rz(2.2777519) q[0];
rz(1.5732425) q[1];
sx q[1];
rz(-1.9116126) q[1];
sx q[1];
rz(0.65043515) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7614131) q[0];
sx q[0];
rz(-2.5963547) q[0];
sx q[0];
rz(-2.7392967) q[0];
rz(-pi) q[1];
rz(2.8312248) q[2];
sx q[2];
rz(-1.5987607) q[2];
sx q[2];
rz(-1.8486277) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0155458) q[1];
sx q[1];
rz(-2.0728557) q[1];
sx q[1];
rz(2.253336) q[1];
x q[2];
rz(0.1273472) q[3];
sx q[3];
rz(-1.6771924) q[3];
sx q[3];
rz(1.1945386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2245366) q[2];
sx q[2];
rz(-1.3755596) q[2];
sx q[2];
rz(-1.3517693) q[2];
rz(-2.3228862) q[3];
sx q[3];
rz(-1.4592146) q[3];
sx q[3];
rz(-1.590135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025892) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(0.57574058) q[0];
rz(-2.2585244) q[1];
sx q[1];
rz(-1.539361) q[1];
sx q[1];
rz(1.8227122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91399804) q[0];
sx q[0];
rz(-2.2994301) q[0];
sx q[0];
rz(-1.2331859) q[0];
rz(-pi) q[1];
rz(-3.1167745) q[2];
sx q[2];
rz(-0.78486004) q[2];
sx q[2];
rz(1.0944686) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79440642) q[1];
sx q[1];
rz(-2.3270895) q[1];
sx q[1];
rz(1.0760572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6448375) q[3];
sx q[3];
rz(-2.7803034) q[3];
sx q[3];
rz(2.3827117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5004702) q[2];
sx q[2];
rz(-1.0176071) q[2];
sx q[2];
rz(-2.9024331) q[2];
rz(-0.33254361) q[3];
sx q[3];
rz(-1.4556689) q[3];
sx q[3];
rz(-1.0913764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34514937) q[0];
sx q[0];
rz(-0.74757663) q[0];
sx q[0];
rz(-2.8969966) q[0];
rz(2.7469475) q[1];
sx q[1];
rz(-0.60401812) q[1];
sx q[1];
rz(-1.2044725) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0255677) q[0];
sx q[0];
rz(-3.126156) q[0];
sx q[0];
rz(-0.77657338) q[0];
x q[1];
rz(1.081624) q[2];
sx q[2];
rz(-1.2252639) q[2];
sx q[2];
rz(-2.6903077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.905758) q[1];
sx q[1];
rz(-1.8595962) q[1];
sx q[1];
rz(1.0314842) q[1];
x q[2];
rz(-1.0706903) q[3];
sx q[3];
rz(-1.2104038) q[3];
sx q[3];
rz(-0.57756027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.34387732) q[2];
sx q[2];
rz(-0.98261967) q[2];
sx q[2];
rz(-0.24125153) q[2];
rz(1.7734211) q[3];
sx q[3];
rz(-1.3894812) q[3];
sx q[3];
rz(-2.1696137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.5695246) q[0];
sx q[0];
rz(-1.854874) q[0];
sx q[0];
rz(1.1014112) q[0];
rz(-1.4934348) q[1];
sx q[1];
rz(-0.28683528) q[1];
sx q[1];
rz(-1.000584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1959609) q[0];
sx q[0];
rz(-1.9871662) q[0];
sx q[0];
rz(-3.065057) q[0];
rz(-pi) q[1];
rz(-1.5550291) q[2];
sx q[2];
rz(-0.55634275) q[2];
sx q[2];
rz(0.76130262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2615255) q[1];
sx q[1];
rz(-1.4698005) q[1];
sx q[1];
rz(1.5965881) q[1];
rz(-pi) q[2];
rz(-1.1492997) q[3];
sx q[3];
rz(-1.7022082) q[3];
sx q[3];
rz(-0.67883713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3558041) q[2];
sx q[2];
rz(-1.4996837) q[2];
sx q[2];
rz(-2.918952) q[2];
rz(1.4786134) q[3];
sx q[3];
rz(-2.7361349) q[3];
sx q[3];
rz(-1.8805898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3123689) q[0];
sx q[0];
rz(-2.0033422) q[0];
sx q[0];
rz(2.9421222) q[0];
rz(-1.624674) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(-0.46474251) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9896444) q[0];
sx q[0];
rz(-1.8590663) q[0];
sx q[0];
rz(3.0040859) q[0];
rz(2.88358) q[2];
sx q[2];
rz(-2.9128053) q[2];
sx q[2];
rz(2.8361671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81632751) q[1];
sx q[1];
rz(-1.8598598) q[1];
sx q[1];
rz(1.2802034) q[1];
x q[2];
rz(1.5681276) q[3];
sx q[3];
rz(-0.87335372) q[3];
sx q[3];
rz(-1.3476882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.042784721) q[2];
sx q[2];
rz(-1.3036737) q[2];
sx q[2];
rz(0.25351563) q[2];
rz(2.2809196) q[3];
sx q[3];
rz(-0.52152514) q[3];
sx q[3];
rz(2.0025939) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0694224) q[0];
sx q[0];
rz(-1.5266029) q[0];
sx q[0];
rz(-1.3602863) q[0];
rz(-2.5379429) q[1];
sx q[1];
rz(-2.344548) q[1];
sx q[1];
rz(0.48702494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2880858) q[0];
sx q[0];
rz(-0.65847337) q[0];
sx q[0];
rz(-0.11866026) q[0];
rz(-pi) q[1];
rz(-0.73871805) q[2];
sx q[2];
rz(-1.6411348) q[2];
sx q[2];
rz(2.3073151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4364061) q[1];
sx q[1];
rz(-1.4049705) q[1];
sx q[1];
rz(-1.2064183) q[1];
rz(2.911024) q[3];
sx q[3];
rz(-1.1784703) q[3];
sx q[3];
rz(0.84298625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46502486) q[2];
sx q[2];
rz(-1.0176696) q[2];
sx q[2];
rz(0.49099311) q[2];
rz(-0.53389126) q[3];
sx q[3];
rz(-0.54755727) q[3];
sx q[3];
rz(0.045031358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.295306) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(-1.2714161) q[0];
rz(-2.2170587) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(2.1163993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7775531) q[0];
sx q[0];
rz(-1.784783) q[0];
sx q[0];
rz(-2.5018033) q[0];
rz(-pi) q[1];
rz(0.94550262) q[2];
sx q[2];
rz(-1.3943814) q[2];
sx q[2];
rz(-2.8636065) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2424836) q[1];
sx q[1];
rz(-2.5216853) q[1];
sx q[1];
rz(-1.0044881) q[1];
x q[2];
rz(1.3253493) q[3];
sx q[3];
rz(-0.87900987) q[3];
sx q[3];
rz(-2.679685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9292494) q[2];
sx q[2];
rz(-0.29444567) q[2];
sx q[2];
rz(-0.92246169) q[2];
rz(-0.38659066) q[3];
sx q[3];
rz(-1.2873193) q[3];
sx q[3];
rz(-1.6283901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59231049) q[0];
sx q[0];
rz(-3.1070502) q[0];
sx q[0];
rz(0.70796815) q[0];
rz(2.6225846) q[1];
sx q[1];
rz(-2.612807) q[1];
sx q[1];
rz(0.66600287) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36668829) q[0];
sx q[0];
rz(-1.6631334) q[0];
sx q[0];
rz(-1.0342741) q[0];
x q[1];
rz(1.9089799) q[2];
sx q[2];
rz(-2.8034503) q[2];
sx q[2];
rz(-1.7902439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4909926) q[1];
sx q[1];
rz(-1.5751936) q[1];
sx q[1];
rz(2.472509) q[1];
rz(-2.186108) q[3];
sx q[3];
rz(-1.59366) q[3];
sx q[3];
rz(2.5290306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4697326) q[2];
sx q[2];
rz(-1.1598776) q[2];
sx q[2];
rz(-1.0603909) q[2];
rz(-1.4276069) q[3];
sx q[3];
rz(-1.5644045) q[3];
sx q[3];
rz(0.76712999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0107182) q[0];
sx q[0];
rz(-2.3227782) q[0];
sx q[0];
rz(-0.70097104) q[0];
rz(2.2531807) q[1];
sx q[1];
rz(-2.7263434) q[1];
sx q[1];
rz(1.7078687) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5007223) q[0];
sx q[0];
rz(-1.6789915) q[0];
sx q[0];
rz(-0.51525028) q[0];
rz(-0.063155576) q[2];
sx q[2];
rz(-0.40663162) q[2];
sx q[2];
rz(-0.26093182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0317185) q[1];
sx q[1];
rz(-1.5496792) q[1];
sx q[1];
rz(1.5892522) q[1];
x q[2];
rz(2.4495226) q[3];
sx q[3];
rz(-1.0502397) q[3];
sx q[3];
rz(0.51560054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99011699) q[2];
sx q[2];
rz(-1.9141804) q[2];
sx q[2];
rz(-1.6647476) q[2];
rz(-0.19674033) q[3];
sx q[3];
rz(-1.1367831) q[3];
sx q[3];
rz(-2.2018532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5340586) q[0];
sx q[0];
rz(-1.5628096) q[0];
sx q[0];
rz(-1.1210572) q[0];
rz(0.47814449) q[1];
sx q[1];
rz(-2.4604764) q[1];
sx q[1];
rz(-2.4268683) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9909458) q[0];
sx q[0];
rz(-2.0523221) q[0];
sx q[0];
rz(1.343438) q[0];
rz(-pi) q[1];
rz(2.7585538) q[2];
sx q[2];
rz(-2.0872119) q[2];
sx q[2];
rz(-0.4412152) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8712232) q[1];
sx q[1];
rz(-2.0107234) q[1];
sx q[1];
rz(1.0031071) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.274005) q[3];
sx q[3];
rz(-1.0699268) q[3];
sx q[3];
rz(1.3102368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5152682) q[2];
sx q[2];
rz(-1.3268665) q[2];
sx q[2];
rz(1.2088306) q[2];
rz(1.6995466) q[3];
sx q[3];
rz(-2.2550826) q[3];
sx q[3];
rz(0.065571688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6616853) q[0];
sx q[0];
rz(-1.3085145) q[0];
sx q[0];
rz(3.0611962) q[0];
rz(-0.54795625) q[1];
sx q[1];
rz(-0.68956551) q[1];
sx q[1];
rz(-1.7908304) q[1];
rz(1.1016079) q[2];
sx q[2];
rz(-1.3528578) q[2];
sx q[2];
rz(2.8982671) q[2];
rz(-2.6790621) q[3];
sx q[3];
rz(-1.1855385) q[3];
sx q[3];
rz(-0.0029467646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
