OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5034135) q[0];
sx q[0];
rz(-1.2794275) q[0];
sx q[0];
rz(0.77603618) q[0];
rz(0.70932055) q[1];
sx q[1];
rz(4.6588916) q[1];
sx q[1];
rz(8.8257596) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4871019) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(2.046665) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9089977) q[2];
sx q[2];
rz(-2.4233492) q[2];
sx q[2];
rz(-0.79923779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30556074) q[1];
sx q[1];
rz(-1.2183883) q[1];
sx q[1];
rz(-0.52959092) q[1];
rz(-pi) q[2];
rz(-2.7077449) q[3];
sx q[3];
rz(-2.3126855) q[3];
sx q[3];
rz(2.4494107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0029995) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(-0.19621672) q[2];
rz(-2.0972706) q[3];
sx q[3];
rz(-0.82572562) q[3];
sx q[3];
rz(2.830982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1021378) q[0];
sx q[0];
rz(-2.4882443) q[0];
sx q[0];
rz(-0.85773221) q[0];
rz(2.6013382) q[1];
sx q[1];
rz(-2.0876355) q[1];
sx q[1];
rz(2.3589755) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7359813) q[0];
sx q[0];
rz(-2.1850039) q[0];
sx q[0];
rz(-0.10615291) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8939607) q[2];
sx q[2];
rz(-0.64715451) q[2];
sx q[2];
rz(-1.1101071) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2007717) q[1];
sx q[1];
rz(-1.7634974) q[1];
sx q[1];
rz(-2.5351943) q[1];
x q[2];
rz(-2.2563062) q[3];
sx q[3];
rz(-1.1210223) q[3];
sx q[3];
rz(-1.3441844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9767849) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(0.59107333) q[2];
rz(1.0278541) q[3];
sx q[3];
rz(-2.5058392) q[3];
sx q[3];
rz(0.13636057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117821) q[0];
sx q[0];
rz(-0.52017838) q[0];
sx q[0];
rz(-2.0853364) q[0];
rz(-2.1504869) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(2.3371005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4933518) q[0];
sx q[0];
rz(-1.3051148) q[0];
sx q[0];
rz(2.4881287) q[0];
x q[1];
rz(1.7011004) q[2];
sx q[2];
rz(-2.5087803) q[2];
sx q[2];
rz(1.9595343) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.04136297) q[1];
sx q[1];
rz(-1.8112665) q[1];
sx q[1];
rz(0.23794707) q[1];
x q[2];
rz(-0.64739703) q[3];
sx q[3];
rz(-1.4167827) q[3];
sx q[3];
rz(-1.9141509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.450401) q[2];
sx q[2];
rz(-2.2914026) q[2];
sx q[2];
rz(-0.56785339) q[2];
rz(1.9495226) q[3];
sx q[3];
rz(-0.53693938) q[3];
sx q[3];
rz(0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464226) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(-1.9544741) q[0];
rz(-2.8861956) q[1];
sx q[1];
rz(-1.4030158) q[1];
sx q[1];
rz(2.752221) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0856649) q[0];
sx q[0];
rz(-1.4367668) q[0];
sx q[0];
rz(2.9197951) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2979173) q[2];
sx q[2];
rz(-1.5483583) q[2];
sx q[2];
rz(-2.0539935) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7070605) q[1];
sx q[1];
rz(-1.3816621) q[1];
sx q[1];
rz(2.2347336) q[1];
rz(-pi) q[2];
rz(-1.7674753) q[3];
sx q[3];
rz(-1.4405319) q[3];
sx q[3];
rz(1.2127339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7369467) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(1.3673937) q[2];
rz(-0.5395475) q[3];
sx q[3];
rz(-1.5522141) q[3];
sx q[3];
rz(1.3246983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0747727) q[0];
sx q[0];
rz(-2.9434151) q[0];
sx q[0];
rz(0.5603801) q[0];
rz(-2.9226774) q[1];
sx q[1];
rz(-2.0603265) q[1];
sx q[1];
rz(-0.17094368) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3363638) q[0];
sx q[0];
rz(-0.85652387) q[0];
sx q[0];
rz(-0.63905604) q[0];
rz(0.60336171) q[2];
sx q[2];
rz(-2.8953279) q[2];
sx q[2];
rz(-2.6880996) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51851398) q[1];
sx q[1];
rz(-2.2840223) q[1];
sx q[1];
rz(0.70116373) q[1];
x q[2];
rz(1.6690977) q[3];
sx q[3];
rz(-0.62255961) q[3];
sx q[3];
rz(-0.05427256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35974744) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(2.183059) q[2];
rz(2.9790699) q[3];
sx q[3];
rz(-2.2992117) q[3];
sx q[3];
rz(-1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.573134) q[0];
sx q[0];
rz(-0.064366654) q[0];
sx q[0];
rz(0.74813133) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-2.7225814) q[1];
sx q[1];
rz(0.31707877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9592764) q[0];
sx q[0];
rz(-1.8166564) q[0];
sx q[0];
rz(1.5005174) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19642475) q[2];
sx q[2];
rz(-1.6187917) q[2];
sx q[2];
rz(2.6065207) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.78817716) q[1];
sx q[1];
rz(-1.4995575) q[1];
sx q[1];
rz(-2.337238) q[1];
rz(-0.040359453) q[3];
sx q[3];
rz(-2.3642614) q[3];
sx q[3];
rz(-3.0632927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3743484) q[2];
sx q[2];
rz(-0.92864645) q[2];
sx q[2];
rz(1.413215) q[2];
rz(-3.0610541) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(-0.22182375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0291075) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(-1.8674194) q[0];
rz(1.796465) q[1];
sx q[1];
rz(-0.74256623) q[1];
sx q[1];
rz(1.3412195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8250834) q[0];
sx q[0];
rz(-2.1232455) q[0];
sx q[0];
rz(-2.1567206) q[0];
rz(2.6244782) q[2];
sx q[2];
rz(-0.54143751) q[2];
sx q[2];
rz(-2.0064029) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6051424) q[1];
sx q[1];
rz(-1.7243694) q[1];
sx q[1];
rz(-1.3773514) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7246095) q[3];
sx q[3];
rz(-2.4857368) q[3];
sx q[3];
rz(0.7747618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93366569) q[2];
sx q[2];
rz(-1.9542481) q[2];
sx q[2];
rz(0.46736091) q[2];
rz(0.76006877) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(-1.8925586) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6811328) q[0];
sx q[0];
rz(-1.0204027) q[0];
sx q[0];
rz(1.4469752) q[0];
rz(0.32577062) q[1];
sx q[1];
rz(-2.9278946) q[1];
sx q[1];
rz(-1.5209341) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64031271) q[0];
sx q[0];
rz(-0.58459832) q[0];
sx q[0];
rz(-2.6452097) q[0];
rz(-pi) q[1];
x q[1];
rz(0.09163945) q[2];
sx q[2];
rz(-1.5704201) q[2];
sx q[2];
rz(0.74972744) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22252488) q[1];
sx q[1];
rz(-1.7456858) q[1];
sx q[1];
rz(-1.9293849) q[1];
rz(2.2748442) q[3];
sx q[3];
rz(-1.4102077) q[3];
sx q[3];
rz(-2.989547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96523607) q[2];
sx q[2];
rz(-2.3945645) q[2];
sx q[2];
rz(2.526324) q[2];
rz(1.8687013) q[3];
sx q[3];
rz(-2.8764909) q[3];
sx q[3];
rz(-0.42241514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0203005) q[0];
sx q[0];
rz(-1.8676119) q[0];
sx q[0];
rz(2.2913388) q[0];
rz(2.0203159) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(0.77176315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4798652) q[0];
sx q[0];
rz(-0.87718946) q[0];
sx q[0];
rz(3.1105009) q[0];
x q[1];
rz(0.36649165) q[2];
sx q[2];
rz(-0.19333177) q[2];
sx q[2];
rz(-1.1761348) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38852019) q[1];
sx q[1];
rz(-2.8917312) q[1];
sx q[1];
rz(1.3461963) q[1];
x q[2];
rz(-0.27098591) q[3];
sx q[3];
rz(-1.2577204) q[3];
sx q[3];
rz(-0.82161689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11543342) q[2];
sx q[2];
rz(-1.5966281) q[2];
sx q[2];
rz(0.8405295) q[2];
rz(3.0606411) q[3];
sx q[3];
rz(-0.5778802) q[3];
sx q[3];
rz(0.24278434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42846546) q[0];
sx q[0];
rz(-0.19793887) q[0];
sx q[0];
rz(1.9848829) q[0];
rz(-1.0276065) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(2.3311232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5549703) q[0];
sx q[0];
rz(-1.7867309) q[0];
sx q[0];
rz(2.3761889) q[0];
x q[1];
rz(-0.54525156) q[2];
sx q[2];
rz(-1.8601189) q[2];
sx q[2];
rz(2.3985753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6582074) q[1];
sx q[1];
rz(-1.2772296) q[1];
sx q[1];
rz(0.29694966) q[1];
rz(-2.9083071) q[3];
sx q[3];
rz(-2.5676227) q[3];
sx q[3];
rz(0.026653224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.345574) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(1.5239117) q[2];
rz(1.1003305) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(-0.17980096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1179467) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(-3.0991411) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(1.4412075) q[2];
sx q[2];
rz(-0.40934206) q[2];
sx q[2];
rz(-0.48624292) q[2];
rz(-0.073033606) q[3];
sx q[3];
rz(-0.45847736) q[3];
sx q[3];
rz(1.7401742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
