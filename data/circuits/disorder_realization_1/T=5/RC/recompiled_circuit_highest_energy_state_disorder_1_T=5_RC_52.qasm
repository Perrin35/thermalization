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
rz(-0.048592351) q[0];
sx q[0];
rz(-1.7162004) q[0];
sx q[0];
rz(-0.26560321) q[0];
rz(2.6939997) q[1];
sx q[1];
rz(-1.5040553) q[1];
sx q[1];
rz(-3.016959) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87089964) q[0];
sx q[0];
rz(-2.2147372) q[0];
sx q[0];
rz(0.55418153) q[0];
rz(1.2536178) q[2];
sx q[2];
rz(-2.7493959) q[2];
sx q[2];
rz(-2.9756096) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9198299) q[1];
sx q[1];
rz(-1.3897093) q[1];
sx q[1];
rz(-2.9403375) q[1];
x q[2];
rz(0.15480583) q[3];
sx q[3];
rz(-1.4034287) q[3];
sx q[3];
rz(0.49661206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0894185) q[2];
sx q[2];
rz(-0.82145059) q[2];
sx q[2];
rz(-0.87665147) q[2];
rz(0.18757251) q[3];
sx q[3];
rz(-2.1552174) q[3];
sx q[3];
rz(-0.11876373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12980421) q[0];
sx q[0];
rz(-0.056948245) q[0];
sx q[0];
rz(2.7563128) q[0];
rz(0.41703364) q[1];
sx q[1];
rz(-0.72154355) q[1];
sx q[1];
rz(2.1293652) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6979264) q[0];
sx q[0];
rz(-1.7886284) q[0];
sx q[0];
rz(-0.99866726) q[0];
x q[1];
rz(2.5778722) q[2];
sx q[2];
rz(-2.4032058) q[2];
sx q[2];
rz(0.2283048) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5972388) q[1];
sx q[1];
rz(-1.8569678) q[1];
sx q[1];
rz(-2.3765122) q[1];
rz(0.89093633) q[3];
sx q[3];
rz(-1.9417986) q[3];
sx q[3];
rz(1.551451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9932844) q[2];
sx q[2];
rz(-0.40996429) q[2];
sx q[2];
rz(-2.6170464) q[2];
rz(-2.2927393) q[3];
sx q[3];
rz(-0.69791228) q[3];
sx q[3];
rz(2.2822288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9808905) q[0];
sx q[0];
rz(-0.41146678) q[0];
sx q[0];
rz(-2.8149862) q[0];
rz(-1.7713361) q[1];
sx q[1];
rz(-2.0854918) q[1];
sx q[1];
rz(-3.1050217) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9243226) q[0];
sx q[0];
rz(-2.0170209) q[0];
sx q[0];
rz(-1.1749595) q[0];
rz(0.70539523) q[2];
sx q[2];
rz(-1.2327895) q[2];
sx q[2];
rz(2.2662804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3044839) q[1];
sx q[1];
rz(-1.3837985) q[1];
sx q[1];
rz(3.0977102) q[1];
rz(-pi) q[2];
rz(2.1362334) q[3];
sx q[3];
rz(-2.812139) q[3];
sx q[3];
rz(-2.0170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6959491) q[2];
sx q[2];
rz(-0.30569884) q[2];
sx q[2];
rz(-1.9574399) q[2];
rz(0.46237692) q[3];
sx q[3];
rz(-2.0591044) q[3];
sx q[3];
rz(-0.46557993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9584123) q[0];
sx q[0];
rz(-1.4875702) q[0];
sx q[0];
rz(2.2271449) q[0];
rz(-2.27521) q[1];
sx q[1];
rz(-1.2976846) q[1];
sx q[1];
rz(-0.31747174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42714092) q[0];
sx q[0];
rz(-1.1981816) q[0];
sx q[0];
rz(1.4535377) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2214041) q[2];
sx q[2];
rz(-1.411806) q[2];
sx q[2];
rz(2.8391389) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13617789) q[1];
sx q[1];
rz(-1.7247611) q[1];
sx q[1];
rz(1.8421768) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0013591493) q[3];
sx q[3];
rz(-1.1906644) q[3];
sx q[3];
rz(1.8111977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0471961) q[2];
sx q[2];
rz(-0.4717584) q[2];
sx q[2];
rz(0.45486927) q[2];
rz(-1.1826078) q[3];
sx q[3];
rz(-0.22795658) q[3];
sx q[3];
rz(-0.14127775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126548) q[0];
sx q[0];
rz(-1.2417355) q[0];
sx q[0];
rz(-0.051359635) q[0];
rz(-2.8142169) q[1];
sx q[1];
rz(-0.86762571) q[1];
sx q[1];
rz(0.059965722) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7284847) q[0];
sx q[0];
rz(-0.93046549) q[0];
sx q[0];
rz(-2.722239) q[0];
x q[1];
rz(-0.65360258) q[2];
sx q[2];
rz(-1.8817668) q[2];
sx q[2];
rz(2.5429436) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9717945) q[1];
sx q[1];
rz(-2.4996539) q[1];
sx q[1];
rz(1.0331912) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7945847) q[3];
sx q[3];
rz(-0.7273376) q[3];
sx q[3];
rz(2.2801496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7472234) q[2];
sx q[2];
rz(-0.9129492) q[2];
sx q[2];
rz(-0.24820122) q[2];
rz(0.40211755) q[3];
sx q[3];
rz(-2.6995903) q[3];
sx q[3];
rz(2.2615652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1345271) q[0];
sx q[0];
rz(-1.2638673) q[0];
sx q[0];
rz(-1.7279351) q[0];
rz(1.9386579) q[1];
sx q[1];
rz(-2.2212432) q[1];
sx q[1];
rz(0.10341067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2884379) q[0];
sx q[0];
rz(-1.5010035) q[0];
sx q[0];
rz(1.5091108) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8273583) q[2];
sx q[2];
rz(-1.7467919) q[2];
sx q[2];
rz(-2.2361148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6573665) q[1];
sx q[1];
rz(-1.5897337) q[1];
sx q[1];
rz(0.12469805) q[1];
rz(0.2528905) q[3];
sx q[3];
rz(-1.1973945) q[3];
sx q[3];
rz(1.4772082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47205135) q[2];
sx q[2];
rz(-2.5504888) q[2];
sx q[2];
rz(-2.9284787) q[2];
rz(1.6040246) q[3];
sx q[3];
rz(-3.0209916) q[3];
sx q[3];
rz(3.0311301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14385001) q[0];
sx q[0];
rz(-2.282833) q[0];
sx q[0];
rz(0.48458883) q[0];
rz(2.6735725) q[1];
sx q[1];
rz(-1.8193529) q[1];
sx q[1];
rz(-3.0048634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28324652) q[0];
sx q[0];
rz(-1.8112881) q[0];
sx q[0];
rz(1.9577615) q[0];
x q[1];
rz(-0.38928208) q[2];
sx q[2];
rz(-1.2503257) q[2];
sx q[2];
rz(2.5581202) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.90367442) q[1];
sx q[1];
rz(-0.72924239) q[1];
sx q[1];
rz(-1.7715471) q[1];
x q[2];
rz(-1.6893295) q[3];
sx q[3];
rz(-1.9876936) q[3];
sx q[3];
rz(-2.0385008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4097269) q[2];
sx q[2];
rz(-3.0146283) q[2];
sx q[2];
rz(-0.48180386) q[2];
rz(1.2305772) q[3];
sx q[3];
rz(-1.9989719) q[3];
sx q[3];
rz(-3.0060911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19732533) q[0];
sx q[0];
rz(-2.8062286) q[0];
sx q[0];
rz(2.770597) q[0];
rz(2.9839997) q[1];
sx q[1];
rz(-1.7540365) q[1];
sx q[1];
rz(2.2612459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77773422) q[0];
sx q[0];
rz(-2.661718) q[0];
sx q[0];
rz(1.721948) q[0];
rz(-pi) q[1];
rz(0.7083986) q[2];
sx q[2];
rz(-1.4428291) q[2];
sx q[2];
rz(-2.5789946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0571361) q[1];
sx q[1];
rz(-1.5669654) q[1];
sx q[1];
rz(-1.6935029) q[1];
rz(-pi) q[2];
rz(2.8572389) q[3];
sx q[3];
rz(-1.8508678) q[3];
sx q[3];
rz(-1.3765311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8627491) q[2];
sx q[2];
rz(-0.66484386) q[2];
sx q[2];
rz(-2.1319907) q[2];
rz(-2.3889715) q[3];
sx q[3];
rz(-1.8372583) q[3];
sx q[3];
rz(2.2032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272334) q[0];
sx q[0];
rz(-0.14150134) q[0];
sx q[0];
rz(-0.046382647) q[0];
rz(-1.4917689) q[1];
sx q[1];
rz(-1.7581419) q[1];
sx q[1];
rz(-0.47328624) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21264874) q[0];
sx q[0];
rz(-1.4601652) q[0];
sx q[0];
rz(-1.7615737) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0436486) q[2];
sx q[2];
rz(-1.2769498) q[2];
sx q[2];
rz(-2.1515255) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9732194) q[1];
sx q[1];
rz(-2.1811317) q[1];
sx q[1];
rz(-1.8944505) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5389082) q[3];
sx q[3];
rz(-1.569935) q[3];
sx q[3];
rz(-2.3501648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9778694) q[2];
sx q[2];
rz(-2.551584) q[2];
sx q[2];
rz(3.1150277) q[2];
rz(0.50655347) q[3];
sx q[3];
rz(-0.81304628) q[3];
sx q[3];
rz(2.626239) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16515054) q[0];
sx q[0];
rz(-2.1858175) q[0];
sx q[0];
rz(2.9859848) q[0];
rz(0.43442976) q[1];
sx q[1];
rz(-2.5948718) q[1];
sx q[1];
rz(0.2880407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018943448) q[0];
sx q[0];
rz(-2.3068006) q[0];
sx q[0];
rz(2.3998866) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6596117) q[2];
sx q[2];
rz(-0.15379158) q[2];
sx q[2];
rz(1.3336934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1120243) q[1];
sx q[1];
rz(-2.3178737) q[1];
sx q[1];
rz(-2.1805641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0218009) q[3];
sx q[3];
rz(-0.70216252) q[3];
sx q[3];
rz(-2.8150866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7274999) q[2];
sx q[2];
rz(-0.45656559) q[2];
sx q[2];
rz(-3.0957733) q[2];
rz(-3.0231061) q[3];
sx q[3];
rz(-2.249497) q[3];
sx q[3];
rz(2.4011325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6427778) q[0];
sx q[0];
rz(-1.5121664) q[0];
sx q[0];
rz(1.2335516) q[0];
rz(-1.9998101) q[1];
sx q[1];
rz(-0.70527609) q[1];
sx q[1];
rz(-1.3377778) q[1];
rz(-2.4579688) q[2];
sx q[2];
rz(-2.0487006) q[2];
sx q[2];
rz(-0.13681199) q[2];
rz(1.508504) q[3];
sx q[3];
rz(-1.5083762) q[3];
sx q[3];
rz(1.4948869) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
