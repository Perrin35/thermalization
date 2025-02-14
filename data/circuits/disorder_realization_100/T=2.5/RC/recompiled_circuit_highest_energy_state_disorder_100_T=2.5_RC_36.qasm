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
rz(-0.21021357) q[0];
sx q[0];
rz(-2.7855594) q[0];
sx q[0];
rz(1.5176679) q[0];
rz(-1.3104982) q[1];
sx q[1];
rz(-0.85135353) q[1];
sx q[1];
rz(0.89827615) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1318645) q[0];
sx q[0];
rz(-2.4180331) q[0];
sx q[0];
rz(-2.1267164) q[0];
rz(-2.078371) q[2];
sx q[2];
rz(-2.5899593) q[2];
sx q[2];
rz(-1.4741613) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15816244) q[1];
sx q[1];
rz(-1.753141) q[1];
sx q[1];
rz(-1.6514062) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0086011767) q[3];
sx q[3];
rz(-1.4066753) q[3];
sx q[3];
rz(-0.65879956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6015168) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(-0.69851056) q[2];
rz(1.6554333) q[3];
sx q[3];
rz(-0.58659068) q[3];
sx q[3];
rz(-1.1521888) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0450714) q[0];
sx q[0];
rz(-2.4424545) q[0];
sx q[0];
rz(-2.844098) q[0];
rz(2.7104764) q[1];
sx q[1];
rz(-0.86687207) q[1];
sx q[1];
rz(1.3131622) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24475141) q[0];
sx q[0];
rz(-2.2479821) q[0];
sx q[0];
rz(2.0170101) q[0];
x q[1];
rz(-1.1651016) q[2];
sx q[2];
rz(-1.5303474) q[2];
sx q[2];
rz(1.8581529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6152108) q[1];
sx q[1];
rz(-1.5783678) q[1];
sx q[1];
rz(-2.3442299) q[1];
rz(-2.8937469) q[3];
sx q[3];
rz(-1.447926) q[3];
sx q[3];
rz(-0.34832277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3678652) q[2];
sx q[2];
rz(-2.6999058) q[2];
sx q[2];
rz(-2.3750677) q[2];
rz(-2.4567228) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(-0.49055704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8399452) q[0];
sx q[0];
rz(-1.8122346) q[0];
sx q[0];
rz(0.30461052) q[0];
rz(1.0003164) q[1];
sx q[1];
rz(-2.0392923) q[1];
sx q[1];
rz(0.89835483) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70750102) q[0];
sx q[0];
rz(-1.487949) q[0];
sx q[0];
rz(-2.104724) q[0];
x q[1];
rz(-2.8091891) q[2];
sx q[2];
rz(-0.89611182) q[2];
sx q[2];
rz(-1.589366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8147874) q[1];
sx q[1];
rz(-1.0070395) q[1];
sx q[1];
rz(-0.98122259) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.091354185) q[3];
sx q[3];
rz(-2.0915789) q[3];
sx q[3];
rz(2.6798673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0072713) q[2];
sx q[2];
rz(-1.4192702) q[2];
sx q[2];
rz(-2.814494) q[2];
rz(1.6359811) q[3];
sx q[3];
rz(-1.4666731) q[3];
sx q[3];
rz(1.8065642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.1148249) q[0];
sx q[0];
rz(-0.46634316) q[0];
sx q[0];
rz(-0.53792167) q[0];
rz(2.3665358) q[1];
sx q[1];
rz(-1.8102976) q[1];
sx q[1];
rz(-2.7937826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3678326) q[0];
sx q[0];
rz(-1.1080386) q[0];
sx q[0];
rz(2.5823043) q[0];
x q[1];
rz(0.78601894) q[2];
sx q[2];
rz(-0.14855222) q[2];
sx q[2];
rz(2.1296322) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2217933) q[1];
sx q[1];
rz(-1.9927653) q[1];
sx q[1];
rz(2.3942562) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1600445) q[3];
sx q[3];
rz(-2.7166945) q[3];
sx q[3];
rz(1.9187669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20370087) q[2];
sx q[2];
rz(-2.4163279) q[2];
sx q[2];
rz(2.4791278) q[2];
rz(1.0558111) q[3];
sx q[3];
rz(-2.8231088) q[3];
sx q[3];
rz(1.8752347) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1556959) q[0];
sx q[0];
rz(-2.6031384) q[0];
sx q[0];
rz(-2.1180617) q[0];
rz(-2.0899978) q[1];
sx q[1];
rz(-1.4902427) q[1];
sx q[1];
rz(0.016062707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9883606) q[0];
sx q[0];
rz(-1.4356336) q[0];
sx q[0];
rz(-2.1188117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93463411) q[2];
sx q[2];
rz(-2.155288) q[2];
sx q[2];
rz(2.5713657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93960436) q[1];
sx q[1];
rz(-1.3496205) q[1];
sx q[1];
rz(-0.53316922) q[1];
rz(-pi) q[2];
rz(-1.8422801) q[3];
sx q[3];
rz(-1.6405655) q[3];
sx q[3];
rz(-1.8673563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2992799) q[2];
sx q[2];
rz(-0.88025847) q[2];
sx q[2];
rz(0.24097815) q[2];
rz(1.109451) q[3];
sx q[3];
rz(-1.2883319) q[3];
sx q[3];
rz(-0.39458767) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0146765) q[0];
sx q[0];
rz(-2.3052445) q[0];
sx q[0];
rz(0.29689223) q[0];
rz(-1.9706005) q[1];
sx q[1];
rz(-1.3208656) q[1];
sx q[1];
rz(-0.5353294) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2123665) q[0];
sx q[0];
rz(-1.7571714) q[0];
sx q[0];
rz(-1.3829689) q[0];
rz(-pi) q[1];
rz(2.3248843) q[2];
sx q[2];
rz(-1.7871408) q[2];
sx q[2];
rz(0.78934455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3228906) q[1];
sx q[1];
rz(-1.358958) q[1];
sx q[1];
rz(-2.0634275) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13958474) q[3];
sx q[3];
rz(-1.4142766) q[3];
sx q[3];
rz(2.8536882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8239596) q[2];
sx q[2];
rz(-0.77482688) q[2];
sx q[2];
rz(-0.63014692) q[2];
rz(-1.0050425) q[3];
sx q[3];
rz(-0.21868394) q[3];
sx q[3];
rz(0.69766831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9610577) q[0];
sx q[0];
rz(-0.068004161) q[0];
sx q[0];
rz(-1.2766174) q[0];
rz(1.0446154) q[1];
sx q[1];
rz(-1.7216543) q[1];
sx q[1];
rz(0.23342625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2058207) q[0];
sx q[0];
rz(-2.290831) q[0];
sx q[0];
rz(0.056179742) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2646003) q[2];
sx q[2];
rz(-2.1900935) q[2];
sx q[2];
rz(-2.3690852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0835159) q[1];
sx q[1];
rz(-0.74878557) q[1];
sx q[1];
rz(-3.0850436) q[1];
rz(-1.5039316) q[3];
sx q[3];
rz(-2.1479448) q[3];
sx q[3];
rz(-1.701783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6650247) q[2];
sx q[2];
rz(-1.8835521) q[2];
sx q[2];
rz(-2.8705719) q[2];
rz(-2.792231) q[3];
sx q[3];
rz(-2.2237399) q[3];
sx q[3];
rz(2.5640986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9409598) q[0];
sx q[0];
rz(-0.9372434) q[0];
sx q[0];
rz(1.8120793) q[0];
rz(0.26893523) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(-1.3806794) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5912) q[0];
sx q[0];
rz(-1.7730646) q[0];
sx q[0];
rz(-1.8254542) q[0];
rz(0.024928851) q[2];
sx q[2];
rz(-2.5453794) q[2];
sx q[2];
rz(2.2298857) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.096752874) q[1];
sx q[1];
rz(-2.1659454) q[1];
sx q[1];
rz(2.9446141) q[1];
rz(-pi) q[2];
rz(-2.9246522) q[3];
sx q[3];
rz(-2.5703781) q[3];
sx q[3];
rz(1.2071213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55714503) q[2];
sx q[2];
rz(-0.86981213) q[2];
sx q[2];
rz(-2.4073041) q[2];
rz(-1.0555438) q[3];
sx q[3];
rz(-1.5922092) q[3];
sx q[3];
rz(-2.1671104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771117) q[0];
sx q[0];
rz(-0.25150126) q[0];
sx q[0];
rz(0.86790458) q[0];
rz(-1.3033298) q[1];
sx q[1];
rz(-1.5497327) q[1];
sx q[1];
rz(-0.23095362) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8036007) q[0];
sx q[0];
rz(-0.70728318) q[0];
sx q[0];
rz(-2.4770205) q[0];
rz(2.7074293) q[2];
sx q[2];
rz(-1.4462398) q[2];
sx q[2];
rz(1.9366154) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7119904) q[1];
sx q[1];
rz(-1.7884521) q[1];
sx q[1];
rz(-2.6764286) q[1];
x q[2];
rz(-2.4695314) q[3];
sx q[3];
rz(-1.5698182) q[3];
sx q[3];
rz(-3.1100215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93297282) q[2];
sx q[2];
rz(-1.4459556) q[2];
sx q[2];
rz(2.7033973) q[2];
rz(0.70133251) q[3];
sx q[3];
rz(-1.067679) q[3];
sx q[3];
rz(2.7856538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0375131) q[0];
sx q[0];
rz(-2.6689745) q[0];
sx q[0];
rz(-1.0802826) q[0];
rz(1.0516306) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(-1.0952605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0290463) q[0];
sx q[0];
rz(-1.284082) q[0];
sx q[0];
rz(2.6081144) q[0];
x q[1];
rz(0.11028408) q[2];
sx q[2];
rz(-0.37026893) q[2];
sx q[2];
rz(-2.8831589) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43021529) q[1];
sx q[1];
rz(-1.681455) q[1];
sx q[1];
rz(0.84215136) q[1];
rz(-pi) q[2];
rz(-1.0686103) q[3];
sx q[3];
rz(-2.1081703) q[3];
sx q[3];
rz(-2.2110232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0206535) q[2];
sx q[2];
rz(-1.8755269) q[2];
sx q[2];
rz(2.2361163) q[2];
rz(2.4261273) q[3];
sx q[3];
rz(-2.4162636) q[3];
sx q[3];
rz(0.65413094) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41659551) q[0];
sx q[0];
rz(-1.7323957) q[0];
sx q[0];
rz(0.61687627) q[0];
rz(2.0116518) q[1];
sx q[1];
rz(-1.8604953) q[1];
sx q[1];
rz(-1.4166191) q[1];
rz(-0.62192179) q[2];
sx q[2];
rz(-0.94246943) q[2];
sx q[2];
rz(3.111919) q[2];
rz(2.9898217) q[3];
sx q[3];
rz(-1.7790773) q[3];
sx q[3];
rz(1.9341999) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
