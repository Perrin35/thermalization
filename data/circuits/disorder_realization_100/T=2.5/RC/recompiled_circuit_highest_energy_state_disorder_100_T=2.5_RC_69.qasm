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
rz(-1.6239248) q[0];
rz(-1.3104982) q[1];
sx q[1];
rz(-0.85135353) q[1];
sx q[1];
rz(-2.2433165) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7018902) q[0];
sx q[0];
rz(-0.97355803) q[0];
sx q[0];
rz(-0.43623121) q[0];
rz(1.0773877) q[2];
sx q[2];
rz(-1.8283683) q[2];
sx q[2];
rz(2.7958564) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5772669) q[1];
sx q[1];
rz(-0.1991867) q[1];
sx q[1];
rz(-0.41173068) q[1];
rz(-pi) q[2];
rz(1.6226852) q[3];
sx q[3];
rz(-2.9772485) q[3];
sx q[3];
rz(0.60620327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6015168) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(2.4430821) q[2];
rz(-1.6554333) q[3];
sx q[3];
rz(-0.58659068) q[3];
sx q[3];
rz(1.1521888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0965213) q[0];
sx q[0];
rz(-0.69913816) q[0];
sx q[0];
rz(2.844098) q[0];
rz(-0.43111626) q[1];
sx q[1];
rz(-0.86687207) q[1];
sx q[1];
rz(-1.8284304) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8968412) q[0];
sx q[0];
rz(-2.2479821) q[0];
sx q[0];
rz(1.1245825) q[0];
x q[1];
rz(-1.4686062) q[2];
sx q[2];
rz(-2.7339978) q[2];
sx q[2];
rz(-0.38122618) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6152108) q[1];
sx q[1];
rz(-1.5783678) q[1];
sx q[1];
rz(-0.79736276) q[1];
x q[2];
rz(2.8937469) q[3];
sx q[3];
rz(-1.447926) q[3];
sx q[3];
rz(0.34832277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3678652) q[2];
sx q[2];
rz(-0.44168681) q[2];
sx q[2];
rz(0.76652491) q[2];
rz(2.4567228) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(-2.6510356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8399452) q[0];
sx q[0];
rz(-1.3293581) q[0];
sx q[0];
rz(-0.30461052) q[0];
rz(-2.1412762) q[1];
sx q[1];
rz(-2.0392923) q[1];
sx q[1];
rz(-2.2432378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72422623) q[0];
sx q[0];
rz(-0.53970102) q[0];
sx q[0];
rz(1.7325355) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2730973) q[2];
sx q[2];
rz(-1.3131333) q[2];
sx q[2];
rz(0.23095498) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8147874) q[1];
sx q[1];
rz(-2.1345532) q[1];
sx q[1];
rz(2.1603701) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0482084) q[3];
sx q[3];
rz(-1.6500123) q[3];
sx q[3];
rz(-2.0780711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0072713) q[2];
sx q[2];
rz(-1.7223225) q[2];
sx q[2];
rz(-2.814494) q[2];
rz(-1.5056115) q[3];
sx q[3];
rz(-1.4666731) q[3];
sx q[3];
rz(-1.3350284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026767749) q[0];
sx q[0];
rz(-0.46634316) q[0];
sx q[0];
rz(-2.603671) q[0];
rz(0.7750569) q[1];
sx q[1];
rz(-1.8102976) q[1];
sx q[1];
rz(2.7937826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93047543) q[0];
sx q[0];
rz(-2.0655091) q[0];
sx q[0];
rz(1.0388264) q[0];
rz(-pi) q[1];
rz(-2.3555737) q[2];
sx q[2];
rz(-0.14855222) q[2];
sx q[2];
rz(-1.0119604) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0671652) q[1];
sx q[1];
rz(-0.8377155) q[1];
sx q[1];
rz(-2.5578686) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1600445) q[3];
sx q[3];
rz(-0.42489811) q[3];
sx q[3];
rz(1.2228257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20370087) q[2];
sx q[2];
rz(-2.4163279) q[2];
sx q[2];
rz(-2.4791278) q[2];
rz(2.0857816) q[3];
sx q[3];
rz(-2.8231088) q[3];
sx q[3];
rz(-1.8752347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1556959) q[0];
sx q[0];
rz(-2.6031384) q[0];
sx q[0];
rz(-1.023531) q[0];
rz(1.0515949) q[1];
sx q[1];
rz(-1.6513499) q[1];
sx q[1];
rz(-0.016062707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6419638) q[0];
sx q[0];
rz(-1.0283386) q[0];
sx q[0];
rz(-2.9835975) q[0];
rz(-pi) q[1];
rz(-2.4532796) q[2];
sx q[2];
rz(-1.0524024) q[2];
sx q[2];
rz(-0.61358085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50243584) q[1];
sx q[1];
rz(-1.0519439) q[1];
sx q[1];
rz(-1.826189) q[1];
rz(-pi) q[2];
rz(-3.06918) q[3];
sx q[3];
rz(-1.2999897) q[3];
sx q[3];
rz(-2.8256302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84231275) q[2];
sx q[2];
rz(-2.2613342) q[2];
sx q[2];
rz(2.9006145) q[2];
rz(2.0321417) q[3];
sx q[3];
rz(-1.2883319) q[3];
sx q[3];
rz(0.39458767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-1.0146765) q[0];
sx q[0];
rz(-0.83634818) q[0];
sx q[0];
rz(-0.29689223) q[0];
rz(-1.1709921) q[1];
sx q[1];
rz(-1.3208656) q[1];
sx q[1];
rz(0.5353294) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7479582) q[0];
sx q[0];
rz(-1.3862594) q[0];
sx q[0];
rz(-2.9519597) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3248843) q[2];
sx q[2];
rz(-1.3544519) q[2];
sx q[2];
rz(-0.78934455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0161617) q[1];
sx q[1];
rz(-0.53278538) q[1];
sx q[1];
rz(1.997581) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4127647) q[3];
sx q[3];
rz(-1.7086638) q[3];
sx q[3];
rz(1.3047895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.31763306) q[2];
sx q[2];
rz(-2.3667658) q[2];
sx q[2];
rz(2.5114457) q[2];
rz(2.1365502) q[3];
sx q[3];
rz(-0.21868394) q[3];
sx q[3];
rz(-2.4439243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9610577) q[0];
sx q[0];
rz(-0.068004161) q[0];
sx q[0];
rz(1.2766174) q[0];
rz(-2.0969773) q[1];
sx q[1];
rz(-1.7216543) q[1];
sx q[1];
rz(0.23342625) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1207399) q[0];
sx q[0];
rz(-0.7218315) q[0];
sx q[0];
rz(-1.634725) q[0];
x q[1];
rz(-0.64197783) q[2];
sx q[2];
rz(-1.8187838) q[2];
sx q[2];
rz(2.5247912) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0580768) q[1];
sx q[1];
rz(-2.3928071) q[1];
sx q[1];
rz(3.0850436) q[1];
rz(-pi) q[2];
rz(1.6376611) q[3];
sx q[3];
rz(-2.1479448) q[3];
sx q[3];
rz(-1.701783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.476568) q[2];
sx q[2];
rz(-1.2580405) q[2];
sx q[2];
rz(-0.27102077) q[2];
rz(-2.792231) q[3];
sx q[3];
rz(-2.2237399) q[3];
sx q[3];
rz(-0.57749403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.2006328) q[0];
sx q[0];
rz(-2.2043493) q[0];
sx q[0];
rz(-1.8120793) q[0];
rz(-2.8726574) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(1.7609133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5503926) q[0];
sx q[0];
rz(-1.7730646) q[0];
sx q[0];
rz(1.8254542) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59606885) q[2];
sx q[2];
rz(-1.5847932) q[2];
sx q[2];
rz(0.67971855) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5854634) q[1];
sx q[1];
rz(-1.4080202) q[1];
sx q[1];
rz(2.1750431) q[1];
x q[2];
rz(-0.56048067) q[3];
sx q[3];
rz(-1.4541601) q[3];
sx q[3];
rz(-0.18033914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5844476) q[2];
sx q[2];
rz(-0.86981213) q[2];
sx q[2];
rz(2.4073041) q[2];
rz(2.0860489) q[3];
sx q[3];
rz(-1.5493834) q[3];
sx q[3];
rz(-0.9744823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.445914) q[0];
sx q[0];
rz(-1.983108) q[0];
sx q[0];
rz(-0.59230174) q[0];
rz(-1.7079389) q[2];
sx q[2];
rz(-2.0013705) q[2];
sx q[2];
rz(0.42335934) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24915299) q[1];
sx q[1];
rz(-1.1174392) q[1];
sx q[1];
rz(-1.3282177) q[1];
rz(0.6720613) q[3];
sx q[3];
rz(-1.5698182) q[3];
sx q[3];
rz(0.031571183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93297282) q[2];
sx q[2];
rz(-1.695637) q[2];
sx q[2];
rz(-0.43819532) q[2];
rz(0.70133251) q[3];
sx q[3];
rz(-2.0739136) q[3];
sx q[3];
rz(0.35593885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0375131) q[0];
sx q[0];
rz(-2.6689745) q[0];
sx q[0];
rz(2.0613101) q[0];
rz(-2.089962) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(2.0463321) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0290463) q[0];
sx q[0];
rz(-1.8575107) q[0];
sx q[0];
rz(-2.6081144) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5280998) q[2];
sx q[2];
rz(-1.9387081) q[2];
sx q[2];
rz(-0.14020444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.124123) q[1];
sx q[1];
rz(-2.4061205) q[1];
sx q[1];
rz(1.7361438) q[1];
rz(-2.5445537) q[3];
sx q[3];
rz(-1.9970915) q[3];
sx q[3];
rz(-2.7753914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1209391) q[2];
sx q[2];
rz(-1.2660657) q[2];
sx q[2];
rz(0.90547639) q[2];
rz(0.71546537) q[3];
sx q[3];
rz(-2.4162636) q[3];
sx q[3];
rz(2.4874617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.62192179) q[2];
sx q[2];
rz(-2.1991232) q[2];
sx q[2];
rz(-0.029673619) q[2];
rz(-1.7814287) q[3];
sx q[3];
rz(-1.4223301) q[3];
sx q[3];
rz(-2.8098047) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
