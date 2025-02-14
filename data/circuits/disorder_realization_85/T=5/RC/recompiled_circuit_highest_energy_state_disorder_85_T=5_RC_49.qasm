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
rz(-0.83952251) q[0];
sx q[0];
rz(-1.5877725) q[0];
sx q[0];
rz(-2.177218) q[0];
rz(-2.7388465) q[1];
sx q[1];
rz(-1.4637113) q[1];
sx q[1];
rz(0.95199624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453934) q[0];
sx q[0];
rz(-2.0598167) q[0];
sx q[0];
rz(-0.10639356) q[0];
rz(1.7875207) q[2];
sx q[2];
rz(-0.43102362) q[2];
sx q[2];
rz(-0.9123411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0178899) q[1];
sx q[1];
rz(-2.5933939) q[1];
sx q[1];
rz(0.75843673) q[1];
rz(-pi) q[2];
rz(0.49523109) q[3];
sx q[3];
rz(-0.81988813) q[3];
sx q[3];
rz(-0.26671975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.52557785) q[2];
sx q[2];
rz(-1.3495477) q[2];
sx q[2];
rz(-2.7363321) q[2];
rz(-1.1303834) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(1.9270886) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2363215) q[0];
sx q[0];
rz(-0.032624809) q[0];
sx q[0];
rz(-0.76500934) q[0];
rz(-2.8969823) q[1];
sx q[1];
rz(-1.8965992) q[1];
sx q[1];
rz(-0.6388706) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7981804) q[0];
sx q[0];
rz(-1.596178) q[0];
sx q[0];
rz(2.4973386) q[0];
rz(2.5193754) q[2];
sx q[2];
rz(-1.8832616) q[2];
sx q[2];
rz(-2.0553596) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42145211) q[1];
sx q[1];
rz(-1.1126592) q[1];
sx q[1];
rz(1.9625743) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8833771) q[3];
sx q[3];
rz(-0.85047837) q[3];
sx q[3];
rz(-2.9434266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8581533) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-0.24432527) q[2];
rz(-1.1089995) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(-0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7636488) q[0];
sx q[0];
rz(-1.0841333) q[0];
sx q[0];
rz(1.9042683) q[0];
rz(-1.411865) q[1];
sx q[1];
rz(-0.4487764) q[1];
sx q[1];
rz(1.3317187) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6520572) q[0];
sx q[0];
rz(-0.79990989) q[0];
sx q[0];
rz(-1.6933939) q[0];
x q[1];
rz(-1.9218512) q[2];
sx q[2];
rz(-1.7760385) q[2];
sx q[2];
rz(0.83265162) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0739286) q[1];
sx q[1];
rz(-0.064753115) q[1];
sx q[1];
rz(-2.6005437) q[1];
x q[2];
rz(0.79257719) q[3];
sx q[3];
rz(-0.78781414) q[3];
sx q[3];
rz(0.28258309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9776913) q[2];
sx q[2];
rz(-2.0745664) q[2];
sx q[2];
rz(2.4816371) q[2];
rz(-1.4466977) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(2.0383294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834943) q[0];
sx q[0];
rz(-2.3907008) q[0];
sx q[0];
rz(1.1335491) q[0];
rz(-2.6253888) q[1];
sx q[1];
rz(-1.3092923) q[1];
sx q[1];
rz(-2.5925327) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20401317) q[0];
sx q[0];
rz(-1.4268865) q[0];
sx q[0];
rz(-1.7372485) q[0];
rz(-pi) q[1];
rz(1.334515) q[2];
sx q[2];
rz(-2.5388814) q[2];
sx q[2];
rz(-0.38379764) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0564894) q[1];
sx q[1];
rz(-1.5055377) q[1];
sx q[1];
rz(0.27537055) q[1];
x q[2];
rz(2.8589949) q[3];
sx q[3];
rz(-1.9535086) q[3];
sx q[3];
rz(-0.0569009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.05222008) q[2];
sx q[2];
rz(-1.0129656) q[2];
sx q[2];
rz(-0.29903665) q[2];
rz(0.61797577) q[3];
sx q[3];
rz(-1.0733913) q[3];
sx q[3];
rz(1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1842136) q[0];
sx q[0];
rz(-0.03802499) q[0];
sx q[0];
rz(2.6755565) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-0.51141089) q[1];
sx q[1];
rz(2.3066511) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3059495) q[0];
sx q[0];
rz(-1.5422652) q[0];
sx q[0];
rz(-3.0898407) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72795002) q[2];
sx q[2];
rz(-3.0228428) q[2];
sx q[2];
rz(2.438022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7606733) q[1];
sx q[1];
rz(-1.6242508) q[1];
sx q[1];
rz(2.0510813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4598875) q[3];
sx q[3];
rz(-0.85452628) q[3];
sx q[3];
rz(-1.1151506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9338108) q[2];
sx q[2];
rz(-0.89263478) q[2];
sx q[2];
rz(-3.1100698) q[2];
rz(2.3837714) q[3];
sx q[3];
rz(-2.0407245) q[3];
sx q[3];
rz(2.5789564) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23444489) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(3.1100519) q[0];
rz(3.0517598) q[1];
sx q[1];
rz(-0.49982163) q[1];
sx q[1];
rz(0.33581844) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0425371) q[0];
sx q[0];
rz(-2.6495948) q[0];
sx q[0];
rz(-0.74721952) q[0];
x q[1];
rz(1.0511984) q[2];
sx q[2];
rz(-0.61695398) q[2];
sx q[2];
rz(1.761957) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0829741) q[1];
sx q[1];
rz(-1.0320283) q[1];
sx q[1];
rz(-2.8584216) q[1];
rz(1.1169711) q[3];
sx q[3];
rz(-1.8403075) q[3];
sx q[3];
rz(-2.5170642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1699367) q[2];
sx q[2];
rz(-1.1145096) q[2];
sx q[2];
rz(-1.8709095) q[2];
rz(3.0830248) q[3];
sx q[3];
rz(-1.4437557) q[3];
sx q[3];
rz(2.7931255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5295277) q[0];
sx q[0];
rz(-2.6030774) q[0];
sx q[0];
rz(-2.8756323) q[0];
rz(0.82414857) q[1];
sx q[1];
rz(-1.3105086) q[1];
sx q[1];
rz(1.6110274) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53479311) q[0];
sx q[0];
rz(-1.3086119) q[0];
sx q[0];
rz(-1.3110089) q[0];
rz(-pi) q[1];
rz(0.058892756) q[2];
sx q[2];
rz(-0.73187553) q[2];
sx q[2];
rz(0.99725396) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2648523) q[1];
sx q[1];
rz(-1.0050259) q[1];
sx q[1];
rz(1.2575722) q[1];
x q[2];
rz(-1.2208185) q[3];
sx q[3];
rz(-0.82982291) q[3];
sx q[3];
rz(-0.32311026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.019235762) q[2];
sx q[2];
rz(-2.0255721) q[2];
sx q[2];
rz(-0.96987152) q[2];
rz(-0.14361778) q[3];
sx q[3];
rz(-2.2031281) q[3];
sx q[3];
rz(2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.403991) q[0];
sx q[0];
rz(-0.25594512) q[0];
sx q[0];
rz(3.0373489) q[0];
rz(0.741611) q[1];
sx q[1];
rz(-2.9545018) q[1];
sx q[1];
rz(0.43824497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7359228) q[0];
sx q[0];
rz(-0.68016648) q[0];
sx q[0];
rz(-2.9783572) q[0];
rz(-pi) q[1];
rz(-0.97317071) q[2];
sx q[2];
rz(-1.3750374) q[2];
sx q[2];
rz(-0.51908609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5129152) q[1];
sx q[1];
rz(-1.6048429) q[1];
sx q[1];
rz(-2.176575) q[1];
x q[2];
rz(0.24298553) q[3];
sx q[3];
rz(-0.54854092) q[3];
sx q[3];
rz(-2.0627956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40615842) q[2];
sx q[2];
rz(-2.8696852) q[2];
sx q[2];
rz(2.820365) q[2];
rz(2.1494703) q[3];
sx q[3];
rz(-2.0480053) q[3];
sx q[3];
rz(1.4120302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.853249) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(2.993592) q[0];
rz(1.1389698) q[1];
sx q[1];
rz(-2.482246) q[1];
sx q[1];
rz(-0.99050561) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873224) q[0];
sx q[0];
rz(-1.089671) q[0];
sx q[0];
rz(-1.4863798) q[0];
x q[1];
rz(-2.4100346) q[2];
sx q[2];
rz(-1.2148982) q[2];
sx q[2];
rz(-1.6936091) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3504178) q[1];
sx q[1];
rz(-1.5759849) q[1];
sx q[1];
rz(1.5395916) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0664987) q[3];
sx q[3];
rz(-1.9073581) q[3];
sx q[3];
rz(2.9838205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.972435) q[2];
sx q[2];
rz(-2.3507698) q[2];
sx q[2];
rz(2.5992744) q[2];
rz(1.4646685) q[3];
sx q[3];
rz(-1.9138391) q[3];
sx q[3];
rz(-2.1577468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2391613) q[0];
sx q[0];
rz(-2.3434174) q[0];
sx q[0];
rz(0.27957988) q[0];
rz(-0.78508776) q[1];
sx q[1];
rz(-2.3468192) q[1];
sx q[1];
rz(2.7395172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1917443) q[0];
sx q[0];
rz(-1.254891) q[0];
sx q[0];
rz(2.6788968) q[0];
rz(-pi) q[1];
rz(-0.78530757) q[2];
sx q[2];
rz(-0.97427893) q[2];
sx q[2];
rz(-1.1889088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8895032) q[1];
sx q[1];
rz(-0.91633233) q[1];
sx q[1];
rz(2.8988125) q[1];
rz(-2.8567186) q[3];
sx q[3];
rz(-1.1732297) q[3];
sx q[3];
rz(-1.8101102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1606007) q[2];
sx q[2];
rz(-0.6610142) q[2];
sx q[2];
rz(0.92657363) q[2];
rz(-2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(-1.2800823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4704623) q[0];
sx q[0];
rz(-2.0736546) q[0];
sx q[0];
rz(0.14508844) q[0];
rz(-3.0628224) q[1];
sx q[1];
rz(-1.7369743) q[1];
sx q[1];
rz(-1.8440934) q[1];
rz(2.3433122) q[2];
sx q[2];
rz(-1.5816806) q[2];
sx q[2];
rz(0.47368373) q[2];
rz(-0.13194247) q[3];
sx q[3];
rz(-2.0880825) q[3];
sx q[3];
rz(-1.1121398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
