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
rz(1.3067955) q[0];
sx q[0];
rz(5.4352221) q[0];
sx q[0];
rz(9.3873831) q[0];
rz(-1.0132064) q[1];
sx q[1];
rz(-0.4816882) q[1];
sx q[1];
rz(-1.6964635) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045528) q[0];
sx q[0];
rz(-1.3890966) q[0];
sx q[0];
rz(-2.7054525) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6736534) q[2];
sx q[2];
rz(-1.5005996) q[2];
sx q[2];
rz(1.6818893) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36276517) q[1];
sx q[1];
rz(-0.76218513) q[1];
sx q[1];
rz(-0.83517167) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3919034) q[3];
sx q[3];
rz(-0.57235347) q[3];
sx q[3];
rz(-1.6088886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41363132) q[2];
sx q[2];
rz(-0.093158826) q[2];
sx q[2];
rz(-1.5067345) q[2];
rz(2.9480751) q[3];
sx q[3];
rz(-0.7842803) q[3];
sx q[3];
rz(0.94582742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042260878) q[0];
sx q[0];
rz(-1.7758545) q[0];
sx q[0];
rz(-2.9625764) q[0];
rz(-1.5366813) q[1];
sx q[1];
rz(-2.8004526) q[1];
sx q[1];
rz(1.590033) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.478708) q[0];
sx q[0];
rz(-0.30618069) q[0];
sx q[0];
rz(1.8992437) q[0];
rz(-pi) q[1];
rz(0.70945246) q[2];
sx q[2];
rz(-1.7691028) q[2];
sx q[2];
rz(2.545216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7187923) q[1];
sx q[1];
rz(-1.8343636) q[1];
sx q[1];
rz(-1.9843319) q[1];
rz(-pi) q[2];
rz(-1.0256944) q[3];
sx q[3];
rz(-1.1602931) q[3];
sx q[3];
rz(-2.5481567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67148525) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(-3.1191077) q[2];
rz(0.89618987) q[3];
sx q[3];
rz(-0.78138566) q[3];
sx q[3];
rz(-1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80980587) q[0];
sx q[0];
rz(-2.9349194) q[0];
sx q[0];
rz(1.608954) q[0];
rz(1.1398075) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(0.87361139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7067817) q[0];
sx q[0];
rz(-2.0764838) q[0];
sx q[0];
rz(-2.3825453) q[0];
rz(-pi) q[1];
rz(2.1223091) q[2];
sx q[2];
rz(-2.1423369) q[2];
sx q[2];
rz(0.71046605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.552678) q[1];
sx q[1];
rz(-0.50929087) q[1];
sx q[1];
rz(0.31917787) q[1];
rz(-pi) q[2];
rz(2.8779712) q[3];
sx q[3];
rz(-0.28214165) q[3];
sx q[3];
rz(-0.96708114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56796271) q[2];
sx q[2];
rz(-2.6831388) q[2];
sx q[2];
rz(-0.47631329) q[2];
rz(1.025398) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(0.29388139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3697701) q[0];
sx q[0];
rz(-0.85404587) q[0];
sx q[0];
rz(0.59637946) q[0];
rz(-2.3884933) q[1];
sx q[1];
rz(-1.785708) q[1];
sx q[1];
rz(0.3515884) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7542587) q[0];
sx q[0];
rz(-1.046858) q[0];
sx q[0];
rz(-1.7225527) q[0];
x q[1];
rz(-2.2542721) q[2];
sx q[2];
rz(-1.9565689) q[2];
sx q[2];
rz(1.9172457) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56573136) q[1];
sx q[1];
rz(-1.8003776) q[1];
sx q[1];
rz(-1.1665535) q[1];
rz(-2.5078689) q[3];
sx q[3];
rz(-1.0175287) q[3];
sx q[3];
rz(2.3330757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4190462) q[2];
sx q[2];
rz(-1.6435577) q[2];
sx q[2];
rz(-2.2124115) q[2];
rz(1.9123214) q[3];
sx q[3];
rz(-0.036673948) q[3];
sx q[3];
rz(-1.2693955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1790328) q[0];
sx q[0];
rz(-2.5293009) q[0];
sx q[0];
rz(0.52781934) q[0];
rz(-2.5714696) q[1];
sx q[1];
rz(-0.56990439) q[1];
sx q[1];
rz(0.39400563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4657426) q[0];
sx q[0];
rz(-2.619474) q[0];
sx q[0];
rz(2.7444849) q[0];
rz(-pi) q[1];
rz(2.6232004) q[2];
sx q[2];
rz(-2.4817991) q[2];
sx q[2];
rz(-2.348748) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3574419) q[1];
sx q[1];
rz(-1.3489745) q[1];
sx q[1];
rz(0.12089575) q[1];
x q[2];
rz(2.4179055) q[3];
sx q[3];
rz(-2.1877648) q[3];
sx q[3];
rz(0.2816412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3923308) q[2];
sx q[2];
rz(-1.8529842) q[2];
sx q[2];
rz(-2.0070845) q[2];
rz(0.025731651) q[3];
sx q[3];
rz(-1.0545571) q[3];
sx q[3];
rz(2.499685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-2.5662017) q[0];
sx q[0];
rz(-1.8113149) q[0];
sx q[0];
rz(0.45912418) q[0];
rz(0.8210012) q[1];
sx q[1];
rz(-2.2586925) q[1];
sx q[1];
rz(-0.51148907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3200681) q[0];
sx q[0];
rz(-1.2787158) q[0];
sx q[0];
rz(-1.0758557) q[0];
rz(-1.3364001) q[2];
sx q[2];
rz(-2.2210741) q[2];
sx q[2];
rz(-2.5087207) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9238115) q[1];
sx q[1];
rz(-2.1379323) q[1];
sx q[1];
rz(1.3918716) q[1];
rz(-pi) q[2];
rz(-1.4019764) q[3];
sx q[3];
rz(-0.69069117) q[3];
sx q[3];
rz(-2.9258779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56200999) q[2];
sx q[2];
rz(-1.1082114) q[2];
sx q[2];
rz(-2.3902334) q[2];
rz(2.5747418) q[3];
sx q[3];
rz(-1.5375429) q[3];
sx q[3];
rz(1.8142627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71996561) q[0];
sx q[0];
rz(-2.9599074) q[0];
sx q[0];
rz(-0.0075465329) q[0];
rz(-2.201572) q[1];
sx q[1];
rz(-2.4766141) q[1];
sx q[1];
rz(3.0090581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36245334) q[0];
sx q[0];
rz(-0.7374239) q[0];
sx q[0];
rz(2.2064184) q[0];
rz(-pi) q[1];
rz(1.7039677) q[2];
sx q[2];
rz(-1.9308763) q[2];
sx q[2];
rz(0.11225004) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.346437) q[1];
sx q[1];
rz(-2.7043685) q[1];
sx q[1];
rz(2.6889685) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1121944) q[3];
sx q[3];
rz(-2.762714) q[3];
sx q[3];
rz(-2.1357812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14564766) q[2];
sx q[2];
rz(-1.6463248) q[2];
sx q[2];
rz(3.0403467) q[2];
rz(0.64949399) q[3];
sx q[3];
rz(-0.5972623) q[3];
sx q[3];
rz(2.7818642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.41530171) q[0];
sx q[0];
rz(-1.2969718) q[0];
sx q[0];
rz(-1.8743961) q[0];
rz(1.2665117) q[1];
sx q[1];
rz(-2.9837065) q[1];
sx q[1];
rz(-2.3983541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38108724) q[0];
sx q[0];
rz(-0.60956565) q[0];
sx q[0];
rz(2.7838092) q[0];
x q[1];
rz(2.350528) q[2];
sx q[2];
rz(-1.3692489) q[2];
sx q[2];
rz(0.1669875) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0580045) q[1];
sx q[1];
rz(-2.8561855) q[1];
sx q[1];
rz(-2.1607375) q[1];
rz(-pi) q[2];
rz(-2.2970639) q[3];
sx q[3];
rz(-2.2457079) q[3];
sx q[3];
rz(-2.9425987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0172952) q[2];
sx q[2];
rz(-2.9399019) q[2];
sx q[2];
rz(-1.4264433) q[2];
rz(1.0157478) q[3];
sx q[3];
rz(-1.7044715) q[3];
sx q[3];
rz(0.82971853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5947241) q[0];
sx q[0];
rz(-0.65383738) q[0];
sx q[0];
rz(-3.1256909) q[0];
rz(-1.6390027) q[1];
sx q[1];
rz(-2.465261) q[1];
sx q[1];
rz(-2.1471088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8581945) q[0];
sx q[0];
rz(-0.97351979) q[0];
sx q[0];
rz(-0.47391639) q[0];
rz(-pi) q[1];
rz(-2.9839462) q[2];
sx q[2];
rz(-0.38262832) q[2];
sx q[2];
rz(-1.1747109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9084103) q[1];
sx q[1];
rz(-0.841978) q[1];
sx q[1];
rz(-0.89836095) q[1];
x q[2];
rz(0.1511622) q[3];
sx q[3];
rz(-2.0796607) q[3];
sx q[3];
rz(0.14948949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1556306) q[2];
sx q[2];
rz(-1.4302284) q[2];
sx q[2];
rz(0.17781167) q[2];
rz(-1.6832247) q[3];
sx q[3];
rz(-0.42176133) q[3];
sx q[3];
rz(-1.873707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630702) q[0];
sx q[0];
rz(-1.2489742) q[0];
sx q[0];
rz(-1.2836237) q[0];
rz(-2.3578857) q[1];
sx q[1];
rz(-2.3678534) q[1];
sx q[1];
rz(1.8342038) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15128216) q[0];
sx q[0];
rz(-1.6626004) q[0];
sx q[0];
rz(0.43938322) q[0];
rz(-3.1112017) q[2];
sx q[2];
rz(-1.9585397) q[2];
sx q[2];
rz(-3.0770242) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5518903) q[1];
sx q[1];
rz(-1.74472) q[1];
sx q[1];
rz(-3.0822192) q[1];
rz(-pi) q[2];
rz(-2.1940007) q[3];
sx q[3];
rz(-0.64631185) q[3];
sx q[3];
rz(1.0912967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54138294) q[2];
sx q[2];
rz(-1.7317438) q[2];
sx q[2];
rz(-0.48995885) q[2];
rz(-1.1146891) q[3];
sx q[3];
rz(-1.1763108) q[3];
sx q[3];
rz(2.8118242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020323044) q[0];
sx q[0];
rz(-1.9135495) q[0];
sx q[0];
rz(1.3187153) q[0];
rz(1.8439138) q[1];
sx q[1];
rz(-0.73328016) q[1];
sx q[1];
rz(2.7136623) q[1];
rz(-2.8823356) q[2];
sx q[2];
rz(-1.8198063) q[2];
sx q[2];
rz(1.168269) q[2];
rz(-1.8036203) q[3];
sx q[3];
rz(-2.4216087) q[3];
sx q[3];
rz(-2.5523228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
