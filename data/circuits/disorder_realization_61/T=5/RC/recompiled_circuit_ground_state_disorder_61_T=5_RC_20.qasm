OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.073590241) q[0];
sx q[0];
rz(4.6103067) q[0];
sx q[0];
rz(8.8972916) q[0];
rz(-2.5975851) q[1];
sx q[1];
rz(-0.44096947) q[1];
sx q[1];
rz(0.0013594065) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6327335) q[0];
sx q[0];
rz(-1.5799045) q[0];
sx q[0];
rz(-1.5004116) q[0];
rz(0.33250471) q[2];
sx q[2];
rz(-0.63897479) q[2];
sx q[2];
rz(-2.9457991) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4365523) q[1];
sx q[1];
rz(-2.4715354) q[1];
sx q[1];
rz(0.64924182) q[1];
rz(1.5415927) q[3];
sx q[3];
rz(-1.7701704) q[3];
sx q[3];
rz(0.77413156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9881543) q[2];
sx q[2];
rz(-2.4997288) q[2];
sx q[2];
rz(1.362907) q[2];
rz(0.07987944) q[3];
sx q[3];
rz(-0.86013836) q[3];
sx q[3];
rz(0.0063272198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36413559) q[0];
sx q[0];
rz(-1.2292925) q[0];
sx q[0];
rz(0.51032132) q[0];
rz(-1.8136884) q[1];
sx q[1];
rz(-2.4102305) q[1];
sx q[1];
rz(2.7046943) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6104159) q[0];
sx q[0];
rz(-2.057909) q[0];
sx q[0];
rz(-0.25419828) q[0];
rz(-pi) q[1];
rz(-1.8067752) q[2];
sx q[2];
rz(-1.1369822) q[2];
sx q[2];
rz(0.12820511) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.092097923) q[1];
sx q[1];
rz(-1.5391401) q[1];
sx q[1];
rz(-2.3718216) q[1];
rz(-1.1673959) q[3];
sx q[3];
rz(-1.511949) q[3];
sx q[3];
rz(1.7047391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83866155) q[2];
sx q[2];
rz(-0.51753664) q[2];
sx q[2];
rz(-0.17520629) q[2];
rz(-3.0025499) q[3];
sx q[3];
rz(-1.4929205) q[3];
sx q[3];
rz(-0.89865941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43056968) q[0];
sx q[0];
rz(-0.88453203) q[0];
sx q[0];
rz(-0.36054605) q[0];
rz(-1.644246) q[1];
sx q[1];
rz(-1.2253573) q[1];
sx q[1];
rz(0.59251934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07537341) q[0];
sx q[0];
rz(-2.1751715) q[0];
sx q[0];
rz(-0.79232596) q[0];
rz(-2.5056867) q[2];
sx q[2];
rz(-1.4472258) q[2];
sx q[2];
rz(1.6963144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.097774769) q[1];
sx q[1];
rz(-1.9182529) q[1];
sx q[1];
rz(-3.0128885) q[1];
rz(-pi) q[2];
rz(-1.1441889) q[3];
sx q[3];
rz(-2.0022829) q[3];
sx q[3];
rz(-0.28080597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7854008) q[2];
sx q[2];
rz(-0.56739002) q[2];
sx q[2];
rz(0.94245911) q[2];
rz(0.59400195) q[3];
sx q[3];
rz(-2.3022251) q[3];
sx q[3];
rz(-0.12282898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9451611) q[0];
sx q[0];
rz(-2.5158947) q[0];
sx q[0];
rz(2.4932267) q[0];
rz(-0.87573373) q[1];
sx q[1];
rz(-1.5114732) q[1];
sx q[1];
rz(-1.3513563) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1296788) q[0];
sx q[0];
rz(-1.5317763) q[0];
sx q[0];
rz(2.3277548) q[0];
x q[1];
rz(2.3301131) q[2];
sx q[2];
rz(-1.7042024) q[2];
sx q[2];
rz(1.8638944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42235628) q[1];
sx q[1];
rz(-1.5420776) q[1];
sx q[1];
rz(2.2442596) q[1];
rz(-pi) q[2];
rz(0.38299527) q[3];
sx q[3];
rz(-1.621709) q[3];
sx q[3];
rz(1.6671772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15328345) q[2];
sx q[2];
rz(-0.51258665) q[2];
sx q[2];
rz(1.2549866) q[2];
rz(2.5399688) q[3];
sx q[3];
rz(-2.2646077) q[3];
sx q[3];
rz(-1.6833444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772142) q[0];
sx q[0];
rz(-2.0266396) q[0];
sx q[0];
rz(0.13378046) q[0];
rz(-1.6556219) q[1];
sx q[1];
rz(-2.8120815) q[1];
sx q[1];
rz(-0.1718743) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3456589) q[0];
sx q[0];
rz(-0.62218124) q[0];
sx q[0];
rz(2.2337929) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3277289) q[2];
sx q[2];
rz(-1.8447734) q[2];
sx q[2];
rz(-2.8629288) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6407224) q[1];
sx q[1];
rz(-1.419104) q[1];
sx q[1];
rz(-1.9239363) q[1];
x q[2];
rz(-0.47953812) q[3];
sx q[3];
rz(-1.1894636) q[3];
sx q[3];
rz(1.6552192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6493426) q[2];
sx q[2];
rz(-0.52521962) q[2];
sx q[2];
rz(1.4225175) q[2];
rz(2.9366142) q[3];
sx q[3];
rz(-0.96115464) q[3];
sx q[3];
rz(0.80406308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64254665) q[0];
sx q[0];
rz(-0.63669425) q[0];
sx q[0];
rz(-2.9472886) q[0];
rz(0.64584541) q[1];
sx q[1];
rz(-1.9648809) q[1];
sx q[1];
rz(-2.7913854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98726942) q[0];
sx q[0];
rz(-1.435017) q[0];
sx q[0];
rz(-2.6063347) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5123213) q[2];
sx q[2];
rz(-2.2221178) q[2];
sx q[2];
rz(-0.68683147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3694075) q[1];
sx q[1];
rz(-0.52474743) q[1];
sx q[1];
rz(-0.16903846) q[1];
rz(-3.0549235) q[3];
sx q[3];
rz(-0.67207591) q[3];
sx q[3];
rz(1.2487401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0575867) q[2];
sx q[2];
rz(-2.7454594) q[2];
sx q[2];
rz(2.8193889) q[2];
rz(-0.75994879) q[3];
sx q[3];
rz(-0.51637572) q[3];
sx q[3];
rz(-3.0907536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26611227) q[0];
sx q[0];
rz(-0.66829824) q[0];
sx q[0];
rz(2.42591) q[0];
rz(3.074805) q[1];
sx q[1];
rz(-2.5804434) q[1];
sx q[1];
rz(-2.7776048) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20446312) q[0];
sx q[0];
rz(-1.9162906) q[0];
sx q[0];
rz(0.95204592) q[0];
x q[1];
rz(2.5189713) q[2];
sx q[2];
rz(-1.004815) q[2];
sx q[2];
rz(-2.7634668) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5426483) q[1];
sx q[1];
rz(-2.3348044) q[1];
sx q[1];
rz(1.548195) q[1];
rz(-pi) q[2];
rz(-2.4956483) q[3];
sx q[3];
rz(-1.6872354) q[3];
sx q[3];
rz(-1.8692335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62465847) q[2];
sx q[2];
rz(-2.2397569) q[2];
sx q[2];
rz(1.9067524) q[2];
rz(-0.45461795) q[3];
sx q[3];
rz(-1.8574628) q[3];
sx q[3];
rz(-3.0668018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55998498) q[0];
sx q[0];
rz(-3.0261664) q[0];
sx q[0];
rz(-0.64062947) q[0];
rz(-0.36264125) q[1];
sx q[1];
rz(-2.8022712) q[1];
sx q[1];
rz(1.6222662) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6453199) q[0];
sx q[0];
rz(-0.64457601) q[0];
sx q[0];
rz(-0.98981895) q[0];
x q[1];
rz(2.8192384) q[2];
sx q[2];
rz(-0.75705069) q[2];
sx q[2];
rz(1.978394) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4942651) q[1];
sx q[1];
rz(-1.4785452) q[1];
sx q[1];
rz(0.89362545) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10102102) q[3];
sx q[3];
rz(-1.1005792) q[3];
sx q[3];
rz(0.322099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92056876) q[2];
sx q[2];
rz(-1.8683044) q[2];
sx q[2];
rz(-1.1576687) q[2];
rz(-2.8122592) q[3];
sx q[3];
rz(-2.9521827) q[3];
sx q[3];
rz(0.91387373) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2374903) q[0];
sx q[0];
rz(-0.79039031) q[0];
sx q[0];
rz(-0.41123408) q[0];
rz(2.5143738) q[1];
sx q[1];
rz(-0.5694446) q[1];
sx q[1];
rz(1.8597182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3498308) q[0];
sx q[0];
rz(-1.786953) q[0];
sx q[0];
rz(2.7695275) q[0];
rz(-0.87780805) q[2];
sx q[2];
rz(-0.77251947) q[2];
sx q[2];
rz(2.7038717) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5465849) q[1];
sx q[1];
rz(-0.91239415) q[1];
sx q[1];
rz(0.14283544) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19880812) q[3];
sx q[3];
rz(-1.947177) q[3];
sx q[3];
rz(2.876296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7810479) q[2];
sx q[2];
rz(-1.5881528) q[2];
sx q[2];
rz(2.1629199) q[2];
rz(0.22748889) q[3];
sx q[3];
rz(-0.74930185) q[3];
sx q[3];
rz(2.6849875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6351673) q[0];
sx q[0];
rz(-3.083482) q[0];
sx q[0];
rz(2.3458922) q[0];
rz(-2.6129163) q[1];
sx q[1];
rz(-0.987459) q[1];
sx q[1];
rz(-0.19435571) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34204179) q[0];
sx q[0];
rz(-0.38368762) q[0];
sx q[0];
rz(1.3253401) q[0];
x q[1];
rz(1.8647381) q[2];
sx q[2];
rz(-1.1938098) q[2];
sx q[2];
rz(-1.3177208) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0463341) q[1];
sx q[1];
rz(-1.4459064) q[1];
sx q[1];
rz(2.6948638) q[1];
rz(-pi) q[2];
rz(-1.0448669) q[3];
sx q[3];
rz(-1.1218438) q[3];
sx q[3];
rz(0.26490213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2928894) q[2];
sx q[2];
rz(-1.6699426) q[2];
sx q[2];
rz(2.8422152) q[2];
rz(2.7401466) q[3];
sx q[3];
rz(-2.5555988) q[3];
sx q[3];
rz(2.6011023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837579) q[0];
sx q[0];
rz(-0.83136375) q[0];
sx q[0];
rz(-1.2663483) q[0];
rz(0.10706317) q[1];
sx q[1];
rz(-0.77580794) q[1];
sx q[1];
rz(-1.2931783) q[1];
rz(-0.54218311) q[2];
sx q[2];
rz(-1.6169006) q[2];
sx q[2];
rz(-2.9666881) q[2];
rz(-0.43184256) q[3];
sx q[3];
rz(-2.8819537) q[3];
sx q[3];
rz(1.0859539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
