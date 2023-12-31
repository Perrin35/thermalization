OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(2.6161939) q[1];
sx q[1];
rz(8.4470196) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5241961) q[0];
sx q[0];
rz(-2.8688736) q[0];
sx q[0];
rz(0.32193907) q[0];
rz(-pi) q[1];
rz(-1.0563645) q[2];
sx q[2];
rz(-1.8407514) q[2];
sx q[2];
rz(-1.1411238) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2797151) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(-3.0008297) q[1];
x q[2];
rz(2.9631859) q[3];
sx q[3];
rz(-0.79531407) q[3];
sx q[3];
rz(-0.91245302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67291659) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(-2.9428234) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7070049) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(-0.36112753) q[0];
rz(-1.7065642) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(0.8180058) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0571787) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(1.8684698) q[0];
x q[1];
rz(3.0188574) q[2];
sx q[2];
rz(-1.9397748) q[2];
sx q[2];
rz(2.4462647) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4600196) q[1];
sx q[1];
rz(-2.1732554) q[1];
sx q[1];
rz(-2.6006992) q[1];
rz(-pi) q[2];
rz(1.4217671) q[3];
sx q[3];
rz(-2.7741198) q[3];
sx q[3];
rz(-0.056882337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.25257418) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(-1.7209631) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(2.7533598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7398359) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(0.31618205) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(-0.2972163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5706144) q[0];
sx q[0];
rz(-1.1665205) q[0];
sx q[0];
rz(1.6398318) q[0];
rz(-2.3702413) q[2];
sx q[2];
rz(-1.0372247) q[2];
sx q[2];
rz(0.16650621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19993648) q[1];
sx q[1];
rz(-1.3841277) q[1];
sx q[1];
rz(3.093064) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66629569) q[3];
sx q[3];
rz(-0.94216457) q[3];
sx q[3];
rz(-2.6578238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-0.42759839) q[2];
rz(1.9528495) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(-0.71290839) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(0.68177044) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101689) q[0];
sx q[0];
rz(-1.3623326) q[0];
sx q[0];
rz(1.3039949) q[0];
rz(-pi) q[1];
rz(2.9245124) q[2];
sx q[2];
rz(-0.97765572) q[2];
sx q[2];
rz(0.3074239) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3817953) q[1];
sx q[1];
rz(-1.1361546) q[1];
sx q[1];
rz(-0.97310193) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81687974) q[3];
sx q[3];
rz(-0.48135346) q[3];
sx q[3];
rz(2.4933185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.008808) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(2.183389) q[2];
rz(1.0664252) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(-2.5812896) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(1.5195742) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9902089) q[0];
sx q[0];
rz(-2.4575893) q[0];
sx q[0];
rz(-0.74599501) q[0];
x q[1];
rz(-1.0137453) q[2];
sx q[2];
rz(-2.1079014) q[2];
sx q[2];
rz(0.86442664) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2677742) q[1];
sx q[1];
rz(-2.3582637) q[1];
sx q[1];
rz(2.8591213) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.394746) q[3];
sx q[3];
rz(-2.5213443) q[3];
sx q[3];
rz(0.32905096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6891629) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(0.2229283) q[2];
rz(0.034742268) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(0.071578659) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(2.0715332) q[0];
rz(1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(1.8575352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.703008) q[0];
sx q[0];
rz(-2.2692338) q[0];
sx q[0];
rz(0.75011487) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9610923) q[2];
sx q[2];
rz(-1.3109866) q[2];
sx q[2];
rz(2.5026472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2388873) q[1];
sx q[1];
rz(-1.1088015) q[1];
sx q[1];
rz(1.8932896) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1070126) q[3];
sx q[3];
rz(-1.3048733) q[3];
sx q[3];
rz(-1.769161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(0.83667886) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(-1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56629431) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(-0.42770806) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(-2.2033851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3283549) q[0];
sx q[0];
rz(-0.52044808) q[0];
sx q[0];
rz(0.92207272) q[0];
x q[1];
rz(1.0023408) q[2];
sx q[2];
rz(-0.92407862) q[2];
sx q[2];
rz(-1.2127753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2141014) q[1];
sx q[1];
rz(-2.2117105) q[1];
sx q[1];
rz(0.44968857) q[1];
rz(-1.8754962) q[3];
sx q[3];
rz(-0.88677553) q[3];
sx q[3];
rz(1.5414433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87970916) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(-1.7555457) q[2];
rz(-1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(1.7077131) q[0];
rz(-1.2738312) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(2.3103255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87330504) q[0];
sx q[0];
rz(-1.9718861) q[0];
sx q[0];
rz(1.7307161) q[0];
rz(-pi) q[1];
rz(-0.93584658) q[2];
sx q[2];
rz(-1.4668462) q[2];
sx q[2];
rz(-2.7452591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89325209) q[1];
sx q[1];
rz(-0.89683956) q[1];
sx q[1];
rz(-2.2158951) q[1];
rz(-2.3723699) q[3];
sx q[3];
rz(-2.8260494) q[3];
sx q[3];
rz(-0.10570082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(1.9780805) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7678541) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(-1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(2.0358553) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39767299) q[0];
sx q[0];
rz(-2.0154675) q[0];
sx q[0];
rz(1.7001274) q[0];
x q[1];
rz(-2.2374723) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(-2.3711575) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3016262) q[1];
sx q[1];
rz(-1.3857538) q[1];
sx q[1];
rz(1.3355096) q[1];
rz(2.4615199) q[3];
sx q[3];
rz(-2.701093) q[3];
sx q[3];
rz(2.3624453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5332807) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(0.14979714) q[2];
rz(1.7685361) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(-1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(2.9737934) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4807178) q[0];
sx q[0];
rz(-1.6010451) q[0];
sx q[0];
rz(1.589993) q[0];
rz(2.5334353) q[2];
sx q[2];
rz(-1.6981352) q[2];
sx q[2];
rz(1.5437768) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2044636) q[1];
sx q[1];
rz(-0.58848721) q[1];
sx q[1];
rz(-0.012339331) q[1];
rz(-0.76370244) q[3];
sx q[3];
rz(-0.16994952) q[3];
sx q[3];
rz(-2.9156239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(0.76114571) q[2];
rz(3.051565) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54820838) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(0.38800115) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.3315196) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(-0.17100632) q[3];
sx q[3];
rz(-2.1033559) q[3];
sx q[3];
rz(-0.81567473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
