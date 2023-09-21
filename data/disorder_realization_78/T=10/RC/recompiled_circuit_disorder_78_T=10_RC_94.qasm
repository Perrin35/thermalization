OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(2.2489927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4757724) q[0];
sx q[0];
rz(-1.2736397) q[0];
sx q[0];
rz(-2.1509403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66944389) q[2];
sx q[2];
rz(-1.1602243) q[2];
sx q[2];
rz(0.66561156) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8436369) q[1];
sx q[1];
rz(-1.2601818) q[1];
sx q[1];
rz(0.85308869) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6139908) q[3];
sx q[3];
rz(-0.82026635) q[3];
sx q[3];
rz(2.299472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(-1.1606476) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(1.2469762) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(1.1670246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8407699) q[0];
sx q[0];
rz(-1.6292028) q[0];
sx q[0];
rz(1.3129243) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8950047) q[2];
sx q[2];
rz(-1.0435259) q[2];
sx q[2];
rz(-1.2621244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0868623) q[1];
sx q[1];
rz(-2.3332101) q[1];
sx q[1];
rz(-3.0562835) q[1];
x q[2];
rz(-0.67006536) q[3];
sx q[3];
rz(-1.9668005) q[3];
sx q[3];
rz(-2.0585287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(0.21437422) q[2];
rz(-0.073444627) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4784933) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(1.0148369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44714156) q[0];
sx q[0];
rz(-1.6269636) q[0];
sx q[0];
rz(0.64356128) q[0];
rz(-1.048676) q[2];
sx q[2];
rz(-2.6819957) q[2];
sx q[2];
rz(-2.2082579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9008357) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(-2.1686173) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2567743) q[3];
sx q[3];
rz(-1.775145) q[3];
sx q[3];
rz(0.94061461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0959452) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(2.4441161) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(-2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(0.43193257) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(2.5057709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63319249) q[0];
sx q[0];
rz(-1.0754555) q[0];
sx q[0];
rz(-0.20423996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3859149) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(-0.76279574) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3969877) q[1];
sx q[1];
rz(-2.3481391) q[1];
sx q[1];
rz(-2.8810487) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5084247) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(-1.4299973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9277966) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-0.14373246) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(2.0671663) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(2.863046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700096) q[0];
sx q[0];
rz(-2.2402813) q[0];
sx q[0];
rz(-0.55218009) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61168806) q[2];
sx q[2];
rz(-1.4809161) q[2];
sx q[2];
rz(2.9181366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46036938) q[1];
sx q[1];
rz(-1.9378098) q[1];
sx q[1];
rz(-0.91194921) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81551084) q[3];
sx q[3];
rz(-1.1482571) q[3];
sx q[3];
rz(2.5583207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(0.29423514) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1086403) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-3.1325353) q[0];
rz(-2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(-0.10805282) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51443433) q[0];
sx q[0];
rz(-1.2540199) q[0];
sx q[0];
rz(0.47382521) q[0];
rz(-pi) q[1];
rz(1.4322386) q[2];
sx q[2];
rz(-1.9493305) q[2];
sx q[2];
rz(-2.0356503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6335771) q[1];
sx q[1];
rz(-2.0532236) q[1];
sx q[1];
rz(-1.1303933) q[1];
rz(-pi) q[2];
rz(1.1507387) q[3];
sx q[3];
rz(-0.58745158) q[3];
sx q[3];
rz(-2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1487427) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(-3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909797) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-0.65761956) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(0.89362842) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4511787) q[0];
sx q[0];
rz(-0.72890857) q[0];
sx q[0];
rz(-2.6364987) q[0];
x q[1];
rz(-2.8473179) q[2];
sx q[2];
rz(-1.1106967) q[2];
sx q[2];
rz(2.0725046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9295846) q[1];
sx q[1];
rz(-1.3699023) q[1];
sx q[1];
rz(2.4561873) q[1];
rz(-0.75974792) q[3];
sx q[3];
rz(-1.500251) q[3];
sx q[3];
rz(-2.7618046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7509193) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(-0.52948362) q[2];
rz(-2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3787518) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(1.09028) q[0];
rz(0.11225637) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.1539248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40076462) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(-3.0853737) q[0];
rz(-0.23837337) q[2];
sx q[2];
rz(-0.56406883) q[2];
sx q[2];
rz(-0.43018815) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6653319) q[1];
sx q[1];
rz(-1.6325103) q[1];
sx q[1];
rz(-2.3803821) q[1];
x q[2];
rz(-0.38903799) q[3];
sx q[3];
rz(-0.93718796) q[3];
sx q[3];
rz(2.6694359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.170084) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78478783) q[0];
sx q[0];
rz(-0.49236449) q[0];
sx q[0];
rz(0.82226336) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27530382) q[2];
sx q[2];
rz(-2.346056) q[2];
sx q[2];
rz(-0.77992935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46718405) q[1];
sx q[1];
rz(-1.7519752) q[1];
sx q[1];
rz(0.49023899) q[1];
x q[2];
rz(-0.13799237) q[3];
sx q[3];
rz(-2.237461) q[3];
sx q[3];
rz(0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3141979) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-0.38273746) q[2];
rz(-2.2132204) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9160354) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(0.25892192) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-2.6616667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7326638) q[0];
sx q[0];
rz(-0.73584475) q[0];
sx q[0];
rz(-1.0622513) q[0];
rz(-pi) q[1];
rz(0.23649044) q[2];
sx q[2];
rz(-0.91954008) q[2];
sx q[2];
rz(-2.1361534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2952134) q[1];
sx q[1];
rz(-1.0615674) q[1];
sx q[1];
rz(-1.9766115) q[1];
rz(2.6142526) q[3];
sx q[3];
rz(-2.1747327) q[3];
sx q[3];
rz(-0.37249836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2991128) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823572) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(1.7469035) q[2];
sx q[2];
rz(-2.0460143) q[2];
sx q[2];
rz(2.569414) q[2];
rz(2.3748623) q[3];
sx q[3];
rz(-0.37692108) q[3];
sx q[3];
rz(2.2212096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
