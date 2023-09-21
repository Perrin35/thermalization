OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(4.6306643) q[0];
sx q[0];
rz(10.319933) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(-1.4562343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.668894) q[0];
sx q[0];
rz(-1.5510674) q[0];
sx q[0];
rz(-1.4399477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0257172) q[2];
sx q[2];
rz(-1.8946049) q[2];
sx q[2];
rz(2.2562502) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40741062) q[1];
sx q[1];
rz(-1.4200746) q[1];
sx q[1];
rz(-0.68024866) q[1];
rz(-pi) q[2];
rz(-2.9772894) q[3];
sx q[3];
rz(-2.8176753) q[3];
sx q[3];
rz(1.8620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(2.6888729) q[2];
rz(2.9833941) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(1.989495) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-2.6706085) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5841056) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(-1.4529865) q[0];
x q[1];
rz(2.9071964) q[2];
sx q[2];
rz(-2.1550551) q[2];
sx q[2];
rz(-0.61569475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60103154) q[1];
sx q[1];
rz(-0.49113501) q[1];
sx q[1];
rz(-1.4547552) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4460375) q[3];
sx q[3];
rz(-1.6323946) q[3];
sx q[3];
rz(1.0775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26329041) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(-0.73633206) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(-2.7022865) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8894316) q[0];
sx q[0];
rz(-1.6003506) q[0];
sx q[0];
rz(1.6512524) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0053824) q[2];
sx q[2];
rz(-0.36362193) q[2];
sx q[2];
rz(-1.9759535) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6046391) q[1];
sx q[1];
rz(-2.1770658) q[1];
sx q[1];
rz(0.21142516) q[1];
rz(-pi) q[2];
rz(-1.0147694) q[3];
sx q[3];
rz(-1.945567) q[3];
sx q[3];
rz(0.32885636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(2.9411194) q[2];
rz(-0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(-1.8970998) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(-2.2669852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133102) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(-1.8391795) q[0];
rz(-pi) q[1];
rz(0.34747296) q[2];
sx q[2];
rz(-2.7345737) q[2];
sx q[2];
rz(2.4328872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9676535) q[1];
sx q[1];
rz(-1.9441009) q[1];
sx q[1];
rz(-1.3268382) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6455669) q[3];
sx q[3];
rz(-1.6755591) q[3];
sx q[3];
rz(1.2391702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.4271663) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(0.57058913) q[0];
rz(2.5866306) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(2.3983009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5071881) q[0];
sx q[0];
rz(-0.92538639) q[0];
sx q[0];
rz(-0.27642823) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18851738) q[2];
sx q[2];
rz(-2.6326615) q[2];
sx q[2];
rz(0.9045507) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6965461) q[1];
sx q[1];
rz(-1.7940709) q[1];
sx q[1];
rz(0.22534196) q[1];
rz(-1.284243) q[3];
sx q[3];
rz(-2.3544469) q[3];
sx q[3];
rz(-0.033165008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-2.8732079) q[2];
rz(-1.0466446) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-1.5140117) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(0.2535893) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(2.0862897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9297766) q[0];
sx q[0];
rz(-1.7535216) q[0];
sx q[0];
rz(2.1972448) q[0];
rz(2.2248473) q[2];
sx q[2];
rz(-1.142821) q[2];
sx q[2];
rz(-0.59631729) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.222059) q[1];
sx q[1];
rz(-0.93898458) q[1];
sx q[1];
rz(0.49948378) q[1];
x q[2];
rz(-0.25552337) q[3];
sx q[3];
rz(-2.3821085) q[3];
sx q[3];
rz(1.4042735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(-2.3690467) q[0];
rz(1.4121217) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-0.57377446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44874292) q[0];
sx q[0];
rz(-2.9208555) q[0];
sx q[0];
rz(2.0953395) q[0];
rz(0.022106604) q[2];
sx q[2];
rz(-1.7709641) q[2];
sx q[2];
rz(2.9526763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3473914) q[1];
sx q[1];
rz(-0.83354356) q[1];
sx q[1];
rz(-2.3458523) q[1];
rz(-pi) q[2];
rz(1.8049559) q[3];
sx q[3];
rz(-1.4244716) q[3];
sx q[3];
rz(0.22406604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1022169) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(-2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-1.8687826) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(-1.2957113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0062795) q[0];
sx q[0];
rz(-1.2426002) q[0];
sx q[0];
rz(-2.6562064) q[0];
rz(0.92408085) q[2];
sx q[2];
rz(-2.7087822) q[2];
sx q[2];
rz(2.2518287) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1295373) q[1];
sx q[1];
rz(-2.4730198) q[1];
sx q[1];
rz(1.936391) q[1];
x q[2];
rz(3.1136884) q[3];
sx q[3];
rz(-1.1301665) q[3];
sx q[3];
rz(2.2925216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(-2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(-0.67767674) q[0];
rz(0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3121376) q[0];
sx q[0];
rz(-1.9337618) q[0];
sx q[0];
rz(-1.3296933) q[0];
rz(0.067473472) q[2];
sx q[2];
rz(-2.3258492) q[2];
sx q[2];
rz(0.93630723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2182525) q[1];
sx q[1];
rz(-0.76949161) q[1];
sx q[1];
rz(-1.8583276) q[1];
rz(-pi) q[2];
rz(2.8832199) q[3];
sx q[3];
rz(-1.6515454) q[3];
sx q[3];
rz(-2.0852058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.2667123) q[2];
rz(-2.1789815) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4203913) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(-2.904073) q[0];
rz(2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(0.231803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7451413) q[0];
sx q[0];
rz(-1.475435) q[0];
sx q[0];
rz(-2.0800637) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3583057) q[2];
sx q[2];
rz(-0.33595339) q[2];
sx q[2];
rz(0.17620262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9634562) q[1];
sx q[1];
rz(-1.8354715) q[1];
sx q[1];
rz(-1.4241649) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3123355) q[3];
sx q[3];
rz(-0.70445326) q[3];
sx q[3];
rz(2.7750912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-2.0533662) q[2];
rz(-0.38816372) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.3303403) q[2];
sx q[2];
rz(-1.2546872) q[2];
sx q[2];
rz(1.6457002) q[2];
rz(1.5512636) q[3];
sx q[3];
rz(-1.7587147) q[3];
sx q[3];
rz(1.3295909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];