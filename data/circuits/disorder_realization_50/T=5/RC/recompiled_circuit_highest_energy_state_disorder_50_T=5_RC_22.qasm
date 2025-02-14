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
rz(0.69831508) q[0];
sx q[0];
rz(2.7628216) q[0];
sx q[0];
rz(8.2674352) q[0];
rz(10.209822) q[1];
sx q[1];
rz(0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136381) q[0];
sx q[0];
rz(-2.1305973) q[0];
sx q[0];
rz(-2.7877121) q[0];
rz(-pi) q[1];
rz(-0.73424642) q[2];
sx q[2];
rz(-0.99319211) q[2];
sx q[2];
rz(-0.89370382) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4545645) q[1];
sx q[1];
rz(-0.44378456) q[1];
sx q[1];
rz(-3.0566178) q[1];
rz(1.5230701) q[3];
sx q[3];
rz(-2.5807947) q[3];
sx q[3];
rz(2.5530397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.779125) q[2];
sx q[2];
rz(-0.96131009) q[2];
sx q[2];
rz(0.15164068) q[2];
rz(-2.9996297) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-2.0040472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4769984) q[0];
sx q[0];
rz(-0.73200309) q[0];
sx q[0];
rz(-0.81749302) q[0];
rz(-3.0290161) q[1];
sx q[1];
rz(-1.45603) q[1];
sx q[1];
rz(2.2419825) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3594204) q[0];
sx q[0];
rz(-0.36672584) q[0];
sx q[0];
rz(-0.72215778) q[0];
rz(-pi) q[1];
rz(-0.67018761) q[2];
sx q[2];
rz(-2.0815583) q[2];
sx q[2];
rz(2.0005884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4790597) q[1];
sx q[1];
rz(-0.823349) q[1];
sx q[1];
rz(-1.4940628) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6255304) q[3];
sx q[3];
rz(-2.1806697) q[3];
sx q[3];
rz(-0.42760951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67109913) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(2.1136843) q[2];
rz(0.73728621) q[3];
sx q[3];
rz(-1.4772011) q[3];
sx q[3];
rz(0.024287311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0055493) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(-1.1935724) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(1.6568291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.68736) q[0];
sx q[0];
rz(-0.052536162) q[0];
sx q[0];
rz(-1.0750858) q[0];
x q[1];
rz(1.7822927) q[2];
sx q[2];
rz(-1.4471608) q[2];
sx q[2];
rz(0.57963003) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3419108) q[1];
sx q[1];
rz(-1.7060301) q[1];
sx q[1];
rz(-2.7970201) q[1];
rz(-pi) q[2];
rz(1.7853224) q[3];
sx q[3];
rz(-2.7709922) q[3];
sx q[3];
rz(1.4166946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9803311) q[2];
sx q[2];
rz(-1.223246) q[2];
sx q[2];
rz(0.70183357) q[2];
rz(-0.069132239) q[3];
sx q[3];
rz(-0.88874236) q[3];
sx q[3];
rz(2.4338636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.055701) q[0];
sx q[0];
rz(-1.1393071) q[0];
sx q[0];
rz(-0.62537801) q[0];
rz(-1.0944132) q[1];
sx q[1];
rz(-1.5233327) q[1];
sx q[1];
rz(-1.3166924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863281) q[0];
sx q[0];
rz(-2.1127955) q[0];
sx q[0];
rz(2.9744451) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28202348) q[2];
sx q[2];
rz(-0.33581734) q[2];
sx q[2];
rz(1.6086325) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5881476) q[1];
sx q[1];
rz(-0.7281174) q[1];
sx q[1];
rz(-0.27256984) q[1];
rz(-2.1573477) q[3];
sx q[3];
rz(-0.86185019) q[3];
sx q[3];
rz(1.9714485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7350498) q[2];
sx q[2];
rz(-0.64261618) q[2];
sx q[2];
rz(3.1114846) q[2];
rz(-0.41664577) q[3];
sx q[3];
rz(-1.7130339) q[3];
sx q[3];
rz(-0.055195181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77867126) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(1.9238506) q[0];
rz(-2.9171004) q[1];
sx q[1];
rz(-1.6480564) q[1];
sx q[1];
rz(-2.7405558) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56479708) q[0];
sx q[0];
rz(-2.5677997) q[0];
sx q[0];
rz(2.3032715) q[0];
x q[1];
rz(-0.46868344) q[2];
sx q[2];
rz(-0.71063738) q[2];
sx q[2];
rz(2.8993894) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89304533) q[1];
sx q[1];
rz(-1.5233938) q[1];
sx q[1];
rz(2.9045312) q[1];
rz(-2.5878536) q[3];
sx q[3];
rz(-0.19036346) q[3];
sx q[3];
rz(0.37356627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2991221) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(-0.50745884) q[2];
rz(-0.92042813) q[3];
sx q[3];
rz(-2.7046552) q[3];
sx q[3];
rz(0.18032716) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4949263) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(-1.0194417) q[0];
rz(1.9693718) q[1];
sx q[1];
rz(-1.9688789) q[1];
sx q[1];
rz(1.9196462) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9393765) q[0];
sx q[0];
rz(-0.76986137) q[0];
sx q[0];
rz(-2.0815996) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.035313531) q[2];
sx q[2];
rz(-2.0853634) q[2];
sx q[2];
rz(-3.0701612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1091369) q[1];
sx q[1];
rz(-1.0278406) q[1];
sx q[1];
rz(1.2486267) q[1];
x q[2];
rz(1.1282519) q[3];
sx q[3];
rz(-2.5199515) q[3];
sx q[3];
rz(0.78237247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4814066) q[2];
sx q[2];
rz(-2.47611) q[2];
sx q[2];
rz(2.0484203) q[2];
rz(-0.53705755) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(-2.4721036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.2375803) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(1.681666) q[0];
rz(-1.6319252) q[1];
sx q[1];
rz(-2.5131112) q[1];
sx q[1];
rz(-0.07930886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58188841) q[0];
sx q[0];
rz(-1.4900582) q[0];
sx q[0];
rz(1.8084007) q[0];
x q[1];
rz(-0.58987381) q[2];
sx q[2];
rz(-2.8980719) q[2];
sx q[2];
rz(2.4517454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7056337) q[1];
sx q[1];
rz(-0.50928947) q[1];
sx q[1];
rz(-0.58756709) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15593021) q[3];
sx q[3];
rz(-0.72375384) q[3];
sx q[3];
rz(2.2035905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0038393) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(2.7286781) q[2];
rz(-1.2952992) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(2.7220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15692784) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(2.8562163) q[0];
rz(-0.05263075) q[1];
sx q[1];
rz(-1.4080518) q[1];
sx q[1];
rz(2.9579128) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7010371) q[0];
sx q[0];
rz(-1.4593048) q[0];
sx q[0];
rz(-3.1058342) q[0];
x q[1];
rz(-2.8561228) q[2];
sx q[2];
rz(-1.7186224) q[2];
sx q[2];
rz(-1.2814652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11631509) q[1];
sx q[1];
rz(-0.59981385) q[1];
sx q[1];
rz(0.13529899) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7643572) q[3];
sx q[3];
rz(-0.68359112) q[3];
sx q[3];
rz(2.205276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24667428) q[2];
sx q[2];
rz(-1.4072714) q[2];
sx q[2];
rz(2.0764009) q[2];
rz(-1.6861606) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(2.5559032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.10373779) q[0];
sx q[0];
rz(-0.81473628) q[0];
sx q[0];
rz(1.6023585) q[0];
rz(0.94929758) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(-1.7112214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1174536) q[0];
sx q[0];
rz(-2.040359) q[0];
sx q[0];
rz(-0.59654327) q[0];
rz(-pi) q[1];
rz(0.075242356) q[2];
sx q[2];
rz(-2.2981567) q[2];
sx q[2];
rz(0.6211578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9897223) q[1];
sx q[1];
rz(-2.8443326) q[1];
sx q[1];
rz(-0.10135915) q[1];
rz(-pi) q[2];
rz(1.1422906) q[3];
sx q[3];
rz(-1.3125889) q[3];
sx q[3];
rz(-2.9203897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0091693) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(0.12651786) q[2];
rz(-2.8070519) q[3];
sx q[3];
rz(-2.8821475) q[3];
sx q[3];
rz(-2.2693966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3009406) q[0];
sx q[0];
rz(-0.88467389) q[0];
sx q[0];
rz(0.51710039) q[0];
rz(-0.05489796) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(-0.089769207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9835998) q[0];
sx q[0];
rz(-1.9741892) q[0];
sx q[0];
rz(0.97121111) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0099118) q[2];
sx q[2];
rz(-1.5176306) q[2];
sx q[2];
rz(-1.9694984) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0193737) q[1];
sx q[1];
rz(-1.8226588) q[1];
sx q[1];
rz(1.4068961) q[1];
x q[2];
rz(-0.2281727) q[3];
sx q[3];
rz(-1.2163278) q[3];
sx q[3];
rz(1.4972403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6101997) q[2];
sx q[2];
rz(-2.0072082) q[2];
sx q[2];
rz(-0.72719491) q[2];
rz(-2.3860892) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(0.21272794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94201921) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(1.749281) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(-0.10689312) q[2];
sx q[2];
rz(-0.86915599) q[2];
sx q[2];
rz(-0.55021777) q[2];
rz(-0.82928113) q[3];
sx q[3];
rz(-2.3622475) q[3];
sx q[3];
rz(-2.0105863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
