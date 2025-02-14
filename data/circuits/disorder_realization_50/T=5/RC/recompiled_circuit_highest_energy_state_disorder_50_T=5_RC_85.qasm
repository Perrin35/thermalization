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
rz(-2.4432776) q[0];
sx q[0];
rz(-2.7628216) q[0];
sx q[0];
rz(1.9842499) q[0];
rz(0.78504374) q[1];
sx q[1];
rz(-0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049126712) q[0];
sx q[0];
rz(-1.2727588) q[0];
sx q[0];
rz(-0.98182337) q[0];
rz(-pi) q[1];
rz(0.73424642) q[2];
sx q[2];
rz(-2.1484005) q[2];
sx q[2];
rz(2.2478888) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96054582) q[1];
sx q[1];
rz(-1.6072453) q[1];
sx q[1];
rz(2.6992082) q[1];
rz(-2.1310813) q[3];
sx q[3];
rz(-1.5454195) q[3];
sx q[3];
rz(-0.94181864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36246768) q[2];
sx q[2];
rz(-0.96131009) q[2];
sx q[2];
rz(-0.15164068) q[2];
rz(-0.14196299) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-1.1375455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(-0.89961019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6647659) q[0];
sx q[0];
rz(-1.3315086) q[0];
sx q[0];
rz(-0.28061687) q[0];
rz(-pi) q[1];
rz(2.1915273) q[2];
sx q[2];
rz(-0.99814864) q[2];
sx q[2];
rz(3.0814296) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.591685) q[1];
sx q[1];
rz(-2.3909759) q[1];
sx q[1];
rz(3.0590712) q[1];
x q[2];
rz(-2.2475776) q[3];
sx q[3];
rz(-1.9872267) q[3];
sx q[3];
rz(2.3125966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4704935) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(2.1136843) q[2];
rz(2.4043064) q[3];
sx q[3];
rz(-1.4772011) q[3];
sx q[3];
rz(3.1173053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1360433) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(-1.9480202) q[0];
rz(-1.2981124) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(-1.4847635) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95052108) q[0];
sx q[0];
rz(-1.6170039) q[0];
sx q[0];
rz(3.1165857) q[0];
x q[1];
rz(0.12642352) q[2];
sx q[2];
rz(-1.3609386) q[2];
sx q[2];
rz(1.017638) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3419108) q[1];
sx q[1];
rz(-1.7060301) q[1];
sx q[1];
rz(2.7970201) q[1];
rz(1.3562702) q[3];
sx q[3];
rz(-2.7709922) q[3];
sx q[3];
rz(-1.4166946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16126157) q[2];
sx q[2];
rz(-1.223246) q[2];
sx q[2];
rz(-2.4397591) q[2];
rz(0.069132239) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(2.4338636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-2.055701) q[0];
sx q[0];
rz(-1.1393071) q[0];
sx q[0];
rz(-0.62537801) q[0];
rz(1.0944132) q[1];
sx q[1];
rz(-1.5233327) q[1];
sx q[1];
rz(1.3166924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024361) q[0];
sx q[0];
rz(-2.5768752) q[0];
sx q[0];
rz(-1.8403649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32347347) q[2];
sx q[2];
rz(-1.4789561) q[2];
sx q[2];
rz(0.30486456) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91176468) q[1];
sx q[1];
rz(-0.87512866) q[1];
sx q[1];
rz(1.3352446) q[1];
rz(0.98424498) q[3];
sx q[3];
rz(-0.86185019) q[3];
sx q[3];
rz(-1.1701442) q[3];
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
rz(-1.4285587) q[3];
sx q[3];
rz(-3.0863975) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3629214) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(1.2177421) q[0];
rz(0.22449224) q[1];
sx q[1];
rz(-1.6480564) q[1];
sx q[1];
rz(0.40103689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5767956) q[0];
sx q[0];
rz(-2.5677997) q[0];
sx q[0];
rz(-2.3032715) q[0];
rz(-pi) q[1];
rz(2.6729092) q[2];
sx q[2];
rz(-0.71063738) q[2];
sx q[2];
rz(-0.24220322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89304533) q[1];
sx q[1];
rz(-1.6181989) q[1];
sx q[1];
rz(-0.23706146) q[1];
rz(1.4698074) q[3];
sx q[3];
rz(-1.4091531) q[3];
sx q[3];
rz(-0.18835959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84247056) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(2.6341338) q[2];
rz(-0.92042813) q[3];
sx q[3];
rz(-2.7046552) q[3];
sx q[3];
rz(0.18032716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64666635) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(1.0194417) q[0];
rz(1.1722209) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(-1.2219465) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98602277) q[0];
sx q[0];
rz(-1.2235852) q[0];
sx q[0];
rz(-0.86882868) q[0];
x q[1];
rz(2.0856306) q[2];
sx q[2];
rz(-1.6015354) q[2];
sx q[2];
rz(1.6248425) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5355365) q[1];
sx q[1];
rz(-2.5185985) q[1];
sx q[1];
rz(-0.48320233) q[1];
rz(1.1282519) q[3];
sx q[3];
rz(-2.5199515) q[3];
sx q[3];
rz(0.78237247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6601861) q[2];
sx q[2];
rz(-2.47611) q[2];
sx q[2];
rz(-2.0484203) q[2];
rz(2.6045351) q[3];
sx q[3];
rz(-1.4635181) q[3];
sx q[3];
rz(2.4721036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2375803) q[0];
sx q[0];
rz(-2.8222988) q[0];
sx q[0];
rz(1.681666) q[0];
rz(1.6319252) q[1];
sx q[1];
rz(-0.62848148) q[1];
sx q[1];
rz(-0.07930886) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96937856) q[0];
sx q[0];
rz(-1.3339808) q[0];
sx q[0];
rz(-0.083061465) q[0];
rz(-pi) q[1];
rz(-0.20360501) q[2];
sx q[2];
rz(-1.7053268) q[2];
sx q[2];
rz(0.30483887) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0876336) q[1];
sx q[1];
rz(-1.1529536) q[1];
sx q[1];
rz(1.2705951) q[1];
rz(2.9856624) q[3];
sx q[3];
rz(-2.4178388) q[3];
sx q[3];
rz(-2.2035905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0038393) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(-0.4129146) q[2];
rz(-1.8462935) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(0.41954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15692784) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(-2.8562163) q[0];
rz(-0.05263075) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(-2.9579128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3899468) q[0];
sx q[0];
rz(-0.11706287) q[0];
sx q[0];
rz(1.2617017) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48641522) q[2];
sx q[2];
rz(-2.8210495) q[2];
sx q[2];
rz(2.965791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27976263) q[1];
sx q[1];
rz(-0.97721902) q[1];
sx q[1];
rz(-1.4788126) q[1];
x q[2];
rz(2.4934108) q[3];
sx q[3];
rz(-1.8055918) q[3];
sx q[3];
rz(0.33644331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24667428) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(1.0651917) q[2];
rz(-1.6861606) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(2.5559032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10373779) q[0];
sx q[0];
rz(-0.81473628) q[0];
sx q[0];
rz(1.5392342) q[0];
rz(2.1922951) q[1];
sx q[1];
rz(-2.0930591) q[1];
sx q[1];
rz(1.4303713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8447782) q[0];
sx q[0];
rz(-2.0955293) q[0];
sx q[0];
rz(-2.120976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.075242356) q[2];
sx q[2];
rz(-0.84343592) q[2];
sx q[2];
rz(0.6211578) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6257244) q[1];
sx q[1];
rz(-1.5411545) q[1];
sx q[1];
rz(-0.29582204) q[1];
x q[2];
rz(-1.999302) q[3];
sx q[3];
rz(-1.3125889) q[3];
sx q[3];
rz(0.22120295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1324233) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(-0.12651786) q[2];
rz(-2.8070519) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(-0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3009406) q[0];
sx q[0];
rz(-2.2569188) q[0];
sx q[0];
rz(-0.51710039) q[0];
rz(-0.05489796) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(3.0518234) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9341016) q[0];
sx q[0];
rz(-2.4330407) q[0];
sx q[0];
rz(-0.92325489) q[0];
rz(-pi) q[1];
rz(-1.6705065) q[2];
sx q[2];
rz(-2.5784628) q[2];
sx q[2];
rz(2.827284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5489784) q[1];
sx q[1];
rz(-1.4121118) q[1];
sx q[1];
rz(0.25513809) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9338984) q[3];
sx q[3];
rz(-1.3570367) q[3];
sx q[3];
rz(0.0068706415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53139293) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(0.72719491) q[2];
rz(0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(0.21272794) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94201921) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(-1.3923116) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(-0.86633273) q[2];
sx q[2];
rz(-1.489218) q[2];
sx q[2];
rz(1.0897286) q[2];
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
