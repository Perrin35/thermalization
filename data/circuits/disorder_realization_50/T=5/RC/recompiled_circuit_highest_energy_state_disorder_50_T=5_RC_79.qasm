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
rz(-2.3565489) q[1];
sx q[1];
rz(-2.4554689) q[1];
sx q[1];
rz(0.52678984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9357932) q[0];
sx q[0];
rz(-2.4895634) q[0];
sx q[0];
rz(-1.0656641) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85028591) q[2];
sx q[2];
rz(-0.97480259) q[2];
sx q[2];
rz(1.13499) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96054582) q[1];
sx q[1];
rz(-1.5343473) q[1];
sx q[1];
rz(-2.6992082) q[1];
rz(-pi) q[2];
rz(-1.6185226) q[3];
sx q[3];
rz(-0.56079799) q[3];
sx q[3];
rz(0.58855295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.779125) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(-2.989952) q[2];
rz(2.9996297) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-1.1375455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4769984) q[0];
sx q[0];
rz(-0.73200309) q[0];
sx q[0];
rz(-2.3240996) q[0];
rz(0.11257653) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(0.89961019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6647659) q[0];
sx q[0];
rz(-1.3315086) q[0];
sx q[0];
rz(2.8609758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1915273) q[2];
sx q[2];
rz(-2.143444) q[2];
sx q[2];
rz(-3.0814296) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1020722) q[1];
sx q[1];
rz(-1.6270492) q[1];
sx q[1];
rz(0.7489167) q[1];
rz(-pi) q[2];
rz(0.95590653) q[3];
sx q[3];
rz(-2.3645176) q[3];
sx q[3];
rz(1.2082548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4704935) q[2];
sx q[2];
rz(-1.6464536) q[2];
sx q[2];
rz(1.0279083) q[2];
rz(-0.73728621) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(0.024287311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1360433) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(-1.9480202) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(-1.4847635) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5201621) q[0];
sx q[0];
rz(-1.545816) q[0];
sx q[0];
rz(-1.6170184) q[0];
x q[1];
rz(-1.7822927) q[2];
sx q[2];
rz(-1.4471608) q[2];
sx q[2];
rz(-0.57963003) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4188267) q[1];
sx q[1];
rz(-1.2294985) q[1];
sx q[1];
rz(-1.7143634) q[1];
x q[2];
rz(-1.2079575) q[3];
sx q[3];
rz(-1.4936183) q[3];
sx q[3];
rz(3.0953323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16126157) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(0.70183357) q[2];
rz(-3.0724604) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(-0.70772901) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0858916) q[0];
sx q[0];
rz(-1.1393071) q[0];
sx q[0];
rz(-0.62537801) q[0];
rz(2.0471795) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(1.3166924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33915658) q[0];
sx q[0];
rz(-0.56471741) q[0];
sx q[0];
rz(-1.3012278) q[0];
rz(-0.32347347) q[2];
sx q[2];
rz(-1.6626366) q[2];
sx q[2];
rz(-0.30486456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91176468) q[1];
sx q[1];
rz(-0.87512866) q[1];
sx q[1];
rz(-1.8063481) q[1];
x q[2];
rz(-0.80009256) q[3];
sx q[3];
rz(-2.0043819) q[3];
sx q[3];
rz(-2.3326002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40654287) q[2];
sx q[2];
rz(-2.4989765) q[2];
sx q[2];
rz(0.030108062) q[2];
rz(-0.41664577) q[3];
sx q[3];
rz(-1.4285587) q[3];
sx q[3];
rz(0.055195181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3629214) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(-1.9238506) q[0];
rz(-2.9171004) q[1];
sx q[1];
rz(-1.4935363) q[1];
sx q[1];
rz(-0.40103689) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3844135) q[0];
sx q[0];
rz(-1.9862439) q[0];
sx q[0];
rz(2.7336304) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.941576) q[2];
sx q[2];
rz(-2.1919498) q[2];
sx q[2];
rz(0.3467243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2485473) q[1];
sx q[1];
rz(-1.6181989) q[1];
sx q[1];
rz(2.9045312) q[1];
x q[2];
rz(-2.5878536) q[3];
sx q[3];
rz(-0.19036346) q[3];
sx q[3];
rz(0.37356627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.84247056) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(0.50745884) q[2];
rz(2.2211645) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(-0.18032716) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4949263) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(1.0194417) q[0];
rz(1.9693718) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(1.2219465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86496124) q[0];
sx q[0];
rz(-2.2233133) q[0];
sx q[0];
rz(2.6990273) q[0];
x q[1];
rz(1.6331633) q[2];
sx q[2];
rz(-2.6259239) q[2];
sx q[2];
rz(3.1413648) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7740357) q[1];
sx q[1];
rz(-1.8453215) q[1];
sx q[1];
rz(-2.5749899) q[1];
x q[2];
rz(-1.1282519) q[3];
sx q[3];
rz(-0.62164111) q[3];
sx q[3];
rz(0.78237247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4814066) q[2];
sx q[2];
rz(-2.47611) q[2];
sx q[2];
rz(2.0484203) q[2];
rz(-0.53705755) q[3];
sx q[3];
rz(-1.4635181) q[3];
sx q[3];
rz(2.4721036) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2375803) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(1.4599266) q[0];
rz(-1.6319252) q[1];
sx q[1];
rz(-0.62848148) q[1];
sx q[1];
rz(0.07930886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722141) q[0];
sx q[0];
rz(-1.3339808) q[0];
sx q[0];
rz(-0.083061465) q[0];
x q[1];
rz(-0.20360501) q[2];
sx q[2];
rz(-1.4362659) q[2];
sx q[2];
rz(2.8367538) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.053959) q[1];
sx q[1];
rz(-1.988639) q[1];
sx q[1];
rz(-1.8709976) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7071869) q[3];
sx q[3];
rz(-0.85771424) q[3];
sx q[3];
rz(-2.4103885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1377533) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(-0.4129146) q[2];
rz(-1.8462935) q[3];
sx q[3];
rz(-0.92419878) q[3];
sx q[3];
rz(-0.41954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9846648) q[0];
sx q[0];
rz(-0.65021896) q[0];
sx q[0];
rz(2.8562163) q[0];
rz(-0.05263075) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(-2.9579128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7010371) q[0];
sx q[0];
rz(-1.4593048) q[0];
sx q[0];
rz(-0.035758408) q[0];
rz(-pi) q[1];
rz(1.4168315) q[2];
sx q[2];
rz(-1.2885258) q[2];
sx q[2];
rz(2.8090614) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27976263) q[1];
sx q[1];
rz(-0.97721902) q[1];
sx q[1];
rz(-1.4788126) q[1];
rz(0.64818188) q[3];
sx q[3];
rz(-1.3360008) q[3];
sx q[3];
rz(0.33644331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8949184) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(-2.0764009) q[2];
rz(1.4554321) q[3];
sx q[3];
rz(-1.8543517) q[3];
sx q[3];
rz(-2.5559032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.10373779) q[0];
sx q[0];
rz(-0.81473628) q[0];
sx q[0];
rz(-1.6023585) q[0];
rz(-0.94929758) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(-1.4303713) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1827138) q[0];
sx q[0];
rz(-2.4005167) q[0];
sx q[0];
rz(2.4070021) q[0];
rz(-pi) q[1];
rz(1.4865506) q[2];
sx q[2];
rz(-0.73053321) q[2];
sx q[2];
rz(2.4075395) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0956968) q[1];
sx q[1];
rz(-1.8664845) q[1];
sx q[1];
rz(1.6017833) q[1];
rz(-1.1422906) q[3];
sx q[3];
rz(-1.3125889) q[3];
sx q[3];
rz(-0.22120295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1324233) q[2];
sx q[2];
rz(-1.8370266) q[2];
sx q[2];
rz(-0.12651786) q[2];
rz(-0.33454076) q[3];
sx q[3];
rz(-2.8821475) q[3];
sx q[3];
rz(-0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3009406) q[0];
sx q[0];
rz(-0.88467389) q[0];
sx q[0];
rz(2.6244923) q[0];
rz(-3.0866947) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(-3.0518234) q[1];
x q[2];
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
rz(0.062762063) q[2];
sx q[2];
rz(-1.0107991) q[2];
sx q[2];
rz(0.43207174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.59261428) q[1];
sx q[1];
rz(-1.7294809) q[1];
sx q[1];
rz(-2.8864546) q[1];
rz(-pi) q[2];
rz(0.2281727) q[3];
sx q[3];
rz(-1.2163278) q[3];
sx q[3];
rz(-1.4972403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6101997) q[2];
sx q[2];
rz(-2.0072082) q[2];
sx q[2];
rz(-2.4143977) q[2];
rz(-0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(-0.21272794) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94201921) q[0];
sx q[0];
rz(-1.6006391) q[0];
sx q[0];
rz(-2.0508456) q[0];
rz(-1.749281) q[1];
sx q[1];
rz(-1.6284457) q[1];
sx q[1];
rz(1.5096691) q[1];
rz(2.2752599) q[2];
sx q[2];
rz(-1.489218) q[2];
sx q[2];
rz(1.0897286) q[2];
rz(2.3123115) q[3];
sx q[3];
rz(-2.3622475) q[3];
sx q[3];
rz(-2.0105863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
