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
rz(0.81646252) q[0];
sx q[0];
rz(3.2433885) q[0];
sx q[0];
rz(9.9643702) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(0.53909477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4350454) q[0];
sx q[0];
rz(-0.4420949) q[0];
sx q[0];
rz(-3.0821783) q[0];
rz(-pi) q[1];
rz(-2.7819958) q[2];
sx q[2];
rz(-1.4289549) q[2];
sx q[2];
rz(-1.3027018) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7254759) q[1];
sx q[1];
rz(-1.892964) q[1];
sx q[1];
rz(0.14741082) q[1];
x q[2];
rz(-1.4629389) q[3];
sx q[3];
rz(-0.65179208) q[3];
sx q[3];
rz(-2.6848313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76217905) q[2];
sx q[2];
rz(-2.2654686) q[2];
sx q[2];
rz(-1.1890821) q[2];
rz(1.9573697) q[3];
sx q[3];
rz(-2.2370179) q[3];
sx q[3];
rz(1.5949465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9963843) q[0];
sx q[0];
rz(-1.5518016) q[0];
sx q[0];
rz(1.3231963) q[0];
rz(0.48201758) q[1];
sx q[1];
rz(-2.2299485) q[1];
sx q[1];
rz(0.97420305) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27069399) q[0];
sx q[0];
rz(-2.2091731) q[0];
sx q[0];
rz(-2.9186072) q[0];
rz(0.09928273) q[2];
sx q[2];
rz(-1.6834604) q[2];
sx q[2];
rz(0.4140062) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6398826) q[1];
sx q[1];
rz(-2.5474265) q[1];
sx q[1];
rz(0.42826786) q[1];
rz(-pi) q[2];
rz(2.5752546) q[3];
sx q[3];
rz(-1.7741331) q[3];
sx q[3];
rz(-1.4161863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1366068) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(-0.74964398) q[2];
rz(-0.43198112) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(1.3031134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62464803) q[0];
sx q[0];
rz(-2.8693146) q[0];
sx q[0];
rz(2.5352449) q[0];
rz(-1.8079405) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(2.8820754) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4496499) q[0];
sx q[0];
rz(-0.29352531) q[0];
sx q[0];
rz(-2.13563) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9664842) q[2];
sx q[2];
rz(-1.3710183) q[2];
sx q[2];
rz(-1.7410884) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9898228) q[1];
sx q[1];
rz(-2.5603316) q[1];
sx q[1];
rz(0.84110028) q[1];
rz(-pi) q[2];
rz(-1.0084148) q[3];
sx q[3];
rz(-1.9451688) q[3];
sx q[3];
rz(1.238747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2567265) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(1.5474896) q[2];
rz(2.3571842) q[3];
sx q[3];
rz(-1.514785) q[3];
sx q[3];
rz(1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6467658) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(-0.25949091) q[0];
rz(-1.9208113) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(-1.6993914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4239765) q[0];
sx q[0];
rz(-2.0956796) q[0];
sx q[0];
rz(0.54846008) q[0];
x q[1];
rz(2.1659508) q[2];
sx q[2];
rz(-1.6791428) q[2];
sx q[2];
rz(2.8841022) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7573141) q[1];
sx q[1];
rz(-2.7065249) q[1];
sx q[1];
rz(2.5237094) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.070017858) q[3];
sx q[3];
rz(-0.88980955) q[3];
sx q[3];
rz(-0.79047608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0356902) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(-0.40994677) q[2];
rz(-1.2116872) q[3];
sx q[3];
rz(-0.68710059) q[3];
sx q[3];
rz(1.80779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425151) q[0];
sx q[0];
rz(-2.4254159) q[0];
sx q[0];
rz(-2.8675365) q[0];
rz(0.43824276) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(2.4058707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5925795) q[0];
sx q[0];
rz(-1.6578478) q[0];
sx q[0];
rz(0.75496952) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9922036) q[2];
sx q[2];
rz(-1.427703) q[2];
sx q[2];
rz(-0.18994513) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.055701697) q[1];
sx q[1];
rz(-0.38905479) q[1];
sx q[1];
rz(0.86430734) q[1];
rz(-0.92145958) q[3];
sx q[3];
rz(-1.9917352) q[3];
sx q[3];
rz(-0.73900797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2030486) q[2];
sx q[2];
rz(-1.4578578) q[2];
sx q[2];
rz(-0.61384821) q[2];
rz(-2.7326873) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(-0.74234211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3530389) q[0];
sx q[0];
rz(-0.63816324) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(1.6963814) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(0.17527418) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69578457) q[0];
sx q[0];
rz(-0.17477594) q[0];
sx q[0];
rz(1.4265027) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66368033) q[2];
sx q[2];
rz(-2.366407) q[2];
sx q[2];
rz(-0.41942393) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8126909) q[1];
sx q[1];
rz(-1.1393424) q[1];
sx q[1];
rz(-0.56568362) q[1];
rz(-pi) q[2];
rz(-2.5644184) q[3];
sx q[3];
rz(-2.1041311) q[3];
sx q[3];
rz(-1.8624901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96482977) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(2.1691587) q[2];
rz(-1.6644299) q[3];
sx q[3];
rz(-1.9921314) q[3];
sx q[3];
rz(0.76178637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85457388) q[0];
sx q[0];
rz(-0.022495689) q[0];
sx q[0];
rz(-0.35312411) q[0];
rz(-2.1137721) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(-2.1305398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38336223) q[0];
sx q[0];
rz(-2.7840021) q[0];
sx q[0];
rz(-0.041186853) q[0];
rz(-pi) q[1];
rz(-0.86507823) q[2];
sx q[2];
rz(-2.0524244) q[2];
sx q[2];
rz(-0.31246802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87349975) q[1];
sx q[1];
rz(-0.38485369) q[1];
sx q[1];
rz(1.9995688) q[1];
rz(-1.8854463) q[3];
sx q[3];
rz(-1.8487747) q[3];
sx q[3];
rz(1.4643471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8280243) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(-1.774452) q[2];
rz(-1.2683055) q[3];
sx q[3];
rz(-1.4788078) q[3];
sx q[3];
rz(-0.84793276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0276412) q[0];
sx q[0];
rz(-0.60751644) q[0];
sx q[0];
rz(-1.129958) q[0];
rz(-0.21895151) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(-0.91046441) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509678) q[0];
sx q[0];
rz(-0.70075894) q[0];
sx q[0];
rz(2.6118082) q[0];
x q[1];
rz(2.6633419) q[2];
sx q[2];
rz(-1.646981) q[2];
sx q[2];
rz(0.26631276) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7607019) q[1];
sx q[1];
rz(-0.58769757) q[1];
sx q[1];
rz(-2.0668849) q[1];
rz(2.4164356) q[3];
sx q[3];
rz(-1.8817543) q[3];
sx q[3];
rz(-0.5288611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8255446) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(-0.46142203) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(2.6958444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6398741) q[0];
sx q[0];
rz(-1.8032782) q[0];
sx q[0];
rz(-2.3097532) q[0];
rz(-2.4141451) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(0.096253455) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.953234) q[0];
sx q[0];
rz(-1.112845) q[0];
sx q[0];
rz(-1.0012549) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8873994) q[2];
sx q[2];
rz(-1.6302178) q[2];
sx q[2];
rz(0.78320247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2029496) q[1];
sx q[1];
rz(-1.4803107) q[1];
sx q[1];
rz(1.7825148) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7670125) q[3];
sx q[3];
rz(-1.726578) q[3];
sx q[3];
rz(0.84583144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7353797) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(-1.5489138) q[2];
rz(-2.5088572) q[3];
sx q[3];
rz(-1.2318719) q[3];
sx q[3];
rz(-0.24937853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525472) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(0.35300514) q[0];
rz(-0.32265916) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(2.7696612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73814017) q[0];
sx q[0];
rz(-2.7785025) q[0];
sx q[0];
rz(-2.1573261) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.811889) q[2];
sx q[2];
rz(-1.7404544) q[2];
sx q[2];
rz(-0.28126954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0363716) q[1];
sx q[1];
rz(-2.7694355) q[1];
sx q[1];
rz(0.45366617) q[1];
rz(-pi) q[2];
rz(3.102469) q[3];
sx q[3];
rz(-2.2949766) q[3];
sx q[3];
rz(-0.20871221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17452621) q[2];
sx q[2];
rz(-1.9517978) q[2];
sx q[2];
rz(1.8444427) q[2];
rz(2.2938812) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90947718) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(-1.3879981) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(-2.3222011) q[2];
sx q[2];
rz(-2.1234305) q[2];
sx q[2];
rz(2.7101868) q[2];
rz(2.7305215) q[3];
sx q[3];
rz(-0.84079327) q[3];
sx q[3];
rz(-0.83415915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
