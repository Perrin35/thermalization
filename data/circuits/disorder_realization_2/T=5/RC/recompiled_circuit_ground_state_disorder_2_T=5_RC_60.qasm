OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6096697) q[0];
sx q[0];
rz(-1.6165531) q[0];
sx q[0];
rz(-1.2582231) q[0];
rz(-1.582107) q[1];
sx q[1];
rz(-2.6546302) q[1];
sx q[1];
rz(2.7064986) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4420051) q[0];
sx q[0];
rz(-2.3908092) q[0];
sx q[0];
rz(0.72608279) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54122178) q[2];
sx q[2];
rz(-0.60784423) q[2];
sx q[2];
rz(1.4262878) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8074602) q[1];
sx q[1];
rz(-1.8402446) q[1];
sx q[1];
rz(-0.079748591) q[1];
rz(-pi) q[2];
rz(1.0886816) q[3];
sx q[3];
rz(-1.0728683) q[3];
sx q[3];
rz(-0.21592797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9293999) q[2];
sx q[2];
rz(-1.8636899) q[2];
sx q[2];
rz(-1.0431935) q[2];
rz(-2.7405401) q[3];
sx q[3];
rz(-1.4570823) q[3];
sx q[3];
rz(-2.5636087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9143518) q[0];
sx q[0];
rz(-0.14509097) q[0];
sx q[0];
rz(-2.5703854) q[0];
rz(0.97194833) q[1];
sx q[1];
rz(-0.74058878) q[1];
sx q[1];
rz(1.1245701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1069431) q[0];
sx q[0];
rz(-0.9616344) q[0];
sx q[0];
rz(0.27584313) q[0];
rz(1.9306514) q[2];
sx q[2];
rz(-0.95155638) q[2];
sx q[2];
rz(-0.93114668) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.827885) q[1];
sx q[1];
rz(-1.6646302) q[1];
sx q[1];
rz(0.83637823) q[1];
x q[2];
rz(0.092512802) q[3];
sx q[3];
rz(-0.20167758) q[3];
sx q[3];
rz(2.8254333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9198415) q[2];
sx q[2];
rz(-1.5566748) q[2];
sx q[2];
rz(0.43608967) q[2];
rz(2.3927205) q[3];
sx q[3];
rz(-2.336453) q[3];
sx q[3];
rz(0.82912904) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0332758) q[0];
sx q[0];
rz(-1.746614) q[0];
sx q[0];
rz(-3.0535611) q[0];
rz(-2.2875359) q[1];
sx q[1];
rz(-2.4805534) q[1];
sx q[1];
rz(1.0850272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4861761) q[0];
sx q[0];
rz(-1.2987505) q[0];
sx q[0];
rz(3.0718417) q[0];
x q[1];
rz(-1.3543812) q[2];
sx q[2];
rz(-1.4206352) q[2];
sx q[2];
rz(-0.79126287) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4702556) q[1];
sx q[1];
rz(-0.47154676) q[1];
sx q[1];
rz(1.9370473) q[1];
x q[2];
rz(-0.38934598) q[3];
sx q[3];
rz(-1.5028477) q[3];
sx q[3];
rz(0.24711497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4714841) q[2];
sx q[2];
rz(-1.7210311) q[2];
sx q[2];
rz(-0.74618435) q[2];
rz(2.9255195) q[3];
sx q[3];
rz(-1.2553299) q[3];
sx q[3];
rz(-2.9358673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7699319) q[0];
sx q[0];
rz(-0.07337229) q[0];
sx q[0];
rz(-1.7727456) q[0];
rz(-2.5313077) q[1];
sx q[1];
rz(-1.7758324) q[1];
sx q[1];
rz(-2.761421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4044426) q[0];
sx q[0];
rz(-1.2594873) q[0];
sx q[0];
rz(-0.4510757) q[0];
rz(-0.45079239) q[2];
sx q[2];
rz(-1.3935699) q[2];
sx q[2];
rz(-1.705436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.026555) q[1];
sx q[1];
rz(-1.7174182) q[1];
sx q[1];
rz(2.0762334) q[1];
x q[2];
rz(1.1690456) q[3];
sx q[3];
rz(-1.6309079) q[3];
sx q[3];
rz(-1.3952122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88468203) q[2];
sx q[2];
rz(-0.94748777) q[2];
sx q[2];
rz(-1.0379637) q[2];
rz(2.7773618) q[3];
sx q[3];
rz(-0.7887775) q[3];
sx q[3];
rz(-0.68575931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6895741) q[0];
sx q[0];
rz(-1.7114534) q[0];
sx q[0];
rz(-0.83056393) q[0];
rz(-0.05016249) q[1];
sx q[1];
rz(-0.12718931) q[1];
sx q[1];
rz(-0.62279472) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1019331) q[0];
sx q[0];
rz(-2.9103511) q[0];
sx q[0];
rz(-0.34666611) q[0];
rz(0.41983126) q[2];
sx q[2];
rz(-1.3908885) q[2];
sx q[2];
rz(-0.77984389) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.924979) q[1];
sx q[1];
rz(-0.93624026) q[1];
sx q[1];
rz(0.44559176) q[1];
rz(-pi) q[2];
rz(-0.0029124459) q[3];
sx q[3];
rz(-3.1191428) q[3];
sx q[3];
rz(-2.8685304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.043776) q[2];
sx q[2];
rz(-2.706683) q[2];
sx q[2];
rz(2.3338976) q[2];
rz(-1.5378599) q[3];
sx q[3];
rz(-1.5065498) q[3];
sx q[3];
rz(-0.22411552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1767126) q[0];
sx q[0];
rz(-1.3310615) q[0];
sx q[0];
rz(1.4759195) q[0];
rz(1.2337947) q[1];
sx q[1];
rz(-1.7314311) q[1];
sx q[1];
rz(2.9880611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3672402) q[0];
sx q[0];
rz(-1.4611547) q[0];
sx q[0];
rz(1.5745274) q[0];
rz(2.6903908) q[2];
sx q[2];
rz(-1.1179194) q[2];
sx q[2];
rz(1.7141327) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0316089) q[1];
sx q[1];
rz(-2.4436207) q[1];
sx q[1];
rz(-0.72781659) q[1];
x q[2];
rz(0.63871049) q[3];
sx q[3];
rz(-1.6473466) q[3];
sx q[3];
rz(-0.18639263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70225707) q[2];
sx q[2];
rz(-2.3372237) q[2];
sx q[2];
rz(-0.55629936) q[2];
rz(-3.124681) q[3];
sx q[3];
rz(-2.1664679) q[3];
sx q[3];
rz(-2.4771966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21338129) q[0];
sx q[0];
rz(-3.0688372) q[0];
sx q[0];
rz(2.4640006) q[0];
rz(-1.821359) q[1];
sx q[1];
rz(-1.0878539) q[1];
sx q[1];
rz(0.56603146) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4844263) q[0];
sx q[0];
rz(-0.56492794) q[0];
sx q[0];
rz(0.66175263) q[0];
rz(-1.5185616) q[2];
sx q[2];
rz(-1.5855316) q[2];
sx q[2];
rz(0.39817223) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2471602) q[1];
sx q[1];
rz(-1.4196779) q[1];
sx q[1];
rz(-0.86455785) q[1];
rz(-pi) q[2];
rz(-1.5291847) q[3];
sx q[3];
rz(-2.0609612) q[3];
sx q[3];
rz(3.0050638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2322106) q[2];
sx q[2];
rz(-0.38572329) q[2];
sx q[2];
rz(0.56987008) q[2];
rz(1.0424403) q[3];
sx q[3];
rz(-1.180155) q[3];
sx q[3];
rz(-1.0038143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7637699) q[0];
sx q[0];
rz(-2.7144987) q[0];
sx q[0];
rz(-2.1183993) q[0];
rz(2.2238253) q[1];
sx q[1];
rz(-2.387391) q[1];
sx q[1];
rz(-0.86407152) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3765434) q[0];
sx q[0];
rz(-1.6934388) q[0];
sx q[0];
rz(-1.3255694) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1064498) q[2];
sx q[2];
rz(-1.9111782) q[2];
sx q[2];
rz(0.92591531) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1004677) q[1];
sx q[1];
rz(-1.0881313) q[1];
sx q[1];
rz(2.3339999) q[1];
rz(-pi) q[2];
rz(0.72743639) q[3];
sx q[3];
rz(-1.2890649) q[3];
sx q[3];
rz(1.4294595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0495905) q[2];
sx q[2];
rz(-1.8747753) q[2];
sx q[2];
rz(-0.76466307) q[2];
rz(-0.90096724) q[3];
sx q[3];
rz(-0.67470208) q[3];
sx q[3];
rz(-2.0425792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.51729274) q[0];
sx q[0];
rz(-2.7462672) q[0];
sx q[0];
rz(2.1584391) q[0];
rz(-0.82955018) q[1];
sx q[1];
rz(-1.7770276) q[1];
sx q[1];
rz(-0.57873908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0163703) q[0];
sx q[0];
rz(-1.1770605) q[0];
sx q[0];
rz(3.1278174) q[0];
rz(-2.5972967) q[2];
sx q[2];
rz(-1.0981264) q[2];
sx q[2];
rz(-0.67677697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1956801) q[1];
sx q[1];
rz(-2.0766313) q[1];
sx q[1];
rz(2.6786184) q[1];
rz(0.10161256) q[3];
sx q[3];
rz(-2.7946473) q[3];
sx q[3];
rz(-2.8941151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4911554) q[2];
sx q[2];
rz(-1.2885685) q[2];
sx q[2];
rz(-0.11162652) q[2];
rz(-2.1857183) q[3];
sx q[3];
rz(-2.1501232) q[3];
sx q[3];
rz(1.8808232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2535506) q[0];
sx q[0];
rz(-2.084806) q[0];
sx q[0];
rz(-3.1274617) q[0];
rz(-1.1125394) q[1];
sx q[1];
rz(-2.2281149) q[1];
sx q[1];
rz(1.5412451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3378572) q[0];
sx q[0];
rz(-1.7183665) q[0];
sx q[0];
rz(-1.7002559) q[0];
rz(-pi) q[1];
rz(-2.2152973) q[2];
sx q[2];
rz(-1.4699114) q[2];
sx q[2];
rz(-2.6580194) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4804876) q[1];
sx q[1];
rz(-0.36927642) q[1];
sx q[1];
rz(-0.1657471) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9716206) q[3];
sx q[3];
rz(-1.408443) q[3];
sx q[3];
rz(1.270783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9666226) q[2];
sx q[2];
rz(-1.9113767) q[2];
sx q[2];
rz(-0.54565412) q[2];
rz(-3.0555861) q[3];
sx q[3];
rz(-1.6274446) q[3];
sx q[3];
rz(0.29451323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1957112) q[0];
sx q[0];
rz(-1.042689) q[0];
sx q[0];
rz(-1.9324017) q[0];
rz(-2.53881) q[1];
sx q[1];
rz(-2.5310015) q[1];
sx q[1];
rz(-2.3878154) q[1];
rz(-1.3134342) q[2];
sx q[2];
rz(-2.4429697) q[2];
sx q[2];
rz(2.2945678) q[2];
rz(2.9908441) q[3];
sx q[3];
rz(-2.5133532) q[3];
sx q[3];
rz(-2.8379074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
