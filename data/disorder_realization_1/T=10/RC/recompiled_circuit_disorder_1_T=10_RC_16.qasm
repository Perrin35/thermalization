OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(5.1685652) q[0];
sx q[0];
rz(9.4246372) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407335) q[0];
sx q[0];
rz(-1.8574323) q[0];
sx q[0];
rz(2.9564234) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5932181) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(0.89474364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5692917) q[1];
sx q[1];
rz(-2.3094059) q[1];
sx q[1];
rz(-0.65170793) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0317806) q[3];
sx q[3];
rz(-1.7870652) q[3];
sx q[3];
rz(-3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(2.0181296) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9155974) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(1.1094088) q[0];
x q[1];
rz(-0.77983071) q[2];
sx q[2];
rz(-1.6155598) q[2];
sx q[2];
rz(-2.0335399) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1883205) q[1];
sx q[1];
rz(-0.76474944) q[1];
sx q[1];
rz(0.79337593) q[1];
x q[2];
rz(1.0081604) q[3];
sx q[3];
rz(-2.1481272) q[3];
sx q[3];
rz(-1.2853704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(2.222555) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-1.1330053) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5223761) q[0];
sx q[0];
rz(-1.4285061) q[0];
sx q[0];
rz(0.0033358047) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3511231) q[2];
sx q[2];
rz(-1.9445473) q[2];
sx q[2];
rz(-2.8369396) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9085711) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(0.33645333) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77550019) q[3];
sx q[3];
rz(-0.79458487) q[3];
sx q[3];
rz(2.9426108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.8481002) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.9807293) q[0];
rz(0.89598957) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4281222) q[0];
sx q[0];
rz(-2.2566416) q[0];
sx q[0];
rz(-1.1629521) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0669078) q[2];
sx q[2];
rz(-1.1876145) q[2];
sx q[2];
rz(0.89786868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96326522) q[1];
sx q[1];
rz(-1.3849003) q[1];
sx q[1];
rz(-2.9796897) q[1];
rz(-pi) q[2];
rz(2.9837708) q[3];
sx q[3];
rz(-1.2994248) q[3];
sx q[3];
rz(-1.8872758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(-2.2616852) q[2];
rz(-3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039466) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(2.1283545) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(2.0577046) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.097682) q[0];
sx q[0];
rz(-1.327276) q[0];
sx q[0];
rz(-0.18885141) q[0];
rz(-pi) q[1];
rz(-2.8903557) q[2];
sx q[2];
rz(-1.305797) q[2];
sx q[2];
rz(-0.83515177) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.52884) q[1];
sx q[1];
rz(-1.6788947) q[1];
sx q[1];
rz(2.0987064) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5203939) q[3];
sx q[3];
rz(-2.084123) q[3];
sx q[3];
rz(2.7798775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(-0.43236732) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(0.17177467) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6796414) q[0];
sx q[0];
rz(-1.2286751) q[0];
sx q[0];
rz(-1.530594) q[0];
rz(-pi) q[1];
rz(0.32161153) q[2];
sx q[2];
rz(-2.4827637) q[2];
sx q[2];
rz(-2.4668601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9858866) q[1];
sx q[1];
rz(-0.17427467) q[1];
sx q[1];
rz(-1.8272912) q[1];
x q[2];
rz(-1.5395245) q[3];
sx q[3];
rz(-0.074723738) q[3];
sx q[3];
rz(-2.5089695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(-2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.2566459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6709267) q[0];
sx q[0];
rz(-1.6111272) q[0];
sx q[0];
rz(-1.3192024) q[0];
rz(0.2101375) q[2];
sx q[2];
rz(-2.0268831) q[2];
sx q[2];
rz(1.2880585) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9659121) q[1];
sx q[1];
rz(-1.6351846) q[1];
sx q[1];
rz(0.28190159) q[1];
x q[2];
rz(2.4323746) q[3];
sx q[3];
rz(-1.3243444) q[3];
sx q[3];
rz(0.49530464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(-1.5301269) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664739) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(-0.30817729) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(-2.7546308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9403801) q[0];
sx q[0];
rz(-1.3450755) q[0];
sx q[0];
rz(0.97198995) q[0];
rz(0.87256356) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(0.97908212) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3124233) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(-1.6664684) q[1];
rz(1.79027) q[3];
sx q[3];
rz(-2.05728) q[3];
sx q[3];
rz(0.51318491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(-0.72775841) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4816149) q[0];
sx q[0];
rz(-2.1351295) q[0];
sx q[0];
rz(0.83156395) q[0];
rz(-pi) q[1];
rz(2.6463037) q[2];
sx q[2];
rz(-2.1914748) q[2];
sx q[2];
rz(-2.4923313) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2652492) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(-2.6190119) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4008972) q[3];
sx q[3];
rz(-1.1974088) q[3];
sx q[3];
rz(-1.5965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(1.059277) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(1.4019029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83090529) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(0.12916818) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14342587) q[2];
sx q[2];
rz(-2.9529245) q[2];
sx q[2];
rz(2.8667237) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5857196) q[1];
sx q[1];
rz(-1.8750637) q[1];
sx q[1];
rz(2.2833706) q[1];
rz(1.7781156) q[3];
sx q[3];
rz(-2.1917686) q[3];
sx q[3];
rz(-0.45351115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4828651) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.520291) q[2];
rz(2.5907717) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14810066) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(2.2254754) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-0.79402906) q[2];
sx q[2];
rz(-0.90894884) q[2];
sx q[2];
rz(-0.85469645) q[2];
rz(-1.6021452) q[3];
sx q[3];
rz(-0.51732224) q[3];
sx q[3];
rz(2.8696685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];