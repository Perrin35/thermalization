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
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(-1.1934086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0573187) q[0];
sx q[0];
rz(-2.8017375) q[0];
sx q[0];
rz(-2.1291332) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5483746) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(2.246849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27543435) q[1];
sx q[1];
rz(-2.1992116) q[1];
sx q[1];
rz(-0.98316146) q[1];
rz(3.0317806) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(-0.11085489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(1.9127282) q[2];
rz(-1.7284283) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(-0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-2.0181296) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9155974) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(1.1094088) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0779755) q[2];
sx q[2];
rz(-0.78084313) q[2];
sx q[2];
rz(-2.7240679) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1264798) q[1];
sx q[1];
rz(-1.0547332) q[1];
sx q[1];
rz(-2.5491787) q[1];
x q[2];
rz(-0.68617679) q[3];
sx q[3];
rz(-0.78305972) q[3];
sx q[3];
rz(-2.1427597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79364395) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(-2.222555) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640901) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-1.1330053) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64273524) q[0];
sx q[0];
rz(-0.14232902) q[0];
sx q[0];
rz(-1.5940773) q[0];
rz(-pi) q[1];
rz(2.6373367) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(2.2222663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7743467) q[1];
sx q[1];
rz(-2.5173752) q[1];
sx q[1];
rz(2.0778076) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5127701) q[3];
sx q[3];
rz(-2.0938794) q[3];
sx q[3];
rz(2.3716225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(1.8815276) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(3.0060351) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82942688) q[0];
sx q[0];
rz(-0.78071785) q[0];
sx q[0];
rz(-2.6902945) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1866456) q[2];
sx q[2];
rz(-1.5015366) q[2];
sx q[2];
rz(-0.70089507) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.902066) q[1];
sx q[1];
rz(-2.8956928) q[1];
sx q[1];
rz(-2.2794101) q[1];
rz(0.1578219) q[3];
sx q[3];
rz(-1.2994248) q[3];
sx q[3];
rz(-1.2543169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(-0.044163477) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039466) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(2.1283545) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-1.0838881) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7679374) q[0];
sx q[0];
rz(-2.8345788) q[0];
sx q[0];
rz(0.92371773) q[0];
rz(-pi) q[1];
rz(1.8439699) q[2];
sx q[2];
rz(-1.8130842) q[2];
sx q[2];
rz(0.66852409) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0005972) q[1];
sx q[1];
rz(-2.6037569) q[1];
sx q[1];
rz(-1.3586033) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51387042) q[3];
sx q[3];
rz(-1.6146982) q[3];
sx q[3];
rz(-1.1843137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
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
rz(-2.005902) q[1];
sx q[1];
rz(2.24618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58127922) q[0];
sx q[0];
rz(-2.7972097) q[0];
sx q[0];
rz(-0.11238213) q[0];
rz(-pi) q[1];
rz(-2.8199811) q[2];
sx q[2];
rz(-0.65882896) q[2];
sx q[2];
rz(-0.67473251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0370334) q[1];
sx q[1];
rz(-1.7393141) q[1];
sx q[1];
rz(-0.044635459) q[1];
rz(-pi) q[2];
rz(-1.496109) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(-1.1903654) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(-1.2566459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089766895) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(3.0999523) q[0];
rz(2.9314552) q[2];
sx q[2];
rz(-2.0268831) q[2];
sx q[2];
rz(1.8535341) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9650967) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(0.22775905) q[1];
rz(1.890896) q[3];
sx q[3];
rz(-2.2543636) q[3];
sx q[3];
rz(-2.2724831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(1.777565) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751188) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(-2.7546308) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2012126) q[0];
sx q[0];
rz(-1.3450755) q[0];
sx q[0];
rz(0.97198995) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1808124) q[2];
sx q[2];
rz(-0.79723251) q[2];
sx q[2];
rz(2.0955992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-2.4657001) q[1];
sx q[1];
rz(-0.11996108) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49658637) q[3];
sx q[3];
rz(-1.7644617) q[3];
sx q[3];
rz(0.95369875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62362921) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(-0.44000885) q[2];
rz(2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6599777) q[0];
sx q[0];
rz(-1.0064631) q[0];
sx q[0];
rz(-0.83156395) q[0];
rz(0.49528893) q[2];
sx q[2];
rz(-0.9501179) q[2];
sx q[2];
rz(0.6492614) q[2];
rz(-pi) q[3];
x q[3];
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
rz(0.52258073) q[1];
rz(-pi) q[2];
rz(0.52600577) q[3];
sx q[3];
rz(-2.3283539) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.4019029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.429075) q[0];
sx q[0];
rz(-1.6970465) q[0];
sx q[0];
rz(1.35668) q[0];
x q[1];
rz(1.5980814) q[2];
sx q[2];
rz(-1.7575043) q[2];
sx q[2];
rz(-0.4208496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2682174) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(-2.7482277) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7781156) q[3];
sx q[3];
rz(-0.94982409) q[3];
sx q[3];
rz(-0.45351115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4828651) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(-0.55082095) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(2.4089087) q[2];
sx q[2];
rz(-0.97326836) q[2];
sx q[2];
rz(1.2748981) q[2];
rz(2.0879073) q[3];
sx q[3];
rz(-1.5862982) q[3];
sx q[3];
rz(1.3261212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];