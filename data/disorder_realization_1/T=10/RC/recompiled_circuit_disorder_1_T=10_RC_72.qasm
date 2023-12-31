OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1228468) q[0];
sx q[0];
rz(-1.393264) q[0];
sx q[0];
rz(-1.2794504) q[0];
rz(-pi) q[1];
rz(2.675406) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(0.28238645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8661583) q[1];
sx q[1];
rz(-2.1992116) q[1];
sx q[1];
rz(2.1584312) q[1];
rz(-pi) q[2];
rz(1.7883349) q[3];
sx q[3];
rz(-1.4635524) q[3];
sx q[3];
rz(-1.6579962) q[3];
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
rz(-2.0404405) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(-0.027659841) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(-1.123463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9155974) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(-2.0321839) q[0];
rz(-pi) q[1];
rz(1.5078817) q[2];
sx q[2];
rz(-2.3496369) q[2];
sx q[2];
rz(-2.6346249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9065735) q[1];
sx q[1];
rz(-1.0636914) q[1];
sx q[1];
rz(-2.1706235) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65621891) q[3];
sx q[3];
rz(-1.1074293) q[3];
sx q[3];
rz(-3.0955293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(-0.91903764) q[2];
rz(2.4675026) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640901) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-1.1330053) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5223761) q[0];
sx q[0];
rz(-1.7130865) q[0];
sx q[0];
rz(3.1382568) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50425597) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(0.91932636) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9085711) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(-0.33645333) q[1];
x q[2];
rz(2.3660925) q[3];
sx q[3];
rz(-2.3470078) q[3];
sx q[3];
rz(2.9426108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(-3.1022762) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(-3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121658) q[0];
sx q[0];
rz(-0.78071785) q[0];
sx q[0];
rz(-2.6902945) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3877669) q[2];
sx q[2];
rz(-0.39003885) q[2];
sx q[2];
rz(2.4412465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.23952661) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(2.2794101) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2961779) q[3];
sx q[3];
rz(-1.4187958) q[3];
sx q[3];
rz(0.27384588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9049412) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(-0.044163477) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(-2.0577046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3736553) q[0];
sx q[0];
rz(-2.8345788) q[0];
sx q[0];
rz(-2.2178749) q[0];
x q[1];
rz(2.8903557) q[2];
sx q[2];
rz(-1.8357956) q[2];
sx q[2];
rz(-0.83515177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.52884) q[1];
sx q[1];
rz(-1.4626979) q[1];
sx q[1];
rz(1.0428863) q[1];
x q[2];
rz(1.5203939) q[3];
sx q[3];
rz(-1.0574697) q[3];
sx q[3];
rz(0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844834) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(2.24618) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58127922) q[0];
sx q[0];
rz(-2.7972097) q[0];
sx q[0];
rz(3.0292105) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8108098) q[2];
sx q[2];
rz(-2.1905106) q[2];
sx q[2];
rz(2.0680708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10455924) q[1];
sx q[1];
rz(-1.4022786) q[1];
sx q[1];
rz(-0.044635459) q[1];
x q[2];
rz(1.6020681) q[3];
sx q[3];
rz(-0.074723738) q[3];
sx q[3];
rz(0.63262313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.133693) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(2.3383979) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068709277) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.8849467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089766895) q[0];
sx q[0];
rz(-1.3194114) q[0];
sx q[0];
rz(0.041640394) q[0];
rz(-pi) q[1];
rz(-1.16876) q[2];
sx q[2];
rz(-0.49905825) q[2];
sx q[2];
rz(1.7390342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7278442) q[1];
sx q[1];
rz(-1.8520978) q[1];
sx q[1];
rz(-1.5037699) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36862936) q[3];
sx q[3];
rz(-2.3978524) q[3];
sx q[3];
rz(-1.3524692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98823035) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(1.777565) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(2.7546308) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-0.63502705) q[0];
sx q[0];
rz(-1.1839266) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96078028) q[2];
sx q[2];
rz(-2.3443601) q[2];
sx q[2];
rz(1.0459935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8291694) q[1];
sx q[1];
rz(-0.90065354) q[1];
sx q[1];
rz(1.4751242) q[1];
rz(2.7510795) q[3];
sx q[3];
rz(-0.53005855) q[3];
sx q[3];
rz(-0.9583677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(-0.7157588) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(-0.72775841) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54287275) q[0];
sx q[0];
rz(-0.96519404) q[0];
sx q[0];
rz(-0.70830317) q[0];
rz(-0.49528893) q[2];
sx q[2];
rz(-2.1914748) q[2];
sx q[2];
rz(-2.4923313) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8763435) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(-2.6190119) q[1];
rz(-0.52600577) q[3];
sx q[3];
rz(-0.81323871) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(-0.72193974) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(-0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(-2.06185) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.4019029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83090529) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(3.0124245) q[0];
rz(-1.5980814) q[2];
sx q[2];
rz(-1.7575043) q[2];
sx q[2];
rz(-2.7207431) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2682174) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(-2.7482277) q[1];
x q[2];
rz(2.8614282) q[3];
sx q[3];
rz(-0.6503085) q[3];
sx q[3];
rz(0.10661099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4828651) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.520291) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(-0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(-2.2254754) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(2.312071) q[2];
sx q[2];
rz(-0.98486949) q[2];
sx q[2];
rz(-2.9688901) q[2];
rz(0.017833088) q[3];
sx q[3];
rz(-1.0537536) q[3];
sx q[3];
rz(-0.23585933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
