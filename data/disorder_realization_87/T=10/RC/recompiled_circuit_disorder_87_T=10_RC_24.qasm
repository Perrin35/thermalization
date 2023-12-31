OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(5.5570931) q[0];
sx q[0];
rz(9.2232016) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(2.2044866) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21298458) q[0];
sx q[0];
rz(-2.1323418) q[0];
sx q[0];
rz(-1.0667849) q[0];
rz(1.0785157) q[2];
sx q[2];
rz(-1.0168076) q[2];
sx q[2];
rz(-0.89821494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9846676) q[1];
sx q[1];
rz(-2.942454) q[1];
sx q[1];
rz(-1.1652633) q[1];
rz(1.9107781) q[3];
sx q[3];
rz(-2.8570606) q[3];
sx q[3];
rz(-1.1839379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(-3.0554331) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(1.8151059) q[0];
rz(1.2558698) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(-0.27145162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6730246) q[0];
sx q[0];
rz(-0.66267555) q[0];
sx q[0];
rz(2.2492692) q[0];
rz(-pi) q[1];
rz(-0.51606744) q[2];
sx q[2];
rz(-1.3000254) q[2];
sx q[2];
rz(0.64112907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2084864) q[1];
sx q[1];
rz(-0.83175627) q[1];
sx q[1];
rz(2.0887124) q[1];
x q[2];
rz(1.9430964) q[3];
sx q[3];
rz(-2.3447403) q[3];
sx q[3];
rz(-0.70355319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0945956) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(-2.5644152) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(-1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24401027) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(0.14257167) q[0];
rz(-1.3525195) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(0.20908633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7861048) q[0];
sx q[0];
rz(-2.3007563) q[0];
sx q[0];
rz(0.79746042) q[0];
x q[1];
rz(-1.9833343) q[2];
sx q[2];
rz(-2.2670855) q[2];
sx q[2];
rz(-2.11889) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28029728) q[1];
sx q[1];
rz(-2.0164844) q[1];
sx q[1];
rz(-2.6764826) q[1];
rz(-pi) q[2];
rz(-0.3029996) q[3];
sx q[3];
rz(-1.3649366) q[3];
sx q[3];
rz(2.1752165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(-1.8557619) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(-0.16734853) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(2.1787815) q[0];
rz(-2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22165933) q[0];
sx q[0];
rz(-2.1533305) q[0];
sx q[0];
rz(-2.4525053) q[0];
rz(-pi) q[1];
rz(0.0041299303) q[2];
sx q[2];
rz(-3.0175856) q[2];
sx q[2];
rz(-2.4069402) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0469989) q[1];
sx q[1];
rz(-1.7013036) q[1];
sx q[1];
rz(0.99884896) q[1];
rz(-0.56941454) q[3];
sx q[3];
rz(-1.6462407) q[3];
sx q[3];
rz(-0.30573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0756388) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(0.11165079) q[2];
rz(2.3305437) q[3];
sx q[3];
rz(-0.45447293) q[3];
sx q[3];
rz(-0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951293) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(0.35650373) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(2.8809663) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286752) q[0];
sx q[0];
rz(-2.9936491) q[0];
sx q[0];
rz(0.28767985) q[0];
rz(-2.8904301) q[2];
sx q[2];
rz(-1.3248982) q[2];
sx q[2];
rz(-2.9159301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.534429) q[1];
sx q[1];
rz(-2.5669614) q[1];
sx q[1];
rz(3.1092005) q[1];
rz(-0.14560933) q[3];
sx q[3];
rz(-1.429261) q[3];
sx q[3];
rz(-1.4841929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2300718) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(2.9122706) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(1.0361766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(1.7274436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8261674) q[0];
sx q[0];
rz(-0.27914473) q[0];
sx q[0];
rz(2.9690353) q[0];
rz(-1.3041777) q[2];
sx q[2];
rz(-2.0732217) q[2];
sx q[2];
rz(0.38602877) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9767178) q[1];
sx q[1];
rz(-2.4396982) q[1];
sx q[1];
rz(0.81275069) q[1];
rz(1.9023444) q[3];
sx q[3];
rz(-2.9132531) q[3];
sx q[3];
rz(-2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3843) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(3.0199155) q[0];
rz(-1.1514459) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-0.40245232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2936195) q[0];
sx q[0];
rz(-1.8403887) q[0];
sx q[0];
rz(-0.86975354) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5791113) q[2];
sx q[2];
rz(-1.5247716) q[2];
sx q[2];
rz(0.41551057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3869785) q[1];
sx q[1];
rz(-1.7474035) q[1];
sx q[1];
rz(-2.8547924) q[1];
rz(-pi) q[2];
rz(3.1137755) q[3];
sx q[3];
rz(-1.2226612) q[3];
sx q[3];
rz(0.58084014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(-0.86223117) q[2];
rz(0.47752738) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(2.4086337) q[0];
rz(3.0015302) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(2.1070811) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1326133) q[0];
sx q[0];
rz(-1.4011369) q[0];
sx q[0];
rz(1.8011814) q[0];
x q[1];
rz(2.5402252) q[2];
sx q[2];
rz(-3.0172918) q[2];
sx q[2];
rz(0.79324317) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58971436) q[1];
sx q[1];
rz(-2.1363746) q[1];
sx q[1];
rz(-3.0977071) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23969527) q[3];
sx q[3];
rz(-2.7125159) q[3];
sx q[3];
rz(2.808188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(2.365716) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4416606) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(-3.124776) q[0];
rz(-3.1230714) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-0.7787849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6803857) q[0];
sx q[0];
rz(-1.9068423) q[0];
sx q[0];
rz(-1.3075605) q[0];
rz(-pi) q[1];
rz(2.0428033) q[2];
sx q[2];
rz(-1.325377) q[2];
sx q[2];
rz(0.62308842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89214954) q[1];
sx q[1];
rz(-1.8037233) q[1];
sx q[1];
rz(1.9473142) q[1];
rz(0.064568297) q[3];
sx q[3];
rz(-1.7913006) q[3];
sx q[3];
rz(-1.8622423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4497711) q[2];
sx q[2];
rz(-1.3042973) q[2];
sx q[2];
rz(2.4251535) q[2];
rz(-1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(-1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.107782) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-2.9901436) q[0];
rz(0.17627136) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(2.418628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.003119) q[0];
sx q[0];
rz(-1.4644633) q[0];
sx q[0];
rz(1.7220201) q[0];
rz(-pi) q[1];
rz(-1.1638219) q[2];
sx q[2];
rz(-0.40728912) q[2];
sx q[2];
rz(-0.92330698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55262676) q[1];
sx q[1];
rz(-2.6347661) q[1];
sx q[1];
rz(-1.456702) q[1];
x q[2];
rz(1.9728174) q[3];
sx q[3];
rz(-2.3178181) q[3];
sx q[3];
rz(1.0812425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67939776) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(-2.6160713) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59239607) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(2.2407871) q[2];
sx q[2];
rz(-1.1849891) q[2];
sx q[2];
rz(1.6580788) q[2];
rz(0.59393926) q[3];
sx q[3];
rz(-2.754302) q[3];
sx q[3];
rz(2.1828628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
