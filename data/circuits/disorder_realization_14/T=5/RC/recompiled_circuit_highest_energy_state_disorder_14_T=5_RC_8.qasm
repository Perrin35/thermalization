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
rz(-1.5934264) q[0];
sx q[0];
rz(-0.46841535) q[0];
sx q[0];
rz(-3.0031437) q[0];
rz(-1.2568714) q[1];
sx q[1];
rz(-1.0705907) q[1];
sx q[1];
rz(-3.1201153) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3250752) q[0];
sx q[0];
rz(-1.3680352) q[0];
sx q[0];
rz(1.2120423) q[0];
rz(0.8201222) q[2];
sx q[2];
rz(-2.409365) q[2];
sx q[2];
rz(2.2158465) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5925046) q[1];
sx q[1];
rz(-0.23356423) q[1];
sx q[1];
rz(-2.2930868) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2091694) q[3];
sx q[3];
rz(-1.0054133) q[3];
sx q[3];
rz(-2.9361182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5452177) q[2];
sx q[2];
rz(-0.45520982) q[2];
sx q[2];
rz(-1.630416) q[2];
rz(-0.049985416) q[3];
sx q[3];
rz(-1.9129916) q[3];
sx q[3];
rz(-0.25240067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2216457) q[0];
sx q[0];
rz(-1.1476465) q[0];
sx q[0];
rz(3.1044712) q[0];
rz(-0.014017398) q[1];
sx q[1];
rz(-0.57828301) q[1];
sx q[1];
rz(0.23606539) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136617) q[0];
sx q[0];
rz(-1.4959916) q[0];
sx q[0];
rz(1.5563677) q[0];
rz(-pi) q[1];
rz(-1.0843956) q[2];
sx q[2];
rz(-0.91569967) q[2];
sx q[2];
rz(-2.4156605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29480793) q[1];
sx q[1];
rz(-0.67344147) q[1];
sx q[1];
rz(0.93017857) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1680683) q[3];
sx q[3];
rz(-0.42713886) q[3];
sx q[3];
rz(1.1534302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1290805) q[2];
sx q[2];
rz(-1.3690288) q[2];
sx q[2];
rz(3.1079666) q[2];
rz(-2.8513837) q[3];
sx q[3];
rz(-0.84607327) q[3];
sx q[3];
rz(-2.7510551) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6212807) q[0];
sx q[0];
rz(-2.1618167) q[0];
sx q[0];
rz(-2.8424971) q[0];
rz(-1.6200804) q[1];
sx q[1];
rz(-0.027093096) q[1];
sx q[1];
rz(-3.1153968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55518245) q[0];
sx q[0];
rz(-1.5497218) q[0];
sx q[0];
rz(-0.0064959244) q[0];
x q[1];
rz(-0.14339672) q[2];
sx q[2];
rz(-1.2931648) q[2];
sx q[2];
rz(1.5609225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67945665) q[1];
sx q[1];
rz(-2.1894881) q[1];
sx q[1];
rz(2.2652744) q[1];
rz(2.6869925) q[3];
sx q[3];
rz(-2.6516783) q[3];
sx q[3];
rz(-2.5471806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4574778) q[2];
sx q[2];
rz(-2.699615) q[2];
sx q[2];
rz(-0.6074062) q[2];
rz(0.36951798) q[3];
sx q[3];
rz(-1.4990467) q[3];
sx q[3];
rz(0.0079689715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5615416) q[0];
sx q[0];
rz(-0.20195584) q[0];
sx q[0];
rz(-1.1308905) q[0];
rz(-2.789403) q[1];
sx q[1];
rz(-2.6512841) q[1];
sx q[1];
rz(-0.21603781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9176005) q[0];
sx q[0];
rz(-0.52187863) q[0];
sx q[0];
rz(-2.6386369) q[0];
rz(-pi) q[1];
rz(-0.44682002) q[2];
sx q[2];
rz(-1.6761314) q[2];
sx q[2];
rz(-3.1329495) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32111383) q[1];
sx q[1];
rz(-2.0778837) q[1];
sx q[1];
rz(0.45407461) q[1];
rz(-pi) q[2];
rz(2.0497672) q[3];
sx q[3];
rz(-2.9505725) q[3];
sx q[3];
rz(-2.614792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.617368) q[2];
sx q[2];
rz(-0.91297954) q[2];
sx q[2];
rz(-0.53140223) q[2];
rz(1.6898539) q[3];
sx q[3];
rz(-0.51133358) q[3];
sx q[3];
rz(2.1112198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.3568929) q[0];
sx q[0];
rz(-2.0578616) q[0];
sx q[0];
rz(2.8392131) q[0];
rz(-2.0284292) q[1];
sx q[1];
rz(-2.1402363) q[1];
sx q[1];
rz(-1.0611634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9386834) q[0];
sx q[0];
rz(-1.6309982) q[0];
sx q[0];
rz(3.0767482) q[0];
x q[1];
rz(-1.6689014) q[2];
sx q[2];
rz(-1.3472234) q[2];
sx q[2];
rz(-1.6983216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1969152) q[1];
sx q[1];
rz(-2.3507171) q[1];
sx q[1];
rz(-2.9612035) q[1];
rz(-pi) q[2];
rz(-0.29814675) q[3];
sx q[3];
rz(-1.8121208) q[3];
sx q[3];
rz(-2.4467447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82659668) q[2];
sx q[2];
rz(-1.7642517) q[2];
sx q[2];
rz(-0.19485168) q[2];
rz(-2.3991614) q[3];
sx q[3];
rz(-0.73115474) q[3];
sx q[3];
rz(-2.8661695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9575397) q[0];
sx q[0];
rz(-2.4803211) q[0];
sx q[0];
rz(-2.0273965) q[0];
rz(2.8320352) q[1];
sx q[1];
rz(-2.20859) q[1];
sx q[1];
rz(1.385744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838309) q[0];
sx q[0];
rz(-1.2886192) q[0];
sx q[0];
rz(0.18309219) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0631752) q[2];
sx q[2];
rz(-1.4709899) q[2];
sx q[2];
rz(2.5863079) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42658994) q[1];
sx q[1];
rz(-2.4414483) q[1];
sx q[1];
rz(-2.7310074) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65151229) q[3];
sx q[3];
rz(-0.88884547) q[3];
sx q[3];
rz(-1.3523952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4625357) q[2];
sx q[2];
rz(-2.6947196) q[2];
sx q[2];
rz(-0.40979579) q[2];
rz(0.077824079) q[3];
sx q[3];
rz(-1.2160559) q[3];
sx q[3];
rz(2.3791544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78793144) q[0];
sx q[0];
rz(-2.9910112) q[0];
sx q[0];
rz(1.6702363) q[0];
rz(-2.3279066) q[1];
sx q[1];
rz(-2.4206471) q[1];
sx q[1];
rz(-2.241316) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532458) q[0];
sx q[0];
rz(-1.4524967) q[0];
sx q[0];
rz(1.6315881) q[0];
x q[1];
rz(1.2834107) q[2];
sx q[2];
rz(-1.6528439) q[2];
sx q[2];
rz(2.1211185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1301816) q[1];
sx q[1];
rz(-0.26618845) q[1];
sx q[1];
rz(1.601851) q[1];
rz(-pi) q[2];
rz(-2.8741292) q[3];
sx q[3];
rz(-0.24175343) q[3];
sx q[3];
rz(0.072168298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6206996) q[2];
sx q[2];
rz(-0.40105477) q[2];
sx q[2];
rz(2.9996784) q[2];
rz(0.17169954) q[3];
sx q[3];
rz(-1.2407691) q[3];
sx q[3];
rz(-0.54952526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5080344) q[0];
sx q[0];
rz(-1.1548076) q[0];
sx q[0];
rz(-1.5800193) q[0];
rz(-0.70478565) q[1];
sx q[1];
rz(-0.90122688) q[1];
sx q[1];
rz(0.27552342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412236) q[0];
sx q[0];
rz(-1.4432206) q[0];
sx q[0];
rz(-3.0795549) q[0];
rz(2.9887096) q[2];
sx q[2];
rz(-2.8934921) q[2];
sx q[2];
rz(2.5923592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0184014) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(1.0276612) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4750214) q[3];
sx q[3];
rz(-1.2351994) q[3];
sx q[3];
rz(-2.1977212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0602818) q[2];
sx q[2];
rz(-1.5014481) q[2];
sx q[2];
rz(-0.29067972) q[2];
rz(-0.79689133) q[3];
sx q[3];
rz(-0.47178888) q[3];
sx q[3];
rz(-2.1577788) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31442916) q[0];
sx q[0];
rz(-1.6570579) q[0];
sx q[0];
rz(2.2737801) q[0];
rz(-0.21226352) q[1];
sx q[1];
rz(-2.2607195) q[1];
sx q[1];
rz(-2.8689522) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6295687) q[0];
sx q[0];
rz(-1.6457412) q[0];
sx q[0];
rz(-1.7864947) q[0];
rz(3.0402203) q[2];
sx q[2];
rz(-2.3194312) q[2];
sx q[2];
rz(-3.0790975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8114879) q[1];
sx q[1];
rz(-1.5327785) q[1];
sx q[1];
rz(2.2056715) q[1];
rz(-pi) q[2];
rz(-2.7271184) q[3];
sx q[3];
rz(-1.8467086) q[3];
sx q[3];
rz(-1.3434354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0923882) q[2];
sx q[2];
rz(-0.35085756) q[2];
sx q[2];
rz(1.1445047) q[2];
rz(2.5149964) q[3];
sx q[3];
rz(-2.383039) q[3];
sx q[3];
rz(-2.8265317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3939683) q[0];
sx q[0];
rz(-0.67248857) q[0];
sx q[0];
rz(-2.5482063) q[0];
rz(-0.74140948) q[1];
sx q[1];
rz(-0.64677042) q[1];
sx q[1];
rz(-1.5470362) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0295059) q[0];
sx q[0];
rz(-1.7884147) q[0];
sx q[0];
rz(3.0299216) q[0];
x q[1];
rz(-0.41489652) q[2];
sx q[2];
rz(-1.7977568) q[2];
sx q[2];
rz(1.0912053) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.64338291) q[1];
sx q[1];
rz(-2.3003182) q[1];
sx q[1];
rz(1.4773083) q[1];
x q[2];
rz(1.7150899) q[3];
sx q[3];
rz(-2.9375667) q[3];
sx q[3];
rz(-2.229759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0934304) q[2];
sx q[2];
rz(-0.54485816) q[2];
sx q[2];
rz(-2.1989934) q[2];
rz(1.7132828) q[3];
sx q[3];
rz(-1.8918248) q[3];
sx q[3];
rz(1.1552756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9485332) q[0];
sx q[0];
rz(-1.5848703) q[0];
sx q[0];
rz(1.7901044) q[0];
rz(-1.1763186) q[1];
sx q[1];
rz(-0.90570025) q[1];
sx q[1];
rz(-0.66014231) q[1];
rz(-0.17778601) q[2];
sx q[2];
rz(-1.8443731) q[2];
sx q[2];
rz(2.2603879) q[2];
rz(2.6082718) q[3];
sx q[3];
rz(-1.7989895) q[3];
sx q[3];
rz(-1.2839272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
