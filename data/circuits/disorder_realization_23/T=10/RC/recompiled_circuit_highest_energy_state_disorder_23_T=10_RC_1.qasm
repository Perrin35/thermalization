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
rz(0.5956369) q[0];
sx q[0];
rz(3.5516153) q[0];
sx q[0];
rz(8.6300996) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(5.3310634) q[1];
sx q[1];
rz(7.4330243) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73939378) q[0];
sx q[0];
rz(-1.5350685) q[0];
sx q[0];
rz(1.2767919) q[0];
x q[1];
rz(1.8914521) q[2];
sx q[2];
rz(-2.3806664) q[2];
sx q[2];
rz(-0.93283949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2268035) q[1];
sx q[1];
rz(-2.7046596) q[1];
sx q[1];
rz(-2.7257257) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47969476) q[3];
sx q[3];
rz(-0.51691662) q[3];
sx q[3];
rz(1.2933047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6940234) q[2];
sx q[2];
rz(-2.1489216) q[2];
sx q[2];
rz(0.64219323) q[2];
rz(1.0688952) q[3];
sx q[3];
rz(-1.5690208) q[3];
sx q[3];
rz(-1.6525035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25817961) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(-0.51134837) q[0];
rz(-1.5072352) q[1];
sx q[1];
rz(-1.2751445) q[1];
sx q[1];
rz(-1.7587597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16973142) q[0];
sx q[0];
rz(-1.3644553) q[0];
sx q[0];
rz(-1.1792762) q[0];
rz(-1.1452008) q[2];
sx q[2];
rz(-0.71768242) q[2];
sx q[2];
rz(1.8459783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1369776) q[1];
sx q[1];
rz(-1.2619234) q[1];
sx q[1];
rz(1.1501794) q[1];
rz(-pi) q[2];
rz(0.52329833) q[3];
sx q[3];
rz(-1.4804258) q[3];
sx q[3];
rz(-3.0449611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4649268) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(-1.0627221) q[3];
sx q[3];
rz(-1.9461742) q[3];
sx q[3];
rz(0.11500558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677143) q[0];
sx q[0];
rz(-2.102484) q[0];
sx q[0];
rz(-0.1097196) q[0];
rz(2.8166855) q[1];
sx q[1];
rz(-1.1512681) q[1];
sx q[1];
rz(2.6085764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1905381) q[0];
sx q[0];
rz(-1.2742654) q[0];
sx q[0];
rz(1.7151136) q[0];
x q[1];
rz(-1.5948086) q[2];
sx q[2];
rz(-1.5188076) q[2];
sx q[2];
rz(-0.94357027) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17688454) q[1];
sx q[1];
rz(-2.1358651) q[1];
sx q[1];
rz(0.76660291) q[1];
x q[2];
rz(1.1148861) q[3];
sx q[3];
rz(-0.81509903) q[3];
sx q[3];
rz(-0.86642735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.816232) q[2];
sx q[2];
rz(-1.3282789) q[2];
sx q[2];
rz(-1.4556966) q[2];
rz(3.0560737) q[3];
sx q[3];
rz(-1.0695846) q[3];
sx q[3];
rz(0.94676179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7607255) q[0];
sx q[0];
rz(-2.4873698) q[0];
sx q[0];
rz(0.72232676) q[0];
rz(1.5708615) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(-0.042472366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4054779) q[0];
sx q[0];
rz(-1.7543989) q[0];
sx q[0];
rz(2.6915789) q[0];
rz(-pi) q[1];
rz(-1.7444403) q[2];
sx q[2];
rz(-3.0499027) q[2];
sx q[2];
rz(-1.1806115) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.153923) q[1];
sx q[1];
rz(-0.81036416) q[1];
sx q[1];
rz(2.8778879) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23353429) q[3];
sx q[3];
rz(-0.81538768) q[3];
sx q[3];
rz(2.5366207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.2461569) q[2];
sx q[2];
rz(-1.3763206) q[2];
sx q[2];
rz(2.8975272) q[2];
rz(2.3501588) q[3];
sx q[3];
rz(-1.8297198) q[3];
sx q[3];
rz(-2.4401276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0303845) q[0];
sx q[0];
rz(-0.1143488) q[0];
sx q[0];
rz(2.9682888) q[0];
rz(0.93470848) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-2.8876143) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.105071) q[0];
sx q[0];
rz(-1.6833841) q[0];
sx q[0];
rz(-2.3238411) q[0];
rz(0.51176519) q[2];
sx q[2];
rz(-0.6590656) q[2];
sx q[2];
rz(-2.2816254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7813999) q[1];
sx q[1];
rz(-1.5915697) q[1];
sx q[1];
rz(2.4482577) q[1];
rz(-pi) q[2];
x q[2];
rz(0.06490223) q[3];
sx q[3];
rz(-1.1838473) q[3];
sx q[3];
rz(-2.153819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0737334) q[2];
sx q[2];
rz(-2.8974055) q[2];
sx q[2];
rz(-1.7013288) q[2];
rz(-2.4654147) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(-0.079782709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22236958) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(2.4134912) q[0];
rz(-1.1657731) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(0.56070915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9084204) q[0];
sx q[0];
rz(-1.6357372) q[0];
sx q[0];
rz(3.0149595) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6892904) q[2];
sx q[2];
rz(-1.0241707) q[2];
sx q[2];
rz(-0.34260633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6627359) q[1];
sx q[1];
rz(-0.4534035) q[1];
sx q[1];
rz(1.0591566) q[1];
rz(3.0440935) q[3];
sx q[3];
rz(-1.6336771) q[3];
sx q[3];
rz(1.1600375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59549436) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(1.072849) q[2];
rz(-0.32235518) q[3];
sx q[3];
rz(-0.41283804) q[3];
sx q[3];
rz(-2.7505007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3748465) q[0];
sx q[0];
rz(-1.5006737) q[0];
sx q[0];
rz(0.40819502) q[0];
rz(0.32632581) q[1];
sx q[1];
rz(-0.34714454) q[1];
sx q[1];
rz(-1.6654642) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6107619) q[0];
sx q[0];
rz(-0.179804) q[0];
sx q[0];
rz(-1.7925949) q[0];
rz(-1.9397975) q[2];
sx q[2];
rz(-0.59518669) q[2];
sx q[2];
rz(3.0734504) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1154815) q[1];
sx q[1];
rz(-0.12899765) q[1];
sx q[1];
rz(-1.1327101) q[1];
rz(-2.9093698) q[3];
sx q[3];
rz(-1.7360064) q[3];
sx q[3];
rz(-2.8286407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8866715) q[2];
sx q[2];
rz(-1.485606) q[2];
sx q[2];
rz(2.2754748) q[2];
rz(0.73650375) q[3];
sx q[3];
rz(-2.0564506) q[3];
sx q[3];
rz(0.88732639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12164584) q[0];
sx q[0];
rz(-0.045504657) q[0];
sx q[0];
rz(1.2787) q[0];
rz(-1.0026898) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(-0.444828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0507752) q[0];
sx q[0];
rz(-1.4231761) q[0];
sx q[0];
rz(-1.0423129) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72250267) q[2];
sx q[2];
rz(-1.926427) q[2];
sx q[2];
rz(-2.1788545) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0623183) q[1];
sx q[1];
rz(-1.4517112) q[1];
sx q[1];
rz(-0.80715413) q[1];
rz(2.7642207) q[3];
sx q[3];
rz(-1.3283973) q[3];
sx q[3];
rz(0.19211543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43850809) q[2];
sx q[2];
rz(-1.727641) q[2];
sx q[2];
rz(0.50602305) q[2];
rz(-2.1972726) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(-0.73616141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.849702) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(3.131409) q[0];
rz(-0.61095515) q[1];
sx q[1];
rz(-0.98038951) q[1];
sx q[1];
rz(-3.0593061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1707538) q[0];
sx q[0];
rz(-2.0505095) q[0];
sx q[0];
rz(-1.9895577) q[0];
rz(-1.8533704) q[2];
sx q[2];
rz(-2.4420441) q[2];
sx q[2];
rz(0.69751537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7499867) q[1];
sx q[1];
rz(-2.9175903) q[1];
sx q[1];
rz(0.27952607) q[1];
rz(-2.4769267) q[3];
sx q[3];
rz(-1.3997404) q[3];
sx q[3];
rz(1.4881575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33132195) q[2];
sx q[2];
rz(-1.4484582) q[2];
sx q[2];
rz(-2.4324379) q[2];
rz(-1.4823312) q[3];
sx q[3];
rz(-0.63153657) q[3];
sx q[3];
rz(-2.6813996) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3191147) q[0];
sx q[0];
rz(-2.3111486) q[0];
sx q[0];
rz(0.6293695) q[0];
rz(-1.7259701) q[1];
sx q[1];
rz(-1.4868088) q[1];
sx q[1];
rz(1.1782882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30692682) q[0];
sx q[0];
rz(-2.3884058) q[0];
sx q[0];
rz(-2.1801394) q[0];
x q[1];
rz(0.35188108) q[2];
sx q[2];
rz(-1.0946254) q[2];
sx q[2];
rz(1.1426403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8161623) q[1];
sx q[1];
rz(-0.79968444) q[1];
sx q[1];
rz(-2.9887448) q[1];
rz(-pi) q[2];
rz(0.62507665) q[3];
sx q[3];
rz(-2.3155594) q[3];
sx q[3];
rz(-2.1476098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2668931) q[2];
sx q[2];
rz(-1.5692254) q[2];
sx q[2];
rz(0.81864041) q[2];
rz(1.5107752) q[3];
sx q[3];
rz(-1.8418334) q[3];
sx q[3];
rz(2.7394133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4928987) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(1.2661487) q[1];
sx q[1];
rz(-1.1504953) q[1];
sx q[1];
rz(-1.6082416) q[1];
rz(1.8902334) q[2];
sx q[2];
rz(-0.70079679) q[2];
sx q[2];
rz(0.23637017) q[2];
rz(-1.6615909) q[3];
sx q[3];
rz(-1.7547073) q[3];
sx q[3];
rz(1.2622866) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
