OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(1.8068846) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(0.97595739) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53544331) q[0];
sx q[0];
rz(-2.1436999) q[0];
sx q[0];
rz(-1.2826305) q[0];
x q[1];
rz(2.8085254) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(-1.2531467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9589899) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(-1.3401003) q[1];
x q[2];
rz(-2.7482412) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(2.2497183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41539899) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(-2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48288229) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9316677) q[0];
sx q[0];
rz(-1.4405662) q[0];
sx q[0];
rz(1.3757214) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86654051) q[2];
sx q[2];
rz(-0.73359493) q[2];
sx q[2];
rz(-1.7689592) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22390631) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(-2.9427337) q[1];
x q[2];
rz(-1.3199602) q[3];
sx q[3];
rz(-1.0227961) q[3];
sx q[3];
rz(0.26720023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(-0.87810278) q[2];
rz(2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(-0.01097824) q[0];
rz(-0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(0.095741622) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11287963) q[0];
sx q[0];
rz(-2.9071147) q[0];
sx q[0];
rz(-2.2857091) q[0];
rz(-1.1126306) q[2];
sx q[2];
rz(-1.1095699) q[2];
sx q[2];
rz(1.065965) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15624554) q[1];
sx q[1];
rz(-0.74838446) q[1];
sx q[1];
rz(2.2845539) q[1];
rz(-0.057283244) q[3];
sx q[3];
rz(-2.5896642) q[3];
sx q[3];
rz(2.3416167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0720955) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(-0.52571458) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(1.0338763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(-2.4497776) q[0];
rz(-0.83861645) q[2];
sx q[2];
rz(-2.4258483) q[2];
sx q[2];
rz(1.0775623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.911072) q[1];
sx q[1];
rz(-1.5057179) q[1];
sx q[1];
rz(-1.7622403) q[1];
x q[2];
rz(0.28169607) q[3];
sx q[3];
rz(-1.8857737) q[3];
sx q[3];
rz(-1.7236934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(1.3467849) q[2];
rz(2.7205617) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088257) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(-2.5674852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78050437) q[0];
sx q[0];
rz(-1.442712) q[0];
sx q[0];
rz(2.1696521) q[0];
x q[1];
rz(-1.5482043) q[2];
sx q[2];
rz(-2.2859757) q[2];
sx q[2];
rz(1.9695645) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1787604) q[1];
sx q[1];
rz(-1.292359) q[1];
sx q[1];
rz(-0.36414418) q[1];
rz(0.19458171) q[3];
sx q[3];
rz(-0.91472018) q[3];
sx q[3];
rz(1.3692828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(1.627702) q[2];
rz(3.1001575) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-0.69818991) q[0];
rz(-2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(-2.4093157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2039316) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(2.926814) q[0];
rz(-0.68533021) q[2];
sx q[2];
rz(-2.2517859) q[2];
sx q[2];
rz(2.6076917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9975035) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(1.287582) q[1];
rz(1.3527649) q[3];
sx q[3];
rz(-2.1302967) q[3];
sx q[3];
rz(2.7263209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7531062) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.1614655) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(-0.15643315) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(0.77004534) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2609445) q[0];
sx q[0];
rz(-2.5136607) q[0];
sx q[0];
rz(2.8448366) q[0];
x q[1];
rz(0.98203512) q[2];
sx q[2];
rz(-1.9141478) q[2];
sx q[2];
rz(-2.8729168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5321977) q[1];
sx q[1];
rz(-0.97071338) q[1];
sx q[1];
rz(0.9432015) q[1];
rz(-pi) q[2];
rz(-1.2392063) q[3];
sx q[3];
rz(-2.050403) q[3];
sx q[3];
rz(0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(1.0151803) q[2];
rz(1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(-2.4027951) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(0.36639211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893433) q[0];
sx q[0];
rz(-2.0848668) q[0];
sx q[0];
rz(-0.71787562) q[0];
rz(-pi) q[1];
rz(1.1364469) q[2];
sx q[2];
rz(-2.1546116) q[2];
sx q[2];
rz(-2.7455612) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3185127) q[1];
sx q[1];
rz(-1.7365343) q[1];
sx q[1];
rz(2.4351099) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1226235) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(-3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-0.88395399) q[0];
rz(0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(0.79137897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7991379) q[0];
sx q[0];
rz(-0.90478169) q[0];
sx q[0];
rz(-2.7333583) q[0];
rz(-pi) q[1];
rz(2.0558526) q[2];
sx q[2];
rz(-1.1527449) q[2];
sx q[2];
rz(-0.56318356) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2453354) q[1];
sx q[1];
rz(-2.5108813) q[1];
sx q[1];
rz(-2.3920822) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16718849) q[3];
sx q[3];
rz(-2.1033923) q[3];
sx q[3];
rz(2.5933468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(-2.7630473) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-2.5337906) q[0];
rz(2.9027477) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.7609319) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35683435) q[0];
sx q[0];
rz(-2.1958302) q[0];
sx q[0];
rz(-2.9805095) q[0];
x q[1];
rz(-1.626255) q[2];
sx q[2];
rz(-2.7055253) q[2];
sx q[2];
rz(0.01576327) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7271991) q[1];
sx q[1];
rz(-0.52800035) q[1];
sx q[1];
rz(-2.2382733) q[1];
rz(-1.4019743) q[3];
sx q[3];
rz(-1.7572548) q[3];
sx q[3];
rz(-0.25690119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7878788) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(0.33622462) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(0.54429383) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(1.9359679) q[2];
sx q[2];
rz(-2.361459) q[2];
sx q[2];
rz(-0.38103719) q[2];
rz(1.0014793) q[3];
sx q[3];
rz(-1.3712728) q[3];
sx q[3];
rz(-1.2029592) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];