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
rz(-1.0567226) q[0];
sx q[0];
rz(-1.280008) q[0];
sx q[0];
rz(-0.7260538) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(-1.0322303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2054133) q[0];
sx q[0];
rz(-1.6717311) q[0];
sx q[0];
rz(-2.4839782) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1025508) q[2];
sx q[2];
rz(-2.0882975) q[2];
sx q[2];
rz(1.3618748) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5121442) q[1];
sx q[1];
rz(-1.2595065) q[1];
sx q[1];
rz(-1.4427079) q[1];
x q[2];
rz(1.0606403) q[3];
sx q[3];
rz(-2.3054696) q[3];
sx q[3];
rz(-1.5165129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44077474) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(-0.087892858) q[2];
rz(-1.4676189) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(2.1427593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217594) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(-1.4922967) q[0];
rz(-2.3536033) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(-0.96510395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5934629) q[0];
sx q[0];
rz(-0.735518) q[0];
sx q[0];
rz(-0.1372565) q[0];
rz(-0.73678686) q[2];
sx q[2];
rz(-1.4649142) q[2];
sx q[2];
rz(2.4491765) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31307855) q[1];
sx q[1];
rz(-1.113849) q[1];
sx q[1];
rz(2.9345153) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9645637) q[3];
sx q[3];
rz(-2.6654589) q[3];
sx q[3];
rz(-3.0974922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1284156) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(-1.3236807) q[2];
rz(0.88661083) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(-2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(2.4460728) q[0];
rz(-0.8356525) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(-2.4998891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1276961) q[0];
sx q[0];
rz(-2.3900716) q[0];
sx q[0];
rz(1.7836545) q[0];
rz(0.71126781) q[2];
sx q[2];
rz(-1.6984273) q[2];
sx q[2];
rz(-0.39943275) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56395951) q[1];
sx q[1];
rz(-2.0490987) q[1];
sx q[1];
rz(2.7863414) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1459177) q[3];
sx q[3];
rz(-2.174587) q[3];
sx q[3];
rz(3.1310981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(1.3345435) q[2];
rz(-1.7415907) q[3];
sx q[3];
rz(-1.497523) q[3];
sx q[3];
rz(-1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982933) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(-0.77348462) q[0];
rz(-3.0553014) q[1];
sx q[1];
rz(-0.50420612) q[1];
sx q[1];
rz(-0.48826826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9241656) q[0];
sx q[0];
rz(-1.1091976) q[0];
sx q[0];
rz(-1.6469) q[0];
rz(-0.69741285) q[2];
sx q[2];
rz(-2.0048755) q[2];
sx q[2];
rz(-0.12034872) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1459256) q[1];
sx q[1];
rz(-2.2856376) q[1];
sx q[1];
rz(2.8413474) q[1];
rz(2.3523496) q[3];
sx q[3];
rz(-2.5960138) q[3];
sx q[3];
rz(-0.18462791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6733072) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(-0.16806531) q[2];
rz(2.7189861) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(-2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.1905097) q[0];
rz(-2.6137784) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(0.0091008069) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9267488) q[0];
sx q[0];
rz(-2.0159449) q[0];
sx q[0];
rz(2.4657004) q[0];
rz(1.509769) q[2];
sx q[2];
rz(-1.6313071) q[2];
sx q[2];
rz(0.64693816) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.121351) q[1];
sx q[1];
rz(-0.96881908) q[1];
sx q[1];
rz(0.51687981) q[1];
rz(0.96076699) q[3];
sx q[3];
rz(-2.4159263) q[3];
sx q[3];
rz(2.8270367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3147543) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(-0.66519386) q[2];
rz(-1.0264171) q[3];
sx q[3];
rz(-2.0668991) q[3];
sx q[3];
rz(-2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196395) q[0];
sx q[0];
rz(-0.53549796) q[0];
sx q[0];
rz(0.39805472) q[0];
rz(-1.4705426) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(-2.3614531) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3624293) q[0];
sx q[0];
rz(-1.0991968) q[0];
sx q[0];
rz(-2.7410024) q[0];
x q[1];
rz(-1.7731116) q[2];
sx q[2];
rz(-1.8087808) q[2];
sx q[2];
rz(-1.1954824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33360219) q[1];
sx q[1];
rz(-1.2871075) q[1];
sx q[1];
rz(-0.70244838) q[1];
x q[2];
rz(2.2557115) q[3];
sx q[3];
rz(-1.8198593) q[3];
sx q[3];
rz(-2.0307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.857343) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(-2.9798689) q[2];
rz(0.11180793) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(-0.73981729) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(0.69851843) q[0];
rz(-1.0922208) q[1];
sx q[1];
rz(-2.3869546) q[1];
sx q[1];
rz(3.0879367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.221401) q[0];
sx q[0];
rz(-2.2760253) q[0];
sx q[0];
rz(0.34366519) q[0];
rz(-pi) q[1];
rz(-2.7611465) q[2];
sx q[2];
rz(-1.9815677) q[2];
sx q[2];
rz(3.041317) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21660575) q[1];
sx q[1];
rz(-1.1007855) q[1];
sx q[1];
rz(-1.6307249) q[1];
rz(-pi) q[2];
rz(-0.92412432) q[3];
sx q[3];
rz(-1.4507626) q[3];
sx q[3];
rz(-1.8450575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0003537) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(-0.21256438) q[2];
rz(-2.8138748) q[3];
sx q[3];
rz(-0.98723427) q[3];
sx q[3];
rz(1.7879965) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.245529) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(-0.32972202) q[0];
rz(0.7715191) q[1];
sx q[1];
rz(-2.3604269) q[1];
sx q[1];
rz(-0.82569295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86733183) q[0];
sx q[0];
rz(-1.4837588) q[0];
sx q[0];
rz(-2.4768193) q[0];
x q[1];
rz(0.078455047) q[2];
sx q[2];
rz(-0.9758853) q[2];
sx q[2];
rz(1.9878146) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9553884) q[1];
sx q[1];
rz(-2.3729747) q[1];
sx q[1];
rz(-0.15665084) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7686033) q[3];
sx q[3];
rz(-1.8231492) q[3];
sx q[3];
rz(2.7276602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9416435) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(-1.6597718) q[2];
rz(3.0485349) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(-0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411497) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(0.39733091) q[0];
rz(0.98995248) q[1];
sx q[1];
rz(-1.7918469) q[1];
sx q[1];
rz(2.1053402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91401228) q[0];
sx q[0];
rz(-1.5591994) q[0];
sx q[0];
rz(0.40475233) q[0];
rz(1.431688) q[2];
sx q[2];
rz(-1.9204307) q[2];
sx q[2];
rz(-2.1099427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81779382) q[1];
sx q[1];
rz(-1.3518055) q[1];
sx q[1];
rz(1.6621672) q[1];
rz(-1.051722) q[3];
sx q[3];
rz(-1.241893) q[3];
sx q[3];
rz(2.2954706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1194666) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(1.8252581) q[2];
rz(-0.25257603) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0040358) q[0];
sx q[0];
rz(-0.79020774) q[0];
sx q[0];
rz(-0.23605119) q[0];
rz(0.099418489) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(1.258446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374219) q[0];
sx q[0];
rz(-1.6485456) q[0];
sx q[0];
rz(0.56044062) q[0];
x q[1];
rz(1.0706606) q[2];
sx q[2];
rz(-0.96958435) q[2];
sx q[2];
rz(-2.380033) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1061069) q[1];
sx q[1];
rz(-2.4508339) q[1];
sx q[1];
rz(-0.18045998) q[1];
rz(-0.40835898) q[3];
sx q[3];
rz(-1.9658739) q[3];
sx q[3];
rz(1.4407106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1666169) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(-2.100259) q[2];
rz(2.5552022) q[3];
sx q[3];
rz(-2.6643463) q[3];
sx q[3];
rz(-2.3311116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7753684) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(2.1495023) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(0.61983776) q[2];
sx q[2];
rz(-0.38417338) q[2];
sx q[2];
rz(-1.2585121) q[2];
rz(-2.548389) q[3];
sx q[3];
rz(-1.2013669) q[3];
sx q[3];
rz(-0.47040924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
