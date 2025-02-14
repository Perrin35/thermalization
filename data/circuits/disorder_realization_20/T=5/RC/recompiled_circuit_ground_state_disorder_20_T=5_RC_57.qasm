OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5372758) q[0];
sx q[0];
rz(-0.24157) q[0];
sx q[0];
rz(-2.8085652) q[0];
rz(-1.1355407) q[1];
sx q[1];
rz(-2.314664) q[1];
sx q[1];
rz(0.64396042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0892093) q[0];
sx q[0];
rz(-2.174447) q[0];
sx q[0];
rz(0.33120819) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30006914) q[2];
sx q[2];
rz(-1.4208394) q[2];
sx q[2];
rz(1.4407002) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17929303) q[1];
sx q[1];
rz(-2.6897618) q[1];
sx q[1];
rz(2.5746114) q[1];
rz(-pi) q[2];
rz(1.1584362) q[3];
sx q[3];
rz(-1.2500016) q[3];
sx q[3];
rz(0.12925805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6551299) q[2];
sx q[2];
rz(-2.6772406) q[2];
sx q[2];
rz(2.4008524) q[2];
rz(-1.5898534) q[3];
sx q[3];
rz(-0.71151763) q[3];
sx q[3];
rz(0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7084259) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(-0.11696996) q[0];
rz(2.6372657) q[1];
sx q[1];
rz(-1.3386936) q[1];
sx q[1];
rz(2.8574944) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4000716) q[0];
sx q[0];
rz(-2.1353545) q[0];
sx q[0];
rz(-0.098640504) q[0];
rz(2.3440222) q[2];
sx q[2];
rz(-2.0013381) q[2];
sx q[2];
rz(-1.0532465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.423838) q[1];
sx q[1];
rz(-1.536213) q[1];
sx q[1];
rz(-1.4486905) q[1];
rz(1.8485214) q[3];
sx q[3];
rz(-1.3284995) q[3];
sx q[3];
rz(1.5577858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.32610193) q[2];
sx q[2];
rz(-2.2559866) q[2];
sx q[2];
rz(-2.3808114) q[2];
rz(-0.86959362) q[3];
sx q[3];
rz(-1.875501) q[3];
sx q[3];
rz(2.6884955) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075439) q[0];
sx q[0];
rz(-1.1203082) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(1.699327) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(-2.7853277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8460444) q[0];
sx q[0];
rz(-1.6303807) q[0];
sx q[0];
rz(-3.1353967) q[0];
x q[1];
rz(-2.5222657) q[2];
sx q[2];
rz(-0.64420036) q[2];
sx q[2];
rz(1.8875811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9911451) q[1];
sx q[1];
rz(-0.65288645) q[1];
sx q[1];
rz(-1.3351403) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8271133) q[3];
sx q[3];
rz(-1.5726834) q[3];
sx q[3];
rz(0.48088117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1189271) q[2];
sx q[2];
rz(-1.3902731) q[2];
sx q[2];
rz(-1.6290132) q[2];
rz(-0.0096983612) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(-2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.4951204) q[0];
sx q[0];
rz(-0.54238129) q[0];
sx q[0];
rz(-2.8908308) q[0];
rz(-1.3465025) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(1.6983324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8478717) q[0];
sx q[0];
rz(-1.3547055) q[0];
sx q[0];
rz(-2.3225075) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94130959) q[2];
sx q[2];
rz(-0.62773365) q[2];
sx q[2];
rz(2.8986487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52407284) q[1];
sx q[1];
rz(-2.1012602) q[1];
sx q[1];
rz(2.4071715) q[1];
x q[2];
rz(0.26224995) q[3];
sx q[3];
rz(-1.3104386) q[3];
sx q[3];
rz(-0.32338705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5859588) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(0.24331681) q[2];
rz(2.5569052) q[3];
sx q[3];
rz(-2.5644315) q[3];
sx q[3];
rz(0.028133597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1481767) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(-0.17955968) q[0];
rz(1.2252294) q[1];
sx q[1];
rz(-1.4045249) q[1];
sx q[1];
rz(-0.45825759) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4438547) q[0];
sx q[0];
rz(-0.30618822) q[0];
sx q[0];
rz(0.26060391) q[0];
rz(-pi) q[1];
rz(-2.413373) q[2];
sx q[2];
rz(-1.3698973) q[2];
sx q[2];
rz(0.52853675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0441664) q[1];
sx q[1];
rz(-1.4961745) q[1];
sx q[1];
rz(1.7057555) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86991258) q[3];
sx q[3];
rz(-1.8509838) q[3];
sx q[3];
rz(-2.3126569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3852343) q[2];
sx q[2];
rz(-2.7172654) q[2];
sx q[2];
rz(-0.9217841) q[2];
rz(-1.2515757) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(-0.061669953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09403041) q[0];
sx q[0];
rz(-0.91218364) q[0];
sx q[0];
rz(-2.9929274) q[0];
rz(-2.8246763) q[1];
sx q[1];
rz(-2.6628559) q[1];
sx q[1];
rz(2.5247578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4267093) q[0];
sx q[0];
rz(-1.5093818) q[0];
sx q[0];
rz(-1.6493357) q[0];
rz(-1.1304752) q[2];
sx q[2];
rz(-1.8991578) q[2];
sx q[2];
rz(0.87000123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65766636) q[1];
sx q[1];
rz(-2.1352508) q[1];
sx q[1];
rz(2.8441843) q[1];
x q[2];
rz(-1.0471763) q[3];
sx q[3];
rz(-1.8111808) q[3];
sx q[3];
rz(2.2464858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2769015) q[2];
sx q[2];
rz(-1.7128877) q[2];
sx q[2];
rz(-0.0083262715) q[2];
rz(2.5152123) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(-0.00055073784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(-0.48208153) q[0];
rz(-3.112402) q[1];
sx q[1];
rz(-1.8455285) q[1];
sx q[1];
rz(-0.76404244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8523435) q[0];
sx q[0];
rz(-1.9411534) q[0];
sx q[0];
rz(-0.072288805) q[0];
rz(1.5381728) q[2];
sx q[2];
rz(-0.66002405) q[2];
sx q[2];
rz(-0.78885022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.205207) q[1];
sx q[1];
rz(-0.97711414) q[1];
sx q[1];
rz(-2.9365345) q[1];
x q[2];
rz(-0.0079844012) q[3];
sx q[3];
rz(-1.1698876) q[3];
sx q[3];
rz(0.29618357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46334106) q[2];
sx q[2];
rz(-1.3193139) q[2];
sx q[2];
rz(2.6574262) q[2];
rz(-0.68228996) q[3];
sx q[3];
rz(-0.65410084) q[3];
sx q[3];
rz(-3.0068523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084455647) q[0];
sx q[0];
rz(-2.3637922) q[0];
sx q[0];
rz(-1.2917668) q[0];
rz(-0.46503398) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(-2.8344287) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7114176) q[0];
sx q[0];
rz(-1.1491551) q[0];
sx q[0];
rz(1.0118075) q[0];
rz(0.91395821) q[2];
sx q[2];
rz(-1.1595396) q[2];
sx q[2];
rz(0.62818254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41514709) q[1];
sx q[1];
rz(-1.7273121) q[1];
sx q[1];
rz(-0.99516258) q[1];
x q[2];
rz(0.028178111) q[3];
sx q[3];
rz(-2.5163076) q[3];
sx q[3];
rz(0.5042838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18053599) q[2];
sx q[2];
rz(-0.44027105) q[2];
sx q[2];
rz(-0.72009909) q[2];
rz(-2.2911206) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(-1.4115964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.20486031) q[0];
sx q[0];
rz(-0.16848773) q[0];
sx q[0];
rz(2.4334461) q[0];
rz(-1.5180961) q[1];
sx q[1];
rz(-0.93815362) q[1];
sx q[1];
rz(2.785397) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1461648) q[0];
sx q[0];
rz(-2.0588027) q[0];
sx q[0];
rz(-2.0515217) q[0];
rz(-3.0999423) q[2];
sx q[2];
rz(-0.61829582) q[2];
sx q[2];
rz(-0.30660812) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6345919) q[1];
sx q[1];
rz(-0.92220491) q[1];
sx q[1];
rz(1.533094) q[1];
rz(2.279874) q[3];
sx q[3];
rz(-0.58363014) q[3];
sx q[3];
rz(0.82635307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0922962) q[2];
sx q[2];
rz(-2.3342817) q[2];
sx q[2];
rz(-0.88579196) q[2];
rz(2.2057335) q[3];
sx q[3];
rz(-1.1805308) q[3];
sx q[3];
rz(-0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7037999) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(1.4495151) q[0];
rz(-2.119078) q[1];
sx q[1];
rz(-2.6578564) q[1];
sx q[1];
rz(1.449301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67774862) q[0];
sx q[0];
rz(-0.41526702) q[0];
sx q[0];
rz(1.4474611) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7141391) q[2];
sx q[2];
rz(-1.9505672) q[2];
sx q[2];
rz(-0.7537656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3811581) q[1];
sx q[1];
rz(-2.1179924) q[1];
sx q[1];
rz(-0.24914279) q[1];
rz(2.9935097) q[3];
sx q[3];
rz(-1.3080773) q[3];
sx q[3];
rz(0.99276357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.25541043) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(-0.6748684) q[2];
rz(0.97149649) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(-1.491588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7991199) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(-2.9091861) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(-0.41396285) q[2];
sx q[2];
rz(-2.064075) q[2];
sx q[2];
rz(1.0309564) q[2];
rz(-0.84216778) q[3];
sx q[3];
rz(-1.5123868) q[3];
sx q[3];
rz(2.2342891) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
