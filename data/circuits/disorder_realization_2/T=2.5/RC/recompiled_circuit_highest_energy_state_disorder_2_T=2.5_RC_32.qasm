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
rz(-0.19718117) q[0];
sx q[0];
rz(4.2137094) q[0];
sx q[0];
rz(11.022047) q[0];
rz(0.3577258) q[1];
sx q[1];
rz(-1.3074713) q[1];
sx q[1];
rz(0.37582418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.887325) q[0];
sx q[0];
rz(-1.4987117) q[0];
sx q[0];
rz(-2.3743126) q[0];
x q[1];
rz(-2.69015) q[2];
sx q[2];
rz(-0.57381781) q[2];
sx q[2];
rz(0.30893477) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5615422) q[1];
sx q[1];
rz(-2.1218178) q[1];
sx q[1];
rz(-2.2660794) q[1];
rz(-pi) q[2];
rz(2.8055111) q[3];
sx q[3];
rz(-2.1770526) q[3];
sx q[3];
rz(-1.1949902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69349217) q[2];
sx q[2];
rz(-2.7585612) q[2];
sx q[2];
rz(0.4105655) q[2];
rz(0.52452123) q[3];
sx q[3];
rz(-2.9295242) q[3];
sx q[3];
rz(-1.1493523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9583424) q[0];
sx q[0];
rz(-2.9999833) q[0];
sx q[0];
rz(2.435922) q[0];
rz(-1.8535829) q[1];
sx q[1];
rz(-0.57669222) q[1];
sx q[1];
rz(-1.1340595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41578226) q[0];
sx q[0];
rz(-2.232612) q[0];
sx q[0];
rz(-0.51915581) q[0];
rz(-pi) q[1];
rz(-0.4749078) q[2];
sx q[2];
rz(-2.7596052) q[2];
sx q[2];
rz(-1.8595921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73423656) q[1];
sx q[1];
rz(-2.3132965) q[1];
sx q[1];
rz(-2.616859) q[1];
rz(-pi) q[2];
rz(2.5763578) q[3];
sx q[3];
rz(-1.4197403) q[3];
sx q[3];
rz(-0.24542576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80295339) q[2];
sx q[2];
rz(-2.6889375) q[2];
sx q[2];
rz(-0.53042859) q[2];
rz(0.28702921) q[3];
sx q[3];
rz(-0.48873264) q[3];
sx q[3];
rz(-2.4052525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0158143) q[0];
sx q[0];
rz(-2.0398572) q[0];
sx q[0];
rz(-1.0544448) q[0];
rz(-1.6619445) q[1];
sx q[1];
rz(-2.2190861) q[1];
sx q[1];
rz(2.0747917) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9196292) q[0];
sx q[0];
rz(-2.9502075) q[0];
sx q[0];
rz(1.3758977) q[0];
rz(-pi) q[1];
rz(1.4786903) q[2];
sx q[2];
rz(-2.2450441) q[2];
sx q[2];
rz(0.2052274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.57745614) q[1];
sx q[1];
rz(-1.3423663) q[1];
sx q[1];
rz(-1.0962568) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9122804) q[3];
sx q[3];
rz(-0.74401281) q[3];
sx q[3];
rz(0.0014425576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58831424) q[2];
sx q[2];
rz(-1.8411627) q[2];
sx q[2];
rz(1.1021357) q[2];
rz(0.45840248) q[3];
sx q[3];
rz(-1.6130684) q[3];
sx q[3];
rz(0.63157356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44963968) q[0];
sx q[0];
rz(-2.3821558) q[0];
sx q[0];
rz(0.48680437) q[0];
rz(-1.5454769) q[1];
sx q[1];
rz(-2.7524452) q[1];
sx q[1];
rz(0.62852377) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31159831) q[0];
sx q[0];
rz(-1.4382558) q[0];
sx q[0];
rz(-1.9674942) q[0];
rz(-pi) q[1];
rz(-0.99765473) q[2];
sx q[2];
rz(-1.4467244) q[2];
sx q[2];
rz(2.5981552) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6316136) q[1];
sx q[1];
rz(-1.9603131) q[1];
sx q[1];
rz(-0.20304598) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0699203) q[3];
sx q[3];
rz(-1.023277) q[3];
sx q[3];
rz(-1.5231093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57807606) q[2];
sx q[2];
rz(-0.11512828) q[2];
sx q[2];
rz(-0.76187491) q[2];
rz(0.21841194) q[3];
sx q[3];
rz(-0.29774791) q[3];
sx q[3];
rz(-1.0485605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8258284) q[0];
sx q[0];
rz(-0.35218069) q[0];
sx q[0];
rz(0.6045565) q[0];
rz(1.8479337) q[1];
sx q[1];
rz(-0.81984055) q[1];
sx q[1];
rz(0.34436071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47681706) q[0];
sx q[0];
rz(-1.6533543) q[0];
sx q[0];
rz(3.1212016) q[0];
x q[1];
rz(1.2372962) q[2];
sx q[2];
rz(-0.97003675) q[2];
sx q[2];
rz(1.5432026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1923848) q[1];
sx q[1];
rz(-0.65509641) q[1];
sx q[1];
rz(-0.75836436) q[1];
rz(-pi) q[2];
rz(-0.46724288) q[3];
sx q[3];
rz(-2.1911216) q[3];
sx q[3];
rz(-1.0091227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0209501) q[2];
sx q[2];
rz(-0.63434333) q[2];
sx q[2];
rz(2.9316588) q[2];
rz(2.0338992) q[3];
sx q[3];
rz(-1.6898797) q[3];
sx q[3];
rz(-2.8797188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8629859) q[0];
sx q[0];
rz(-0.59026533) q[0];
sx q[0];
rz(-0.069409542) q[0];
rz(-3.1076) q[1];
sx q[1];
rz(-2.4885204) q[1];
sx q[1];
rz(-0.45190826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6959608) q[0];
sx q[0];
rz(-1.9602141) q[0];
sx q[0];
rz(-1.3971863) q[0];
rz(-pi) q[1];
x q[1];
rz(2.185427) q[2];
sx q[2];
rz(-2.1611804) q[2];
sx q[2];
rz(-2.0942795) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76594162) q[1];
sx q[1];
rz(-0.022863511) q[1];
sx q[1];
rz(1.2506865) q[1];
x q[2];
rz(2.8699401) q[3];
sx q[3];
rz(-1.8166887) q[3];
sx q[3];
rz(3.0761443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41428337) q[2];
sx q[2];
rz(-2.3392129) q[2];
sx q[2];
rz(-1.1451954) q[2];
rz(-0.9582054) q[3];
sx q[3];
rz(-1.2088135) q[3];
sx q[3];
rz(0.27829471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77022839) q[0];
sx q[0];
rz(-0.48786491) q[0];
sx q[0];
rz(0.86139739) q[0];
rz(0.96249145) q[1];
sx q[1];
rz(-0.95314127) q[1];
sx q[1];
rz(-2.9846587) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033534) q[0];
sx q[0];
rz(-0.80488759) q[0];
sx q[0];
rz(2.0867086) q[0];
rz(-pi) q[1];
rz(0.10617039) q[2];
sx q[2];
rz(-1.2154801) q[2];
sx q[2];
rz(0.67007845) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90060252) q[1];
sx q[1];
rz(-2.1442736) q[1];
sx q[1];
rz(-2.1585224) q[1];
x q[2];
rz(1.7926927) q[3];
sx q[3];
rz(-2.923344) q[3];
sx q[3];
rz(1.143592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.018205) q[2];
sx q[2];
rz(-0.79424477) q[2];
sx q[2];
rz(2.0920853) q[2];
rz(-2.6617995) q[3];
sx q[3];
rz(-0.19194651) q[3];
sx q[3];
rz(0.16201365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4341693) q[0];
sx q[0];
rz(-2.8987085) q[0];
sx q[0];
rz(2.6450787) q[0];
rz(-1.6560582) q[1];
sx q[1];
rz(-0.86692202) q[1];
sx q[1];
rz(-0.022580126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88779792) q[0];
sx q[0];
rz(-0.27844515) q[0];
sx q[0];
rz(0.87775567) q[0];
rz(2.0099925) q[2];
sx q[2];
rz(-2.9832057) q[2];
sx q[2];
rz(-1.2513551) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1068017) q[1];
sx q[1];
rz(-0.80134922) q[1];
sx q[1];
rz(-1.31557) q[1];
rz(-pi) q[2];
rz(-1.0362421) q[3];
sx q[3];
rz(-1.4009985) q[3];
sx q[3];
rz(1.9056729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5217343) q[2];
sx q[2];
rz(-1.1405742) q[2];
sx q[2];
rz(-0.015259585) q[2];
rz(0.36733019) q[3];
sx q[3];
rz(-0.44706774) q[3];
sx q[3];
rz(-0.36146155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5527282) q[0];
sx q[0];
rz(-2.9155154) q[0];
sx q[0];
rz(-3.0638301) q[0];
rz(-0.11847682) q[1];
sx q[1];
rz(-0.34093726) q[1];
sx q[1];
rz(-2.4854787) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3518302) q[0];
sx q[0];
rz(-1.2675722) q[0];
sx q[0];
rz(3.0089761) q[0];
rz(0.17123789) q[2];
sx q[2];
rz(-1.3067553) q[2];
sx q[2];
rz(-2.03351) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39086884) q[1];
sx q[1];
rz(-1.0926529) q[1];
sx q[1];
rz(0.45420809) q[1];
x q[2];
rz(-0.55226294) q[3];
sx q[3];
rz(-1.1202462) q[3];
sx q[3];
rz(2.0678064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6180827) q[2];
sx q[2];
rz(-1.9205576) q[2];
sx q[2];
rz(2.6005884) q[2];
rz(-2.6909761) q[3];
sx q[3];
rz(-1.6559947) q[3];
sx q[3];
rz(-2.165152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4377874) q[0];
sx q[0];
rz(-1.7923651) q[0];
sx q[0];
rz(3.034814) q[0];
rz(1.5497426) q[1];
sx q[1];
rz(-0.50250643) q[1];
sx q[1];
rz(1.7639311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4431245) q[0];
sx q[0];
rz(-1.5266935) q[0];
sx q[0];
rz(-1.1209784) q[0];
rz(-pi) q[1];
rz(-1.7936034) q[2];
sx q[2];
rz(-1.2651332) q[2];
sx q[2];
rz(-0.7204186) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71860524) q[1];
sx q[1];
rz(-1.3235301) q[1];
sx q[1];
rz(2.894677) q[1];
rz(0.53394188) q[3];
sx q[3];
rz(-2.0755872) q[3];
sx q[3];
rz(-1.8478888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0489073) q[2];
sx q[2];
rz(-1.4212298) q[2];
sx q[2];
rz(-2.8468724) q[2];
rz(-0.60485351) q[3];
sx q[3];
rz(-1.9113144) q[3];
sx q[3];
rz(0.10943432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7034364) q[0];
sx q[0];
rz(-1.823517) q[0];
sx q[0];
rz(2.1448268) q[0];
rz(-3.1238212) q[1];
sx q[1];
rz(-1.8780864) q[1];
sx q[1];
rz(-1.3395739) q[1];
rz(0.80452917) q[2];
sx q[2];
rz(-2.011841) q[2];
sx q[2];
rz(2.9622215) q[2];
rz(0.83497077) q[3];
sx q[3];
rz(-2.2899262) q[3];
sx q[3];
rz(-0.038577608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
