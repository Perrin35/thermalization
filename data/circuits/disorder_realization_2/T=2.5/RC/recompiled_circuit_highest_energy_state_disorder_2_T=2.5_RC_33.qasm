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
rz(2.9444115) q[0];
sx q[0];
rz(-1.0721167) q[0];
sx q[0];
rz(1.5443235) q[0];
rz(-2.7838669) q[1];
sx q[1];
rz(-1.8341213) q[1];
sx q[1];
rz(2.7657685) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25426769) q[0];
sx q[0];
rz(-1.642881) q[0];
sx q[0];
rz(2.3743126) q[0];
rz(-pi) q[1];
rz(-0.52679072) q[2];
sx q[2];
rz(-1.3317022) q[2];
sx q[2];
rz(-1.4930489) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5800505) q[1];
sx q[1];
rz(-2.1218178) q[1];
sx q[1];
rz(-0.8755133) q[1];
rz(0.9373215) q[3];
sx q[3];
rz(-1.8452462) q[3];
sx q[3];
rz(2.9622674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4481005) q[2];
sx q[2];
rz(-2.7585612) q[2];
sx q[2];
rz(0.4105655) q[2];
rz(2.6170714) q[3];
sx q[3];
rz(-2.9295242) q[3];
sx q[3];
rz(-1.9922403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.18325026) q[0];
sx q[0];
rz(-2.9999833) q[0];
sx q[0];
rz(0.70567065) q[0];
rz(1.8535829) q[1];
sx q[1];
rz(-2.5649004) q[1];
sx q[1];
rz(-1.1340595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7258104) q[0];
sx q[0];
rz(-2.232612) q[0];
sx q[0];
rz(0.51915581) q[0];
x q[1];
rz(-0.34313029) q[2];
sx q[2];
rz(-1.742082) q[2];
sx q[2];
rz(2.4076381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6993076) q[1];
sx q[1];
rz(-0.87934994) q[1];
sx q[1];
rz(-2.0705018) q[1];
rz(1.7491327) q[3];
sx q[3];
rz(-1.0127676) q[3];
sx q[3];
rz(-1.7210646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80295339) q[2];
sx q[2];
rz(-0.45265517) q[2];
sx q[2];
rz(-2.6111641) q[2];
rz(0.28702921) q[3];
sx q[3];
rz(-0.48873264) q[3];
sx q[3];
rz(0.73634017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257783) q[0];
sx q[0];
rz(-1.1017355) q[0];
sx q[0];
rz(-1.0544448) q[0];
rz(1.6619445) q[1];
sx q[1];
rz(-2.2190861) q[1];
sx q[1];
rz(1.066801) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9196292) q[0];
sx q[0];
rz(-2.9502075) q[0];
sx q[0];
rz(1.765695) q[0];
rz(1.4786903) q[2];
sx q[2];
rz(-0.89654857) q[2];
sx q[2];
rz(-0.2052274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7330526) q[1];
sx q[1];
rz(-2.618737) q[1];
sx q[1];
rz(-2.0414635) q[1];
rz(-pi) q[2];
rz(-1.2293123) q[3];
sx q[3];
rz(-0.74401281) q[3];
sx q[3];
rz(3.1401501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5532784) q[2];
sx q[2];
rz(-1.8411627) q[2];
sx q[2];
rz(2.039457) q[2];
rz(2.6831902) q[3];
sx q[3];
rz(-1.6130684) q[3];
sx q[3];
rz(-0.63157356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.691953) q[0];
sx q[0];
rz(-0.75943685) q[0];
sx q[0];
rz(-2.6547883) q[0];
rz(1.5961157) q[1];
sx q[1];
rz(-0.38914746) q[1];
sx q[1];
rz(2.5130689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31159831) q[0];
sx q[0];
rz(-1.7033368) q[0];
sx q[0];
rz(-1.9674942) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99765473) q[2];
sx q[2];
rz(-1.6948683) q[2];
sx q[2];
rz(-2.5981552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0027517) q[1];
sx q[1];
rz(-1.758444) q[1];
sx q[1];
rz(-1.9676588) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0221338) q[3];
sx q[3];
rz(-1.5096153) q[3];
sx q[3];
rz(3.1312628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57807606) q[2];
sx q[2];
rz(-3.0264644) q[2];
sx q[2];
rz(0.76187491) q[2];
rz(2.9231807) q[3];
sx q[3];
rz(-0.29774791) q[3];
sx q[3];
rz(-2.0930321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4223347) q[0];
sx q[0];
rz(-3.0565593) q[0];
sx q[0];
rz(-1.3291977) q[0];
rz(-1.2372962) q[2];
sx q[2];
rz(-2.1715559) q[2];
sx q[2];
rz(-1.5983901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3185259) q[1];
sx q[1];
rz(-1.112655) q[1];
sx q[1];
rz(1.0847102) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8959534) q[3];
sx q[3];
rz(-1.1955441) q[3];
sx q[3];
rz(0.2764052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0209501) q[2];
sx q[2];
rz(-2.5072493) q[2];
sx q[2];
rz(-0.20993385) q[2];
rz(2.0338992) q[3];
sx q[3];
rz(-1.6898797) q[3];
sx q[3];
rz(0.2618739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27860677) q[0];
sx q[0];
rz(-0.59026533) q[0];
sx q[0];
rz(-0.069409542) q[0];
rz(3.1076) q[1];
sx q[1];
rz(-0.65307224) q[1];
sx q[1];
rz(2.6896844) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87837093) q[0];
sx q[0];
rz(-2.7170406) q[0];
sx q[0];
rz(2.7431737) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4546409) q[2];
sx q[2];
rz(-1.0712306) q[2];
sx q[2];
rz(2.9925008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44575366) q[1];
sx q[1];
rz(-1.5924982) q[1];
sx q[1];
rz(3.134397) q[1];
rz(0.27165256) q[3];
sx q[3];
rz(-1.324904) q[3];
sx q[3];
rz(3.0761443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7273093) q[2];
sx q[2];
rz(-2.3392129) q[2];
sx q[2];
rz(-1.9963973) q[2];
rz(2.1833873) q[3];
sx q[3];
rz(-1.2088135) q[3];
sx q[3];
rz(0.27829471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3713643) q[0];
sx q[0];
rz(-2.6537277) q[0];
sx q[0];
rz(0.86139739) q[0];
rz(-2.1791012) q[1];
sx q[1];
rz(-0.95314127) q[1];
sx q[1];
rz(0.15693396) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033534) q[0];
sx q[0];
rz(-2.3367051) q[0];
sx q[0];
rz(1.0548841) q[0];
rz(-1.2136344) q[2];
sx q[2];
rz(-1.6703122) q[2];
sx q[2];
rz(-0.9377756) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2409901) q[1];
sx q[1];
rz(-2.1442736) q[1];
sx q[1];
rz(-2.1585224) q[1];
x q[2];
rz(1.3488999) q[3];
sx q[3];
rz(-2.923344) q[3];
sx q[3];
rz(1.9980007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.018205) q[2];
sx q[2];
rz(-2.3473479) q[2];
sx q[2];
rz(-1.0495074) q[2];
rz(2.6617995) q[3];
sx q[3];
rz(-2.9496461) q[3];
sx q[3];
rz(-2.979579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70742339) q[0];
sx q[0];
rz(-0.24288414) q[0];
sx q[0];
rz(2.6450787) q[0];
rz(-1.6560582) q[1];
sx q[1];
rz(-2.2746706) q[1];
sx q[1];
rz(0.022580126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.784784) q[0];
sx q[0];
rz(-1.7473146) q[0];
sx q[0];
rz(-1.7872732) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4272251) q[2];
sx q[2];
rz(-1.6379135) q[2];
sx q[2];
rz(0.75380933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39358006) q[1];
sx q[1];
rz(-0.8023387) q[1];
sx q[1];
rz(-0.2549766) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2461698) q[3];
sx q[3];
rz(-0.55836779) q[3];
sx q[3];
rz(-0.61290395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5217343) q[2];
sx q[2];
rz(-1.1405742) q[2];
sx q[2];
rz(3.1263331) q[2];
rz(-0.36733019) q[3];
sx q[3];
rz(-2.6945249) q[3];
sx q[3];
rz(2.7801311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5527282) q[0];
sx q[0];
rz(-2.9155154) q[0];
sx q[0];
rz(3.0638301) q[0];
rz(3.0231158) q[1];
sx q[1];
rz(-0.34093726) q[1];
sx q[1];
rz(-2.4854787) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25877798) q[0];
sx q[0];
rz(-1.4442632) q[0];
sx q[0];
rz(-1.8765429) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3030197) q[2];
sx q[2];
rz(-1.736044) q[2];
sx q[2];
rz(2.6337773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9356797) q[1];
sx q[1];
rz(-2.4944843) q[1];
sx q[1];
rz(-2.2733105) q[1];
rz(-pi) q[2];
rz(-2.3967152) q[3];
sx q[3];
rz(-0.69760579) q[3];
sx q[3];
rz(-3.0234203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52350998) q[2];
sx q[2];
rz(-1.221035) q[2];
sx q[2];
rz(0.54100424) q[2];
rz(2.6909761) q[3];
sx q[3];
rz(-1.485598) q[3];
sx q[3];
rz(-2.165152) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4377874) q[0];
sx q[0];
rz(-1.7923651) q[0];
sx q[0];
rz(-3.034814) q[0];
rz(1.5497426) q[1];
sx q[1];
rz(-0.50250643) q[1];
sx q[1];
rz(-1.3776616) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6984682) q[0];
sx q[0];
rz(-1.6148991) q[0];
sx q[0];
rz(2.0206142) q[0];
rz(-0.3129199) q[2];
sx q[2];
rz(-1.358479) q[2];
sx q[2];
rz(2.2231399) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.351016) q[1];
sx q[1];
rz(-1.3315397) q[1];
sx q[1];
rz(1.3161216) q[1];
rz(-pi) q[2];
rz(2.1414457) q[3];
sx q[3];
rz(-1.1091145) q[3];
sx q[3];
rz(-3.1401855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0926853) q[2];
sx q[2];
rz(-1.4212298) q[2];
sx q[2];
rz(2.8468724) q[2];
rz(-0.60485351) q[3];
sx q[3];
rz(-1.9113144) q[3];
sx q[3];
rz(0.10943432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43815628) q[0];
sx q[0];
rz(-1.3180757) q[0];
sx q[0];
rz(-0.99676589) q[0];
rz(-3.1238212) q[1];
sx q[1];
rz(-1.8780864) q[1];
sx q[1];
rz(-1.3395739) q[1];
rz(-2.1684767) q[2];
sx q[2];
rz(-2.2804307) q[2];
sx q[2];
rz(-1.3326999) q[2];
rz(0.65405398) q[3];
sx q[3];
rz(-2.1624343) q[3];
sx q[3];
rz(2.1611675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
