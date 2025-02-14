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
rz(1.1178782) q[0];
sx q[0];
rz(-2.0562545) q[0];
sx q[0];
rz(0.35882741) q[0];
rz(-1.8817512) q[1];
sx q[1];
rz(-2.0460195) q[1];
sx q[1];
rz(-2.1020558) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22408501) q[0];
sx q[0];
rz(-1.9244902) q[0];
sx q[0];
rz(-0.68266305) q[0];
rz(1.5284003) q[2];
sx q[2];
rz(-1.5982407) q[2];
sx q[2];
rz(-1.8191486) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7233091) q[1];
sx q[1];
rz(-1.2110577) q[1];
sx q[1];
rz(-0.32029256) q[1];
rz(-1.667439) q[3];
sx q[3];
rz(-2.2386239) q[3];
sx q[3];
rz(-2.5092595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.871668) q[2];
sx q[2];
rz(-0.85447997) q[2];
sx q[2];
rz(-0.16647767) q[2];
rz(-3.0050333) q[3];
sx q[3];
rz(-0.86524335) q[3];
sx q[3];
rz(1.4599919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0110375) q[0];
sx q[0];
rz(-0.32259536) q[0];
sx q[0];
rz(0.36001298) q[0];
rz(0.75045466) q[1];
sx q[1];
rz(-2.2513335) q[1];
sx q[1];
rz(3.0976565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.982807) q[0];
sx q[0];
rz(-1.4318529) q[0];
sx q[0];
rz(1.7278633) q[0];
rz(-pi) q[1];
rz(-3.0925748) q[2];
sx q[2];
rz(-0.92713461) q[2];
sx q[2];
rz(1.1903534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39831623) q[1];
sx q[1];
rz(-2.3651402) q[1];
sx q[1];
rz(-2.8321596) q[1];
rz(-pi) q[2];
x q[2];
rz(2.206541) q[3];
sx q[3];
rz(-0.32554829) q[3];
sx q[3];
rz(-2.9071484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7875849) q[2];
sx q[2];
rz(-0.44469357) q[2];
sx q[2];
rz(-1.0178052) q[2];
rz(-2.9338037) q[3];
sx q[3];
rz(-2.5354247) q[3];
sx q[3];
rz(0.35473216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6924234) q[0];
sx q[0];
rz(-2.848931) q[0];
sx q[0];
rz(-1.0798651) q[0];
rz(1.0049459) q[1];
sx q[1];
rz(-2.4506073) q[1];
sx q[1];
rz(1.2452004) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5678788) q[0];
sx q[0];
rz(-1.2612997) q[0];
sx q[0];
rz(1.482016) q[0];
rz(1.7342308) q[2];
sx q[2];
rz(-1.15112) q[2];
sx q[2];
rz(0.070022665) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1221262) q[1];
sx q[1];
rz(-1.2858741) q[1];
sx q[1];
rz(2.1287005) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9726282) q[3];
sx q[3];
rz(-0.95674435) q[3];
sx q[3];
rz(-1.799063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24091844) q[2];
sx q[2];
rz(-1.2718879) q[2];
sx q[2];
rz(2.7980878) q[2];
rz(3.1268696) q[3];
sx q[3];
rz(-0.69535178) q[3];
sx q[3];
rz(-1.2408313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13046509) q[0];
sx q[0];
rz(-1.0437597) q[0];
sx q[0];
rz(-2.4567132) q[0];
rz(-0.20826805) q[1];
sx q[1];
rz(-1.3083649) q[1];
sx q[1];
rz(0.79316717) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356157) q[0];
sx q[0];
rz(-0.55287939) q[0];
sx q[0];
rz(1.7589491) q[0];
x q[1];
rz(2.2440282) q[2];
sx q[2];
rz(-2.8688768) q[2];
sx q[2];
rz(0.9630335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85911688) q[1];
sx q[1];
rz(-2.1819677) q[1];
sx q[1];
rz(3.1219367) q[1];
rz(1.9408405) q[3];
sx q[3];
rz(-1.368282) q[3];
sx q[3];
rz(0.56087342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3785279) q[2];
sx q[2];
rz(-1.6918162) q[2];
sx q[2];
rz(2.6587963) q[2];
rz(1.5040846) q[3];
sx q[3];
rz(-0.62862325) q[3];
sx q[3];
rz(2.8211236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046431635) q[0];
sx q[0];
rz(-0.23155364) q[0];
sx q[0];
rz(0.10145536) q[0];
rz(2.9265112) q[1];
sx q[1];
rz(-1.2098034) q[1];
sx q[1];
rz(-1.2717517) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5520855) q[0];
sx q[0];
rz(-0.55098039) q[0];
sx q[0];
rz(1.0066342) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4319806) q[2];
sx q[2];
rz(-0.73118787) q[2];
sx q[2];
rz(-3.0572011) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78465652) q[1];
sx q[1];
rz(-1.1788834) q[1];
sx q[1];
rz(2.980113) q[1];
rz(2.1121096) q[3];
sx q[3];
rz(-0.39697638) q[3];
sx q[3];
rz(3.027806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7868906) q[2];
sx q[2];
rz(-0.44822794) q[2];
sx q[2];
rz(0.1498214) q[2];
rz(0.24770728) q[3];
sx q[3];
rz(-0.37023538) q[3];
sx q[3];
rz(0.80055922) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.596452) q[0];
sx q[0];
rz(-0.62262028) q[0];
sx q[0];
rz(-1.2860292) q[0];
rz(-2.3226358) q[1];
sx q[1];
rz(-0.30636925) q[1];
sx q[1];
rz(0.18316306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4990627) q[0];
sx q[0];
rz(-2.1234491) q[0];
sx q[0];
rz(-1.8091317) q[0];
rz(-pi) q[1];
rz(0.18503616) q[2];
sx q[2];
rz(-1.5079947) q[2];
sx q[2];
rz(2.7887995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0350041) q[1];
sx q[1];
rz(-0.87751203) q[1];
sx q[1];
rz(-2.6412852) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8882033) q[3];
sx q[3];
rz(-1.0930233) q[3];
sx q[3];
rz(-1.0047411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5554819) q[2];
sx q[2];
rz(-1.7038944) q[2];
sx q[2];
rz(-0.5371896) q[2];
rz(1.4026583) q[3];
sx q[3];
rz(-1.2638998) q[3];
sx q[3];
rz(-0.62854832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1527767) q[0];
sx q[0];
rz(-2.7429136) q[0];
sx q[0];
rz(0.11845778) q[0];
rz(1.1064103) q[1];
sx q[1];
rz(-1.2245919) q[1];
sx q[1];
rz(-2.4647663) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22960104) q[0];
sx q[0];
rz(-1.4957447) q[0];
sx q[0];
rz(1.6261392) q[0];
rz(-2.9666128) q[2];
sx q[2];
rz(-2.4799441) q[2];
sx q[2];
rz(1.511919) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0829777) q[1];
sx q[1];
rz(-1.737002) q[1];
sx q[1];
rz(0.18402255) q[1];
rz(1.6366129) q[3];
sx q[3];
rz(-3.1269347) q[3];
sx q[3];
rz(2.3848118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38830882) q[2];
sx q[2];
rz(-0.94677418) q[2];
sx q[2];
rz(3.0906265) q[2];
rz(0.1554337) q[3];
sx q[3];
rz(-1.4916462) q[3];
sx q[3];
rz(1.0203993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8841356) q[0];
sx q[0];
rz(-1.9866885) q[0];
sx q[0];
rz(-2.1575523) q[0];
rz(-2.6718196) q[1];
sx q[1];
rz(-1.677294) q[1];
sx q[1];
rz(-2.9205172) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4656847) q[0];
sx q[0];
rz(-3.0455112) q[0];
sx q[0];
rz(-0.19739993) q[0];
x q[1];
rz(1.9479772) q[2];
sx q[2];
rz(-2.1862891) q[2];
sx q[2];
rz(-2.7831558) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.03121491) q[1];
sx q[1];
rz(-1.3450012) q[1];
sx q[1];
rz(-1.8757191) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4793977) q[3];
sx q[3];
rz(-2.0948296) q[3];
sx q[3];
rz(-2.9063922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85349637) q[2];
sx q[2];
rz(-0.6820389) q[2];
sx q[2];
rz(1.811553) q[2];
rz(-2.5174777) q[3];
sx q[3];
rz(-0.95804405) q[3];
sx q[3];
rz(-2.7522855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57988155) q[0];
sx q[0];
rz(-1.4917253) q[0];
sx q[0];
rz(1.0505744) q[0];
rz(-2.1613278) q[1];
sx q[1];
rz(-2.0591718) q[1];
sx q[1];
rz(1.6260653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78895348) q[0];
sx q[0];
rz(-1.727956) q[0];
sx q[0];
rz(0.53257863) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8063861) q[2];
sx q[2];
rz(-2.0353999) q[2];
sx q[2];
rz(-0.36201477) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5055464) q[1];
sx q[1];
rz(-2.6719173) q[1];
sx q[1];
rz(-2.4081025) q[1];
rz(-1.100176) q[3];
sx q[3];
rz(-1.017316) q[3];
sx q[3];
rz(0.054892232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76422292) q[2];
sx q[2];
rz(-1.919596) q[2];
sx q[2];
rz(-0.73966217) q[2];
rz(-0.1554648) q[3];
sx q[3];
rz(-0.37853095) q[3];
sx q[3];
rz(2.9445061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71843201) q[0];
sx q[0];
rz(-0.53373706) q[0];
sx q[0];
rz(-2.0045795) q[0];
rz(-0.63277376) q[1];
sx q[1];
rz(-2.4686765) q[1];
sx q[1];
rz(-0.64579642) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.150313) q[0];
sx q[0];
rz(-1.8984573) q[0];
sx q[0];
rz(1.3929266) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77984758) q[2];
sx q[2];
rz(-1.4044242) q[2];
sx q[2];
rz(-0.8188687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.551522) q[1];
sx q[1];
rz(-2.3209718) q[1];
sx q[1];
rz(-0.094610081) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8237173) q[3];
sx q[3];
rz(-1.0680305) q[3];
sx q[3];
rz(-1.1989145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48468963) q[2];
sx q[2];
rz(-0.36269665) q[2];
sx q[2];
rz(3.1032069) q[2];
rz(3.0302327) q[3];
sx q[3];
rz(-0.42698082) q[3];
sx q[3];
rz(-0.67489433) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4511694) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.3896598) q[0];
rz(1.2870862) q[1];
sx q[1];
rz(-0.99936395) q[1];
sx q[1];
rz(-2.0157464) q[1];
rz(2.1846175) q[2];
sx q[2];
rz(-2.0718832) q[2];
sx q[2];
rz(2.8399443) q[2];
rz(0.75931924) q[3];
sx q[3];
rz(-1.6990468) q[3];
sx q[3];
rz(2.852462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
