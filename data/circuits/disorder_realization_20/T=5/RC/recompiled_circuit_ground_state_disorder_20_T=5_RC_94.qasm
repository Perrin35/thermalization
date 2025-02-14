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
rz(0.33302745) q[0];
rz(2.0060519) q[1];
sx q[1];
rz(-0.82692868) q[1];
sx q[1];
rz(2.4976322) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0892093) q[0];
sx q[0];
rz(-2.174447) q[0];
sx q[0];
rz(2.8103845) q[0];
x q[1];
rz(-2.6690588) q[2];
sx q[2];
rz(-0.3344377) q[2];
sx q[2];
rz(0.57991934) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17929303) q[1];
sx q[1];
rz(-0.4518309) q[1];
sx q[1];
rz(-2.5746114) q[1];
x q[2];
rz(1.1584362) q[3];
sx q[3];
rz(-1.2500016) q[3];
sx q[3];
rz(-3.0123346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48646271) q[2];
sx q[2];
rz(-2.6772406) q[2];
sx q[2];
rz(-0.74074024) q[2];
rz(-1.5898534) q[3];
sx q[3];
rz(-0.71151763) q[3];
sx q[3];
rz(0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084259) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(-3.0246227) q[0];
rz(-0.50432694) q[1];
sx q[1];
rz(-1.802899) q[1];
sx q[1];
rz(0.28409827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11782538) q[0];
sx q[0];
rz(-1.4875011) q[0];
sx q[0];
rz(-2.1375594) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79757046) q[2];
sx q[2];
rz(-2.0013381) q[2];
sx q[2];
rz(-2.0883462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.423838) q[1];
sx q[1];
rz(-1.6053797) q[1];
sx q[1];
rz(-1.4486905) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8900364) q[3];
sx q[3];
rz(-1.3013892) q[3];
sx q[3];
rz(0.081307383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8154907) q[2];
sx q[2];
rz(-0.88560605) q[2];
sx q[2];
rz(-2.3808114) q[2];
rz(-2.271999) q[3];
sx q[3];
rz(-1.875501) q[3];
sx q[3];
rz(0.45309711) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075439) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(-2.889422) q[0];
rz(-1.699327) q[1];
sx q[1];
rz(-0.95528722) q[1];
sx q[1];
rz(-2.7853277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39922455) q[0];
sx q[0];
rz(-3.0816874) q[0];
sx q[0];
rz(1.674288) q[0];
rz(-pi) q[1];
rz(-0.54889955) q[2];
sx q[2];
rz(-1.9268914) q[2];
sx q[2];
rz(2.9402972) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1504476) q[1];
sx q[1];
rz(-2.4887062) q[1];
sx q[1];
rz(-1.3351403) q[1];
rz(-pi) q[2];
rz(-0.0019508501) q[3];
sx q[3];
rz(-1.3144799) q[3];
sx q[3];
rz(2.0511829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1189271) q[2];
sx q[2];
rz(-1.3902731) q[2];
sx q[2];
rz(1.6290132) q[2];
rz(-0.0096983612) q[3];
sx q[3];
rz(-1.2303979) q[3];
sx q[3];
rz(2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951204) q[0];
sx q[0];
rz(-2.5992114) q[0];
sx q[0];
rz(-2.8908308) q[0];
rz(1.7950902) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(-1.4432602) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6667056) q[0];
sx q[0];
rz(-2.3009818) q[0];
sx q[0];
rz(-0.29191125) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0403149) q[2];
sx q[2];
rz(-1.217739) q[2];
sx q[2];
rz(-0.7952035) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5233744) q[1];
sx q[1];
rz(-0.95444767) q[1];
sx q[1];
rz(2.2395654) q[1];
rz(-pi) q[2];
rz(2.8793427) q[3];
sx q[3];
rz(-1.3104386) q[3];
sx q[3];
rz(0.32338705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.5556339) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(2.8982758) q[2];
rz(-2.5569052) q[3];
sx q[3];
rz(-2.5644315) q[3];
sx q[3];
rz(-0.028133597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1481767) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(-2.962033) q[0];
rz(1.2252294) q[1];
sx q[1];
rz(-1.4045249) q[1];
sx q[1];
rz(-0.45825759) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1711463) q[0];
sx q[0];
rz(-1.2752646) q[0];
sx q[0];
rz(-1.4895205) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72821961) q[2];
sx q[2];
rz(-1.3698973) q[2];
sx q[2];
rz(-0.52853675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0974262) q[1];
sx q[1];
rz(-1.6454182) q[1];
sx q[1];
rz(1.4358372) q[1];
rz(-pi) q[2];
rz(-1.1510994) q[3];
sx q[3];
rz(-2.3956798) q[3];
sx q[3];
rz(-1.0585566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7563584) q[2];
sx q[2];
rz(-2.7172654) q[2];
sx q[2];
rz(2.2198086) q[2];
rz(-1.2515757) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-3.0475622) q[0];
sx q[0];
rz(-0.91218364) q[0];
sx q[0];
rz(0.14866522) q[0];
rz(2.8246763) q[1];
sx q[1];
rz(-2.6628559) q[1];
sx q[1];
rz(-2.5247578) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7148833) q[0];
sx q[0];
rz(-1.6322109) q[0];
sx q[0];
rz(1.4922569) q[0];
rz(0.89646879) q[2];
sx q[2];
rz(-0.54276641) q[2];
sx q[2];
rz(-1.3010058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4839263) q[1];
sx q[1];
rz(-2.1352508) q[1];
sx q[1];
rz(0.29740833) q[1];
x q[2];
rz(-2.0944164) q[3];
sx q[3];
rz(-1.3304119) q[3];
sx q[3];
rz(2.2464858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2769015) q[2];
sx q[2];
rz(-1.428705) q[2];
sx q[2];
rz(-3.1332664) q[2];
rz(2.5152123) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(-0.00055073784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(2.6595111) q[0];
rz(-3.112402) q[1];
sx q[1];
rz(-1.2960641) q[1];
sx q[1];
rz(-2.3775502) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2892492) q[0];
sx q[0];
rz(-1.2004392) q[0];
sx q[0];
rz(-0.072288805) q[0];
rz(0.025310658) q[2];
sx q[2];
rz(-0.9111852) q[2];
sx q[2];
rz(2.3940304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58043874) q[1];
sx q[1];
rz(-2.517546) q[1];
sx q[1];
rz(1.2777722) q[1];
rz(-pi) q[2];
rz(-0.0079844012) q[3];
sx q[3];
rz(-1.1698876) q[3];
sx q[3];
rz(0.29618357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46334106) q[2];
sx q[2];
rz(-1.3193139) q[2];
sx q[2];
rz(0.48416644) q[2];
rz(-0.68228996) q[3];
sx q[3];
rz(-2.4874918) q[3];
sx q[3];
rz(3.0068523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084455647) q[0];
sx q[0];
rz(-0.77780044) q[0];
sx q[0];
rz(-1.8498259) q[0];
rz(-0.46503398) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(-2.8344287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43873337) q[0];
sx q[0];
rz(-0.68638681) q[0];
sx q[0];
rz(2.2728361) q[0];
rz(-2.6382006) q[2];
sx q[2];
rz(-0.97676986) q[2];
sx q[2];
rz(2.4979748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.41514709) q[1];
sx q[1];
rz(-1.4142805) q[1];
sx q[1];
rz(-2.1464301) q[1];
rz(0.62509663) q[3];
sx q[3];
rz(-1.5543043) q[3];
sx q[3];
rz(1.0436637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9610567) q[2];
sx q[2];
rz(-0.44027105) q[2];
sx q[2];
rz(-0.72009909) q[2];
rz(2.2911206) q[3];
sx q[3];
rz(-1.974778) q[3];
sx q[3];
rz(-1.4115964) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9367323) q[0];
sx q[0];
rz(-2.9731049) q[0];
sx q[0];
rz(2.4334461) q[0];
rz(-1.6234966) q[1];
sx q[1];
rz(-0.93815362) q[1];
sx q[1];
rz(-2.785397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8339539) q[0];
sx q[0];
rz(-2.4705268) q[0];
sx q[0];
rz(2.4249581) q[0];
rz(-2.5237066) q[2];
sx q[2];
rz(-1.5949342) q[2];
sx q[2];
rz(1.9113505) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.041014) q[1];
sx q[1];
rz(-1.5407526) q[1];
sx q[1];
rz(0.64893367) q[1];
rz(-pi) q[2];
rz(2.735504) q[3];
sx q[3];
rz(-2.0022939) q[3];
sx q[3];
rz(3.114498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0922962) q[2];
sx q[2];
rz(-0.80731097) q[2];
sx q[2];
rz(0.88579196) q[2];
rz(0.93585912) q[3];
sx q[3];
rz(-1.9610619) q[3];
sx q[3];
rz(-0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43779272) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(1.6920775) q[0];
rz(-1.0225147) q[1];
sx q[1];
rz(-0.48373628) q[1];
sx q[1];
rz(-1.6922916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67774862) q[0];
sx q[0];
rz(-2.7263256) q[0];
sx q[0];
rz(1.4474611) q[0];
rz(0.42745356) q[2];
sx q[2];
rz(-1.9505672) q[2];
sx q[2];
rz(-2.3878271) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30561799) q[1];
sx q[1];
rz(-0.59594369) q[1];
sx q[1];
rz(1.1862331) q[1];
rz(-pi) q[2];
rz(-0.14808296) q[3];
sx q[3];
rz(-1.8335153) q[3];
sx q[3];
rz(-0.99276357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.25541043) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(-0.6748684) q[2];
rz(2.1700962) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(-1.6500047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(0.2324066) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(-2.7276298) q[2];
sx q[2];
rz(-1.0775177) q[2];
sx q[2];
rz(-2.1106363) q[2];
rz(0.078217004) q[3];
sx q[3];
rz(-2.297904) q[3];
sx q[3];
rz(-2.5301508) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
