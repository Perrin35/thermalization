OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(2.0895045) q[0];
sx q[0];
rz(10.917561) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(2.066943) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91905347) q[0];
sx q[0];
rz(-0.018964501) q[0];
sx q[0];
rz(-0.28199621) q[0];
rz(-pi) q[1];
rz(-0.15901788) q[2];
sx q[2];
rz(-0.79610014) q[2];
sx q[2];
rz(-0.12479347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5584471) q[1];
sx q[1];
rz(-1.4902643) q[1];
sx q[1];
rz(-3.0196378) q[1];
x q[2];
rz(-2.6082637) q[3];
sx q[3];
rz(-1.1733857) q[3];
sx q[3];
rz(-1.0226137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(-1.1734022) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(0.2127969) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56869498) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(-2.9808295) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(-1.6404023) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6340856) q[0];
sx q[0];
rz(-1.5459783) q[0];
sx q[0];
rz(0.12269845) q[0];
x q[1];
rz(-2.8295849) q[2];
sx q[2];
rz(-1.1216251) q[2];
sx q[2];
rz(-2.1829407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7306819) q[1];
sx q[1];
rz(-2.4446179) q[1];
sx q[1];
rz(-2.8721786) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5853911) q[3];
sx q[3];
rz(-0.39108927) q[3];
sx q[3];
rz(0.99816834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(-1.5141053) q[2];
rz(-1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.5037781) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30341879) q[0];
sx q[0];
rz(-1.6359328) q[0];
sx q[0];
rz(1.6584048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90942803) q[2];
sx q[2];
rz(-0.69790188) q[2];
sx q[2];
rz(0.94240377) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9864195) q[1];
sx q[1];
rz(-1.6134312) q[1];
sx q[1];
rz(2.8405872) q[1];
rz(0.13111968) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(-2.2281447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(0.25742325) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96823111) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(2.0898576) q[0];
rz(2.5412718) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7950492) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(1.6853149) q[0];
rz(-pi) q[1];
rz(2.7096268) q[2];
sx q[2];
rz(-1.0216121) q[2];
sx q[2];
rz(2.2941022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60050868) q[1];
sx q[1];
rz(-2.2693172) q[1];
sx q[1];
rz(-0.38703106) q[1];
rz(-pi) q[2];
rz(1.2265008) q[3];
sx q[3];
rz(-0.85840423) q[3];
sx q[3];
rz(-2.815849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2465308) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(-2.4728298) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(-3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(-2.0892129) q[0];
rz(-1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(-0.88422424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28750557) q[0];
sx q[0];
rz(-2.8794718) q[0];
sx q[0];
rz(0.8192807) q[0];
x q[1];
rz(-2.0074559) q[2];
sx q[2];
rz(-1.1982802) q[2];
sx q[2];
rz(2.8357752) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4864038) q[1];
sx q[1];
rz(-2.9826394) q[1];
sx q[1];
rz(-2.7093191) q[1];
rz(-2.103881) q[3];
sx q[3];
rz(-1.8928877) q[3];
sx q[3];
rz(-1.0853634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23652442) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(-1.0926251) q[2];
rz(2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(0.75105914) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(1.3177692) q[0];
rz(-0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(2.172519) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58020619) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(-1.6798621) q[0];
rz(-0.71408748) q[2];
sx q[2];
rz(-1.6220777) q[2];
sx q[2];
rz(-0.72061348) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5321188) q[1];
sx q[1];
rz(-1.9874959) q[1];
sx q[1];
rz(2.6344224) q[1];
rz(-pi) q[2];
x q[2];
rz(-8*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(-1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.6102012) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(-1.4490022) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(1.4356027) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.766151) q[0];
sx q[0];
rz(-2.0631844) q[0];
sx q[0];
rz(-0.96813162) q[0];
rz(-pi) q[1];
rz(0.82201634) q[2];
sx q[2];
rz(-0.026958131) q[2];
sx q[2];
rz(2.6798623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83511931) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(-1.9496099) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2613951) q[3];
sx q[3];
rz(-2.5185563) q[3];
sx q[3];
rz(-0.33720371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.7763604) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2030905) q[0];
sx q[0];
rz(-1.7068958) q[0];
sx q[0];
rz(0.076119856) q[0];
rz(-pi) q[1];
rz(1.442996) q[2];
sx q[2];
rz(-1.5926077) q[2];
sx q[2];
rz(-2.4335361) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33278388) q[1];
sx q[1];
rz(-2.60751) q[1];
sx q[1];
rz(0.10314718) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.385231) q[3];
sx q[3];
rz(-1.984333) q[3];
sx q[3];
rz(2.6698649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0662971) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(-0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(2.2176567) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3850708) q[0];
sx q[0];
rz(-1.1413045) q[0];
sx q[0];
rz(0.71682741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6529979) q[2];
sx q[2];
rz(-1.5602419) q[2];
sx q[2];
rz(2.7985364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20272217) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(0.90941456) q[1];
x q[2];
rz(2.1903673) q[3];
sx q[3];
rz(-2.1853672) q[3];
sx q[3];
rz(0.28704498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.3941992) q[2];
rz(1.9723069) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(-0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37344638) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.753153) q[0];
rz(0.0069847981) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(3.1034234) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.105135) q[0];
sx q[0];
rz(-1.5974701) q[0];
sx q[0];
rz(2.1795576) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1881234) q[2];
sx q[2];
rz(-1.6560935) q[2];
sx q[2];
rz(-1.6558937) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.1306418) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(1.6846659) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90922728) q[3];
sx q[3];
rz(-2.0709166) q[3];
sx q[3];
rz(1.7928894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.4981131) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-0.11001982) q[2];
sx q[2];
rz(-1.5129871) q[2];
sx q[2];
rz(-1.5469698) q[2];
rz(2.8244143) q[3];
sx q[3];
rz(-0.89430292) q[3];
sx q[3];
rz(3.111307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
