OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5858894) q[0];
sx q[0];
rz(-1.1893505) q[0];
sx q[0];
rz(-1.1638292) q[0];
rz(2.3948506) q[1];
sx q[1];
rz(-2.876694) q[1];
sx q[1];
rz(2.9704111) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8590916) q[0];
sx q[0];
rz(-1.4177979) q[0];
sx q[0];
rz(1.7923844) q[0];
rz(-pi) q[1];
rz(-0.88340448) q[2];
sx q[2];
rz(-2.2210178) q[2];
sx q[2];
rz(1.6028849) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65217742) q[1];
sx q[1];
rz(-2.846183) q[1];
sx q[1];
rz(1.5687464) q[1];
rz(3.1247507) q[3];
sx q[3];
rz(-1.5653658) q[3];
sx q[3];
rz(3.0433486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.091398679) q[2];
sx q[2];
rz(-1.5768496) q[2];
sx q[2];
rz(1.0696627) q[2];
rz(1.7012677) q[3];
sx q[3];
rz(-2.1170719) q[3];
sx q[3];
rz(1.1721771) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1759724) q[0];
sx q[0];
rz(-1.6930641) q[0];
sx q[0];
rz(0.58797055) q[0];
rz(1.8486456) q[1];
sx q[1];
rz(-0.99423948) q[1];
sx q[1];
rz(-2.9867244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42255369) q[0];
sx q[0];
rz(-2.6309359) q[0];
sx q[0];
rz(-2.9209748) q[0];
x q[1];
rz(2.9005342) q[2];
sx q[2];
rz(-1.4285285) q[2];
sx q[2];
rz(0.011474284) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9714635) q[1];
sx q[1];
rz(-1.9057353) q[1];
sx q[1];
rz(-1.0748802) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4343461) q[3];
sx q[3];
rz(-1.3605355) q[3];
sx q[3];
rz(-2.7790359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4984442) q[2];
sx q[2];
rz(-2.4109106) q[2];
sx q[2];
rz(3.0291962) q[2];
rz(-0.82320881) q[3];
sx q[3];
rz(-1.9197437) q[3];
sx q[3];
rz(-2.6196151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68558973) q[0];
sx q[0];
rz(-2.0528448) q[0];
sx q[0];
rz(-2.1734557) q[0];
rz(1.118842) q[1];
sx q[1];
rz(-1.6066931) q[1];
sx q[1];
rz(1.8082089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1407993) q[0];
sx q[0];
rz(-0.62441545) q[0];
sx q[0];
rz(-1.8483759) q[0];
x q[1];
rz(2.5381375) q[2];
sx q[2];
rz(-1.1759182) q[2];
sx q[2];
rz(1.4970571) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1044429) q[1];
sx q[1];
rz(-2.5698476) q[1];
sx q[1];
rz(-1.7556095) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2127688) q[3];
sx q[3];
rz(-1.4910526) q[3];
sx q[3];
rz(2.2984576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94990388) q[2];
sx q[2];
rz(-2.1046941) q[2];
sx q[2];
rz(-2.8766768) q[2];
rz(0.31629899) q[3];
sx q[3];
rz(-2.8079872) q[3];
sx q[3];
rz(1.7501638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3259657) q[0];
sx q[0];
rz(-2.9209324) q[0];
sx q[0];
rz(-1.6478446) q[0];
rz(-1.7296467) q[1];
sx q[1];
rz(-1.0271881) q[1];
sx q[1];
rz(0.6573917) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9566204) q[0];
sx q[0];
rz(-1.4669384) q[0];
sx q[0];
rz(1.4968027) q[0];
x q[1];
rz(-2.3573586) q[2];
sx q[2];
rz(-1.1335229) q[2];
sx q[2];
rz(1.0904877) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58739793) q[1];
sx q[1];
rz(-0.19853354) q[1];
sx q[1];
rz(1.3819329) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74243706) q[3];
sx q[3];
rz(-2.1898139) q[3];
sx q[3];
rz(0.05832626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9221981) q[2];
sx q[2];
rz(-1.9191091) q[2];
sx q[2];
rz(-0.40327367) q[2];
rz(-1.9379617) q[3];
sx q[3];
rz(-1.5091242) q[3];
sx q[3];
rz(0.50886124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68172425) q[0];
sx q[0];
rz(-0.42400703) q[0];
sx q[0];
rz(-2.5176609) q[0];
rz(-2.4195747) q[1];
sx q[1];
rz(-1.3474418) q[1];
sx q[1];
rz(3.0375979) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0777587) q[0];
sx q[0];
rz(-0.88853271) q[0];
sx q[0];
rz(1.3880678) q[0];
rz(0.76461069) q[2];
sx q[2];
rz(-0.48073623) q[2];
sx q[2];
rz(1.802747) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38253575) q[1];
sx q[1];
rz(-2.0811354) q[1];
sx q[1];
rz(-2.3063763) q[1];
x q[2];
rz(0.89678371) q[3];
sx q[3];
rz(-2.1851903) q[3];
sx q[3];
rz(-2.3562285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7907052) q[2];
sx q[2];
rz(-0.71983379) q[2];
sx q[2];
rz(-2.3525815) q[2];
rz(-1.8770494) q[3];
sx q[3];
rz(-1.0435373) q[3];
sx q[3];
rz(-2.1544382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8378545) q[0];
sx q[0];
rz(-0.721295) q[0];
sx q[0];
rz(-2.9413132) q[0];
rz(1.6750977) q[1];
sx q[1];
rz(-1.8148345) q[1];
sx q[1];
rz(-1.5178348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2643833) q[0];
sx q[0];
rz(-2.0743999) q[0];
sx q[0];
rz(-1.1491386) q[0];
rz(-pi) q[1];
rz(2.0744223) q[2];
sx q[2];
rz(-1.4563515) q[2];
sx q[2];
rz(1.1888072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1788097) q[1];
sx q[1];
rz(-0.97985044) q[1];
sx q[1];
rz(1.5750118) q[1];
rz(2.3430941) q[3];
sx q[3];
rz(-1.3991465) q[3];
sx q[3];
rz(0.6397748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62310654) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(-2.8002807) q[2];
rz(-2.2845279) q[3];
sx q[3];
rz(-1.9144446) q[3];
sx q[3];
rz(2.6763693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.051006) q[0];
sx q[0];
rz(-1.0336646) q[0];
sx q[0];
rz(-1.8940014) q[0];
rz(-0.11820758) q[1];
sx q[1];
rz(-1.8659614) q[1];
sx q[1];
rz(-3.1059713) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6929631) q[0];
sx q[0];
rz(-1.5933196) q[0];
sx q[0];
rz(1.6883786) q[0];
rz(1.1664671) q[2];
sx q[2];
rz(-1.7448062) q[2];
sx q[2];
rz(-1.6650944) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.25917398) q[1];
sx q[1];
rz(-2.4738564) q[1];
sx q[1];
rz(0.40294934) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12322085) q[3];
sx q[3];
rz(-2.356592) q[3];
sx q[3];
rz(2.8383925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9901765) q[2];
sx q[2];
rz(-1.6665062) q[2];
sx q[2];
rz(-2.0659633) q[2];
rz(0.52078024) q[3];
sx q[3];
rz(-2.5213089) q[3];
sx q[3];
rz(-1.552593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.90466475) q[0];
sx q[0];
rz(-2.1187145) q[0];
sx q[0];
rz(-2.2383595) q[0];
rz(1.5029933) q[1];
sx q[1];
rz(-1.4475977) q[1];
sx q[1];
rz(1.2618056) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55106581) q[0];
sx q[0];
rz(-1.3414978) q[0];
sx q[0];
rz(-0.14044827) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4440358) q[2];
sx q[2];
rz(-2.4930525) q[2];
sx q[2];
rz(1.4319789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8324229) q[1];
sx q[1];
rz(-1.6214402) q[1];
sx q[1];
rz(-3.0502351) q[1];
x q[2];
rz(-0.66916211) q[3];
sx q[3];
rz(-2.4186196) q[3];
sx q[3];
rz(-2.1694136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5624076) q[2];
sx q[2];
rz(-0.6520485) q[2];
sx q[2];
rz(-0.11422608) q[2];
rz(1.3769897) q[3];
sx q[3];
rz(-2.0012794) q[3];
sx q[3];
rz(-0.50447869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0196411) q[0];
sx q[0];
rz(-2.1421102) q[0];
sx q[0];
rz(-1.9695388) q[0];
rz(-1.3837586) q[1];
sx q[1];
rz(-1.1446605) q[1];
sx q[1];
rz(-0.11464548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5929554) q[0];
sx q[0];
rz(-2.1229738) q[0];
sx q[0];
rz(2.1594285) q[0];
rz(-1.6646321) q[2];
sx q[2];
rz(-2.0619446) q[2];
sx q[2];
rz(-2.0744086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.35797) q[1];
sx q[1];
rz(-1.3232035) q[1];
sx q[1];
rz(-3.0686343) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10520868) q[3];
sx q[3];
rz(-2.3423839) q[3];
sx q[3];
rz(1.5952283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5961479) q[2];
sx q[2];
rz(-0.20716509) q[2];
sx q[2];
rz(-2.2607415) q[2];
rz(-2.7018069) q[3];
sx q[3];
rz(-1.9763016) q[3];
sx q[3];
rz(1.1254651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3038444) q[0];
sx q[0];
rz(-1.6551908) q[0];
sx q[0];
rz(3.0434171) q[0];
rz(2.5671666) q[1];
sx q[1];
rz(-1.6822633) q[1];
sx q[1];
rz(-2.548545) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0837729) q[0];
sx q[0];
rz(-1.2811986) q[0];
sx q[0];
rz(0.82478351) q[0];
rz(-pi) q[1];
rz(-2.2111528) q[2];
sx q[2];
rz(-0.64210063) q[2];
sx q[2];
rz(0.053089945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1758842) q[1];
sx q[1];
rz(-1.5820326) q[1];
sx q[1];
rz(-1.5574993) q[1];
rz(0.63238588) q[3];
sx q[3];
rz(-1.2999855) q[3];
sx q[3];
rz(0.89030823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9888931) q[2];
sx q[2];
rz(-1.9579192) q[2];
sx q[2];
rz(2.2031671) q[2];
rz(-1.3931795) q[3];
sx q[3];
rz(-1.8311071) q[3];
sx q[3];
rz(-0.2383298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68168454) q[0];
sx q[0];
rz(-0.60612283) q[0];
sx q[0];
rz(-1.7932307) q[0];
rz(-1.0472736) q[1];
sx q[1];
rz(-2.2127163) q[1];
sx q[1];
rz(-0.97102078) q[1];
rz(-0.045108724) q[2];
sx q[2];
rz(-1.5061629) q[2];
sx q[2];
rz(-0.20284222) q[2];
rz(-0.88225928) q[3];
sx q[3];
rz(-2.8591446) q[3];
sx q[3];
rz(2.1459116) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
