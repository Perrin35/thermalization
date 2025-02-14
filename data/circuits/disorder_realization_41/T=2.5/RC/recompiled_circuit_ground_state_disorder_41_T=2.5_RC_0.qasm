OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(1.289225) q[0];
sx q[0];
rz(9.3318648) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(-0.73268259) q[1];
sx q[1];
rz(-1.2394152) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6042955) q[0];
sx q[0];
rz(-1.822246) q[0];
sx q[0];
rz(2.9049113) q[0];
rz(-pi) q[1];
rz(-2.5074717) q[2];
sx q[2];
rz(-0.32358746) q[2];
sx q[2];
rz(2.7860118) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8653633) q[1];
sx q[1];
rz(-1.5412093) q[1];
sx q[1];
rz(-0.33025708) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4437583) q[3];
sx q[3];
rz(-0.69785944) q[3];
sx q[3];
rz(2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6713082) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(1.7926463) q[2];
rz(2.3979483) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(0.17717895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.0002366) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(-0.29362383) q[0];
rz(-1.0307505) q[1];
sx q[1];
rz(-1.3615969) q[1];
sx q[1];
rz(-0.34293276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2156649) q[0];
sx q[0];
rz(-0.56680381) q[0];
sx q[0];
rz(1.9606007) q[0];
rz(-pi) q[1];
rz(-1.8937102) q[2];
sx q[2];
rz(-2.7512449) q[2];
sx q[2];
rz(2.5935612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3174177) q[1];
sx q[1];
rz(-1.3567827) q[1];
sx q[1];
rz(-2.432968) q[1];
rz(0.94128709) q[3];
sx q[3];
rz(-1.7751179) q[3];
sx q[3];
rz(-2.2245537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0127516) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(0.020708474) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6193806) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(-0.0044599175) q[0];
rz(2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-2.3562145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4606708) q[0];
sx q[0];
rz(-0.17853949) q[0];
sx q[0];
rz(-2.6723249) q[0];
rz(-pi) q[1];
rz(-0.17526971) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(-0.64495211) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2783616) q[1];
sx q[1];
rz(-1.0309217) q[1];
sx q[1];
rz(0.079915483) q[1];
x q[2];
rz(0.80180577) q[3];
sx q[3];
rz(-0.84322819) q[3];
sx q[3];
rz(-1.8319195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-0.9202756) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(2.8010098) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(-2.3497439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444645) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(2.2564364) q[0];
rz(-2.3947233) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(-1.0964099) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.261026) q[0];
sx q[0];
rz(-1.9184904) q[0];
sx q[0];
rz(0.33190042) q[0];
x q[1];
rz(3.1109516) q[2];
sx q[2];
rz(-1.6302135) q[2];
sx q[2];
rz(0.63842809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29624048) q[1];
sx q[1];
rz(-1.3120323) q[1];
sx q[1];
rz(0.47611632) q[1];
rz(0.73170264) q[3];
sx q[3];
rz(-1.9644794) q[3];
sx q[3];
rz(1.3418152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5235644) q[2];
sx q[2];
rz(-0.21054331) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(-2.5060182) q[3];
sx q[3];
rz(-0.98819757) q[3];
sx q[3];
rz(-2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4516975) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(-2.6746993) q[0];
rz(3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(2.0707524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0127481) q[0];
sx q[0];
rz(-1.8086047) q[0];
sx q[0];
rz(1.1046011) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51009615) q[2];
sx q[2];
rz(-2.7802711) q[2];
sx q[2];
rz(2.9053743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5731949) q[1];
sx q[1];
rz(-2.0802976) q[1];
sx q[1];
rz(-2.0054818) q[1];
rz(2.4399906) q[3];
sx q[3];
rz(-1.0402586) q[3];
sx q[3];
rz(0.70487937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0232627) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(-3.1039589) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0090050176) q[0];
sx q[0];
rz(-1.9629033) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(-1.5124849) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(1.9992794) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5286251) q[0];
sx q[0];
rz(-0.36699793) q[0];
sx q[0];
rz(-0.27530833) q[0];
x q[1];
rz(1.7785759) q[2];
sx q[2];
rz(-1.5099974) q[2];
sx q[2];
rz(-0.5375934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9276322) q[1];
sx q[1];
rz(-2.7134656) q[1];
sx q[1];
rz(2.5090748) q[1];
rz(-pi) q[2];
rz(0.31123881) q[3];
sx q[3];
rz(-1.7390307) q[3];
sx q[3];
rz(0.5420891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(3.0211871) q[2];
rz(-1.1558007) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(-2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8026546) q[0];
sx q[0];
rz(-1.3405565) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(-1.5441719) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(-0.45101756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9266201) q[0];
sx q[0];
rz(-0.30719137) q[0];
sx q[0];
rz(0.19812576) q[0];
rz(-pi) q[1];
rz(-1.0181997) q[2];
sx q[2];
rz(-0.42358735) q[2];
sx q[2];
rz(-2.6747963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1731733) q[1];
sx q[1];
rz(-2.0035158) q[1];
sx q[1];
rz(1.7041901) q[1];
x q[2];
rz(3.1322828) q[3];
sx q[3];
rz(-0.81039372) q[3];
sx q[3];
rz(1.9805465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94577998) q[2];
sx q[2];
rz(-1.7273629) q[2];
sx q[2];
rz(1.3607402) q[2];
rz(2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(-1.7520693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239531) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-2.4355167) q[0];
rz(-1.5977244) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(-0.85618883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9129974) q[0];
sx q[0];
rz(-1.5513591) q[0];
sx q[0];
rz(-3.134208) q[0];
rz(-pi) q[1];
rz(-2.0308308) q[2];
sx q[2];
rz(-0.80321124) q[2];
sx q[2];
rz(-2.3873461) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4256712) q[1];
sx q[1];
rz(-1.7737264) q[1];
sx q[1];
rz(-0.58260609) q[1];
rz(-pi) q[2];
rz(3.007902) q[3];
sx q[3];
rz(-2.055759) q[3];
sx q[3];
rz(-0.07459379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84264821) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(-0.16743463) q[2];
rz(1.5873448) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7635968) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(1.7077712) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(0.30002123) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8269129) q[0];
sx q[0];
rz(-1.2859435) q[0];
sx q[0];
rz(-1.6919943) q[0];
rz(-pi) q[1];
rz(-0.062252684) q[2];
sx q[2];
rz(-0.57800284) q[2];
sx q[2];
rz(0.5731155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91050657) q[1];
sx q[1];
rz(-0.66201895) q[1];
sx q[1];
rz(-0.20259095) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2737502) q[3];
sx q[3];
rz(-1.6658837) q[3];
sx q[3];
rz(2.8114708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47664777) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(2.7139968) q[2];
rz(1.5852196) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(-2.3868886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72117358) q[0];
sx q[0];
rz(-2.6009646) q[0];
sx q[0];
rz(-0.46947259) q[0];
rz(-2.2162614) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(-3.1184149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7849451) q[0];
sx q[0];
rz(-1.0482422) q[0];
sx q[0];
rz(-1.7607259) q[0];
rz(-pi) q[1];
rz(2.2397348) q[2];
sx q[2];
rz(-1.0791856) q[2];
sx q[2];
rz(2.7265446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1286298) q[1];
sx q[1];
rz(-0.2405162) q[1];
sx q[1];
rz(0.81324767) q[1];
rz(-1.7816824) q[3];
sx q[3];
rz(-2.4004705) q[3];
sx q[3];
rz(-1.1531342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8699441) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(0.49087697) q[2];
rz(1.0902181) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(0.31392613) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28521095) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(0.27676997) q[1];
sx q[1];
rz(-2.3634187) q[1];
sx q[1];
rz(1.707911) q[1];
rz(-2.7748952) q[2];
sx q[2];
rz(-2.4187805) q[2];
sx q[2];
rz(0.4772966) q[2];
rz(-0.033944081) q[3];
sx q[3];
rz(-0.53474075) q[3];
sx q[3];
rz(-0.026215601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
