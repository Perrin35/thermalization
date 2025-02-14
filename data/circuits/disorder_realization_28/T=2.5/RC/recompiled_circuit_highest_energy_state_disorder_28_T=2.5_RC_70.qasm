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
rz(-2.3856491) q[0];
sx q[0];
rz(-0.33060253) q[0];
sx q[0];
rz(-2.8688353) q[0];
rz(-0.35562149) q[1];
sx q[1];
rz(-1.0300809) q[1];
sx q[1];
rz(1.570809) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2251303) q[0];
sx q[0];
rz(-1.4989841) q[0];
sx q[0];
rz(-1.4738002) q[0];
rz(-pi) q[1];
rz(1.9000969) q[2];
sx q[2];
rz(-2.0138171) q[2];
sx q[2];
rz(1.4141413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94704365) q[1];
sx q[1];
rz(-1.5284555) q[1];
sx q[1];
rz(2.4621099) q[1];
x q[2];
rz(0.2491643) q[3];
sx q[3];
rz(-0.37940413) q[3];
sx q[3];
rz(0.5642341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1640799) q[2];
sx q[2];
rz(-1.432632) q[2];
sx q[2];
rz(2.254503) q[2];
rz(-1.2648434) q[3];
sx q[3];
rz(-1.0586459) q[3];
sx q[3];
rz(3.1177055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.3905268) q[0];
sx q[0];
rz(-0.85705119) q[0];
sx q[0];
rz(-0.29399011) q[0];
rz(1.7901621) q[1];
sx q[1];
rz(-1.0963115) q[1];
sx q[1];
rz(-1.1658123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.074075) q[0];
sx q[0];
rz(-2.7505028) q[0];
sx q[0];
rz(-2.142201) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1140536) q[2];
sx q[2];
rz(-1.6438494) q[2];
sx q[2];
rz(-2.6932584) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9820751) q[1];
sx q[1];
rz(-2.0802092) q[1];
sx q[1];
rz(2.0768836) q[1];
x q[2];
rz(2.7044917) q[3];
sx q[3];
rz(-2.5270998) q[3];
sx q[3];
rz(1.3853427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68842781) q[2];
sx q[2];
rz(-2.5025949) q[2];
sx q[2];
rz(1.2028018) q[2];
rz(-1.5996251) q[3];
sx q[3];
rz(-1.3381693) q[3];
sx q[3];
rz(-0.28551027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315392) q[0];
sx q[0];
rz(-0.078131214) q[0];
sx q[0];
rz(1.1745289) q[0];
rz(0.9749167) q[1];
sx q[1];
rz(-2.2623623) q[1];
sx q[1];
rz(2.6851795) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5147091) q[0];
sx q[0];
rz(-2.868342) q[0];
sx q[0];
rz(1.4529626) q[0];
x q[1];
rz(-2.1295616) q[2];
sx q[2];
rz(-0.44545275) q[2];
sx q[2];
rz(-0.73778668) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2043201) q[1];
sx q[1];
rz(-1.0144925) q[1];
sx q[1];
rz(-0.67436237) q[1];
rz(-pi) q[2];
rz(2.1022845) q[3];
sx q[3];
rz(-2.3520497) q[3];
sx q[3];
rz(1.9122151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5259033) q[2];
sx q[2];
rz(-2.758785) q[2];
sx q[2];
rz(2.3846386) q[2];
rz(-3.1225539) q[3];
sx q[3];
rz(-1.8502539) q[3];
sx q[3];
rz(-2.3058057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2794936) q[0];
sx q[0];
rz(-0.71074301) q[0];
sx q[0];
rz(2.1918462) q[0];
rz(2.9277335) q[1];
sx q[1];
rz(-0.93997926) q[1];
sx q[1];
rz(-0.59923879) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1480162) q[0];
sx q[0];
rz(-2.2207119) q[0];
sx q[0];
rz(-2.9666569) q[0];
x q[1];
rz(-2.9733359) q[2];
sx q[2];
rz(-1.0609385) q[2];
sx q[2];
rz(2.5349701) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.487452) q[1];
sx q[1];
rz(-1.0448487) q[1];
sx q[1];
rz(2.9000645) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20885373) q[3];
sx q[3];
rz(-2.5503272) q[3];
sx q[3];
rz(2.8610736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5102101) q[2];
sx q[2];
rz(-1.5082794) q[2];
sx q[2];
rz(-0.82376662) q[2];
rz(-2.0508164) q[3];
sx q[3];
rz(-0.80086509) q[3];
sx q[3];
rz(2.4763079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7432231) q[0];
sx q[0];
rz(-0.38341612) q[0];
sx q[0];
rz(-2.1552591) q[0];
rz(-0.24364722) q[1];
sx q[1];
rz(-1.2657093) q[1];
sx q[1];
rz(-2.9820014) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361317) q[0];
sx q[0];
rz(-1.8842995) q[0];
sx q[0];
rz(0.82833293) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99981941) q[2];
sx q[2];
rz(-2.2717064) q[2];
sx q[2];
rz(0.69055218) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9265598) q[1];
sx q[1];
rz(-1.4825645) q[1];
sx q[1];
rz(-2.0162321) q[1];
rz(0.87017228) q[3];
sx q[3];
rz(-1.7990094) q[3];
sx q[3];
rz(2.1765577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39449447) q[2];
sx q[2];
rz(-1.2622086) q[2];
sx q[2];
rz(2.0060284) q[2];
rz(-2.6731532) q[3];
sx q[3];
rz(-1.2218385) q[3];
sx q[3];
rz(-2.2170317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0202476) q[0];
sx q[0];
rz(-1.0806885) q[0];
sx q[0];
rz(-0.24170804) q[0];
rz(1.3555591) q[1];
sx q[1];
rz(-0.87239289) q[1];
sx q[1];
rz(0.88368574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51118219) q[0];
sx q[0];
rz(-2.1616461) q[0];
sx q[0];
rz(-2.1686337) q[0];
rz(-pi) q[1];
rz(1.8275798) q[2];
sx q[2];
rz(-1.2681539) q[2];
sx q[2];
rz(-3.1336477) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3436529) q[1];
sx q[1];
rz(-2.6561497) q[1];
sx q[1];
rz(1.14231) q[1];
x q[2];
rz(-1.9817686) q[3];
sx q[3];
rz(-1.5710063) q[3];
sx q[3];
rz(2.2157912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9200865) q[2];
sx q[2];
rz(-2.2111427) q[2];
sx q[2];
rz(-0.42631701) q[2];
rz(0.91655556) q[3];
sx q[3];
rz(-1.5896268) q[3];
sx q[3];
rz(2.0385273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.351848) q[0];
sx q[0];
rz(-1.9639356) q[0];
sx q[0];
rz(1.1258997) q[0];
rz(3.0806091) q[1];
sx q[1];
rz(-2.1911502) q[1];
sx q[1];
rz(-1.0152063) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39097745) q[0];
sx q[0];
rz(-2.1732507) q[0];
sx q[0];
rz(-2.0842805) q[0];
x q[1];
rz(-1.9092375) q[2];
sx q[2];
rz(-1.6403753) q[2];
sx q[2];
rz(-1.4078946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8440666) q[1];
sx q[1];
rz(-1.6324784) q[1];
sx q[1];
rz(1.9703377) q[1];
x q[2];
rz(0.33608243) q[3];
sx q[3];
rz(-3.0214632) q[3];
sx q[3];
rz(1.1230113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31275493) q[2];
sx q[2];
rz(-1.9373031) q[2];
sx q[2];
rz(-0.22605669) q[2];
rz(-0.95212805) q[3];
sx q[3];
rz(-0.13855562) q[3];
sx q[3];
rz(2.3329195) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1299745) q[0];
sx q[0];
rz(-1.8267153) q[0];
sx q[0];
rz(-0.34039482) q[0];
rz(0.099523425) q[1];
sx q[1];
rz(-1.1025905) q[1];
sx q[1];
rz(-1.5460825) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9700546) q[0];
sx q[0];
rz(-0.70101372) q[0];
sx q[0];
rz(-2.6215977) q[0];
rz(-1.298102) q[2];
sx q[2];
rz(-1.6249949) q[2];
sx q[2];
rz(1.7555838) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73879646) q[1];
sx q[1];
rz(-1.5288884) q[1];
sx q[1];
rz(0.099385029) q[1];
rz(-2.9603954) q[3];
sx q[3];
rz(-2.6799723) q[3];
sx q[3];
rz(-2.4571854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.978329) q[2];
sx q[2];
rz(-1.1922057) q[2];
sx q[2];
rz(1.0015798) q[2];
rz(-1.7892276) q[3];
sx q[3];
rz(-2.6632301) q[3];
sx q[3];
rz(-1.9549687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773862) q[0];
sx q[0];
rz(-1.0564251) q[0];
sx q[0];
rz(3.1412083) q[0];
rz(2.7035825) q[1];
sx q[1];
rz(-1.6338467) q[1];
sx q[1];
rz(-2.6843574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3246888) q[0];
sx q[0];
rz(-1.4131372) q[0];
sx q[0];
rz(-0.05013843) q[0];
x q[1];
rz(2.1887378) q[2];
sx q[2];
rz(-1.0818726) q[2];
sx q[2];
rz(3.0143154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.286736) q[1];
sx q[1];
rz(-1.348362) q[1];
sx q[1];
rz(0.35448928) q[1];
x q[2];
rz(-0.5064241) q[3];
sx q[3];
rz(-1.862414) q[3];
sx q[3];
rz(-2.9675561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3823201) q[2];
sx q[2];
rz(-2.1542532) q[2];
sx q[2];
rz(-1.6299204) q[2];
rz(-0.090653732) q[3];
sx q[3];
rz(-1.0414618) q[3];
sx q[3];
rz(-2.4618728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35336211) q[0];
sx q[0];
rz(-1.7857977) q[0];
sx q[0];
rz(0.28025383) q[0];
rz(-1.7578112) q[1];
sx q[1];
rz(-0.46418142) q[1];
sx q[1];
rz(1.2253449) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31497629) q[0];
sx q[0];
rz(-2.9626863) q[0];
sx q[0];
rz(2.5681483) q[0];
rz(0.88718517) q[2];
sx q[2];
rz(-1.9211743) q[2];
sx q[2];
rz(-2.7840707) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6514488) q[1];
sx q[1];
rz(-1.8134535) q[1];
sx q[1];
rz(2.3745684) q[1];
rz(-1.8605385) q[3];
sx q[3];
rz(-1.8612766) q[3];
sx q[3];
rz(1.2742467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19259351) q[2];
sx q[2];
rz(-1.5263564) q[2];
sx q[2];
rz(0.68142778) q[2];
rz(-2.0415908) q[3];
sx q[3];
rz(-2.2097578) q[3];
sx q[3];
rz(-1.3070235) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.451402) q[0];
sx q[0];
rz(-2.207353) q[0];
sx q[0];
rz(1.9718715) q[0];
rz(2.7611217) q[1];
sx q[1];
rz(-0.66743757) q[1];
sx q[1];
rz(0.22421511) q[1];
rz(-0.26421122) q[2];
sx q[2];
rz(-2.254494) q[2];
sx q[2];
rz(1.2848134) q[2];
rz(-2.1627103) q[3];
sx q[3];
rz(-1.9364193) q[3];
sx q[3];
rz(1.3649552) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
