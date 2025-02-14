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
rz(-2.9450077) q[0];
sx q[0];
rz(-0.44214806) q[0];
sx q[0];
rz(-2.2801939) q[0];
rz(1.5098894) q[1];
sx q[1];
rz(-1.8169401) q[1];
sx q[1];
rz(-0.035592508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0738348) q[0];
sx q[0];
rz(-1.7661716) q[0];
sx q[0];
rz(-2.8003576) q[0];
x q[1];
rz(-2.8350949) q[2];
sx q[2];
rz(-1.674429) q[2];
sx q[2];
rz(0.95970861) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0390748) q[1];
sx q[1];
rz(-0.84474428) q[1];
sx q[1];
rz(-1.5229358) q[1];
rz(-1.5930327) q[3];
sx q[3];
rz(-2.4614868) q[3];
sx q[3];
rz(0.82397991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.55277905) q[2];
sx q[2];
rz(-1.8401044) q[2];
sx q[2];
rz(-2.0471052) q[2];
rz(1.6257446) q[3];
sx q[3];
rz(-2.2690014) q[3];
sx q[3];
rz(-0.14357963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0933519) q[0];
sx q[0];
rz(-0.99026647) q[0];
sx q[0];
rz(2.7675203) q[0];
rz(2.2834942) q[1];
sx q[1];
rz(-0.94081196) q[1];
sx q[1];
rz(-1.0757974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6158524) q[0];
sx q[0];
rz(-1.6739556) q[0];
sx q[0];
rz(3.0391058) q[0];
rz(1.3316989) q[2];
sx q[2];
rz(-2.1054724) q[2];
sx q[2];
rz(0.83132529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0829216) q[1];
sx q[1];
rz(-1.4545014) q[1];
sx q[1];
rz(0.40177675) q[1];
x q[2];
rz(2.3736453) q[3];
sx q[3];
rz(-2.30184) q[3];
sx q[3];
rz(-2.6382382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51521236) q[2];
sx q[2];
rz(-0.43312803) q[2];
sx q[2];
rz(2.058378) q[2];
rz(-1.2927239) q[3];
sx q[3];
rz(-2.0079398) q[3];
sx q[3];
rz(0.92866549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2556297) q[0];
sx q[0];
rz(-2.3704381) q[0];
sx q[0];
rz(2.6166925) q[0];
rz(-1.6820172) q[1];
sx q[1];
rz(-0.73760215) q[1];
sx q[1];
rz(-0.14579138) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0130769) q[0];
sx q[0];
rz(-1.4798963) q[0];
sx q[0];
rz(0.10581067) q[0];
rz(-0.86051382) q[2];
sx q[2];
rz(-2.4053221) q[2];
sx q[2];
rz(0.78303534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5151095) q[1];
sx q[1];
rz(-1.6609539) q[1];
sx q[1];
rz(-0.43403352) q[1];
x q[2];
rz(0.57246142) q[3];
sx q[3];
rz(-2.2160071) q[3];
sx q[3];
rz(-1.087134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24474457) q[2];
sx q[2];
rz(-1.1620099) q[2];
sx q[2];
rz(2.9761918) q[2];
rz(0.8616972) q[3];
sx q[3];
rz(-2.4498037) q[3];
sx q[3];
rz(-2.944788) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2794613) q[0];
sx q[0];
rz(-2.617351) q[0];
sx q[0];
rz(2.4133546) q[0];
rz(-2.8184452) q[1];
sx q[1];
rz(-1.0341945) q[1];
sx q[1];
rz(-2.1023777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6704039) q[0];
sx q[0];
rz(-1.146933) q[0];
sx q[0];
rz(-0.43760646) q[0];
rz(-2.2065998) q[2];
sx q[2];
rz(-1.3424917) q[2];
sx q[2];
rz(1.993597) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2715214) q[1];
sx q[1];
rz(-1.0876552) q[1];
sx q[1];
rz(-1.4663694) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81890727) q[3];
sx q[3];
rz(-1.0207796) q[3];
sx q[3];
rz(-0.14980042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6497961) q[2];
sx q[2];
rz(-1.0504664) q[2];
sx q[2];
rz(0.51898471) q[2];
rz(-1.0745878) q[3];
sx q[3];
rz(-1.2858177) q[3];
sx q[3];
rz(-0.49526596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0786667) q[0];
sx q[0];
rz(-2.2195897) q[0];
sx q[0];
rz(0.17157383) q[0];
rz(2.0473139) q[1];
sx q[1];
rz(-1.3010052) q[1];
sx q[1];
rz(2.9069854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0701092) q[0];
sx q[0];
rz(-0.30710328) q[0];
sx q[0];
rz(-2.0060517) q[0];
rz(-pi) q[1];
rz(0.50723664) q[2];
sx q[2];
rz(-1.6167069) q[2];
sx q[2];
rz(-3.1243589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7947096) q[1];
sx q[1];
rz(-0.88635072) q[1];
sx q[1];
rz(-2.3331932) q[1];
rz(-pi) q[2];
rz(-0.42886491) q[3];
sx q[3];
rz(-0.75733105) q[3];
sx q[3];
rz(2.0067818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91308633) q[2];
sx q[2];
rz(-1.9501016) q[2];
sx q[2];
rz(-1.4027493) q[2];
rz(2.5168822) q[3];
sx q[3];
rz(-2.1731845) q[3];
sx q[3];
rz(1.3431842) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45873555) q[0];
sx q[0];
rz(-0.0077489297) q[0];
sx q[0];
rz(-0.32753456) q[0];
rz(3.1134743) q[1];
sx q[1];
rz(-1.997812) q[1];
sx q[1];
rz(-0.99172529) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018143749) q[0];
sx q[0];
rz(-0.94052343) q[0];
sx q[0];
rz(-0.92395497) q[0];
x q[1];
rz(-2.6985833) q[2];
sx q[2];
rz(-2.8923419) q[2];
sx q[2];
rz(-2.5009048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36054143) q[1];
sx q[1];
rz(-2.326093) q[1];
sx q[1];
rz(-0.057572854) q[1];
rz(-pi) q[2];
rz(-1.5434274) q[3];
sx q[3];
rz(-2.0413766) q[3];
sx q[3];
rz(-1.1740299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.907054) q[2];
sx q[2];
rz(-1.7639284) q[2];
sx q[2];
rz(-1.482359) q[2];
rz(1.0659069) q[3];
sx q[3];
rz(-1.9612471) q[3];
sx q[3];
rz(1.0970241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7557573) q[0];
sx q[0];
rz(-3.0297854) q[0];
sx q[0];
rz(-2.8461611) q[0];
rz(0.23751986) q[1];
sx q[1];
rz(-0.93370456) q[1];
sx q[1];
rz(-0.74737731) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6866561) q[0];
sx q[0];
rz(-2.1224408) q[0];
sx q[0];
rz(2.174189) q[0];
rz(-1.0355936) q[2];
sx q[2];
rz(-2.1280625) q[2];
sx q[2];
rz(-1.1320499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.69773) q[1];
sx q[1];
rz(-1.9404812) q[1];
sx q[1];
rz(-1.0401985) q[1];
rz(-pi) q[2];
rz(-1.9456057) q[3];
sx q[3];
rz(-2.0829933) q[3];
sx q[3];
rz(-2.1846445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6682917) q[2];
sx q[2];
rz(-1.95582) q[2];
sx q[2];
rz(-2.8866923) q[2];
rz(-2.9292246) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(0.00096850639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2905228) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(-3.0176924) q[0];
rz(-2.2663785) q[1];
sx q[1];
rz(-2.3085935) q[1];
sx q[1];
rz(-2.5828054) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66984719) q[0];
sx q[0];
rz(-2.5515243) q[0];
sx q[0];
rz(-2.0315995) q[0];
rz(-pi) q[1];
rz(2.0697849) q[2];
sx q[2];
rz(-1.3200511) q[2];
sx q[2];
rz(2.8607228) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3271823) q[1];
sx q[1];
rz(-2.5375527) q[1];
sx q[1];
rz(-2.6534897) q[1];
rz(-pi) q[2];
rz(0.93620091) q[3];
sx q[3];
rz(-1.2290599) q[3];
sx q[3];
rz(-2.5777011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.045804068) q[2];
sx q[2];
rz(-0.49171058) q[2];
sx q[2];
rz(-2.4208505) q[2];
rz(0.36302429) q[3];
sx q[3];
rz(-1.2237153) q[3];
sx q[3];
rz(1.5117517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95661288) q[0];
sx q[0];
rz(-0.19421254) q[0];
sx q[0];
rz(-0.089056253) q[0];
rz(0.36674276) q[1];
sx q[1];
rz(-2.2373503) q[1];
sx q[1];
rz(1.6820224) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3066481) q[0];
sx q[0];
rz(-2.1900926) q[0];
sx q[0];
rz(-0.28138103) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6834246) q[2];
sx q[2];
rz(-1.0976085) q[2];
sx q[2];
rz(0.13147182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2024398) q[1];
sx q[1];
rz(-1.9382361) q[1];
sx q[1];
rz(-0.82384404) q[1];
x q[2];
rz(2.5427548) q[3];
sx q[3];
rz(-1.2019079) q[3];
sx q[3];
rz(1.9112596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9350932) q[2];
sx q[2];
rz(-1.7743899) q[2];
sx q[2];
rz(-1.8915141) q[2];
rz(1.9153204) q[3];
sx q[3];
rz(-0.63392249) q[3];
sx q[3];
rz(-2.5026076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117689) q[0];
sx q[0];
rz(-1.1331695) q[0];
sx q[0];
rz(-2.0015707) q[0];
rz(-0.58311588) q[1];
sx q[1];
rz(-1.4864328) q[1];
sx q[1];
rz(2.4737632) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4224189) q[0];
sx q[0];
rz(-1.2454709) q[0];
sx q[0];
rz(0.76119411) q[0];
rz(2.2726353) q[2];
sx q[2];
rz(-1.7545106) q[2];
sx q[2];
rz(0.76317235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9676303) q[1];
sx q[1];
rz(-2.6625427) q[1];
sx q[1];
rz(1.1381989) q[1];
rz(-pi) q[2];
rz(2.0290501) q[3];
sx q[3];
rz(-0.83247165) q[3];
sx q[3];
rz(-2.3979681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84393152) q[2];
sx q[2];
rz(-2.6466978) q[2];
sx q[2];
rz(-2.3893791) q[2];
rz(0.080605896) q[3];
sx q[3];
rz(-1.3496496) q[3];
sx q[3];
rz(-2.5940671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3362296) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(2.2930131) q[1];
sx q[1];
rz(-2.4129557) q[1];
sx q[1];
rz(1.6853263) q[1];
rz(-2.3586629) q[2];
sx q[2];
rz(-0.83275262) q[2];
sx q[2];
rz(0.61093753) q[2];
rz(0.6497339) q[3];
sx q[3];
rz(-2.5841373) q[3];
sx q[3];
rz(2.0852603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
