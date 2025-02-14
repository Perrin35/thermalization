OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.10915169) q[0];
sx q[0];
rz(4.2946058) q[0];
sx q[0];
rz(9.499318) q[0];
rz(-1.5300765) q[1];
sx q[1];
rz(-0.059377436) q[1];
sx q[1];
rz(0.52409726) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4573114) q[0];
sx q[0];
rz(-2.1403098) q[0];
sx q[0];
rz(2.7049541) q[0];
rz(-0.33676001) q[2];
sx q[2];
rz(-1.3408061) q[2];
sx q[2];
rz(-0.57631341) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0007038) q[1];
sx q[1];
rz(-0.86727625) q[1];
sx q[1];
rz(2.3402653) q[1];
x q[2];
rz(-1.1537179) q[3];
sx q[3];
rz(-3.0410389) q[3];
sx q[3];
rz(-1.2300958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3268299) q[2];
sx q[2];
rz(-0.49850285) q[2];
sx q[2];
rz(0.89547431) q[2];
rz(0.26252663) q[3];
sx q[3];
rz(-0.34590507) q[3];
sx q[3];
rz(-3.053022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.217591) q[0];
sx q[0];
rz(-1.1955248) q[0];
sx q[0];
rz(0.1161282) q[0];
rz(-2.5092292) q[1];
sx q[1];
rz(-0.26591161) q[1];
sx q[1];
rz(1.6557453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7065379) q[0];
sx q[0];
rz(-3.0710906) q[0];
sx q[0];
rz(0.32891957) q[0];
rz(1.952195) q[2];
sx q[2];
rz(-0.57514319) q[2];
sx q[2];
rz(-0.084897651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8101059) q[1];
sx q[1];
rz(-1.4432194) q[1];
sx q[1];
rz(1.677631) q[1];
x q[2];
rz(2.2712098) q[3];
sx q[3];
rz(-1.4705949) q[3];
sx q[3];
rz(-0.23476511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8302725) q[2];
sx q[2];
rz(-0.98036426) q[2];
sx q[2];
rz(-0.92973989) q[2];
rz(-0.43944198) q[3];
sx q[3];
rz(-1.6012871) q[3];
sx q[3];
rz(0.75696993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0759401) q[0];
sx q[0];
rz(-2.3891698) q[0];
sx q[0];
rz(-0.886985) q[0];
rz(1.2607964) q[1];
sx q[1];
rz(-2.0340684) q[1];
sx q[1];
rz(1.6541803) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58081478) q[0];
sx q[0];
rz(-0.11753035) q[0];
sx q[0];
rz(-1.5663901) q[0];
x q[1];
rz(1.4753129) q[2];
sx q[2];
rz(-1.9621358) q[2];
sx q[2];
rz(-3.0961159) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1765114) q[1];
sx q[1];
rz(-1.9002894) q[1];
sx q[1];
rz(-1.552187) q[1];
rz(-pi) q[2];
rz(-2.6657899) q[3];
sx q[3];
rz(-2.5047205) q[3];
sx q[3];
rz(2.3726557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6416574) q[2];
sx q[2];
rz(-0.96330088) q[2];
sx q[2];
rz(-2.7063043) q[2];
rz(0.38517243) q[3];
sx q[3];
rz(-1.9305072) q[3];
sx q[3];
rz(2.1729573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4061072) q[0];
sx q[0];
rz(-3.0574953) q[0];
sx q[0];
rz(-3.0870068) q[0];
rz(1.6361884) q[1];
sx q[1];
rz(-1.2137698) q[1];
sx q[1];
rz(2.7240567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82334119) q[0];
sx q[0];
rz(-2.7007036) q[0];
sx q[0];
rz(0.13680451) q[0];
x q[1];
rz(2.7052077) q[2];
sx q[2];
rz(-2.4616995) q[2];
sx q[2];
rz(-1.6298378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7723267) q[1];
sx q[1];
rz(-1.7290745) q[1];
sx q[1];
rz(-1.5637828) q[1];
rz(-pi) q[2];
rz(-1.3261111) q[3];
sx q[3];
rz(-0.26991329) q[3];
sx q[3];
rz(-3.1343366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2893082) q[2];
sx q[2];
rz(-2.4384629) q[2];
sx q[2];
rz(-3.1404148) q[2];
rz(-1.3834472) q[3];
sx q[3];
rz(-2.8774084) q[3];
sx q[3];
rz(2.0980825) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99172878) q[0];
sx q[0];
rz(-0.1606476) q[0];
sx q[0];
rz(-0.37586656) q[0];
rz(-2.1916892) q[1];
sx q[1];
rz(-1.4774731) q[1];
sx q[1];
rz(2.4163767) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0501971) q[0];
sx q[0];
rz(-1.689347) q[0];
sx q[0];
rz(-1.4539745) q[0];
x q[1];
rz(-2.2967417) q[2];
sx q[2];
rz(-2.2965659) q[2];
sx q[2];
rz(-1.8595075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8456373) q[1];
sx q[1];
rz(-1.6682677) q[1];
sx q[1];
rz(-1.5766162) q[1];
rz(1.3106724) q[3];
sx q[3];
rz(-2.2571466) q[3];
sx q[3];
rz(-2.7174866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7627365) q[2];
sx q[2];
rz(-1.4521658) q[2];
sx q[2];
rz(2.6689996) q[2];
rz(0.11911123) q[3];
sx q[3];
rz(-2.5178858) q[3];
sx q[3];
rz(0.81010336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.8238207) q[0];
sx q[0];
rz(-0.20109421) q[0];
sx q[0];
rz(-0.066545181) q[0];
rz(-2.9464974) q[1];
sx q[1];
rz(-1.0761484) q[1];
sx q[1];
rz(-1.1544352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46259743) q[0];
sx q[0];
rz(-2.6021311) q[0];
sx q[0];
rz(-2.0360002) q[0];
rz(-pi) q[1];
rz(2.6850291) q[2];
sx q[2];
rz(-1.4924268) q[2];
sx q[2];
rz(-2.5541277) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30224702) q[1];
sx q[1];
rz(-1.6776507) q[1];
sx q[1];
rz(1.9356506) q[1];
rz(1.5510173) q[3];
sx q[3];
rz(-0.58073509) q[3];
sx q[3];
rz(1.494921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2164312) q[2];
sx q[2];
rz(-1.602828) q[2];
sx q[2];
rz(-0.72539854) q[2];
rz(1.409449) q[3];
sx q[3];
rz(-1.9603445) q[3];
sx q[3];
rz(2.5349272) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055939097) q[0];
sx q[0];
rz(-0.5760718) q[0];
sx q[0];
rz(-2.2357909) q[0];
rz(0.55465758) q[1];
sx q[1];
rz(-2.3034818) q[1];
sx q[1];
rz(-0.14403266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54948037) q[0];
sx q[0];
rz(-1.0735816) q[0];
sx q[0];
rz(2.8449932) q[0];
rz(-3.0564791) q[2];
sx q[2];
rz(-2.0817698) q[2];
sx q[2];
rz(2.2764595) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0121978) q[1];
sx q[1];
rz(-1.8499814) q[1];
sx q[1];
rz(-1.9703945) q[1];
rz(-pi) q[2];
rz(-1.3200754) q[3];
sx q[3];
rz(-1.7169239) q[3];
sx q[3];
rz(-3.0878029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6683228) q[2];
sx q[2];
rz(-0.39445764) q[2];
sx q[2];
rz(-0.99297601) q[2];
rz(-2.9571577) q[3];
sx q[3];
rz(-2.1858229) q[3];
sx q[3];
rz(-2.1485645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1297146) q[0];
sx q[0];
rz(-2.5373902) q[0];
sx q[0];
rz(-2.7469444) q[0];
rz(-2.6036085) q[1];
sx q[1];
rz(-0.6901651) q[1];
sx q[1];
rz(-1.1916377) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93617986) q[0];
sx q[0];
rz(-1.8260341) q[0];
sx q[0];
rz(-2.5827239) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3778119) q[2];
sx q[2];
rz(-1.9375083) q[2];
sx q[2];
rz(0.16736469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0816457) q[1];
sx q[1];
rz(-0.73127247) q[1];
sx q[1];
rz(-0.88986963) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1147145) q[3];
sx q[3];
rz(-2.1563765) q[3];
sx q[3];
rz(0.28156137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3523606) q[2];
sx q[2];
rz(-1.5287986) q[2];
sx q[2];
rz(-2.1996876) q[2];
rz(-0.76147979) q[3];
sx q[3];
rz(-2.0165636) q[3];
sx q[3];
rz(3.1010845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86904675) q[0];
sx q[0];
rz(-1.7750374) q[0];
sx q[0];
rz(2.0696562) q[0];
rz(-2.5275224) q[1];
sx q[1];
rz(-2.2449988) q[1];
sx q[1];
rz(-0.44642064) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75335129) q[0];
sx q[0];
rz(-0.89453012) q[0];
sx q[0];
rz(0.36305289) q[0];
x q[1];
rz(0.69009366) q[2];
sx q[2];
rz(-1.141667) q[2];
sx q[2];
rz(-2.5370425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3371967) q[1];
sx q[1];
rz(-0.91584086) q[1];
sx q[1];
rz(1.4319929) q[1];
rz(-2.8458251) q[3];
sx q[3];
rz(-0.52433521) q[3];
sx q[3];
rz(1.7400896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3093695) q[2];
sx q[2];
rz(-2.1840405) q[2];
sx q[2];
rz(1.8831801) q[2];
rz(-1.9084515) q[3];
sx q[3];
rz(-0.62658739) q[3];
sx q[3];
rz(-1.8241749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719139) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.423214) q[0];
rz(0.24249679) q[1];
sx q[1];
rz(-2.6624694) q[1];
sx q[1];
rz(1.8877782) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4187188) q[0];
sx q[0];
rz(-0.44298816) q[0];
sx q[0];
rz(2.0110058) q[0];
x q[1];
rz(-0.10842936) q[2];
sx q[2];
rz(-3.0745818) q[2];
sx q[2];
rz(2.5052414) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36404836) q[1];
sx q[1];
rz(-1.2723039) q[1];
sx q[1];
rz(-1.0320028) q[1];
rz(2.9217072) q[3];
sx q[3];
rz(-1.9337092) q[3];
sx q[3];
rz(-0.21838926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8157876) q[2];
sx q[2];
rz(-0.82745224) q[2];
sx q[2];
rz(3.117756) q[2];
rz(1.8627953) q[3];
sx q[3];
rz(-2.8102504) q[3];
sx q[3];
rz(2.3648025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1635638) q[0];
sx q[0];
rz(-1.4373056) q[0];
sx q[0];
rz(1.9594255) q[0];
rz(-2.4143746) q[1];
sx q[1];
rz(-1.4677508) q[1];
sx q[1];
rz(-0.8263091) q[1];
rz(-3.1221409) q[2];
sx q[2];
rz(-1.7218334) q[2];
sx q[2];
rz(-1.3109372) q[2];
rz(0.047688382) q[3];
sx q[3];
rz(-1.449122) q[3];
sx q[3];
rz(2.4301152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
