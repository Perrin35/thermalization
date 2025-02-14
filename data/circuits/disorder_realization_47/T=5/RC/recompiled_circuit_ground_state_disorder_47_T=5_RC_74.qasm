OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(-1.2824143) q[0];
sx q[0];
rz(3.0106944) q[0];
rz(-2.6137597) q[1];
sx q[1];
rz(-0.37017828) q[1];
sx q[1];
rz(0.17732492) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41577121) q[0];
sx q[0];
rz(-2.3292589) q[0];
sx q[0];
rz(-1.2052631) q[0];
x q[1];
rz(0.22819569) q[2];
sx q[2];
rz(-1.2773809) q[2];
sx q[2];
rz(-2.5876837) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6244391) q[1];
sx q[1];
rz(-0.53688795) q[1];
sx q[1];
rz(-1.0978903) q[1];
rz(-pi) q[2];
rz(0.27324851) q[3];
sx q[3];
rz(-2.1299703) q[3];
sx q[3];
rz(-0.76598734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6191972) q[2];
sx q[2];
rz(-2.0530687) q[2];
sx q[2];
rz(1.8001451) q[2];
rz(-1.3522735) q[3];
sx q[3];
rz(-1.5315346) q[3];
sx q[3];
rz(1.8488098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.751048) q[0];
sx q[0];
rz(-1.9123257) q[0];
sx q[0];
rz(-2.6265662) q[0];
rz(2.7797049) q[1];
sx q[1];
rz(-1.3970951) q[1];
sx q[1];
rz(0.99641689) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1086209) q[0];
sx q[0];
rz(-1.7713799) q[0];
sx q[0];
rz(-1.8114835) q[0];
rz(-2.9662232) q[2];
sx q[2];
rz(-1.9992454) q[2];
sx q[2];
rz(-1.0340978) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60153466) q[1];
sx q[1];
rz(-2.6121579) q[1];
sx q[1];
rz(-1.387818) q[1];
rz(-2.1085694) q[3];
sx q[3];
rz(-1.9578079) q[3];
sx q[3];
rz(-1.8333216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4126052) q[2];
sx q[2];
rz(-0.94634405) q[2];
sx q[2];
rz(-1.4808572) q[2];
rz(-2.6442773) q[3];
sx q[3];
rz(-0.9484843) q[3];
sx q[3];
rz(2.4174387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6052674) q[0];
sx q[0];
rz(-1.8999506) q[0];
sx q[0];
rz(-1.190825) q[0];
rz(-0.39769998) q[1];
sx q[1];
rz(-1.5813446) q[1];
sx q[1];
rz(-0.89888987) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2069103) q[0];
sx q[0];
rz(-0.3461306) q[0];
sx q[0];
rz(-1.7680566) q[0];
x q[1];
rz(1.4853046) q[2];
sx q[2];
rz(-1.3807266) q[2];
sx q[2];
rz(-2.9051733) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11375107) q[1];
sx q[1];
rz(-1.6154624) q[1];
sx q[1];
rz(2.1034481) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8679138) q[3];
sx q[3];
rz(-2.6032902) q[3];
sx q[3];
rz(0.75002128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1206104) q[2];
sx q[2];
rz(-1.1889428) q[2];
sx q[2];
rz(-0.19213842) q[2];
rz(-0.7297248) q[3];
sx q[3];
rz(-2.9967522) q[3];
sx q[3];
rz(-0.63989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9925053) q[0];
sx q[0];
rz(-2.3893116) q[0];
sx q[0];
rz(-1.3080904) q[0];
rz(-0.84450841) q[1];
sx q[1];
rz(-1.6788071) q[1];
sx q[1];
rz(0.62215296) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34595153) q[0];
sx q[0];
rz(-1.1591732) q[0];
sx q[0];
rz(-0.53070416) q[0];
x q[1];
rz(-1.6255369) q[2];
sx q[2];
rz(-1.3856882) q[2];
sx q[2];
rz(0.4601882) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6250861) q[1];
sx q[1];
rz(-1.7343905) q[1];
sx q[1];
rz(1.6229936) q[1];
rz(-pi) q[2];
rz(2.801311) q[3];
sx q[3];
rz(-2.8064686) q[3];
sx q[3];
rz(-1.3763031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1239803) q[2];
sx q[2];
rz(-1.7698741) q[2];
sx q[2];
rz(0.89517108) q[2];
rz(2.7923942) q[3];
sx q[3];
rz(-0.57217351) q[3];
sx q[3];
rz(2.9077742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.8643841) q[0];
sx q[0];
rz(-1.207749) q[0];
sx q[0];
rz(0.54339093) q[0];
rz(-0.1132938) q[1];
sx q[1];
rz(-1.7430867) q[1];
sx q[1];
rz(-1.2921804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77063507) q[0];
sx q[0];
rz(-1.9538625) q[0];
sx q[0];
rz(0.84324734) q[0];
x q[1];
rz(-1.0903484) q[2];
sx q[2];
rz(-0.70864973) q[2];
sx q[2];
rz(2.6330269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1067321) q[1];
sx q[1];
rz(-2.1550064) q[1];
sx q[1];
rz(0.75967914) q[1];
x q[2];
rz(-2.4025147) q[3];
sx q[3];
rz(-1.0024286) q[3];
sx q[3];
rz(-3.0773602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8285404) q[2];
sx q[2];
rz(-1.7834168) q[2];
sx q[2];
rz(0.48946112) q[2];
rz(-2.475907) q[3];
sx q[3];
rz(-0.61167115) q[3];
sx q[3];
rz(0.78589511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996443) q[0];
sx q[0];
rz(-2.2388832) q[0];
sx q[0];
rz(-0.06074252) q[0];
rz(-1.863106) q[1];
sx q[1];
rz(-0.91438952) q[1];
sx q[1];
rz(-1.0260822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36878186) q[0];
sx q[0];
rz(-2.368481) q[0];
sx q[0];
rz(-0.85526534) q[0];
x q[1];
rz(-2.0572971) q[2];
sx q[2];
rz(-0.62556534) q[2];
sx q[2];
rz(-1.6589669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85426165) q[1];
sx q[1];
rz(-1.8214487) q[1];
sx q[1];
rz(-1.7118368) q[1];
rz(-pi) q[2];
rz(1.8965263) q[3];
sx q[3];
rz(-2.0754078) q[3];
sx q[3];
rz(-2.4852999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8817899) q[2];
sx q[2];
rz(-1.5564352) q[2];
sx q[2];
rz(3.0628824) q[2];
rz(-1.6591266) q[3];
sx q[3];
rz(-2.8251298) q[3];
sx q[3];
rz(-0.44221529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.215613) q[0];
sx q[0];
rz(-2.7613566) q[0];
sx q[0];
rz(1.6967787) q[0];
rz(0.80698693) q[1];
sx q[1];
rz(-1.746256) q[1];
sx q[1];
rz(3.1296465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1211469) q[0];
sx q[0];
rz(-2.5021126) q[0];
sx q[0];
rz(-1.6156275) q[0];
rz(-pi) q[1];
rz(-2.9482466) q[2];
sx q[2];
rz(-0.51562658) q[2];
sx q[2];
rz(2.8597997) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1216189) q[1];
sx q[1];
rz(-1.498292) q[1];
sx q[1];
rz(0.95901476) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9868287) q[3];
sx q[3];
rz(-1.1763363) q[3];
sx q[3];
rz(0.607777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5936467) q[2];
sx q[2];
rz(-0.30210364) q[2];
sx q[2];
rz(-2.3054403) q[2];
rz(0.052308403) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(-0.93594319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1702105) q[0];
sx q[0];
rz(-2.2844071) q[0];
sx q[0];
rz(-0.49825391) q[0];
rz(-2.1360629) q[1];
sx q[1];
rz(-1.4509095) q[1];
sx q[1];
rz(0.11140579) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.191649) q[0];
sx q[0];
rz(-0.43131098) q[0];
sx q[0];
rz(1.3383207) q[0];
rz(-pi) q[1];
rz(1.3206215) q[2];
sx q[2];
rz(-1.9873575) q[2];
sx q[2];
rz(0.60384679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.028424358) q[1];
sx q[1];
rz(-2.3481124) q[1];
sx q[1];
rz(-1.9774578) q[1];
rz(-0.34078719) q[3];
sx q[3];
rz(-0.36431815) q[3];
sx q[3];
rz(-2.7076569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9446543) q[2];
sx q[2];
rz(-1.780218) q[2];
sx q[2];
rz(-1.4062175) q[2];
rz(1.1574636) q[3];
sx q[3];
rz(-0.99298733) q[3];
sx q[3];
rz(1.5895313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0208397) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(-1.7027759) q[0];
rz(-2.2976047) q[1];
sx q[1];
rz(-1.8344717) q[1];
sx q[1];
rz(-2.9356975) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015259764) q[0];
sx q[0];
rz(-1.5301989) q[0];
sx q[0];
rz(0.6313398) q[0];
rz(-1.6240142) q[2];
sx q[2];
rz(-2.2578265) q[2];
sx q[2];
rz(1.991697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.89005) q[1];
sx q[1];
rz(-1.492842) q[1];
sx q[1];
rz(-0.13549094) q[1];
x q[2];
rz(-0.67365065) q[3];
sx q[3];
rz(-2.5360772) q[3];
sx q[3];
rz(-0.57693627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0042808) q[2];
sx q[2];
rz(-1.0575123) q[2];
sx q[2];
rz(1.6005969) q[2];
rz(2.6565523) q[3];
sx q[3];
rz(-1.78616) q[3];
sx q[3];
rz(0.11390991) q[3];
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
rz(-1.4892905) q[0];
sx q[0];
rz(-0.052209608) q[0];
sx q[0];
rz(1.2963649) q[0];
rz(-1.0507874) q[1];
sx q[1];
rz(-1.6810345) q[1];
sx q[1];
rz(2.5349862) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8243931) q[0];
sx q[0];
rz(-1.6034295) q[0];
sx q[0];
rz(1.4072818) q[0];
x q[1];
rz(2.9398389) q[2];
sx q[2];
rz(-2.5725908) q[2];
sx q[2];
rz(0.099009939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7672528) q[1];
sx q[1];
rz(-1.8757959) q[1];
sx q[1];
rz(-1.5067596) q[1];
rz(-0.38207558) q[3];
sx q[3];
rz(-0.30508546) q[3];
sx q[3];
rz(-2.8918653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4198833) q[2];
sx q[2];
rz(-2.0100644) q[2];
sx q[2];
rz(-2.4427872) q[2];
rz(1.4612259) q[3];
sx q[3];
rz(-1.0137839) q[3];
sx q[3];
rz(-0.47718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2496495) q[0];
sx q[0];
rz(-1.5105381) q[0];
sx q[0];
rz(1.426209) q[0];
rz(-2.7632948) q[1];
sx q[1];
rz(-2.5143647) q[1];
sx q[1];
rz(-0.73077269) q[1];
rz(-1.0964285) q[2];
sx q[2];
rz(-1.8468936) q[2];
sx q[2];
rz(1.5441128) q[2];
rz(1.6017687) q[3];
sx q[3];
rz(-2.1270558) q[3];
sx q[3];
rz(1.1598641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
