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
rz(-3.1373625) q[0];
sx q[0];
rz(-2.5047996) q[0];
sx q[0];
rz(1.1871583) q[0];
rz(4.123426) q[1];
sx q[1];
rz(0.47458664) q[1];
sx q[1];
rz(8.4835806) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5309187) q[0];
sx q[0];
rz(-2.0661021) q[0];
sx q[0];
rz(-1.5727059) q[0];
x q[1];
rz(-2.1332987) q[2];
sx q[2];
rz(-2.9037711) q[2];
sx q[2];
rz(0.49006762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6189673) q[1];
sx q[1];
rz(-1.9272447) q[1];
sx q[1];
rz(2.5438479) q[1];
rz(-pi) q[2];
rz(3.0676316) q[3];
sx q[3];
rz(-1.1152905) q[3];
sx q[3];
rz(0.8616254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59419751) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(-1.0198062) q[2];
rz(-2.8504573) q[3];
sx q[3];
rz(-2.9499493) q[3];
sx q[3];
rz(1.8118793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6137467) q[0];
sx q[0];
rz(-0.83714038) q[0];
sx q[0];
rz(1.5035195) q[0];
rz(1.232052) q[1];
sx q[1];
rz(-0.82545311) q[1];
sx q[1];
rz(-1.8881316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6604174) q[0];
sx q[0];
rz(-1.7965797) q[0];
sx q[0];
rz(1.1522473) q[0];
rz(-pi) q[1];
rz(-2.6297036) q[2];
sx q[2];
rz(-1.8266668) q[2];
sx q[2];
rz(-0.40718174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4644891) q[1];
sx q[1];
rz(-1.3687412) q[1];
sx q[1];
rz(1.7176164) q[1];
rz(-2.0157865) q[3];
sx q[3];
rz(-2.2644694) q[3];
sx q[3];
rz(1.3563658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2848009) q[2];
sx q[2];
rz(-0.54232001) q[2];
sx q[2];
rz(2.8112603) q[2];
rz(-0.10279837) q[3];
sx q[3];
rz(-1.2448575) q[3];
sx q[3];
rz(-0.61906329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61515808) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(1.0636299) q[0];
rz(-0.10387736) q[1];
sx q[1];
rz(-0.79236284) q[1];
sx q[1];
rz(1.1054543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1310467) q[0];
sx q[0];
rz(-1.5233114) q[0];
sx q[0];
rz(0.032104339) q[0];
x q[1];
rz(2.3850307) q[2];
sx q[2];
rz(-1.3420449) q[2];
sx q[2];
rz(0.33332664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4729541) q[1];
sx q[1];
rz(-2.611043) q[1];
sx q[1];
rz(1.8133241) q[1];
rz(1.6126339) q[3];
sx q[3];
rz(-1.5825019) q[3];
sx q[3];
rz(0.042543471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0217454) q[2];
sx q[2];
rz(-0.80588078) q[2];
sx q[2];
rz(2.9160807) q[2];
rz(1.3269904) q[3];
sx q[3];
rz(-0.83695379) q[3];
sx q[3];
rz(0.64083159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91318146) q[0];
sx q[0];
rz(-1.5786194) q[0];
sx q[0];
rz(-0.6382424) q[0];
rz(-1.111521) q[1];
sx q[1];
rz(-2.1941954) q[1];
sx q[1];
rz(-2.9543455) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53786862) q[0];
sx q[0];
rz(-0.6938254) q[0];
sx q[0];
rz(-2.0286454) q[0];
rz(-pi) q[1];
rz(1.6333605) q[2];
sx q[2];
rz(-0.79088941) q[2];
sx q[2];
rz(-3.084201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61758047) q[1];
sx q[1];
rz(-2.1033546) q[1];
sx q[1];
rz(1.4790227) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71336915) q[3];
sx q[3];
rz(-1.6371048) q[3];
sx q[3];
rz(1.8153657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1049261) q[2];
sx q[2];
rz(-1.5035524) q[2];
sx q[2];
rz(0.46910134) q[2];
rz(0.30096287) q[3];
sx q[3];
rz(-2.8807237) q[3];
sx q[3];
rz(1.4891967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8982573) q[0];
sx q[0];
rz(-1.5168334) q[0];
sx q[0];
rz(-2.6014056) q[0];
rz(-0.46145269) q[1];
sx q[1];
rz(-1.6199473) q[1];
sx q[1];
rz(0.25925055) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.777939) q[0];
sx q[0];
rz(-2.724132) q[0];
sx q[0];
rz(-0.71939845) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44633003) q[2];
sx q[2];
rz(-2.0392373) q[2];
sx q[2];
rz(-2.6233752) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1446602) q[1];
sx q[1];
rz(-1.0029666) q[1];
sx q[1];
rz(-2.9564942) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2397887) q[3];
sx q[3];
rz(-1.3223303) q[3];
sx q[3];
rz(1.0120165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.79869142) q[2];
sx q[2];
rz(-0.68621245) q[2];
sx q[2];
rz(2.2314824) q[2];
rz(-0.56509194) q[3];
sx q[3];
rz(-1.3684973) q[3];
sx q[3];
rz(-2.4250987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8703576) q[0];
sx q[0];
rz(-2.7483181) q[0];
sx q[0];
rz(-1.2275335) q[0];
rz(2.6365989) q[1];
sx q[1];
rz(-2.125449) q[1];
sx q[1];
rz(1.9379157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1211091) q[0];
sx q[0];
rz(-1.6899037) q[0];
sx q[0];
rz(0.24863338) q[0];
rz(-pi) q[1];
rz(-0.9616379) q[2];
sx q[2];
rz(-2.3949902) q[2];
sx q[2];
rz(-0.37140977) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54654361) q[1];
sx q[1];
rz(-2.1814934) q[1];
sx q[1];
rz(0.86169075) q[1];
rz(1.1111383) q[3];
sx q[3];
rz(-1.8098003) q[3];
sx q[3];
rz(-2.8515332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0248854) q[2];
sx q[2];
rz(-1.025277) q[2];
sx q[2];
rz(-2.1825979) q[2];
rz(2.8420908) q[3];
sx q[3];
rz(-0.12226573) q[3];
sx q[3];
rz(1.6925156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2285948) q[0];
sx q[0];
rz(-1.8916425) q[0];
sx q[0];
rz(-2.6477497) q[0];
rz(-0.85810581) q[1];
sx q[1];
rz(-1.9789507) q[1];
sx q[1];
rz(1.3546622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0070902) q[0];
sx q[0];
rz(-1.0225774) q[0];
sx q[0];
rz(-0.28315084) q[0];
rz(-pi) q[1];
rz(1.0787215) q[2];
sx q[2];
rz(-1.4495696) q[2];
sx q[2];
rz(1.4722784) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87397447) q[1];
sx q[1];
rz(-2.3268891) q[1];
sx q[1];
rz(-1.6674629) q[1];
rz(-0.94791651) q[3];
sx q[3];
rz(-0.77548945) q[3];
sx q[3];
rz(0.38748565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6845503) q[2];
sx q[2];
rz(-2.5164618) q[2];
sx q[2];
rz(-0.54614145) q[2];
rz(-2.9438733) q[3];
sx q[3];
rz(-1.0309912) q[3];
sx q[3];
rz(2.8945967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.2208743) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(-2.0350631) q[0];
rz(-0.58440343) q[1];
sx q[1];
rz(-1.9677275) q[1];
sx q[1];
rz(-1.0027142) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2276777) q[0];
sx q[0];
rz(-1.7310237) q[0];
sx q[0];
rz(-1.1315617) q[0];
rz(-2.2539576) q[2];
sx q[2];
rz(-0.37345593) q[2];
sx q[2];
rz(-1.4591726) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7103191) q[1];
sx q[1];
rz(-0.51352507) q[1];
sx q[1];
rz(3.0679614) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2190621) q[3];
sx q[3];
rz(-2.1603843) q[3];
sx q[3];
rz(-2.0812579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57548412) q[2];
sx q[2];
rz(-2.2169952) q[2];
sx q[2];
rz(-1.3893348) q[2];
rz(-0.53650457) q[3];
sx q[3];
rz(-1.6335231) q[3];
sx q[3];
rz(-2.2947521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0480807) q[0];
sx q[0];
rz(-0.41963136) q[0];
sx q[0];
rz(-2.804948) q[0];
rz(3.0632784) q[1];
sx q[1];
rz(-1.2093733) q[1];
sx q[1];
rz(1.9858817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63756779) q[0];
sx q[0];
rz(-2.2462566) q[0];
sx q[0];
rz(-0.43677377) q[0];
rz(-pi) q[1];
rz(-1.9843654) q[2];
sx q[2];
rz(-1.6317657) q[2];
sx q[2];
rz(-2.3496036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.80336055) q[1];
sx q[1];
rz(-0.50932099) q[1];
sx q[1];
rz(1.2912911) q[1];
rz(-pi) q[2];
rz(2.8248057) q[3];
sx q[3];
rz(-1.5419349) q[3];
sx q[3];
rz(-2.3742418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0597421) q[2];
sx q[2];
rz(-1.0826702) q[2];
sx q[2];
rz(0.82320631) q[2];
rz(2.1246223) q[3];
sx q[3];
rz(-0.40139324) q[3];
sx q[3];
rz(3.1034191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.3391649) q[0];
sx q[0];
rz(-0.22295727) q[0];
sx q[0];
rz(1.4270576) q[0];
rz(-1.4786221) q[1];
sx q[1];
rz(-2.069811) q[1];
sx q[1];
rz(0.85085416) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6991147) q[0];
sx q[0];
rz(-0.48865396) q[0];
sx q[0];
rz(-1.3241934) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1673604) q[2];
sx q[2];
rz(-0.15359226) q[2];
sx q[2];
rz(-0.099791137) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46063328) q[1];
sx q[1];
rz(-2.6151028) q[1];
sx q[1];
rz(-0.10641702) q[1];
rz(-pi) q[2];
rz(3.0012742) q[3];
sx q[3];
rz(-1.5068973) q[3];
sx q[3];
rz(-1.4177223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.458638) q[2];
sx q[2];
rz(-1.3502324) q[2];
sx q[2];
rz(-2.3214935) q[2];
rz(1.594515) q[3];
sx q[3];
rz(-1.7087414) q[3];
sx q[3];
rz(-1.9471751) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1632669) q[0];
sx q[0];
rz(-1.6303202) q[0];
sx q[0];
rz(-1.6250961) q[0];
rz(-2.3691879) q[1];
sx q[1];
rz(-0.62320566) q[1];
sx q[1];
rz(-2.1355245) q[1];
rz(-2.8683581) q[2];
sx q[2];
rz(-0.71003503) q[2];
sx q[2];
rz(2.4207921) q[2];
rz(0.77843546) q[3];
sx q[3];
rz(-0.62494559) q[3];
sx q[3];
rz(0.86452053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
