OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(2.6846057) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(3.8222651) q[1];
sx q[1];
rz(8.1487976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9306643) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(-1.7095837) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.090473526) q[2];
sx q[2];
rz(-1.6455368) q[2];
sx q[2];
rz(0.59404101) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(3.0800746) q[1];
rz(-pi) q[2];
rz(-2.071322) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(-0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-0.16513744) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(-0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052658) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(3.0965565) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(-2.3388458) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.734056) q[0];
sx q[0];
rz(-1.7965839) q[0];
sx q[0];
rz(1.9682103) q[0];
rz(2.0741834) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(-2.3561321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3992577) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.6595112) q[1];
rz(-pi) q[2];
rz(2.2400168) q[3];
sx q[3];
rz(-1.672861) q[3];
sx q[3];
rz(-1.7166184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7654483) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(-1.7155898) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95214343) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(2.8564575) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(-1.3396938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0403255) q[0];
sx q[0];
rz(-2.2650532) q[0];
sx q[0];
rz(-0.56157748) q[0];
rz(-pi) q[1];
rz(-1.9402903) q[2];
sx q[2];
rz(-2.4174147) q[2];
sx q[2];
rz(-2.5232814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9078319) q[1];
sx q[1];
rz(-1.7572174) q[1];
sx q[1];
rz(-0.16252653) q[1];
x q[2];
rz(1.1717103) q[3];
sx q[3];
rz(-2.2200924) q[3];
sx q[3];
rz(0.89023512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.408067) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-3.0696707) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(-0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(0.66657153) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65581215) q[0];
sx q[0];
rz(-0.27944767) q[0];
sx q[0];
rz(2.4026472) q[0];
x q[1];
rz(-2.4273657) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(-0.54745882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92723144) q[1];
sx q[1];
rz(-0.9775369) q[1];
sx q[1];
rz(2.4638922) q[1];
rz(-pi) q[2];
rz(0.781226) q[3];
sx q[3];
rz(-0.20877148) q[3];
sx q[3];
rz(2.3695932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-0.56046265) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-2.3045325) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.9821092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54407185) q[0];
sx q[0];
rz(-0.87771767) q[0];
sx q[0];
rz(1.9835299) q[0];
rz(-0.97255623) q[2];
sx q[2];
rz(-1.4417357) q[2];
sx q[2];
rz(-1.8051403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24134484) q[1];
sx q[1];
rz(-0.7175788) q[1];
sx q[1];
rz(-0.89300968) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5043143) q[3];
sx q[3];
rz(-2.0162458) q[3];
sx q[3];
rz(0.74419903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(-0.064855382) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5746675) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(3.0976345) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(-0.10087068) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3374702) q[0];
sx q[0];
rz(-2.1967948) q[0];
sx q[0];
rz(-1.7581598) q[0];
rz(-pi) q[1];
rz(-0.51668824) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(-1.237243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3979891) q[1];
sx q[1];
rz(-1.7730224) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6359667) q[3];
sx q[3];
rz(-2.6856832) q[3];
sx q[3];
rz(1.0491919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-2.2163088) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(-1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15935414) q[0];
sx q[0];
rz(-1.237545) q[0];
sx q[0];
rz(0.23157816) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22600941) q[2];
sx q[2];
rz(-1.2603716) q[2];
sx q[2];
rz(-1.7058536) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5016985) q[1];
sx q[1];
rz(-1.5475029) q[1];
sx q[1];
rz(3.0122019) q[1];
x q[2];
rz(1.5958105) q[3];
sx q[3];
rz(-0.56846148) q[3];
sx q[3];
rz(2.3908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(0.98085105) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(3.0623073) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.1987196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90889701) q[0];
sx q[0];
rz(-0.43982738) q[0];
sx q[0];
rz(2.9798685) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.086026831) q[2];
sx q[2];
rz(-1.9326397) q[2];
sx q[2];
rz(-1.3081683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8504387) q[1];
sx q[1];
rz(-1.5209232) q[1];
sx q[1];
rz(-2.6384141) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49097455) q[3];
sx q[3];
rz(-2.6262865) q[3];
sx q[3];
rz(1.5332424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(-2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-2.5532706) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844922) q[0];
sx q[0];
rz(-1.4643837) q[0];
sx q[0];
rz(-1.6305001) q[0];
x q[1];
rz(-2.80745) q[2];
sx q[2];
rz(-1.9311937) q[2];
sx q[2];
rz(0.44597086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0578627) q[1];
sx q[1];
rz(-1.607556) q[1];
sx q[1];
rz(0.68395331) q[1];
rz(1.2133573) q[3];
sx q[3];
rz(-0.94944984) q[3];
sx q[3];
rz(-0.77222094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(2.771647) q[2];
rz(2.2423819) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(-2.1883011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.1870435) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(0.6533587) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95589721) q[0];
sx q[0];
rz(-0.16793185) q[0];
sx q[0];
rz(0.76783617) q[0];
x q[1];
rz(-2.3751971) q[2];
sx q[2];
rz(-0.99457127) q[2];
sx q[2];
rz(-3.0255084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2059584) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(2.2531855) q[1];
x q[2];
rz(0.023601836) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595173) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(0.72485483) q[2];
sx q[2];
rz(-0.58314322) q[2];
sx q[2];
rz(-2.4287139) q[2];
rz(2.1609859) q[3];
sx q[3];
rz(-1.4242314) q[3];
sx q[3];
rz(-0.84606708) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];