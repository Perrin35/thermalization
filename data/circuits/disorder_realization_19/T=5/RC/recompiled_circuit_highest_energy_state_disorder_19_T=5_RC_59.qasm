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
rz(-2.6416566) q[0];
sx q[0];
rz(-1.9498107) q[0];
sx q[0];
rz(-1.7697822) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(-0.86725441) q[1];
sx q[1];
rz(0.39630085) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7157519) q[0];
sx q[0];
rz(-1.4914728) q[0];
sx q[0];
rz(0.21004814) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1566341) q[2];
sx q[2];
rz(-2.7348619) q[2];
sx q[2];
rz(-1.0309645) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9170678) q[1];
sx q[1];
rz(-0.5366508) q[1];
sx q[1];
rz(1.3858252) q[1];
x q[2];
rz(0.39537786) q[3];
sx q[3];
rz(-1.3002743) q[3];
sx q[3];
rz(-0.45365712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1402011) q[2];
sx q[2];
rz(-0.80696693) q[2];
sx q[2];
rz(-1.3414471) q[2];
rz(2.0724824) q[3];
sx q[3];
rz(-2.9954698) q[3];
sx q[3];
rz(1.3816381) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57217252) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(-0.65895748) q[0];
rz(-3.068889) q[1];
sx q[1];
rz(-1.0560938) q[1];
sx q[1];
rz(1.8812077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2598686) q[0];
sx q[0];
rz(-1.6521075) q[0];
sx q[0];
rz(-2.0242244) q[0];
rz(-pi) q[1];
rz(0.44620338) q[2];
sx q[2];
rz(-2.0409191) q[2];
sx q[2];
rz(-2.0044495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94174332) q[1];
sx q[1];
rz(-1.3056271) q[1];
sx q[1];
rz(-1.6380235) q[1];
rz(-0.96003344) q[3];
sx q[3];
rz(-1.6690147) q[3];
sx q[3];
rz(-1.4032422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8947123) q[2];
sx q[2];
rz(-0.74328819) q[2];
sx q[2];
rz(1.3512208) q[2];
rz(-0.61659914) q[3];
sx q[3];
rz(-2.1117881) q[3];
sx q[3];
rz(-2.8504347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63224822) q[0];
sx q[0];
rz(-0.46724874) q[0];
sx q[0];
rz(0.21009357) q[0];
rz(-3.0585152) q[1];
sx q[1];
rz(-1.2385085) q[1];
sx q[1];
rz(-0.067213623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7930896) q[0];
sx q[0];
rz(-1.5843035) q[0];
sx q[0];
rz(-1.5587139) q[0];
x q[1];
rz(1.018173) q[2];
sx q[2];
rz(-2.6913096) q[2];
sx q[2];
rz(-2.1076815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1695425) q[1];
sx q[1];
rz(-1.2360555) q[1];
sx q[1];
rz(-2.2856345) q[1];
x q[2];
rz(0.10566575) q[3];
sx q[3];
rz(-2.4332231) q[3];
sx q[3];
rz(1.2670704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82789603) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(1.3649887) q[2];
rz(0.40417534) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(-1.9425758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2760524) q[0];
sx q[0];
rz(-1.7409538) q[0];
sx q[0];
rz(2.6633967) q[0];
rz(-2.1887691) q[1];
sx q[1];
rz(-2.9915504) q[1];
sx q[1];
rz(-2.731954) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8119558) q[0];
sx q[0];
rz(-2.4853737) q[0];
sx q[0];
rz(3.0983221) q[0];
rz(2.9443789) q[2];
sx q[2];
rz(-2.6773415) q[2];
sx q[2];
rz(2.0418389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9497021) q[1];
sx q[1];
rz(-2.161945) q[1];
sx q[1];
rz(0.77456919) q[1];
rz(-1.1616108) q[3];
sx q[3];
rz(-2.7936862) q[3];
sx q[3];
rz(0.51300853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5459368) q[2];
sx q[2];
rz(-1.7031534) q[2];
sx q[2];
rz(1.1379918) q[2];
rz(-0.99500895) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(-0.086624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74395162) q[0];
sx q[0];
rz(-1.3084363) q[0];
sx q[0];
rz(2.8635136) q[0];
rz(1.2644348) q[1];
sx q[1];
rz(-1.2828628) q[1];
sx q[1];
rz(-1.5637195) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0959543) q[0];
sx q[0];
rz(-1.6976893) q[0];
sx q[0];
rz(-1.2058099) q[0];
rz(1.9327963) q[2];
sx q[2];
rz(-0.81564192) q[2];
sx q[2];
rz(1.0175616) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3914096) q[1];
sx q[1];
rz(-2.5538008) q[1];
sx q[1];
rz(-0.21601917) q[1];
rz(1.4117354) q[3];
sx q[3];
rz(-1.6624699) q[3];
sx q[3];
rz(-1.5810458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2084674) q[2];
sx q[2];
rz(-1.6105904) q[2];
sx q[2];
rz(-0.36766407) q[2];
rz(1.095088) q[3];
sx q[3];
rz(-2.118066) q[3];
sx q[3];
rz(0.17525214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4472189) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(-1.8023941) q[0];
rz(-0.12204349) q[1];
sx q[1];
rz(-1.3587147) q[1];
sx q[1];
rz(-2.7809714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80194762) q[0];
sx q[0];
rz(-2.3390798) q[0];
sx q[0];
rz(-2.5967477) q[0];
rz(-pi) q[1];
rz(-1.9908268) q[2];
sx q[2];
rz(-1.5345062) q[2];
sx q[2];
rz(-2.0535713) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5721467) q[1];
sx q[1];
rz(-1.1317557) q[1];
sx q[1];
rz(2.8296058) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8887094) q[3];
sx q[3];
rz(-1.5759625) q[3];
sx q[3];
rz(-0.21765003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.7755166) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(-2.0029081) q[2];
rz(-2.8876997) q[3];
sx q[3];
rz(-0.64520276) q[3];
sx q[3];
rz(-0.74786389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.89895407) q[0];
sx q[0];
rz(-2.2518318) q[0];
sx q[0];
rz(-2.1658072) q[0];
rz(1.1294533) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(-1.3536369) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5022875) q[0];
sx q[0];
rz(-2.3717518) q[0];
sx q[0];
rz(2.7147311) q[0];
rz(3.0616272) q[2];
sx q[2];
rz(-1.1478394) q[2];
sx q[2];
rz(-2.1048196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5013988) q[1];
sx q[1];
rz(-0.36194776) q[1];
sx q[1];
rz(-2.7622403) q[1];
x q[2];
rz(-2.4134992) q[3];
sx q[3];
rz(-2.545176) q[3];
sx q[3];
rz(2.8567258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56753105) q[2];
sx q[2];
rz(-0.83063829) q[2];
sx q[2];
rz(2.2243824) q[2];
rz(1.5926825) q[3];
sx q[3];
rz(-0.33532381) q[3];
sx q[3];
rz(-2.3622021) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9354644) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(-3.1167378) q[0];
rz(1.8553597) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(1.53055) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0760374) q[0];
sx q[0];
rz(-0.60638705) q[0];
sx q[0];
rz(-0.72640149) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3127682) q[2];
sx q[2];
rz(-0.90094968) q[2];
sx q[2];
rz(-2.9473398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2678044) q[1];
sx q[1];
rz(-0.28364983) q[1];
sx q[1];
rz(3.0171418) q[1];
rz(-0.1389109) q[3];
sx q[3];
rz(-1.362097) q[3];
sx q[3];
rz(-3.1155966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.020869104) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(-1.2845767) q[2];
rz(2.8203216) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(0.2215213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016841737) q[0];
sx q[0];
rz(-1.1469954) q[0];
sx q[0];
rz(-2.2220213) q[0];
rz(-0.59182566) q[1];
sx q[1];
rz(-1.0706736) q[1];
sx q[1];
rz(0.33669534) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78354708) q[0];
sx q[0];
rz(-1.8382731) q[0];
sx q[0];
rz(-3.0022183) q[0];
rz(-1.6918963) q[2];
sx q[2];
rz(-1.9316439) q[2];
sx q[2];
rz(-1.0570414) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8996657) q[1];
sx q[1];
rz(-1.4414235) q[1];
sx q[1];
rz(-2.632377) q[1];
rz(-pi) q[2];
rz(-2.982364) q[3];
sx q[3];
rz(-1.8223636) q[3];
sx q[3];
rz(2.9081374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0261953) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(-0.78286147) q[2];
rz(0.10910263) q[3];
sx q[3];
rz(-1.0166054) q[3];
sx q[3];
rz(1.8946064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03610177) q[0];
sx q[0];
rz(-1.6721268) q[0];
sx q[0];
rz(2.2176657) q[0];
rz(0.52664122) q[1];
sx q[1];
rz(-1.275332) q[1];
sx q[1];
rz(-0.99871666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0629421) q[0];
sx q[0];
rz(-1.5178158) q[0];
sx q[0];
rz(0.15356252) q[0];
rz(-pi) q[1];
rz(-1.8615103) q[2];
sx q[2];
rz(-1.6502585) q[2];
sx q[2];
rz(-0.13979152) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1371138) q[1];
sx q[1];
rz(-2.5527337) q[1];
sx q[1];
rz(1.5084672) q[1];
rz(-pi) q[2];
rz(-0.67855723) q[3];
sx q[3];
rz(-1.863409) q[3];
sx q[3];
rz(-1.7998526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1363137) q[2];
sx q[2];
rz(-2.2515991) q[2];
sx q[2];
rz(0.33373731) q[2];
rz(-0.8606832) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(-3.1349643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46398791) q[0];
sx q[0];
rz(-1.1900359) q[0];
sx q[0];
rz(0.66587454) q[0];
rz(-2.5869276) q[1];
sx q[1];
rz(-1.8634836) q[1];
sx q[1];
rz(-2.9992933) q[1];
rz(-1.3210422) q[2];
sx q[2];
rz(-1.142923) q[2];
sx q[2];
rz(-0.10071071) q[2];
rz(2.5139204) q[3];
sx q[3];
rz(-1.097642) q[3];
sx q[3];
rz(-0.69246694) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
