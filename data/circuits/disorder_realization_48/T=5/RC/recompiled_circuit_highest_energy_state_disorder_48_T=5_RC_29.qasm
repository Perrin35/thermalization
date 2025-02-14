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
rz(0.91520619) q[0];
sx q[0];
rz(-0.45867607) q[0];
sx q[0];
rz(0.31804481) q[0];
rz(7.6492352) q[1];
sx q[1];
rz(0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2554277) q[0];
sx q[0];
rz(-1.0647827) q[0];
sx q[0];
rz(0.54765986) q[0];
x q[1];
rz(-1.0452966) q[2];
sx q[2];
rz(-2.3901148) q[2];
sx q[2];
rz(2.1906302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1789238) q[1];
sx q[1];
rz(-2.0589863) q[1];
sx q[1];
rz(2.8250384) q[1];
rz(-pi) q[2];
rz(0.3955807) q[3];
sx q[3];
rz(-2.1115344) q[3];
sx q[3];
rz(2.9313425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.723145) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(-2.1212228) q[2];
rz(-2.4127035) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(1.6093904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6272524) q[0];
sx q[0];
rz(-1.1662551) q[0];
sx q[0];
rz(3.1033206) q[0];
rz(-0.12380869) q[1];
sx q[1];
rz(-0.47272155) q[1];
sx q[1];
rz(1.5708057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8077652) q[0];
sx q[0];
rz(-3.1196731) q[0];
sx q[0];
rz(-2.0448565) q[0];
rz(-pi) q[1];
rz(0.77124243) q[2];
sx q[2];
rz(-1.4995769) q[2];
sx q[2];
rz(-2.3598755) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1727708) q[1];
sx q[1];
rz(-1.5712196) q[1];
sx q[1];
rz(1.569773) q[1];
rz(-pi) q[2];
rz(-2.6255076) q[3];
sx q[3];
rz(-1.4654963) q[3];
sx q[3];
rz(0.95296958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52246419) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(-2.5627356) q[2];
rz(-2.0959496) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(2.3845909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66913644) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-2.6876167) q[0];
rz(2.9767735) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(-1.4705315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31904083) q[0];
sx q[0];
rz(-1.7421573) q[0];
sx q[0];
rz(1.1771855) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0963832) q[2];
sx q[2];
rz(-1.8586577) q[2];
sx q[2];
rz(-0.20341104) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0297894) q[1];
sx q[1];
rz(-1.8743427) q[1];
sx q[1];
rz(-2.4412267) q[1];
x q[2];
rz(2.4999714) q[3];
sx q[3];
rz(-2.372916) q[3];
sx q[3];
rz(-2.8856483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37956023) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(-1.7851768) q[2];
rz(2.6421269) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(2.0904026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22719638) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(1.0585349) q[0];
rz(-1.9805485) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(0.11722359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7331084) q[0];
sx q[0];
rz(-0.54382864) q[0];
sx q[0];
rz(3.0658019) q[0];
rz(2.2160596) q[2];
sx q[2];
rz(-1.7286073) q[2];
sx q[2];
rz(2.7881546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30018729) q[1];
sx q[1];
rz(-2.3267728) q[1];
sx q[1];
rz(-0.38736613) q[1];
rz(0.53026188) q[3];
sx q[3];
rz(-2.432309) q[3];
sx q[3];
rz(1.6146631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3566572) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(-1.3680722) q[2];
rz(1.28537) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(-2.4135597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0990937) q[0];
sx q[0];
rz(-2.2137401) q[0];
sx q[0];
rz(1.4121144) q[0];
rz(-2.1389351) q[1];
sx q[1];
rz(-2.899677) q[1];
sx q[1];
rz(1.5589421) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5413645) q[0];
sx q[0];
rz(-1.1523655) q[0];
sx q[0];
rz(0.61516841) q[0];
x q[1];
rz(-3.0465165) q[2];
sx q[2];
rz(-2.4209967) q[2];
sx q[2];
rz(-2.528185) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0274899) q[1];
sx q[1];
rz(-1.5472425) q[1];
sx q[1];
rz(2.1211758) q[1];
x q[2];
rz(0.094314625) q[3];
sx q[3];
rz(-1.8911165) q[3];
sx q[3];
rz(-1.8220779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.58854181) q[2];
sx q[2];
rz(-1.43247) q[2];
sx q[2];
rz(-0.39349619) q[2];
rz(2.1816175) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(2.1737607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680962) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(-0.18336503) q[0];
rz(0.49194899) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(1.9221745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26579612) q[0];
sx q[0];
rz(-1.4816493) q[0];
sx q[0];
rz(1.4899979) q[0];
rz(-pi) q[1];
rz(-1.9398795) q[2];
sx q[2];
rz(-3.0574692) q[2];
sx q[2];
rz(-0.28757986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9907836) q[1];
sx q[1];
rz(-1.1649302) q[1];
sx q[1];
rz(-1.1836786) q[1];
x q[2];
rz(1.4063356) q[3];
sx q[3];
rz(-2.6643848) q[3];
sx q[3];
rz(0.26613126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(-0.08237002) q[2];
rz(2.2149337) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(-1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6019186) q[0];
sx q[0];
rz(-0.32783666) q[0];
sx q[0];
rz(-2.4251921) q[0];
rz(2.7377103) q[1];
sx q[1];
rz(-1.5682181) q[1];
sx q[1];
rz(-1.7049888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8536896) q[0];
sx q[0];
rz(-0.40659621) q[0];
sx q[0];
rz(1.2000183) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7030969) q[2];
sx q[2];
rz(-1.0463593) q[2];
sx q[2];
rz(0.99614267) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0246504) q[1];
sx q[1];
rz(-2.7305909) q[1];
sx q[1];
rz(-2.890439) q[1];
rz(0.8088423) q[3];
sx q[3];
rz(-2.1385962) q[3];
sx q[3];
rz(-2.0140935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3639823) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(3.0762365) q[2];
rz(-2.2743547) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(-1.3245827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48677483) q[0];
sx q[0];
rz(-1.4210533) q[0];
sx q[0];
rz(0.73100334) q[0];
rz(-2.3664318) q[1];
sx q[1];
rz(-0.33439264) q[1];
sx q[1];
rz(2.7269272) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054798445) q[0];
sx q[0];
rz(-2.2761184) q[0];
sx q[0];
rz(2.2948625) q[0];
x q[1];
rz(2.6111905) q[2];
sx q[2];
rz(-1.3866405) q[2];
sx q[2];
rz(-0.14094409) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29293167) q[1];
sx q[1];
rz(-2.7457016) q[1];
sx q[1];
rz(-2.4421874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4067134) q[3];
sx q[3];
rz(-1.7030431) q[3];
sx q[3];
rz(-1.924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9097462) q[2];
sx q[2];
rz(-0.48603386) q[2];
sx q[2];
rz(3.0492142) q[2];
rz(2.3653024) q[3];
sx q[3];
rz(-1.4624566) q[3];
sx q[3];
rz(-0.20364729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6542776) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(-0.017729433) q[0];
rz(1.0461294) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(2.0808992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8413139) q[0];
sx q[0];
rz(-1.5567686) q[0];
sx q[0];
rz(0.068723516) q[0];
rz(-pi) q[1];
rz(-1.9477884) q[2];
sx q[2];
rz(-2.5193938) q[2];
sx q[2];
rz(2.5925328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5405766) q[1];
sx q[1];
rz(-0.81588826) q[1];
sx q[1];
rz(2.9040082) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52379108) q[3];
sx q[3];
rz(-1.6219211) q[3];
sx q[3];
rz(1.8845975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1549418) q[2];
sx q[2];
rz(-2.3255746) q[2];
sx q[2];
rz(2.4972534) q[2];
rz(2.2996357) q[3];
sx q[3];
rz(-1.569845) q[3];
sx q[3];
rz(2.2346066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7613206) q[0];
sx q[0];
rz(-2.705882) q[0];
sx q[0];
rz(-2.6897588) q[0];
rz(-1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(1.570328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752351) q[0];
sx q[0];
rz(-1.5646828) q[0];
sx q[0];
rz(-0.17019043) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0958448) q[2];
sx q[2];
rz(-1.6745076) q[2];
sx q[2];
rz(-0.46670676) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6856076) q[1];
sx q[1];
rz(-0.41049126) q[1];
sx q[1];
rz(-1.8985156) q[1];
rz(0.7223981) q[3];
sx q[3];
rz(-1.5766451) q[3];
sx q[3];
rz(-1.4129782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5181804) q[2];
sx q[2];
rz(-1.0717012) q[2];
sx q[2];
rz(1.5706459) q[2];
rz(3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(0.51641974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.327772) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(-2.9227921) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(1.1280288) q[2];
sx q[2];
rz(-1.3793066) q[2];
sx q[2];
rz(-2.1366091) q[2];
rz(-2.6811213) q[3];
sx q[3];
rz(-2.4302767) q[3];
sx q[3];
rz(-2.5828491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
