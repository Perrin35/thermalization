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
rz(0.49993604) q[0];
sx q[0];
rz(1.9498107) q[0];
sx q[0];
rz(10.796588) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(5.4159309) q[1];
sx q[1];
rz(9.8210788) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16184743) q[0];
sx q[0];
rz(-1.3614185) q[0];
sx q[0];
rz(-1.6518948) q[0];
x q[1];
rz(-1.9849586) q[2];
sx q[2];
rz(-2.7348619) q[2];
sx q[2];
rz(-2.1106281) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.010155023) q[1];
sx q[1];
rz(-1.0442631) q[1];
sx q[1];
rz(-0.10897691) q[1];
x q[2];
rz(1.2788775) q[3];
sx q[3];
rz(-1.9510428) q[3];
sx q[3];
rz(-1.0060665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1402011) q[2];
sx q[2];
rz(-2.3346257) q[2];
sx q[2];
rz(-1.8001455) q[2];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(-1.7599546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57217252) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(-0.65895748) q[0];
rz(3.068889) q[1];
sx q[1];
rz(-1.0560938) q[1];
sx q[1];
rz(1.2603849) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9958268) q[0];
sx q[0];
rz(-0.46016177) q[0];
sx q[0];
rz(1.3868679) q[0];
rz(-pi) q[1];
x q[1];
rz(1.057823) q[2];
sx q[2];
rz(-1.9657081) q[2];
sx q[2];
rz(-0.22024834) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61140984) q[1];
sx q[1];
rz(-1.6356704) q[1];
sx q[1];
rz(2.8758509) q[1];
rz(-0.96003344) q[3];
sx q[3];
rz(-1.472578) q[3];
sx q[3];
rz(1.4032422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2468804) q[2];
sx q[2];
rz(-0.74328819) q[2];
sx q[2];
rz(1.7903719) q[2];
rz(0.61659914) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.63224822) q[0];
sx q[0];
rz(-2.6743439) q[0];
sx q[0];
rz(0.21009357) q[0];
rz(3.0585152) q[1];
sx q[1];
rz(-1.9030842) q[1];
sx q[1];
rz(3.074379) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7602762) q[0];
sx q[0];
rz(-3.1234703) q[0];
sx q[0];
rz(2.4118547) q[0];
rz(-pi) q[1];
rz(-1.018173) q[2];
sx q[2];
rz(-2.6913096) q[2];
sx q[2];
rz(2.1076815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1695425) q[1];
sx q[1];
rz(-1.2360555) q[1];
sx q[1];
rz(0.85595815) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70560734) q[3];
sx q[3];
rz(-1.5021245) q[3];
sx q[3];
rz(-0.22335438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82789603) q[2];
sx q[2];
rz(-1.4915497) q[2];
sx q[2];
rz(1.3649887) q[2];
rz(-2.7374173) q[3];
sx q[3];
rz(-1.8242691) q[3];
sx q[3];
rz(-1.1990168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2760524) q[0];
sx q[0];
rz(-1.4006389) q[0];
sx q[0];
rz(2.6633967) q[0];
rz(-2.1887691) q[1];
sx q[1];
rz(-2.9915504) q[1];
sx q[1];
rz(0.40963867) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3296369) q[0];
sx q[0];
rz(-2.4853737) q[0];
sx q[0];
rz(-0.043270525) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9443789) q[2];
sx q[2];
rz(-2.6773415) q[2];
sx q[2];
rz(1.0997538) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.86120104) q[1];
sx q[1];
rz(-0.93547677) q[1];
sx q[1];
rz(-0.7648356) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1616108) q[3];
sx q[3];
rz(-0.34790643) q[3];
sx q[3];
rz(2.6285841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5459368) q[2];
sx q[2];
rz(-1.4384392) q[2];
sx q[2];
rz(1.1379918) q[2];
rz(-2.1465837) q[3];
sx q[3];
rz(-0.55750877) q[3];
sx q[3];
rz(-0.086624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.397641) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(-2.8635136) q[0];
rz(-1.2644348) q[1];
sx q[1];
rz(-1.8587298) q[1];
sx q[1];
rz(1.5778731) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9362617) q[0];
sx q[0];
rz(-2.7561128) q[0];
sx q[0];
rz(-1.9140713) q[0];
rz(-pi) q[1];
rz(-2.3529542) q[2];
sx q[2];
rz(-1.3099726) q[2];
sx q[2];
rz(-2.3344085) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63999004) q[1];
sx q[1];
rz(-1.6899365) q[1];
sx q[1];
rz(0.57699012) q[1];
rz(1.044908) q[3];
sx q[3];
rz(-2.9581986) q[3];
sx q[3];
rz(-2.6332651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9331253) q[2];
sx q[2];
rz(-1.5310023) q[2];
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
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4472189) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(1.8023941) q[0];
rz(-3.0195492) q[1];
sx q[1];
rz(-1.3587147) q[1];
sx q[1];
rz(2.7809714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0568389) q[0];
sx q[0];
rz(-0.90843102) q[0];
sx q[0];
rz(2.063089) q[0];
rz(-pi) q[1];
rz(1.6595938) q[2];
sx q[2];
rz(-0.42150233) q[2];
sx q[2];
rz(0.56383609) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0619658) q[1];
sx q[1];
rz(-2.6089416) q[1];
sx q[1];
rz(-2.1497186) q[1];
rz(0.020645647) q[3];
sx q[3];
rz(-0.25293487) q[3];
sx q[3];
rz(-1.7684574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3660761) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(1.1386846) q[2];
rz(-2.8876997) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(-2.3937288) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89895407) q[0];
sx q[0];
rz(-2.2518318) q[0];
sx q[0];
rz(2.1658072) q[0];
rz(-2.0121393) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(-1.3536369) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24714805) q[0];
sx q[0];
rz(-1.2784875) q[0];
sx q[0];
rz(-0.72299374) q[0];
rz(-pi) q[1];
rz(1.1466402) q[2];
sx q[2];
rz(-1.643702) q[2];
sx q[2];
rz(-0.50114252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8541214) q[1];
sx q[1];
rz(-1.439289) q[1];
sx q[1];
rz(-2.8033957) q[1];
rz(-pi) q[2];
rz(-0.46903543) q[3];
sx q[3];
rz(-1.9538662) q[3];
sx q[3];
rz(-1.9214326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56753105) q[2];
sx q[2];
rz(-0.83063829) q[2];
sx q[2];
rz(2.2243824) q[2];
rz(1.5489102) q[3];
sx q[3];
rz(-0.33532381) q[3];
sx q[3];
rz(2.3622021) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2061283) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(3.1167378) q[0];
rz(1.8553597) q[1];
sx q[1];
rz(-1.356946) q[1];
sx q[1];
rz(1.6110427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2516305) q[0];
sx q[0];
rz(-1.1306835) q[0];
sx q[0];
rz(2.0024695) q[0];
x q[1];
rz(2.8299061) q[2];
sx q[2];
rz(-2.4309845) q[2];
sx q[2];
rz(0.59618261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87378824) q[1];
sx q[1];
rz(-2.8579428) q[1];
sx q[1];
rz(3.0171418) q[1];
x q[2];
rz(2.1498411) q[3];
sx q[3];
rz(-2.8914521) q[3];
sx q[3];
rz(2.5219805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.020869104) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(-1.2845767) q[2];
rz(-2.8203216) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(2.9200714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1247509) q[0];
sx q[0];
rz(-1.1469954) q[0];
sx q[0];
rz(-0.91957134) q[0];
rz(0.59182566) q[1];
sx q[1];
rz(-2.070919) q[1];
sx q[1];
rz(-2.8048973) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29558173) q[0];
sx q[0];
rz(-2.8407556) q[0];
sx q[0];
rz(-1.1016125) q[0];
rz(-pi) q[1];
rz(0.36328237) q[2];
sx q[2];
rz(-1.4575301) q[2];
sx q[2];
rz(2.670778) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7408161) q[1];
sx q[1];
rz(-1.0662406) q[1];
sx q[1];
rz(-1.4228805) q[1];
rz(-pi) q[2];
rz(1.0180264) q[3];
sx q[3];
rz(-2.8447731) q[3];
sx q[3];
rz(-2.8021013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0261953) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(0.78286147) q[2];
rz(0.10910263) q[3];
sx q[3];
rz(-1.0166054) q[3];
sx q[3];
rz(-1.2469863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03610177) q[0];
sx q[0];
rz(-1.4694659) q[0];
sx q[0];
rz(-0.92392695) q[0];
rz(-0.52664122) q[1];
sx q[1];
rz(-1.275332) q[1];
sx q[1];
rz(-2.142876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0629421) q[0];
sx q[0];
rz(-1.6237769) q[0];
sx q[0];
rz(0.15356252) q[0];
x q[1];
rz(1.8417712) q[2];
sx q[2];
rz(-2.8405115) q[2];
sx q[2];
rz(-1.6903413) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6560581) q[1];
sx q[1];
rz(-1.5361934) q[1];
sx q[1];
rz(-2.1587579) q[1];
rz(2.4630354) q[3];
sx q[3];
rz(-1.863409) q[3];
sx q[3];
rz(1.34174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0052789) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(-2.8078553) q[2];
rz(-0.8606832) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(-3.1349643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6776047) q[0];
sx q[0];
rz(-1.1900359) q[0];
sx q[0];
rz(0.66587454) q[0];
rz(0.55466501) q[1];
sx q[1];
rz(-1.8634836) q[1];
sx q[1];
rz(-2.9992933) q[1];
rz(-1.3210422) q[2];
sx q[2];
rz(-1.142923) q[2];
sx q[2];
rz(-0.10071071) q[2];
rz(-1.0068245) q[3];
sx q[3];
rz(-1.020731) q[3];
sx q[3];
rz(1.1976477) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
