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
rz(-1.191782) q[0];
sx q[0];
rz(-1.3718104) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(5.4159309) q[1];
sx q[1];
rz(9.8210788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21090487) q[0];
sx q[0];
rz(-0.22432029) q[0];
sx q[0];
rz(2.7773662) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9464067) q[2];
sx q[2];
rz(-1.7306788) q[2];
sx q[2];
rz(-0.92354666) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.635879) q[1];
sx q[1];
rz(-1.6649655) q[1];
sx q[1];
rz(1.0416743) q[1];
rz(-pi) q[2];
rz(-2.5175573) q[3];
sx q[3];
rz(-0.47500941) q[3];
sx q[3];
rz(1.6866682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1402011) q[2];
sx q[2];
rz(-0.80696693) q[2];
sx q[2];
rz(-1.3414471) q[2];
rz(-1.0691102) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.7599546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57217252) q[0];
sx q[0];
rz(-0.91745806) q[0];
sx q[0];
rz(-0.65895748) q[0];
rz(-0.072703687) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(-1.2603849) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9958268) q[0];
sx q[0];
rz(-0.46016177) q[0];
sx q[0];
rz(1.7547248) q[0];
x q[1];
rz(-0.44620338) q[2];
sx q[2];
rz(-1.1006736) q[2];
sx q[2];
rz(1.1371431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.61140984) q[1];
sx q[1];
rz(-1.5059222) q[1];
sx q[1];
rz(-0.26574175) q[1];
rz(-pi) q[2];
rz(3.0218868) q[3];
sx q[3];
rz(-2.1781892) q[3];
sx q[3];
rz(0.23609438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8947123) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(-1.3512208) q[2];
rz(-0.61659914) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(2.8504347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5093444) q[0];
sx q[0];
rz(-0.46724874) q[0];
sx q[0];
rz(-2.9314991) q[0];
rz(3.0585152) q[1];
sx q[1];
rz(-1.9030842) q[1];
sx q[1];
rz(3.074379) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7930896) q[0];
sx q[0];
rz(-1.5572892) q[0];
sx q[0];
rz(-1.5828787) q[0];
rz(1.1804586) q[2];
sx q[2];
rz(-1.3403041) q[2];
sx q[2];
rz(1.0437488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3785079) q[1];
sx q[1];
rz(-0.77662599) q[1];
sx q[1];
rz(1.08294) q[1];
rz(-pi) q[2];
rz(3.0359269) q[3];
sx q[3];
rz(-0.70836954) q[3];
sx q[3];
rz(-1.8745223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82789603) q[2];
sx q[2];
rz(-1.4915497) q[2];
sx q[2];
rz(1.776604) q[2];
rz(0.40417534) q[3];
sx q[3];
rz(-1.8242691) q[3];
sx q[3];
rz(-1.1990168) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86554027) q[0];
sx q[0];
rz(-1.4006389) q[0];
sx q[0];
rz(-2.6633967) q[0];
rz(-0.95282355) q[1];
sx q[1];
rz(-0.15004221) q[1];
sx q[1];
rz(0.40963867) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3296369) q[0];
sx q[0];
rz(-2.4853737) q[0];
sx q[0];
rz(0.043270525) q[0];
x q[1];
rz(2.9443789) q[2];
sx q[2];
rz(-0.46425113) q[2];
sx q[2];
rz(-2.0418389) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19189056) q[1];
sx q[1];
rz(-0.97964761) q[1];
sx q[1];
rz(2.3670235) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8919935) q[3];
sx q[3];
rz(-1.4347335) q[3];
sx q[3];
rz(-0.67067671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5459368) q[2];
sx q[2];
rz(-1.4384392) q[2];
sx q[2];
rz(2.0036009) q[2];
rz(-0.99500895) q[3];
sx q[3];
rz(-0.55750877) q[3];
sx q[3];
rz(-3.0549684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74395162) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(-0.27807903) q[0];
rz(1.2644348) q[1];
sx q[1];
rz(-1.2828628) q[1];
sx q[1];
rz(1.5778731) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9362617) q[0];
sx q[0];
rz(-0.38547984) q[0];
sx q[0];
rz(-1.9140713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78863849) q[2];
sx q[2];
rz(-1.3099726) q[2];
sx q[2];
rz(-2.3344085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1335781) q[1];
sx q[1];
rz(-0.99841324) q[1];
sx q[1];
rz(-1.7126669) q[1];
rz(1.044908) q[3];
sx q[3];
rz(-2.9581986) q[3];
sx q[3];
rz(0.50832752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2084674) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4472189) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(1.3391986) q[0];
rz(-3.0195492) q[1];
sx q[1];
rz(-1.3587147) q[1];
sx q[1];
rz(-0.36062127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80194762) q[0];
sx q[0];
rz(-2.3390798) q[0];
sx q[0];
rz(-0.54484493) q[0];
x q[1];
rz(-1.6595938) q[2];
sx q[2];
rz(-0.42150233) q[2];
sx q[2];
rz(2.5777566) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0619658) q[1];
sx q[1];
rz(-2.6089416) q[1];
sx q[1];
rz(-2.1497186) q[1];
rz(-pi) q[2];
rz(0.25288323) q[3];
sx q[3];
rz(-1.5656302) q[3];
sx q[3];
rz(2.9239426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3660761) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(-1.1386846) q[2];
rz(0.25389296) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(0.74786389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89895407) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(2.1658072) q[0];
rz(-1.1294533) q[1];
sx q[1];
rz(-2.6934846) q[1];
sx q[1];
rz(1.7879558) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0746552) q[0];
sx q[0];
rz(-0.88464175) q[0];
sx q[0];
rz(-1.1891436) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7464306) q[2];
sx q[2];
rz(-0.43000107) q[2];
sx q[2];
rz(-1.912009) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8541214) q[1];
sx q[1];
rz(-1.439289) q[1];
sx q[1];
rz(2.8033957) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46903543) q[3];
sx q[3];
rz(-1.9538662) q[3];
sx q[3];
rz(1.22016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5740616) q[2];
sx q[2];
rz(-2.3109544) q[2];
sx q[2];
rz(-0.91721025) q[2];
rz(-1.5926825) q[3];
sx q[3];
rz(-2.8062688) q[3];
sx q[3];
rz(-2.3622021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2061283) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(-3.1167378) q[0];
rz(1.286233) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(-1.53055) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2516305) q[0];
sx q[0];
rz(-1.1306835) q[0];
sx q[0];
rz(1.1391231) q[0];
rz(-pi) q[1];
rz(-0.31168657) q[2];
sx q[2];
rz(-0.71060813) q[2];
sx q[2];
rz(-0.59618261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2678044) q[1];
sx q[1];
rz(-0.28364983) q[1];
sx q[1];
rz(-3.0171418) q[1];
x q[2];
rz(2.1498411) q[3];
sx q[3];
rz(-0.25014057) q[3];
sx q[3];
rz(-2.5219805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.020869104) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(-1.2845767) q[2];
rz(-0.32127109) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(-2.9200714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-0.59182566) q[1];
sx q[1];
rz(-1.0706736) q[1];
sx q[1];
rz(-2.8048973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3580456) q[0];
sx q[0];
rz(-1.8382731) q[0];
sx q[0];
rz(-0.13937431) q[0];
rz(-pi) q[1];
rz(-2.8317802) q[2];
sx q[2];
rz(-0.37978077) q[2];
sx q[2];
rz(1.3889695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.241927) q[1];
sx q[1];
rz(-1.4414235) q[1];
sx q[1];
rz(-0.50921567) q[1];
rz(-2.982364) q[3];
sx q[3];
rz(-1.8223636) q[3];
sx q[3];
rz(-0.23345527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1153974) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(-2.3587312) q[2];
rz(0.10910263) q[3];
sx q[3];
rz(-2.1249873) q[3];
sx q[3];
rz(1.2469863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03610177) q[0];
sx q[0];
rz(-1.6721268) q[0];
sx q[0];
rz(-0.92392695) q[0];
rz(-0.52664122) q[1];
sx q[1];
rz(-1.8662607) q[1];
sx q[1];
rz(2.142876) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.963388) q[0];
sx q[0];
rz(-2.9792157) q[0];
sx q[0];
rz(2.8078662) q[0];
rz(-pi) q[1];
rz(1.2998215) q[2];
sx q[2];
rz(-2.8405115) q[2];
sx q[2];
rz(-1.4512514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0044788) q[1];
sx q[1];
rz(-2.5527337) q[1];
sx q[1];
rz(-1.6331255) q[1];
rz(-pi) q[2];
rz(2.6941006) q[3];
sx q[3];
rz(-2.4119141) q[3];
sx q[3];
rz(0.11451572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0052789) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(-0.33373731) q[2];
rz(-0.8606832) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(-3.1349643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6776047) q[0];
sx q[0];
rz(-1.9515568) q[0];
sx q[0];
rz(-2.4757181) q[0];
rz(0.55466501) q[1];
sx q[1];
rz(-1.8634836) q[1];
sx q[1];
rz(-2.9992933) q[1];
rz(-2.6449345) q[2];
sx q[2];
rz(-2.6500812) q[2];
sx q[2];
rz(-0.65190114) q[2];
rz(2.4246115) q[3];
sx q[3];
rz(-0.76631279) q[3];
sx q[3];
rz(-2.8240639) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
