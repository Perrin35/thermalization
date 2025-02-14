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
rz(1.0072768) q[0];
sx q[0];
rz(2.4957823) q[0];
sx q[0];
rz(9.0584005) q[0];
rz(-2.3925048) q[1];
sx q[1];
rz(-1.6269416) q[1];
sx q[1];
rz(0.14263022) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7932972) q[0];
sx q[0];
rz(-1.0401588) q[0];
sx q[0];
rz(-2.2133618) q[0];
rz(-pi) q[1];
rz(-0.22444522) q[2];
sx q[2];
rz(-1.7068752) q[2];
sx q[2];
rz(-2.4405757) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6444589) q[1];
sx q[1];
rz(-2.6503149) q[1];
sx q[1];
rz(1.3158448) q[1];
rz(-2.5589988) q[3];
sx q[3];
rz(-2.5911463) q[3];
sx q[3];
rz(-0.96874505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35965219) q[2];
sx q[2];
rz(-0.47986978) q[2];
sx q[2];
rz(1.5521607) q[2];
rz(1.747067) q[3];
sx q[3];
rz(-2.7701869) q[3];
sx q[3];
rz(-2.4594405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29779103) q[0];
sx q[0];
rz(-0.60972917) q[0];
sx q[0];
rz(-0.30526701) q[0];
rz(-1.8318532) q[1];
sx q[1];
rz(-0.42116183) q[1];
sx q[1];
rz(3.0895244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0355988) q[0];
sx q[0];
rz(-1.6806423) q[0];
sx q[0];
rz(0.1479827) q[0];
rz(-pi) q[1];
rz(1.270625) q[2];
sx q[2];
rz(-1.9294293) q[2];
sx q[2];
rz(1.7110362) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1250097) q[1];
sx q[1];
rz(-0.31453153) q[1];
sx q[1];
rz(2.1279863) q[1];
rz(-pi) q[2];
rz(1.7379187) q[3];
sx q[3];
rz(-2.4666836) q[3];
sx q[3];
rz(3.1233146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0569968) q[2];
sx q[2];
rz(-1.8085542) q[2];
sx q[2];
rz(1.4196654) q[2];
rz(-2.2325884) q[3];
sx q[3];
rz(-0.82908583) q[3];
sx q[3];
rz(-1.7349294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97515714) q[0];
sx q[0];
rz(-1.9571914) q[0];
sx q[0];
rz(1.3982406) q[0];
rz(0.94683975) q[1];
sx q[1];
rz(-2.0473174) q[1];
sx q[1];
rz(1.2907226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4104267) q[0];
sx q[0];
rz(-1.7058347) q[0];
sx q[0];
rz(0.076083994) q[0];
x q[1];
rz(0.78611908) q[2];
sx q[2];
rz(-1.2183471) q[2];
sx q[2];
rz(0.64482433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5980893) q[1];
sx q[1];
rz(-0.072797983) q[1];
sx q[1];
rz(-1.1612215) q[1];
x q[2];
rz(2.7644058) q[3];
sx q[3];
rz(-2.6335149) q[3];
sx q[3];
rz(2.7722904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1281841) q[2];
sx q[2];
rz(-0.42542294) q[2];
sx q[2];
rz(-0.95995963) q[2];
rz(-0.32402447) q[3];
sx q[3];
rz(-0.5158546) q[3];
sx q[3];
rz(2.7014151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0078761) q[0];
sx q[0];
rz(-0.77819264) q[0];
sx q[0];
rz(0.20009759) q[0];
rz(-2.5307185) q[1];
sx q[1];
rz(-0.69823825) q[1];
sx q[1];
rz(-0.79455882) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.824423) q[0];
sx q[0];
rz(-0.18704913) q[0];
sx q[0];
rz(-0.6424128) q[0];
x q[1];
rz(-1.346454) q[2];
sx q[2];
rz(-0.68549978) q[2];
sx q[2];
rz(-2.3239674) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17758372) q[1];
sx q[1];
rz(-1.7224041) q[1];
sx q[1];
rz(1.8880009) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4157905) q[3];
sx q[3];
rz(-2.0575541) q[3];
sx q[3];
rz(-0.85909708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8841298) q[2];
sx q[2];
rz(-0.56371671) q[2];
sx q[2];
rz(-1.5719315) q[2];
rz(-1.3563159) q[3];
sx q[3];
rz(-1.1607728) q[3];
sx q[3];
rz(0.13264382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2785579) q[0];
sx q[0];
rz(-0.25535169) q[0];
sx q[0];
rz(2.417946) q[0];
rz(-3.0408995) q[1];
sx q[1];
rz(-1.0483402) q[1];
sx q[1];
rz(0.06631276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1785051) q[0];
sx q[0];
rz(-1.6338991) q[0];
sx q[0];
rz(-1.0020688) q[0];
rz(2.5939717) q[2];
sx q[2];
rz(-2.1518618) q[2];
sx q[2];
rz(0.18085322) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1182013) q[1];
sx q[1];
rz(-1.7382335) q[1];
sx q[1];
rz(0.91806478) q[1];
rz(0.39855115) q[3];
sx q[3];
rz(-1.8722187) q[3];
sx q[3];
rz(2.4749746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2012607) q[2];
sx q[2];
rz(-1.2137493) q[2];
sx q[2];
rz(0.21855375) q[2];
rz(2.7471733) q[3];
sx q[3];
rz(-2.4557178) q[3];
sx q[3];
rz(-0.39943892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.047693096) q[0];
sx q[0];
rz(-2.221929) q[0];
sx q[0];
rz(-3.1374875) q[0];
rz(-0.63402367) q[1];
sx q[1];
rz(-1.6396294) q[1];
sx q[1];
rz(-2.1852469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30247453) q[0];
sx q[0];
rz(-1.7540521) q[0];
sx q[0];
rz(-1.6123346) q[0];
x q[1];
rz(-1.9955098) q[2];
sx q[2];
rz(-2.0419952) q[2];
sx q[2];
rz(-1.6933407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9117212) q[1];
sx q[1];
rz(-1.2930065) q[1];
sx q[1];
rz(-1.5650657) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31305571) q[3];
sx q[3];
rz(-2.7767973) q[3];
sx q[3];
rz(-0.49653253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6213106) q[2];
sx q[2];
rz(-1.7063528) q[2];
sx q[2];
rz(0.83462805) q[2];
rz(2.7708715) q[3];
sx q[3];
rz(-2.4407237) q[3];
sx q[3];
rz(2.4247775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653559) q[0];
sx q[0];
rz(-0.014572425) q[0];
sx q[0];
rz(-2.6915349) q[0];
rz(2.693148) q[1];
sx q[1];
rz(-0.83201718) q[1];
sx q[1];
rz(-2.4041596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0423575) q[0];
sx q[0];
rz(-2.6747534) q[0];
sx q[0];
rz(-1.1506266) q[0];
x q[1];
rz(0.040410553) q[2];
sx q[2];
rz(-1.4295661) q[2];
sx q[2];
rz(0.34858957) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0998675) q[1];
sx q[1];
rz(-1.8908354) q[1];
sx q[1];
rz(-2.7766224) q[1];
x q[2];
rz(-1.3025018) q[3];
sx q[3];
rz(-2.3327069) q[3];
sx q[3];
rz(-0.20352645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2927148) q[2];
sx q[2];
rz(-0.60056168) q[2];
sx q[2];
rz(2.4265477) q[2];
rz(2.6333366) q[3];
sx q[3];
rz(-1.6786989) q[3];
sx q[3];
rz(-1.7432632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4051064) q[0];
sx q[0];
rz(-0.54623258) q[0];
sx q[0];
rz(2.7832094) q[0];
rz(0.37250039) q[1];
sx q[1];
rz(-1.3595711) q[1];
sx q[1];
rz(0.43982664) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0459691) q[0];
sx q[0];
rz(-1.6120503) q[0];
sx q[0];
rz(2.3455363) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2829214) q[2];
sx q[2];
rz(-0.57980782) q[2];
sx q[2];
rz(-2.498954) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32268128) q[1];
sx q[1];
rz(-2.8214294) q[1];
sx q[1];
rz(0.59150954) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6251011) q[3];
sx q[3];
rz(-1.2479932) q[3];
sx q[3];
rz(1.0564547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0899352) q[2];
sx q[2];
rz(-3.0445485) q[2];
sx q[2];
rz(0.70866054) q[2];
rz(-0.25157252) q[3];
sx q[3];
rz(-2.1422062) q[3];
sx q[3];
rz(-3.1075509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.2947023) q[0];
sx q[0];
rz(-2.1069694) q[0];
sx q[0];
rz(-2.5501472) q[0];
rz(3.1239608) q[1];
sx q[1];
rz(-1.0523187) q[1];
sx q[1];
rz(3.0406521) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8629193) q[0];
sx q[0];
rz(-1.2569208) q[0];
sx q[0];
rz(0.14213965) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5048961) q[2];
sx q[2];
rz(-1.5392188) q[2];
sx q[2];
rz(2.4898138) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18848914) q[1];
sx q[1];
rz(-1.3425867) q[1];
sx q[1];
rz(-3.0191849) q[1];
x q[2];
rz(-1.798257) q[3];
sx q[3];
rz(-0.69256562) q[3];
sx q[3];
rz(2.8253959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4883604) q[2];
sx q[2];
rz(-2.4182352) q[2];
sx q[2];
rz(-1.4867268) q[2];
rz(-0.15445736) q[3];
sx q[3];
rz(-1.1052701) q[3];
sx q[3];
rz(2.7753593) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3514997) q[0];
sx q[0];
rz(-2.9993151) q[0];
sx q[0];
rz(2.067814) q[0];
rz(2.0692661) q[1];
sx q[1];
rz(-0.25389478) q[1];
sx q[1];
rz(-0.059018746) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85262839) q[0];
sx q[0];
rz(-1.5375332) q[0];
sx q[0];
rz(-2.9423174) q[0];
rz(1.0605766) q[2];
sx q[2];
rz(-1.3395044) q[2];
sx q[2];
rz(2.9289322) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7714951) q[1];
sx q[1];
rz(-2.1019249) q[1];
sx q[1];
rz(-0.70651502) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87484545) q[3];
sx q[3];
rz(-0.84934649) q[3];
sx q[3];
rz(1.9375999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3286288) q[2];
sx q[2];
rz(-1.2920516) q[2];
sx q[2];
rz(2.9126419) q[2];
rz(1.1936584) q[3];
sx q[3];
rz(-0.3242068) q[3];
sx q[3];
rz(1.4540023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1345632) q[0];
sx q[0];
rz(-0.80637359) q[0];
sx q[0];
rz(-0.73643186) q[0];
rz(-1.7560584) q[1];
sx q[1];
rz(-1.7693188) q[1];
sx q[1];
rz(1.7973695) q[1];
rz(-2.110021) q[2];
sx q[2];
rz(-1.5466718) q[2];
sx q[2];
rz(-0.9654733) q[2];
rz(2.8577639) q[3];
sx q[3];
rz(-2.5093439) q[3];
sx q[3];
rz(0.96998246) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
