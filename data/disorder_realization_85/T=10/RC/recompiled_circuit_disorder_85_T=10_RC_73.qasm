OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(-2.5512295) q[0];
sx q[0];
rz(-0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57173079) q[0];
sx q[0];
rz(-0.96643448) q[0];
sx q[0];
rz(-2.5992924) q[0];
x q[1];
rz(-1.0797834) q[2];
sx q[2];
rz(-0.91677374) q[2];
sx q[2];
rz(-0.71066463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7264759) q[1];
sx q[1];
rz(-1.2435902) q[1];
sx q[1];
rz(0.16023689) q[1];
rz(-pi) q[2];
rz(-0.052470603) q[3];
sx q[3];
rz(-0.82004181) q[3];
sx q[3];
rz(-2.6657871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(-2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.20724021) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(-0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-2.4172799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121719) q[0];
sx q[0];
rz(-0.90733007) q[0];
sx q[0];
rz(2.7396766) q[0];
rz(1.6547336) q[2];
sx q[2];
rz(-1.3720023) q[2];
sx q[2];
rz(-2.4059911) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38072941) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(2.1058583) q[1];
rz(2.904326) q[3];
sx q[3];
rz(-1.3788584) q[3];
sx q[3];
rz(2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(3.0959685) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.2190855) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1299382) q[0];
sx q[0];
rz(-1.219795) q[0];
sx q[0];
rz(-2.8610693) q[0];
x q[1];
rz(1.2246386) q[2];
sx q[2];
rz(-2.4097754) q[2];
sx q[2];
rz(-0.099230448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8428426) q[1];
sx q[1];
rz(-1.0423653) q[1];
sx q[1];
rz(0.28038402) q[1];
rz(-1.838802) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(0.55004317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(2.6233327) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(-0.82675654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0543538) q[0];
sx q[0];
rz(-1.0687437) q[0];
sx q[0];
rz(-2.0421844) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3739359) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(1.3354288) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0930867) q[1];
sx q[1];
rz(-0.48492453) q[1];
sx q[1];
rz(-2.4946458) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2590253) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(1.242897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(1.6453843) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(2.856423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171705) q[0];
sx q[0];
rz(-1.6132014) q[0];
sx q[0];
rz(2.0457343) q[0];
rz(-pi) q[1];
rz(-1.1826913) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(-2.2270122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4080216) q[1];
sx q[1];
rz(-1.5712149) q[1];
sx q[1];
rz(-1.8838521) q[1];
rz(-pi) q[2];
rz(-2.7399555) q[3];
sx q[3];
rz(-1.65997) q[3];
sx q[3];
rz(3.0107486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-2.6021393) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(0.15701292) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(1.3670115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.052664) q[0];
sx q[0];
rz(-0.0757218) q[0];
sx q[0];
rz(-1.1195539) q[0];
x q[1];
rz(2.5622257) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(0.52057779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94783084) q[1];
sx q[1];
rz(-1.8596706) q[1];
sx q[1];
rz(3.0572901) q[1];
x q[2];
rz(0.35297024) q[3];
sx q[3];
rz(-1.4735231) q[3];
sx q[3];
rz(0.73247611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(0.40346754) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24213174) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(1.1706932) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255575) q[0];
sx q[0];
rz(-3.1177525) q[0];
sx q[0];
rz(-0.26160474) q[0];
rz(-pi) q[1];
rz(-1.5198176) q[2];
sx q[2];
rz(-1.4547252) q[2];
sx q[2];
rz(-0.34250868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9782941) q[1];
sx q[1];
rz(-0.98493176) q[1];
sx q[1];
rz(-0.75978029) q[1];
rz(-2.6550754) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(-0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(-1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(2.4954605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.811759) q[0];
sx q[0];
rz(-0.23010294) q[0];
sx q[0];
rz(-2.0287201) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6693194) q[2];
sx q[2];
rz(-1.4204645) q[2];
sx q[2];
rz(2.3538102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7568126) q[1];
sx q[1];
rz(-1.2287178) q[1];
sx q[1];
rz(0.53201075) q[1];
rz(-pi) q[2];
rz(2.2046702) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(-0.015451775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-0.45483744) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(-0.79750693) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-3.033175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78128821) q[0];
sx q[0];
rz(-0.66626781) q[0];
sx q[0];
rz(-1.4772619) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3921146) q[2];
sx q[2];
rz(-1.6166501) q[2];
sx q[2];
rz(-0.97464857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3474116) q[1];
sx q[1];
rz(-1.3550183) q[1];
sx q[1];
rz(-0.17098917) q[1];
x q[2];
rz(0.18687825) q[3];
sx q[3];
rz(-0.87516057) q[3];
sx q[3];
rz(2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0124399) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(2.8016395) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(0.72559124) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(2.7774096) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(0.055158786) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63294166) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(2.7642194) q[0];
x q[1];
rz(-0.55969413) q[2];
sx q[2];
rz(-1.392138) q[2];
sx q[2];
rz(2.6924804) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89796472) q[1];
sx q[1];
rz(-2.1568255) q[1];
sx q[1];
rz(0.90378739) q[1];
rz(-2.0268029) q[3];
sx q[3];
rz(-1.6037233) q[3];
sx q[3];
rz(-1.959521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.4162083) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(-2.9329119) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(-1.4738884) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(0.24003868) q[3];
sx q[3];
rz(-1.5272899) q[3];
sx q[3];
rz(1.8361113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];