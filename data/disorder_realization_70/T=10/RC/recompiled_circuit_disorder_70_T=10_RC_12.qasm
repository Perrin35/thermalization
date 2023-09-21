OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2274322) q[0];
sx q[0];
rz(-0.81613805) q[0];
sx q[0];
rz(0.039365191) q[0];
rz(-pi) q[1];
rz(2.2251031) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(0.33545845) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12373911) q[1];
sx q[1];
rz(-2.1542319) q[1];
sx q[1];
rz(-0.53637335) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7414054) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(-0.17822972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.263164) q[2];
rz(-1.8566711) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(0.22110573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9843922) q[0];
sx q[0];
rz(-1.7565787) q[0];
sx q[0];
rz(1.1093344) q[0];
rz(-pi) q[1];
rz(-1.2731304) q[2];
sx q[2];
rz(-0.3634828) q[2];
sx q[2];
rz(-0.01089451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4202538) q[1];
sx q[1];
rz(-1.4495279) q[1];
sx q[1];
rz(-2.6810357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17640555) q[3];
sx q[3];
rz(-2.6325912) q[3];
sx q[3];
rz(-0.41364663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(2.5533) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(-1.0480405) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(2.1411238) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(0.75769889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43514566) q[0];
sx q[0];
rz(-2.5120814) q[0];
sx q[0];
rz(-2.5802617) q[0];
rz(-pi) q[1];
rz(1.7311086) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(-2.2018873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0211027) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(-0.91336577) q[1];
rz(-pi) q[2];
rz(0.14962872) q[3];
sx q[3];
rz(-0.64390874) q[3];
sx q[3];
rz(0.097188918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(-0.83646742) q[2];
rz(-1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(-0.20382717) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214685) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.3285332) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851345) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(-1.7865208) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6287597) q[2];
sx q[2];
rz(-0.83216681) q[2];
sx q[2];
rz(1.8682478) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7718539) q[1];
sx q[1];
rz(-0.26680294) q[1];
sx q[1];
rz(-1.5688483) q[1];
rz(-pi) q[2];
x q[2];
rz(2.02416) q[3];
sx q[3];
rz(-2.8957267) q[3];
sx q[3];
rz(-2.9256431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(0.47248653) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(2.8033946) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11855928) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(-0.27079196) q[0];
rz(-pi) q[1];
rz(2.675823) q[2];
sx q[2];
rz(-1.2301187) q[2];
sx q[2];
rz(-0.53616947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5711172) q[1];
sx q[1];
rz(-1.8091396) q[1];
sx q[1];
rz(1.8829324) q[1];
rz(-pi) q[2];
rz(-0.91315956) q[3];
sx q[3];
rz(-2.1214161) q[3];
sx q[3];
rz(0.0084358128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(-1.6882287) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(-1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-1.0908303) q[0];
rz(-2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(0.18879034) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81297368) q[0];
sx q[0];
rz(-2.6972174) q[0];
sx q[0];
rz(0.61291738) q[0];
rz(1.9300869) q[2];
sx q[2];
rz(-2.852716) q[2];
sx q[2];
rz(-1.8747683) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69953883) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(-0.26775189) q[1];
x q[2];
rz(0.64340274) q[3];
sx q[3];
rz(-0.2989558) q[3];
sx q[3];
rz(-1.9770196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.171689) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(-1.3151273) q[2];
rz(-1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090102) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.3791929) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8934879) q[0];
sx q[0];
rz(-1.6525434) q[0];
sx q[0];
rz(-1.3260613) q[0];
rz(-pi) q[1];
rz(-0.71521476) q[2];
sx q[2];
rz(-2.5589802) q[2];
sx q[2];
rz(-1.6974534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52723253) q[1];
sx q[1];
rz(-2.7990606) q[1];
sx q[1];
rz(-2.637898) q[1];
x q[2];
rz(0.93382436) q[3];
sx q[3];
rz(-2.056042) q[3];
sx q[3];
rz(-2.3683734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(0.76422894) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.2932628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18023597) q[0];
sx q[0];
rz(-2.0016252) q[0];
sx q[0];
rz(0.21813099) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13883491) q[2];
sx q[2];
rz(-0.60563696) q[2];
sx q[2];
rz(-2.5790737) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1845491) q[1];
sx q[1];
rz(-0.42889412) q[1];
sx q[1];
rz(1.2317608) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.771365) q[3];
sx q[3];
rz(-1.1784394) q[3];
sx q[3];
rz(1.4810824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3747037) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(2.5602706) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.8301615) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8563961) q[0];
sx q[0];
rz(-0.49001339) q[0];
sx q[0];
rz(-1.6243837) q[0];
rz(1.7436696) q[2];
sx q[2];
rz(-1.0908974) q[2];
sx q[2];
rz(-2.929504) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3661256) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(-2.8833564) q[1];
x q[2];
rz(-1.5136396) q[3];
sx q[3];
rz(-0.74712979) q[3];
sx q[3];
rz(0.77880083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.320497) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(-2.9750032) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2607516) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-2.0064158) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(-2.7744055) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8242278) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(0.048006417) q[0];
rz(-pi) q[1];
rz(0.62717168) q[2];
sx q[2];
rz(-1.3417202) q[2];
sx q[2];
rz(2.0723745) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6348833) q[1];
sx q[1];
rz(-2.2193529) q[1];
sx q[1];
rz(1.8222005) q[1];
x q[2];
rz(-2.6606584) q[3];
sx q[3];
rz(-0.98610611) q[3];
sx q[3];
rz(-1.1993053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(1.665303) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-3.0565312) q[2];
sx q[2];
rz(-0.75123514) q[2];
sx q[2];
rz(-3.0820465) q[2];
rz(2.9885837) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];