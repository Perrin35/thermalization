OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6025699) q[0];
sx q[0];
rz(-0.56350001) q[0];
sx q[0];
rz(-2.6846057) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(3.8222651) q[1];
sx q[1];
rz(8.1487976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63859343) q[0];
sx q[0];
rz(-2.8832957) q[0];
sx q[0];
rz(0.55708416) q[0];
rz(-pi) q[1];
x q[1];
rz(0.692042) q[2];
sx q[2];
rz(-3.0243052) q[2];
sx q[2];
rz(0.28809822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16986019) q[1];
sx q[1];
rz(-1.6310551) q[1];
sx q[1];
rz(1.7737284) q[1];
x q[2];
rz(2.7492417) q[3];
sx q[3];
rz(-2.1808743) q[3];
sx q[3];
rz(3.0685134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(3.0965565) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(2.3388458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.734056) q[0];
sx q[0];
rz(-1.3450087) q[0];
sx q[0];
rz(1.9682103) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10404189) q[2];
sx q[2];
rz(-2.0718956) q[2];
sx q[2];
rz(0.83545557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3853559) q[1];
sx q[1];
rz(-1.6223575) q[1];
sx q[1];
rz(2.1897584) q[1];
x q[2];
rz(2.2400168) q[3];
sx q[3];
rz(-1.4687317) q[3];
sx q[3];
rz(-1.4249742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37614432) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(0.60570335) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95214343) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(2.8564575) q[0];
rz(2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.8018988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1012672) q[0];
sx q[0];
rz(-2.2650532) q[0];
sx q[0];
rz(-2.5800152) q[0];
rz(-pi) q[1];
rz(2.8324098) q[2];
sx q[2];
rz(-2.2367466) q[2];
sx q[2];
rz(-0.14112976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2337607) q[1];
sx q[1];
rz(-1.7572174) q[1];
sx q[1];
rz(2.9790661) q[1];
rz(-pi) q[2];
rz(-0.47311802) q[3];
sx q[3];
rz(-2.3948673) q[3];
sx q[3];
rz(2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(1.408067) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(-3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(2.9273422) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(-2.4750211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4144856) q[0];
sx q[0];
rz(-1.3654728) q[0];
sx q[0];
rz(-1.7617102) q[0];
x q[1];
rz(-0.33118403) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(2.3664898) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9208593) q[1];
sx q[1];
rz(-2.1174866) q[1];
sx q[1];
rz(2.2842555) q[1];
rz(-1.4227082) q[3];
sx q[3];
rz(-1.7184966) q[3];
sx q[3];
rz(-1.5642017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.6289571) q[0];
sx q[0];
rz(-1.2197138) q[0];
sx q[0];
rz(-2.3045325) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(-1.1594835) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449074) q[0];
sx q[0];
rz(-0.7888182) q[0];
sx q[0];
rz(-0.44992723) q[0];
x q[1];
rz(-0.15578606) q[2];
sx q[2];
rz(-2.163379) q[2];
sx q[2];
rz(0.32183811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24134484) q[1];
sx q[1];
rz(-0.7175788) q[1];
sx q[1];
rz(-0.89300968) q[1];
x q[2];
rz(1.0347576) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(1.1353726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7065457) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(-2.0909677) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(-3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(2.0264453) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(3.040722) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80412241) q[0];
sx q[0];
rz(-2.1967948) q[0];
sx q[0];
rz(1.7581598) q[0];
rz(-pi) q[1];
rz(-0.48730536) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(-3.0234408) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3979891) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(-0.40525873) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(-2.1586771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826542) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(-1.5902279) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(-2.4538453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6763517) q[0];
sx q[0];
rz(-0.40333336) q[0];
sx q[0];
rz(-2.1562735) q[0];
rz(-pi) q[1];
rz(2.9155832) q[2];
sx q[2];
rz(-1.8812211) q[2];
sx q[2];
rz(-1.4357391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2076599) q[1];
sx q[1];
rz(-1.4414409) q[1];
sx q[1];
rz(-1.5473066) q[1];
rz(-1.5457821) q[3];
sx q[3];
rz(-0.56846148) q[3];
sx q[3];
rz(2.3908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(0.079285346) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(1.1987196) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51533651) q[0];
sx q[0];
rz(-1.5021828) q[0];
sx q[0];
rz(-0.4347883) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.347581) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(1.5471293) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5115576) q[1];
sx q[1];
rz(-0.50543284) q[1];
sx q[1];
rz(0.10314421) q[1];
x q[2];
rz(-0.46320398) q[3];
sx q[3];
rz(-1.8052881) q[3];
sx q[3];
rz(2.6686058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(-0.78906995) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.566074) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(0.9334329) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8844922) q[0];
sx q[0];
rz(-1.677209) q[0];
sx q[0];
rz(1.5110925) q[0];
rz(-pi) q[1];
rz(-1.1912212) q[2];
sx q[2];
rz(-1.2588725) q[2];
sx q[2];
rz(2.1385857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.55798462) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(3.0834553) q[1];
rz(0.45456072) q[3];
sx q[3];
rz(-2.436736) q[3];
sx q[3];
rz(-2.9398033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(-0.36994568) q[2];
rz(0.89921078) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(-0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545492) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(-0.18383372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606124) q[0];
sx q[0];
rz(-1.4502589) q[0];
sx q[0];
rz(1.4535849) q[0];
rz(0.76639558) q[2];
sx q[2];
rz(-2.1470214) q[2];
sx q[2];
rz(3.0255084) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2059584) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(-2.2531855) q[1];
rz(3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(-1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(1.1584875) q[2];
sx q[2];
rz(-1.145913) q[2];
sx q[2];
rz(-0.10213012) q[2];
rz(-0.98060676) q[3];
sx q[3];
rz(-1.4242314) q[3];
sx q[3];
rz(-0.84606708) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];