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
rz(-0.4975118) q[0];
sx q[0];
rz(4.4805718) q[0];
sx q[0];
rz(6.5988402) q[0];
rz(-1.9714126) q[1];
sx q[1];
rz(-2.7538731) q[1];
sx q[1];
rz(2.5375836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3919298) q[0];
sx q[0];
rz(-1.6214341) q[0];
sx q[0];
rz(1.873436) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96678218) q[2];
sx q[2];
rz(-0.76672115) q[2];
sx q[2];
rz(-1.7974092) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0993201) q[1];
sx q[1];
rz(-1.0612773) q[1];
sx q[1];
rz(-2.7541942) q[1];
rz(-1.2420869) q[3];
sx q[3];
rz(-0.32530537) q[3];
sx q[3];
rz(-1.6360375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9903367) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(-1.3585496) q[2];
rz(1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(-1.6319298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58594054) q[0];
sx q[0];
rz(-0.63783115) q[0];
sx q[0];
rz(-1.1974539) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(-1.7053568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.013151) q[0];
sx q[0];
rz(-1.801898) q[0];
sx q[0];
rz(0.9856772) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4714437) q[2];
sx q[2];
rz(-0.87026419) q[2];
sx q[2];
rz(-1.6068899) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2982218) q[1];
sx q[1];
rz(-0.22208115) q[1];
sx q[1];
rz(0.48413094) q[1];
x q[2];
rz(-0.45577502) q[3];
sx q[3];
rz(-1.3823798) q[3];
sx q[3];
rz(-2.6236603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88840914) q[2];
sx q[2];
rz(-1.2359572) q[2];
sx q[2];
rz(-0.02296981) q[2];
rz(-2.2394771) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(-2.7149916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4334634) q[0];
sx q[0];
rz(-1.7925649) q[0];
sx q[0];
rz(0.9915114) q[0];
rz(-1.0333215) q[1];
sx q[1];
rz(-2.8585377) q[1];
sx q[1];
rz(1.906377) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69312364) q[0];
sx q[0];
rz(-0.64769324) q[0];
sx q[0];
rz(-0.41843398) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5318145) q[2];
sx q[2];
rz(-1.1705361) q[2];
sx q[2];
rz(-2.5230809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.24656235) q[1];
sx q[1];
rz(-0.8324648) q[1];
sx q[1];
rz(2.8042996) q[1];
x q[2];
rz(-0.83241141) q[3];
sx q[3];
rz(-0.73736546) q[3];
sx q[3];
rz(2.5016145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11423763) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(2.1503964) q[2];
rz(-1.8441955) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(1.4170925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7436413) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(-2.3241296) q[0];
rz(-0.3262597) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(-2.356333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0164675) q[0];
sx q[0];
rz(-0.40832357) q[0];
sx q[0];
rz(2.295408) q[0];
rz(-pi) q[1];
rz(-1.7677069) q[2];
sx q[2];
rz(-2.2047538) q[2];
sx q[2];
rz(1.9206573) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7689159) q[1];
sx q[1];
rz(-1.1296037) q[1];
sx q[1];
rz(0.4695453) q[1];
x q[2];
rz(1.4912075) q[3];
sx q[3];
rz(-1.2168988) q[3];
sx q[3];
rz(2.6864664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0205445) q[2];
sx q[2];
rz(-2.497017) q[2];
sx q[2];
rz(0.56309593) q[2];
rz(0.8935039) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23984443) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(-1.5455986) q[0];
rz(-2.1196938) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(2.9603069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5820532) q[0];
sx q[0];
rz(-1.3783558) q[0];
sx q[0];
rz(2.0700702) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6659545) q[2];
sx q[2];
rz(-1.5283094) q[2];
sx q[2];
rz(-1.795639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.681057) q[1];
sx q[1];
rz(-1.3486282) q[1];
sx q[1];
rz(3.1128037) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6986609) q[3];
sx q[3];
rz(-1.7128403) q[3];
sx q[3];
rz(1.0012817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.01866092) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(-1.2459416) q[2];
rz(-2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(-2.3004801) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(-2.8327508) q[0];
rz(0.32288512) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(-0.13993941) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2353781) q[0];
sx q[0];
rz(-1.8188634) q[0];
sx q[0];
rz(0.95309044) q[0];
x q[1];
rz(0.094801589) q[2];
sx q[2];
rz(-2.0249512) q[2];
sx q[2];
rz(2.4081782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7227712) q[1];
sx q[1];
rz(-2.3496637) q[1];
sx q[1];
rz(1.140711) q[1];
x q[2];
rz(2.2298584) q[3];
sx q[3];
rz(-1.3305802) q[3];
sx q[3];
rz(1.6068589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66190019) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(2.8996186) q[2];
rz(-2.3954929) q[3];
sx q[3];
rz(-1.3662162) q[3];
sx q[3];
rz(1.411875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.11874966) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(1.2710849) q[0];
rz(2.5785043) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(-1.7005327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1313167) q[0];
sx q[0];
rz(-0.97945853) q[0];
sx q[0];
rz(2.813126) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7092429) q[2];
sx q[2];
rz(-1.9488153) q[2];
sx q[2];
rz(-1.0919309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.064441) q[1];
sx q[1];
rz(-1.3485326) q[1];
sx q[1];
rz(0.42776107) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8981179) q[3];
sx q[3];
rz(-2.4015275) q[3];
sx q[3];
rz(-0.26672637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9219804) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(-1.316635) q[2];
rz(2.5804139) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(-1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0372666) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-2.2547146) q[0];
rz(0.2001702) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(-1.0221457) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3322395) q[0];
sx q[0];
rz(-3.0080171) q[0];
sx q[0];
rz(-1.7627526) q[0];
rz(2.3273507) q[2];
sx q[2];
rz(-2.1541284) q[2];
sx q[2];
rz(2.5098206) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8238942) q[1];
sx q[1];
rz(-2.3168457) q[1];
sx q[1];
rz(-3.1391678) q[1];
rz(-0.51428719) q[3];
sx q[3];
rz(-2.5011241) q[3];
sx q[3];
rz(-0.91292229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6518121) q[2];
sx q[2];
rz(-2.688789) q[2];
sx q[2];
rz(0.81857267) q[2];
rz(0.98322785) q[3];
sx q[3];
rz(-0.95663095) q[3];
sx q[3];
rz(-0.53808588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035456903) q[0];
sx q[0];
rz(-1.9945194) q[0];
sx q[0];
rz(-1.5552833) q[0];
rz(0.69681329) q[1];
sx q[1];
rz(-2.9882444) q[1];
sx q[1];
rz(-2.3568025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9826384) q[0];
sx q[0];
rz(-1.0305911) q[0];
sx q[0];
rz(-2.0870952) q[0];
x q[1];
rz(1.3680787) q[2];
sx q[2];
rz(-1.4215111) q[2];
sx q[2];
rz(-1.210807) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5604374) q[1];
sx q[1];
rz(-2.0056917) q[1];
sx q[1];
rz(1.6078609) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9905213) q[3];
sx q[3];
rz(-1.0015206) q[3];
sx q[3];
rz(1.9939533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4697504) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(2.3308241) q[2];
rz(-1.2189216) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(-1.0651275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(0.78275457) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(-1.4240356) q[0];
rz(-2.3639823) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(2.3861859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9616325) q[0];
sx q[0];
rz(-2.9316219) q[0];
sx q[0];
rz(-1.5487973) q[0];
rz(-pi) q[1];
rz(-2.7451186) q[2];
sx q[2];
rz(-1.5546397) q[2];
sx q[2];
rz(1.1262058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0499784) q[1];
sx q[1];
rz(-2.0777933) q[1];
sx q[1];
rz(1.9732287) q[1];
rz(-pi) q[2];
rz(-1.2437263) q[3];
sx q[3];
rz(-0.9216412) q[3];
sx q[3];
rz(-2.7130733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5767673) q[2];
sx q[2];
rz(-2.8456523) q[2];
sx q[2];
rz(2.9887065) q[2];
rz(-1.2711924) q[3];
sx q[3];
rz(-1.8891687) q[3];
sx q[3];
rz(0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85970238) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(-1.2234756) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(-2.0064932) q[2];
sx q[2];
rz(-2.1232599) q[2];
sx q[2];
rz(2.9386414) q[2];
rz(0.85930227) q[3];
sx q[3];
rz(-0.51325428) q[3];
sx q[3];
rz(-0.2556066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
