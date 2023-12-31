OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8653712) q[0];
sx q[0];
rz(-2.2844391) q[0];
sx q[0];
rz(3.0091118) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(5.6981882) q[1];
sx q[1];
rz(11.873801) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8938457) q[0];
sx q[0];
rz(-0.5286628) q[0];
sx q[0];
rz(-1.371944) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64451005) q[2];
sx q[2];
rz(-1.1877726) q[2];
sx q[2];
rz(-0.884998) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3092279) q[1];
sx q[1];
rz(-1.698306) q[1];
sx q[1];
rz(3.068919) q[1];
x q[2];
rz(-2.8928738) q[3];
sx q[3];
rz(-1.1097483) q[3];
sx q[3];
rz(2.98364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.5365323) q[2];
rz(1.6202554) q[3];
sx q[3];
rz(-1.6531569) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3261616) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.2492299) q[0];
rz(-2.5800887) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(-2.5610279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8037572) q[0];
sx q[0];
rz(-0.32807402) q[0];
sx q[0];
rz(0.34123811) q[0];
rz(-1.3524019) q[2];
sx q[2];
rz(-0.94422715) q[2];
sx q[2];
rz(-0.44109694) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3483352) q[1];
sx q[1];
rz(-1.0151334) q[1];
sx q[1];
rz(0.72850119) q[1];
rz(-pi) q[2];
rz(-2.7892116) q[3];
sx q[3];
rz(-2.0003194) q[3];
sx q[3];
rz(0.22892117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(2.7495524) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(2.5382606) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6648401) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(-3.1402918) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(3.1087648) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884018) q[0];
sx q[0];
rz(-1.8817888) q[0];
sx q[0];
rz(-1.8131959) q[0];
rz(-pi) q[1];
rz(0.29067729) q[2];
sx q[2];
rz(-1.0376087) q[2];
sx q[2];
rz(-2.3454587) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9877316) q[1];
sx q[1];
rz(-0.78480936) q[1];
sx q[1];
rz(-1.3373109) q[1];
x q[2];
rz(-1.0266617) q[3];
sx q[3];
rz(-1.2715724) q[3];
sx q[3];
rz(-2.2194089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34439987) q[2];
sx q[2];
rz(-1.025528) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(-2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(-2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(0.26279703) q[0];
rz(-0.2335877) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-2.3707726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642826) q[0];
sx q[0];
rz(-1.5543803) q[0];
sx q[0];
rz(1.55127) q[0];
rz(-pi) q[1];
rz(-3.0371573) q[2];
sx q[2];
rz(-2.0586788) q[2];
sx q[2];
rz(1.7761531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82476393) q[1];
sx q[1];
rz(-1.0526122) q[1];
sx q[1];
rz(1.1512685) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1905258) q[3];
sx q[3];
rz(-1.7912205) q[3];
sx q[3];
rz(0.64173736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(-2.4285994) q[2];
rz(-2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7825496) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(-1.4404526) q[0];
rz(3.0474995) q[1];
sx q[1];
rz(-2.4021939) q[1];
sx q[1];
rz(-2.9715911) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1124135) q[0];
sx q[0];
rz(-2.3110483) q[0];
sx q[0];
rz(-1.2519757) q[0];
rz(2.0001569) q[2];
sx q[2];
rz(-0.51157198) q[2];
sx q[2];
rz(-0.19829743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69668329) q[1];
sx q[1];
rz(-1.0001567) q[1];
sx q[1];
rz(-2.6615104) q[1];
x q[2];
rz(-0.87336297) q[3];
sx q[3];
rz(-2.1057099) q[3];
sx q[3];
rz(-1.3267335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(1.4040995) q[2];
rz(1.5935625) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(-1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(-0.45853841) q[0];
rz(0.25587747) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(-2.4564254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4165669) q[0];
sx q[0];
rz(-1.5013114) q[0];
sx q[0];
rz(2.2216703) q[0];
rz(-pi) q[1];
rz(2.9951727) q[2];
sx q[2];
rz(-0.78260566) q[2];
sx q[2];
rz(-0.27461068) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7428776) q[1];
sx q[1];
rz(-2.0600973) q[1];
sx q[1];
rz(-2.2687885) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5983726) q[3];
sx q[3];
rz(-0.58610361) q[3];
sx q[3];
rz(-2.5930149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(1.4292599) q[2];
rz(1.0990934) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.012506164) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(-2.4122453) q[0];
rz(-0.29306456) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(1.1475295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0070533) q[0];
sx q[0];
rz(-0.42047406) q[0];
sx q[0];
rz(0.18376952) q[0];
rz(-pi) q[1];
rz(-1.0853026) q[2];
sx q[2];
rz(-2.448423) q[2];
sx q[2];
rz(1.8144516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16312576) q[1];
sx q[1];
rz(-1.5297869) q[1];
sx q[1];
rz(0.48467111) q[1];
rz(-pi) q[2];
rz(-2.494874) q[3];
sx q[3];
rz(-0.59520703) q[3];
sx q[3];
rz(1.7531542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3532233) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(-0.74907556) q[2];
rz(2.4979112) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(-2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(2.3102982) q[0];
rz(-1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-2.7430699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9996223) q[0];
sx q[0];
rz(-0.49641434) q[0];
sx q[0];
rz(-1.9103785) q[0];
rz(1.4648828) q[2];
sx q[2];
rz(-0.97230655) q[2];
sx q[2];
rz(-2.523409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1917845) q[1];
sx q[1];
rz(-2.5439918) q[1];
sx q[1];
rz(-1.530184) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6284292) q[3];
sx q[3];
rz(-1.0740105) q[3];
sx q[3];
rz(-2.7548807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1901671) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(-1.8438967) q[2];
rz(1.1249582) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(2.4979533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-3.1289537) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(1.7522316) q[0];
rz(1.5147491) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(-2.0432037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59745715) q[0];
sx q[0];
rz(-1.8034593) q[0];
sx q[0];
rz(-2.068589) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3214925) q[2];
sx q[2];
rz(-1.2250049) q[2];
sx q[2];
rz(1.0320013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64975148) q[1];
sx q[1];
rz(-1.5713308) q[1];
sx q[1];
rz(-2.5623296) q[1];
x q[2];
rz(-1.3994201) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(2.6220654) q[2];
rz(1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(0.46646068) q[0];
rz(-2.9699504) q[1];
sx q[1];
rz(-1.9263093) q[1];
sx q[1];
rz(-2.5126273) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91613149) q[0];
sx q[0];
rz(-0.3779419) q[0];
sx q[0];
rz(1.7037237) q[0];
rz(-pi) q[1];
rz(1.8877108) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(0.72310477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65044636) q[1];
sx q[1];
rz(-1.5728587) q[1];
sx q[1];
rz(-2.6994929) q[1];
rz(-pi) q[2];
rz(-1.2260776) q[3];
sx q[3];
rz(-2.260672) q[3];
sx q[3];
rz(-0.41050875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(-1.0591327) q[2];
rz(-2.4641666) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(-2.2476851) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8582936) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(2.4181096) q[2];
sx q[2];
rz(-1.6037446) q[2];
sx q[2];
rz(-1.4708191) q[2];
rz(1.6737291) q[3];
sx q[3];
rz(-2.5522305) q[3];
sx q[3];
rz(1.1141368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
