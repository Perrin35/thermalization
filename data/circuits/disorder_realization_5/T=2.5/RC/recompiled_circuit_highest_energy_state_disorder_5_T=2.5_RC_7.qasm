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
rz(-0.37762168) q[0];
sx q[0];
rz(-2.7130337) q[0];
sx q[0];
rz(-0.23634401) q[0];
rz(-1.6726681) q[1];
sx q[1];
rz(-1.5863215) q[1];
sx q[1];
rz(-0.16170391) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16888389) q[0];
sx q[0];
rz(-2.5492269) q[0];
sx q[0];
rz(1.374872) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8558883) q[2];
sx q[2];
rz(-0.20390715) q[2];
sx q[2];
rz(-3.0406024) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.573805) q[1];
sx q[1];
rz(-1.5777138) q[1];
sx q[1];
rz(1.597639) q[1];
rz(-pi) q[2];
rz(2.6246965) q[3];
sx q[3];
rz(-0.42498838) q[3];
sx q[3];
rz(-1.6551415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8990495) q[2];
sx q[2];
rz(-2.264302) q[2];
sx q[2];
rz(0.38929942) q[2];
rz(-2.9111275) q[3];
sx q[3];
rz(-0.018229818) q[3];
sx q[3];
rz(-0.61483312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56735754) q[0];
sx q[0];
rz(-2.1901972) q[0];
sx q[0];
rz(-1.4943328) q[0];
rz(-1.5859454) q[1];
sx q[1];
rz(-2.9289398) q[1];
sx q[1];
rz(2.0071323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26274219) q[0];
sx q[0];
rz(-1.3171742) q[0];
sx q[0];
rz(0.72611673) q[0];
rz(-pi) q[1];
rz(2.0087162) q[2];
sx q[2];
rz(-2.2575402) q[2];
sx q[2];
rz(-1.9599635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4006869) q[1];
sx q[1];
rz(-2.4016651) q[1];
sx q[1];
rz(-0.94018971) q[1];
rz(-pi) q[2];
rz(-2.6051571) q[3];
sx q[3];
rz(-0.89563771) q[3];
sx q[3];
rz(-1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.050934164) q[2];
sx q[2];
rz(-0.91973534) q[2];
sx q[2];
rz(-1.301379) q[2];
rz(2.0672412) q[3];
sx q[3];
rz(-0.31000546) q[3];
sx q[3];
rz(-1.5460792) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.018547) q[0];
sx q[0];
rz(-0.30236852) q[0];
sx q[0];
rz(-2.5240335) q[0];
rz(1.1005719) q[1];
sx q[1];
rz(-0.01958422) q[1];
sx q[1];
rz(2.7380131) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7779008) q[0];
sx q[0];
rz(-0.11523031) q[0];
sx q[0];
rz(-1.6621424) q[0];
rz(-pi) q[1];
rz(-0.52578747) q[2];
sx q[2];
rz(-1.4830768) q[2];
sx q[2];
rz(1.45873) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8858635) q[1];
sx q[1];
rz(-1.7049864) q[1];
sx q[1];
rz(1.7291452) q[1];
rz(2.7652214) q[3];
sx q[3];
rz(-0.59089243) q[3];
sx q[3];
rz(-1.2749689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6214211) q[2];
sx q[2];
rz(-1.4474063) q[2];
sx q[2];
rz(-0.78110313) q[2];
rz(2.805294) q[3];
sx q[3];
rz(-1.4399485) q[3];
sx q[3];
rz(2.4380016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.6853365) q[0];
sx q[0];
rz(-2.6254613) q[0];
sx q[0];
rz(-1.2180895) q[0];
rz(1.7958027) q[1];
sx q[1];
rz(-3.1333874) q[1];
sx q[1];
rz(1.2340612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6910373) q[0];
sx q[0];
rz(-2.9515186) q[0];
sx q[0];
rz(0.051817373) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61361112) q[2];
sx q[2];
rz(-1.8119115) q[2];
sx q[2];
rz(2.5054153) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8414536) q[1];
sx q[1];
rz(-0.44504224) q[1];
sx q[1];
rz(0.95325327) q[1];
rz(-pi) q[2];
rz(-0.59438057) q[3];
sx q[3];
rz(-0.012033741) q[3];
sx q[3];
rz(1.8569225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1951083) q[2];
sx q[2];
rz(-2.7184964) q[2];
sx q[2];
rz(-1.9931603) q[2];
rz(-0.75442433) q[3];
sx q[3];
rz(-1.9781338) q[3];
sx q[3];
rz(0.26328009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7719583) q[0];
sx q[0];
rz(-0.088005528) q[0];
sx q[0];
rz(-2.5391286) q[0];
rz(0.98214904) q[1];
sx q[1];
rz(-0.0063449675) q[1];
sx q[1];
rz(2.6314661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23778039) q[0];
sx q[0];
rz(-2.8717715) q[0];
sx q[0];
rz(-2.3348007) q[0];
x q[1];
rz(1.6922166) q[2];
sx q[2];
rz(-1.6573141) q[2];
sx q[2];
rz(1.0457116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1593455) q[1];
sx q[1];
rz(-0.67894113) q[1];
sx q[1];
rz(2.5117876) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.00058393936) q[3];
sx q[3];
rz(-2.3272132) q[3];
sx q[3];
rz(-1.3834805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.043857418) q[2];
sx q[2];
rz(-1.1172224) q[2];
sx q[2];
rz(2.1065693) q[2];
rz(-1.6739316) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(0.58290946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8408836) q[0];
sx q[0];
rz(-2.1568334) q[0];
sx q[0];
rz(1.0915407) q[0];
rz(-0.40065271) q[1];
sx q[1];
rz(-0.0023829208) q[1];
sx q[1];
rz(0.56652743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0877675) q[0];
sx q[0];
rz(-1.7051518) q[0];
sx q[0];
rz(2.3162486) q[0];
rz(-pi) q[1];
rz(-2.3344757) q[2];
sx q[2];
rz(-0.87701488) q[2];
sx q[2];
rz(-1.4556231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.57921806) q[1];
sx q[1];
rz(-1.7953582) q[1];
sx q[1];
rz(0.57266219) q[1];
x q[2];
rz(-0.57265307) q[3];
sx q[3];
rz(-1.6653456) q[3];
sx q[3];
rz(1.270164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94622508) q[2];
sx q[2];
rz(-2.3905498) q[2];
sx q[2];
rz(-1.5283778) q[2];
rz(2.1402806) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(-2.9010469) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840851) q[0];
sx q[0];
rz(-0.79428285) q[0];
sx q[0];
rz(1.9965782) q[0];
rz(1.6192294) q[1];
sx q[1];
rz(-0.018773627) q[1];
sx q[1];
rz(1.9782664) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881108) q[0];
sx q[0];
rz(-1.6170752) q[0];
sx q[0];
rz(1.2812216) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9923963) q[2];
sx q[2];
rz(-2.2393804) q[2];
sx q[2];
rz(-0.72982349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8790354) q[1];
sx q[1];
rz(-1.702629) q[1];
sx q[1];
rz(2.1798032) q[1];
rz(1.6895503) q[3];
sx q[3];
rz(-2.2014479) q[3];
sx q[3];
rz(1.5753559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6254977) q[2];
sx q[2];
rz(-2.1558546) q[2];
sx q[2];
rz(2.7682448) q[2];
rz(0.18800023) q[3];
sx q[3];
rz(-1.5865654) q[3];
sx q[3];
rz(-1.1628304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9173376) q[0];
sx q[0];
rz(-0.52977109) q[0];
sx q[0];
rz(-0.24714558) q[0];
rz(-0.22357926) q[1];
sx q[1];
rz(-3.1378742) q[1];
sx q[1];
rz(-1.6252958) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2151536) q[0];
sx q[0];
rz(-1.5064539) q[0];
sx q[0];
rz(-1.45432) q[0];
x q[1];
rz(0.76659492) q[2];
sx q[2];
rz(-1.8803036) q[2];
sx q[2];
rz(2.7628037) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1268218) q[1];
sx q[1];
rz(-0.39431652) q[1];
sx q[1];
rz(-3.1000351) q[1];
x q[2];
rz(-1.6187632) q[3];
sx q[3];
rz(-2.164451) q[3];
sx q[3];
rz(-1.1831212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94459263) q[2];
sx q[2];
rz(-0.71114117) q[2];
sx q[2];
rz(1.5947343) q[2];
rz(2.366015) q[3];
sx q[3];
rz(-1.8577134) q[3];
sx q[3];
rz(1.4625134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
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
rz(0.32770661) q[0];
sx q[0];
rz(-2.1729108) q[0];
sx q[0];
rz(-2.0342597) q[0];
rz(-0.92512098) q[1];
sx q[1];
rz(-3.1396301) q[1];
sx q[1];
rz(0.75606871) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32996477) q[0];
sx q[0];
rz(-2.1058308) q[0];
sx q[0];
rz(-2.7828548) q[0];
rz(0.0063757141) q[2];
sx q[2];
rz(-1.0399605) q[2];
sx q[2];
rz(2.6311925) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.392627) q[1];
sx q[1];
rz(-0.64425978) q[1];
sx q[1];
rz(2.6436252) q[1];
rz(2.954079) q[3];
sx q[3];
rz(-2.5616259) q[3];
sx q[3];
rz(-0.43646508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6206616) q[2];
sx q[2];
rz(-0.88520092) q[2];
sx q[2];
rz(2.1405641) q[2];
rz(1.3641317) q[3];
sx q[3];
rz(-0.93327156) q[3];
sx q[3];
rz(2.6076243) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022024632) q[0];
sx q[0];
rz(-1.7689995) q[0];
sx q[0];
rz(-2.6764828) q[0];
rz(1.8067092) q[1];
sx q[1];
rz(-2.7692134) q[1];
sx q[1];
rz(-1.575527) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4441285) q[0];
sx q[0];
rz(-2.9078472) q[0];
sx q[0];
rz(-1.759619) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.100824) q[2];
sx q[2];
rz(-1.3000856) q[2];
sx q[2];
rz(1.9418429) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5425092) q[1];
sx q[1];
rz(-1.5712106) q[1];
sx q[1];
rz(-3.1387657) q[1];
rz(-3.0483143) q[3];
sx q[3];
rz(-0.3808379) q[3];
sx q[3];
rz(-2.8990428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89432013) q[2];
sx q[2];
rz(-3.0970116) q[2];
sx q[2];
rz(1.0856005) q[2];
rz(-1.2224489) q[3];
sx q[3];
rz(-0.51210755) q[3];
sx q[3];
rz(-1.2781757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4492252) q[0];
sx q[0];
rz(-1.4341555) q[0];
sx q[0];
rz(1.7644802) q[0];
rz(-1.5583246) q[1];
sx q[1];
rz(-2.227034) q[1];
sx q[1];
rz(-2.9169678) q[1];
rz(3.0408411) q[2];
sx q[2];
rz(-1.5755079) q[2];
sx q[2];
rz(1.8625349) q[2];
rz(0.23595594) q[3];
sx q[3];
rz(-1.1370084) q[3];
sx q[3];
rz(-1.3113547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
