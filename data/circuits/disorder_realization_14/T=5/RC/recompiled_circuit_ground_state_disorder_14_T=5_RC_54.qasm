OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8811964) q[0];
sx q[0];
rz(-1.908778) q[0];
sx q[0];
rz(2.8137299) q[0];
rz(-0.39363632) q[1];
sx q[1];
rz(4.7596158) q[1];
sx q[1];
rz(14.16852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0868139) q[0];
sx q[0];
rz(-2.940455) q[0];
sx q[0];
rz(-3.1166409) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0225564) q[2];
sx q[2];
rz(-1.3371144) q[2];
sx q[2];
rz(-1.7811687) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5757619) q[1];
sx q[1];
rz(-1.1850433) q[1];
sx q[1];
rz(1.7931531) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1957715) q[3];
sx q[3];
rz(-1.2030691) q[3];
sx q[3];
rz(-1.036552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.1800805) q[2];
sx q[2];
rz(-2.1891258) q[2];
sx q[2];
rz(2.3251779) q[2];
rz(0.69711971) q[3];
sx q[3];
rz(-0.37728798) q[3];
sx q[3];
rz(0.94059801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78696752) q[0];
sx q[0];
rz(-1.0579728) q[0];
sx q[0];
rz(2.4355198) q[0];
rz(0.1708897) q[1];
sx q[1];
rz(-2.6306174) q[1];
sx q[1];
rz(2.1461646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4137693) q[0];
sx q[0];
rz(-1.2301155) q[0];
sx q[0];
rz(-1.2321939) q[0];
x q[1];
rz(2.6508743) q[2];
sx q[2];
rz(-0.92923895) q[2];
sx q[2];
rz(-3.0131755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4793206) q[1];
sx q[1];
rz(-2.8779128) q[1];
sx q[1];
rz(2.1946043) q[1];
rz(-pi) q[2];
rz(0.63748116) q[3];
sx q[3];
rz(-1.5439171) q[3];
sx q[3];
rz(2.0826552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.55676111) q[2];
sx q[2];
rz(-0.91412592) q[2];
sx q[2];
rz(-2.8625028) q[2];
rz(-2.4537405) q[3];
sx q[3];
rz(-0.22356859) q[3];
sx q[3];
rz(-2.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85182178) q[0];
sx q[0];
rz(-2.2230447) q[0];
sx q[0];
rz(2.7626792) q[0];
rz(-1.1907578) q[1];
sx q[1];
rz(-2.7749116) q[1];
sx q[1];
rz(2.3653638) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60769865) q[0];
sx q[0];
rz(-1.726681) q[0];
sx q[0];
rz(1.1381696) q[0];
x q[1];
rz(1.7227145) q[2];
sx q[2];
rz(-1.3891274) q[2];
sx q[2];
rz(0.40975299) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4632193) q[1];
sx q[1];
rz(-1.8253981) q[1];
sx q[1];
rz(2.4237952) q[1];
x q[2];
rz(0.046136304) q[3];
sx q[3];
rz(-1.0469808) q[3];
sx q[3];
rz(-2.5906117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1227405) q[2];
sx q[2];
rz(-1.9032225) q[2];
sx q[2];
rz(-0.59993258) q[2];
rz(-1.9020724) q[3];
sx q[3];
rz(-1.4116838) q[3];
sx q[3];
rz(-0.6111353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3032853) q[0];
sx q[0];
rz(-0.67876434) q[0];
sx q[0];
rz(-1.7093866) q[0];
rz(-2.2749061) q[1];
sx q[1];
rz(-1.7984093) q[1];
sx q[1];
rz(1.0923045) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537395) q[0];
sx q[0];
rz(-0.9747552) q[0];
sx q[0];
rz(0.88576646) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1020911) q[2];
sx q[2];
rz(-0.73248011) q[2];
sx q[2];
rz(-1.8341482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6505423) q[1];
sx q[1];
rz(-2.67727) q[1];
sx q[1];
rz(0.004267172) q[1];
rz(0.79444973) q[3];
sx q[3];
rz(-2.2243119) q[3];
sx q[3];
rz(0.78032035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0563353) q[2];
sx q[2];
rz(-1.1577572) q[2];
sx q[2];
rz(-1.866327) q[2];
rz(-1.7665675) q[3];
sx q[3];
rz(-1.2396038) q[3];
sx q[3];
rz(1.0331155) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.105724) q[0];
sx q[0];
rz(-2.2585456) q[0];
sx q[0];
rz(-1.3637654) q[0];
rz(-2.5719602) q[1];
sx q[1];
rz(-0.49915794) q[1];
sx q[1];
rz(0.5035351) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7020126) q[0];
sx q[0];
rz(-1.3786714) q[0];
sx q[0];
rz(-1.4441326) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3634144) q[2];
sx q[2];
rz(-2.2433837) q[2];
sx q[2];
rz(1.961381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9787748) q[1];
sx q[1];
rz(-1.6520588) q[1];
sx q[1];
rz(-0.012465076) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77025099) q[3];
sx q[3];
rz(-1.7903084) q[3];
sx q[3];
rz(1.3971922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86294952) q[2];
sx q[2];
rz(-2.7549665) q[2];
sx q[2];
rz(1.6031727) q[2];
rz(-3.0053511) q[3];
sx q[3];
rz(-1.7377661) q[3];
sx q[3];
rz(-1.2374102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513026) q[0];
sx q[0];
rz(-2.5232115) q[0];
sx q[0];
rz(-0.34143099) q[0];
rz(-2.4104207) q[1];
sx q[1];
rz(-2.5656504) q[1];
sx q[1];
rz(-2.0700571) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6302231) q[0];
sx q[0];
rz(-1.9804738) q[0];
sx q[0];
rz(-1.131534) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5093083) q[2];
sx q[2];
rz(-2.1849423) q[2];
sx q[2];
rz(1.5100069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9042957) q[1];
sx q[1];
rz(-1.132442) q[1];
sx q[1];
rz(-2.3890952) q[1];
rz(-1.5789746) q[3];
sx q[3];
rz(-2.6371752) q[3];
sx q[3];
rz(-0.8917419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7858481) q[2];
sx q[2];
rz(-1.637849) q[2];
sx q[2];
rz(2.6453633) q[2];
rz(-1.548454) q[3];
sx q[3];
rz(-1.691247) q[3];
sx q[3];
rz(2.1350071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9736495) q[0];
sx q[0];
rz(-1.4730467) q[0];
sx q[0];
rz(3.1267401) q[0];
rz(-1.4133833) q[1];
sx q[1];
rz(-1.1057066) q[1];
sx q[1];
rz(2.3766439) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0524372) q[0];
sx q[0];
rz(-0.90298072) q[0];
sx q[0];
rz(-0.74217702) q[0];
rz(-pi) q[1];
rz(0.42301793) q[2];
sx q[2];
rz(-2.4126965) q[2];
sx q[2];
rz(-2.4409082) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97100869) q[1];
sx q[1];
rz(-0.39675823) q[1];
sx q[1];
rz(-3.1301857) q[1];
rz(-2.025795) q[3];
sx q[3];
rz(-0.46543803) q[3];
sx q[3];
rz(2.1707141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0689653) q[2];
sx q[2];
rz(-1.5494538) q[2];
sx q[2];
rz(2.3907982) q[2];
rz(-0.99610656) q[3];
sx q[3];
rz(-0.80544296) q[3];
sx q[3];
rz(0.44710818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4873416) q[0];
sx q[0];
rz(-2.5387006) q[0];
sx q[0];
rz(2.531429) q[0];
rz(2.8918686) q[1];
sx q[1];
rz(-1.6619253) q[1];
sx q[1];
rz(-3.0009559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67540231) q[0];
sx q[0];
rz(-1.1151214) q[0];
sx q[0];
rz(-0.82929096) q[0];
x q[1];
rz(-1.4678367) q[2];
sx q[2];
rz(-1.067458) q[2];
sx q[2];
rz(1.9416941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7761155) q[1];
sx q[1];
rz(-1.8473313) q[1];
sx q[1];
rz(1.4955669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28293149) q[3];
sx q[3];
rz(-2.3608876) q[3];
sx q[3];
rz(1.9819575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3934624) q[2];
sx q[2];
rz(-2.738214) q[2];
sx q[2];
rz(2.051579) q[2];
rz(-1.2202834) q[3];
sx q[3];
rz(-1.0459432) q[3];
sx q[3];
rz(-2.2530341) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1423993) q[0];
sx q[0];
rz(-1.1141454) q[0];
sx q[0];
rz(2.2166369) q[0];
rz(-1.5026745) q[1];
sx q[1];
rz(-1.0097367) q[1];
sx q[1];
rz(0.22770539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8186491) q[0];
sx q[0];
rz(-2.5997351) q[0];
sx q[0];
rz(1.3635554) q[0];
rz(2.5617719) q[2];
sx q[2];
rz(-1.3353773) q[2];
sx q[2];
rz(-1.5832242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.40434228) q[1];
sx q[1];
rz(-1.3793762) q[1];
sx q[1];
rz(0.36083513) q[1];
x q[2];
rz(-1.3761767) q[3];
sx q[3];
rz(-0.96590079) q[3];
sx q[3];
rz(0.61974398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1091252) q[2];
sx q[2];
rz(-0.6296857) q[2];
sx q[2];
rz(-0.49051782) q[2];
rz(2.0971175) q[3];
sx q[3];
rz(-2.677768) q[3];
sx q[3];
rz(-3.131955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57546416) q[0];
sx q[0];
rz(-2.6818891) q[0];
sx q[0];
rz(1.8756961) q[0];
rz(1.5220386) q[1];
sx q[1];
rz(-1.252389) q[1];
sx q[1];
rz(2.6957846) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4981549) q[0];
sx q[0];
rz(-2.0863976) q[0];
sx q[0];
rz(2.7650973) q[0];
rz(-pi) q[1];
rz(-1.5238026) q[2];
sx q[2];
rz(-1.0596152) q[2];
sx q[2];
rz(-1.5906478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83611402) q[1];
sx q[1];
rz(-0.96313253) q[1];
sx q[1];
rz(1.4038248) q[1];
rz(0.021546797) q[3];
sx q[3];
rz(-1.2974707) q[3];
sx q[3];
rz(0.97238982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1294935) q[2];
sx q[2];
rz(-0.98126498) q[2];
sx q[2];
rz(-2.1141466) q[2];
rz(1.3198613) q[3];
sx q[3];
rz(-2.8362995) q[3];
sx q[3];
rz(2.6125438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32578596) q[0];
sx q[0];
rz(-1.7740213) q[0];
sx q[0];
rz(-1.0446145) q[0];
rz(2.2141937) q[1];
sx q[1];
rz(-1.8721885) q[1];
sx q[1];
rz(2.4354557) q[1];
rz(-0.76005475) q[2];
sx q[2];
rz(-2.7031923) q[2];
sx q[2];
rz(-2.2491982) q[2];
rz(-2.7488662) q[3];
sx q[3];
rz(-0.65734335) q[3];
sx q[3];
rz(-1.7398294) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
