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
rz(-3.0109316) q[0];
sx q[0];
rz(-1.8008512) q[0];
sx q[0];
rz(-1.2641579) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(-1.0882508) q[1];
sx q[1];
rz(0.01297125) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2269079) q[0];
sx q[0];
rz(-2.6400262) q[0];
sx q[0];
rz(-2.8109994) q[0];
rz(-0.16275587) q[2];
sx q[2];
rz(-1.6784894) q[2];
sx q[2];
rz(-2.8329605) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88718677) q[1];
sx q[1];
rz(-1.9209083) q[1];
sx q[1];
rz(1.2302047) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5586224) q[3];
sx q[3];
rz(-1.3912001) q[3];
sx q[3];
rz(-2.666134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8791447) q[2];
sx q[2];
rz(-1.4538572) q[2];
sx q[2];
rz(-3.0461779) q[2];
rz(2.6573507) q[3];
sx q[3];
rz(-0.80317322) q[3];
sx q[3];
rz(-1.4580844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5937623) q[0];
sx q[0];
rz(-1.4365124) q[0];
sx q[0];
rz(2.4875212) q[0];
rz(-1.3520799) q[1];
sx q[1];
rz(-2.4857931) q[1];
sx q[1];
rz(-0.10993122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93811518) q[0];
sx q[0];
rz(-0.02349668) q[0];
sx q[0];
rz(-2.1142939) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1651697) q[2];
sx q[2];
rz(-2.4196673) q[2];
sx q[2];
rz(0.41625574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1812352) q[1];
sx q[1];
rz(-1.3981888) q[1];
sx q[1];
rz(0.45876518) q[1];
rz(2.5853048) q[3];
sx q[3];
rz(-2.6835052) q[3];
sx q[3];
rz(-2.213495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8249417) q[2];
sx q[2];
rz(-1.9330838) q[2];
sx q[2];
rz(-1.6744772) q[2];
rz(2.1064827) q[3];
sx q[3];
rz(-0.78101522) q[3];
sx q[3];
rz(1.3181814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255945) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(2.3509534) q[0];
rz(0.74360338) q[1];
sx q[1];
rz(-1.598282) q[1];
sx q[1];
rz(2.5620983) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6181439) q[0];
sx q[0];
rz(-1.0001527) q[0];
sx q[0];
rz(-0.6499001) q[0];
rz(-pi) q[1];
rz(-2.8786009) q[2];
sx q[2];
rz(-1.2858258) q[2];
sx q[2];
rz(-2.6927352) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6404785) q[1];
sx q[1];
rz(-1.6081297) q[1];
sx q[1];
rz(1.3126612) q[1];
rz(-pi) q[2];
rz(0.73852818) q[3];
sx q[3];
rz(-0.66940847) q[3];
sx q[3];
rz(0.51046342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4641331) q[2];
sx q[2];
rz(-0.83659283) q[2];
sx q[2];
rz(-2.7602688) q[2];
rz(-0.85121202) q[3];
sx q[3];
rz(-0.92654735) q[3];
sx q[3];
rz(-2.4066431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320936) q[0];
sx q[0];
rz(-2.5272326) q[0];
sx q[0];
rz(0.36439782) q[0];
rz(-0.20026194) q[1];
sx q[1];
rz(-1.4118782) q[1];
sx q[1];
rz(0.099954896) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6969589) q[0];
sx q[0];
rz(-1.4958188) q[0];
sx q[0];
rz(3.1211389) q[0];
rz(-pi) q[1];
rz(-1.1420629) q[2];
sx q[2];
rz(-1.6661281) q[2];
sx q[2];
rz(-2.7470061) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2435088) q[1];
sx q[1];
rz(-1.5931411) q[1];
sx q[1];
rz(0.0925272) q[1];
rz(-0.95619802) q[3];
sx q[3];
rz(-2.0609988) q[3];
sx q[3];
rz(-2.4407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0692856) q[2];
sx q[2];
rz(-2.0733209) q[2];
sx q[2];
rz(-0.11492534) q[2];
rz(-1.140444) q[3];
sx q[3];
rz(-1.8585669) q[3];
sx q[3];
rz(0.97525245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2045778) q[0];
sx q[0];
rz(-0.95252043) q[0];
sx q[0];
rz(2.9691147) q[0];
rz(-2.1306439) q[1];
sx q[1];
rz(-2.5806249) q[1];
sx q[1];
rz(2.9867461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60733838) q[0];
sx q[0];
rz(-1.7887428) q[0];
sx q[0];
rz(1.8740669) q[0];
x q[1];
rz(-0.24598083) q[2];
sx q[2];
rz(-2.2005251) q[2];
sx q[2];
rz(0.069815947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78945827) q[1];
sx q[1];
rz(-0.43757868) q[1];
sx q[1];
rz(0.34559135) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9633767) q[3];
sx q[3];
rz(-1.6807739) q[3];
sx q[3];
rz(2.4770155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0871206) q[2];
sx q[2];
rz(-2.8643769) q[2];
sx q[2];
rz(0.37230125) q[2];
rz(0.39792684) q[3];
sx q[3];
rz(-1.1436661) q[3];
sx q[3];
rz(0.00042375617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.18846866) q[0];
sx q[0];
rz(-0.57250452) q[0];
sx q[0];
rz(1.5484126) q[0];
rz(-0.36987034) q[1];
sx q[1];
rz(-1.7431755) q[1];
sx q[1];
rz(2.5028548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34015481) q[0];
sx q[0];
rz(-1.7246913) q[0];
sx q[0];
rz(1.811054) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0453141) q[2];
sx q[2];
rz(-1.3945701) q[2];
sx q[2];
rz(-0.8187364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6795599) q[1];
sx q[1];
rz(-2.5278494) q[1];
sx q[1];
rz(-2.2824953) q[1];
rz(-pi) q[2];
rz(2.1356167) q[3];
sx q[3];
rz(-0.86274877) q[3];
sx q[3];
rz(-2.4410409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0729735) q[2];
sx q[2];
rz(-2.0899453) q[2];
sx q[2];
rz(0.2529141) q[2];
rz(0.84826338) q[3];
sx q[3];
rz(-1.4784644) q[3];
sx q[3];
rz(-3.0566791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10097583) q[0];
sx q[0];
rz(-2.9315797) q[0];
sx q[0];
rz(-2.9926391) q[0];
rz(1.8303998) q[1];
sx q[1];
rz(-1.7313749) q[1];
sx q[1];
rz(3.0874918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.220517) q[0];
sx q[0];
rz(-2.069325) q[0];
sx q[0];
rz(1.4709298) q[0];
rz(-pi) q[1];
x q[1];
rz(1.615754) q[2];
sx q[2];
rz(-2.3727594) q[2];
sx q[2];
rz(-2.2368778) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5230519) q[1];
sx q[1];
rz(-2.4344011) q[1];
sx q[1];
rz(2.9721353) q[1];
rz(-1.5255494) q[3];
sx q[3];
rz(-0.47569573) q[3];
sx q[3];
rz(3.0632116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3397303) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(-0.96699634) q[2];
rz(-2.4407834) q[3];
sx q[3];
rz(-0.98027027) q[3];
sx q[3];
rz(0.35082671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12857777) q[0];
sx q[0];
rz(-1.1739434) q[0];
sx q[0];
rz(1.5933734) q[0];
rz(-0.37823996) q[1];
sx q[1];
rz(-0.75235569) q[1];
sx q[1];
rz(2.7517448) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5711576) q[0];
sx q[0];
rz(-1.7661372) q[0];
sx q[0];
rz(0.51049149) q[0];
rz(1.8748449) q[2];
sx q[2];
rz(-1.5866733) q[2];
sx q[2];
rz(2.7823663) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15738641) q[1];
sx q[1];
rz(-2.0573924) q[1];
sx q[1];
rz(1.9686048) q[1];
rz(2.8865783) q[3];
sx q[3];
rz(-1.7863635) q[3];
sx q[3];
rz(0.6620342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0924015) q[2];
sx q[2];
rz(-1.1213877) q[2];
sx q[2];
rz(1.258705) q[2];
rz(2.8973268) q[3];
sx q[3];
rz(-1.786307) q[3];
sx q[3];
rz(-0.99610966) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8213537) q[0];
sx q[0];
rz(-1.2293674) q[0];
sx q[0];
rz(-1.7852596) q[0];
rz(2.9755196) q[1];
sx q[1];
rz(-0.95172721) q[1];
sx q[1];
rz(-2.70539) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1937418) q[0];
sx q[0];
rz(-2.0213038) q[0];
sx q[0];
rz(1.2240901) q[0];
x q[1];
rz(-2.2460552) q[2];
sx q[2];
rz(-0.51273275) q[2];
sx q[2];
rz(-2.2487179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6375745) q[1];
sx q[1];
rz(-2.3406041) q[1];
sx q[1];
rz(-1.2474508) q[1];
x q[2];
rz(1.5028788) q[3];
sx q[3];
rz(-1.1350067) q[3];
sx q[3];
rz(-0.21623789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1205552) q[2];
sx q[2];
rz(-1.977481) q[2];
sx q[2];
rz(0.54350054) q[2];
rz(0.049526878) q[3];
sx q[3];
rz(-1.6560358) q[3];
sx q[3];
rz(1.653695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.6108342) q[0];
sx q[0];
rz(-3.0093091) q[0];
sx q[0];
rz(0.079205967) q[0];
rz(-2.7414956) q[1];
sx q[1];
rz(-2.2875417) q[1];
sx q[1];
rz(-1.0265464) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8456789) q[0];
sx q[0];
rz(-2.2845553) q[0];
sx q[0];
rz(2.8299324) q[0];
x q[1];
rz(1.4756938) q[2];
sx q[2];
rz(-1.9355023) q[2];
sx q[2];
rz(-1.1191747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76128298) q[1];
sx q[1];
rz(-1.4766465) q[1];
sx q[1];
rz(-2.8172819) q[1];
x q[2];
rz(2.9273716) q[3];
sx q[3];
rz(-2.7832) q[3];
sx q[3];
rz(2.551799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2058699) q[2];
sx q[2];
rz(-1.8368072) q[2];
sx q[2];
rz(-1.6892461) q[2];
rz(0.82990372) q[3];
sx q[3];
rz(-2.2645576) q[3];
sx q[3];
rz(-0.95480603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6185388) q[0];
sx q[0];
rz(-1.1445615) q[0];
sx q[0];
rz(-1.6339697) q[0];
rz(1.2670831) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(0.98536678) q[2];
sx q[2];
rz(-1.1078688) q[2];
sx q[2];
rz(-2.3175793) q[2];
rz(2.5821154) q[3];
sx q[3];
rz(-1.8636129) q[3];
sx q[3];
rz(-1.158798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
