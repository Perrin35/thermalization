OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(-0.49931061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12257523) q[0];
sx q[0];
rz(-0.28261533) q[0];
sx q[0];
rz(-0.22960381) q[0];
rz(-pi) q[1];
rz(0.84988611) q[2];
sx q[2];
rz(-1.4404785) q[2];
sx q[2];
rz(2.7352114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7826739) q[1];
sx q[1];
rz(-0.75190699) q[1];
sx q[1];
rz(0.2207252) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1148881) q[3];
sx q[3];
rz(-1.6741802) q[3];
sx q[3];
rz(2.0863233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22380655) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(-2.1271558) q[2];
rz(-2.9075918) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4734128) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(-1.0043253) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-1.0027592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6111345) q[0];
sx q[0];
rz(-2.3728752) q[0];
sx q[0];
rz(-2.4918633) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82680306) q[2];
sx q[2];
rz(-0.25233341) q[2];
sx q[2];
rz(1.2834872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.037307449) q[1];
sx q[1];
rz(-0.98854317) q[1];
sx q[1];
rz(-0.18584713) q[1];
rz(1.9243794) q[3];
sx q[3];
rz(-2.1401322) q[3];
sx q[3];
rz(-2.6667037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(0.41444591) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(-0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931483) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(-0.88358203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012264) q[0];
sx q[0];
rz(-1.4063615) q[0];
sx q[0];
rz(2.7357487) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1951139) q[2];
sx q[2];
rz(-0.4779856) q[2];
sx q[2];
rz(-1.9566655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59144943) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(-1.2858461) q[1];
rz(-pi) q[2];
rz(2.9900842) q[3];
sx q[3];
rz(-2.476859) q[3];
sx q[3];
rz(2.1745149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-2.591419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.195939) q[0];
sx q[0];
rz(-1.494207) q[0];
sx q[0];
rz(0.73777317) q[0];
x q[1];
rz(2.8470464) q[2];
sx q[2];
rz(-1.6606332) q[2];
sx q[2];
rz(1.3603269) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95851129) q[1];
sx q[1];
rz(-1.6861048) q[1];
sx q[1];
rz(2.1469841) q[1];
rz(-1.6226193) q[3];
sx q[3];
rz(-1.341815) q[3];
sx q[3];
rz(-0.94829544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52175534) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(-0.11166212) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-2.2263288) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278496) q[0];
sx q[0];
rz(-0.62999524) q[0];
sx q[0];
rz(1.4907452) q[0];
x q[1];
rz(-1.5259597) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(-2.5068138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1716869) q[1];
sx q[1];
rz(-1.7261506) q[1];
sx q[1];
rz(0.032226493) q[1];
rz(-2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(-1.4828009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9203732) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-2.7450096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30124861) q[0];
sx q[0];
rz(-2.6123971) q[0];
sx q[0];
rz(1.0521786) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55690342) q[2];
sx q[2];
rz(-2.5957426) q[2];
sx q[2];
rz(-0.57005537) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5571594) q[1];
sx q[1];
rz(-1.6646619) q[1];
sx q[1];
rz(-0.90969109) q[1];
rz(-pi) q[2];
rz(1.8144238) q[3];
sx q[3];
rz(-1.0001839) q[3];
sx q[3];
rz(2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(-1.5552103) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(-1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-1.126359) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92857498) q[0];
sx q[0];
rz(-1.2733766) q[0];
sx q[0];
rz(1.4777044) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0790714) q[2];
sx q[2];
rz(-1.4462496) q[2];
sx q[2];
rz(-1.6342271) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4451051) q[1];
sx q[1];
rz(-1.3828053) q[1];
sx q[1];
rz(2.6871215) q[1];
rz(-1.2231636) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(1.4060494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(-3.1088366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(-2.2420922) q[0];
x q[1];
rz(2.1342437) q[2];
sx q[2];
rz(-1.9512366) q[2];
sx q[2];
rz(0.13605875) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.36383648) q[1];
sx q[1];
rz(-1.6297852) q[1];
sx q[1];
rz(-1.9883518) q[1];
x q[2];
rz(-0.0891536) q[3];
sx q[3];
rz(-1.1925863) q[3];
sx q[3];
rz(-1.13812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(0.62210554) q[2];
rz(1.9744251) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099982925) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.6149678) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(0.48535767) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55886666) q[0];
sx q[0];
rz(-0.14478806) q[0];
sx q[0];
rz(-2.7607714) q[0];
rz(-2.543407) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(-0.26192947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7945054) q[1];
sx q[1];
rz(-0.80650389) q[1];
sx q[1];
rz(-2.0054714) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63391179) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-0.67542911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0043250672) q[0];
sx q[0];
rz(-2.8398872) q[0];
sx q[0];
rz(0.76122491) q[0];
rz(-1.422545) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(1.6809747) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64893374) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(1.9365063) q[1];
x q[2];
rz(2.3365799) q[3];
sx q[3];
rz(-0.22634889) q[3];
sx q[3];
rz(1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-0.65199488) q[2];
rz(0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.2364173) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(0.61959735) q[2];
sx q[2];
rz(-1.6710812) q[2];
sx q[2];
rz(-2.8166213) q[2];
rz(-2.7491309) q[3];
sx q[3];
rz(-2.2069208) q[3];
sx q[3];
rz(1.8451286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
