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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36131418) q[0];
sx q[0];
rz(-1.8457992) q[0];
sx q[0];
rz(1.5048024) q[0];
x q[1];
rz(2.2917065) q[2];
sx q[2];
rz(-1.7011142) q[2];
sx q[2];
rz(2.7352114) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0807304) q[1];
sx q[1];
rz(-2.3002491) q[1];
sx q[1];
rz(-1.7727477) q[1];
rz(-1.6742168) q[3];
sx q[3];
rz(-1.5442344) q[3];
sx q[3];
rz(0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(-2.9075918) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(2.1388334) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53045814) q[0];
sx q[0];
rz(-0.76871745) q[0];
sx q[0];
rz(0.64972933) q[0];
rz(-pi) q[1];
rz(-1.7582714) q[2];
sx q[2];
rz(-1.7406929) q[2];
sx q[2];
rz(2.1263009) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7111295) q[1];
sx q[1];
rz(-1.4158447) q[1];
sx q[1];
rz(0.98053812) q[1];
x q[2];
rz(0.49591222) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(-1.0752614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26943794) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(-0.41444591) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(-0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(2.2580106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2471836) q[0];
sx q[0];
rz(-0.43617019) q[0];
sx q[0];
rz(2.7437074) q[0];
rz(-2.8475464) q[2];
sx q[2];
rz(-1.1883192) q[2];
sx q[2];
rz(-1.2750212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5501432) q[1];
sx q[1];
rz(-1.5933697) q[1];
sx q[1];
rz(-1.2858461) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6885353) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(-1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(-0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(0.5501737) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94565369) q[0];
sx q[0];
rz(-1.494207) q[0];
sx q[0];
rz(-2.4038195) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8407211) q[2];
sx q[2];
rz(-0.30756018) q[2];
sx q[2];
rz(0.077066271) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.354047) q[1];
sx q[1];
rz(-0.58632942) q[1];
sx q[1];
rz(-1.7802618) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5189734) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(-2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(-0.91526389) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6633776) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(-2.1992654) q[0];
x q[1];
rz(1.615633) q[2];
sx q[2];
rz(-1.0460639) q[2];
sx q[2];
rz(2.5068138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76449672) q[1];
sx q[1];
rz(-0.15863523) q[1];
sx q[1];
rz(-1.3678958) q[1];
rz(-0.81766537) q[3];
sx q[3];
rz(-0.92901232) q[3];
sx q[3];
rz(1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(2.0625431) q[0];
rz(2.5462529) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(2.7450096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2561589) q[0];
sx q[0];
rz(-2.0246756) q[0];
sx q[0];
rz(-2.8594349) q[0];
rz(2.6655212) q[2];
sx q[2];
rz(-1.2928315) q[2];
sx q[2];
rz(-1.6518041) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5844333) q[1];
sx q[1];
rz(-1.6646619) q[1];
sx q[1];
rz(2.2319016) q[1];
x q[2];
rz(1.8144238) q[3];
sx q[3];
rz(-2.1414087) q[3];
sx q[3];
rz(-2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5863824) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-2.0152337) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267245) q[0];
sx q[0];
rz(-1.4818026) q[0];
sx q[0];
rz(2.8429549) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8225841) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(0.28282794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78325242) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9184291) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(-1.7355433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9324947) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.37384) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2643124) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(-0.89950048) q[0];
x q[1];
rz(2.1342437) q[2];
sx q[2];
rz(-1.9512366) q[2];
sx q[2];
rz(-3.0055339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0667463) q[1];
sx q[1];
rz(-2.7201338) q[1];
sx q[1];
rz(-1.4261817) q[1];
rz(-pi) q[2];
rz(-0.0891536) q[3];
sx q[3];
rz(-1.1925863) q[3];
sx q[3];
rz(2.0034727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-2.7511403) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.6149678) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-2.656235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55886666) q[0];
sx q[0];
rz(-2.9968046) q[0];
sx q[0];
rz(-0.38082122) q[0];
rz(-pi) q[1];
rz(0.59818563) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(-0.26192947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8977412) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(-2.7276911) q[1];
rz(-pi) q[2];
rz(-2.5599307) q[3];
sx q[3];
rz(-1.9546486) q[3];
sx q[3];
rz(1.4130842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(1.4341226) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(0.5982582) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1372676) q[0];
sx q[0];
rz(-2.8398872) q[0];
sx q[0];
rz(0.76122491) q[0];
x q[1];
rz(0.21226378) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-3.0002909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4926589) q[1];
sx q[1];
rz(-1.6655386) q[1];
sx q[1];
rz(1.2050864) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80501276) q[3];
sx q[3];
rz(-0.22634889) q[3];
sx q[3];
rz(-1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.2364173) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-0.17157208) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(2.0486352) q[3];
sx q[3];
rz(-0.7328877) q[3];
sx q[3];
rz(2.4536798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
