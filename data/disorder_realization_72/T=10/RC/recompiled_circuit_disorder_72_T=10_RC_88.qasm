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
rz(1.7680661) q[0];
sx q[0];
rz(11.058523) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7802785) q[0];
sx q[0];
rz(-1.2957934) q[0];
sx q[0];
rz(-1.6367903) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17272858) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(-1.2781065) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0807304) q[1];
sx q[1];
rz(-2.3002491) q[1];
sx q[1];
rz(-1.3688449) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8226837) q[3];
sx q[3];
rz(-3.0348274) q[3];
sx q[3];
rz(1.8330542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22380655) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28264499) q[0];
sx q[0];
rz(-2.1574321) q[0];
sx q[0];
rz(2.100201) q[0];
rz(1.7582714) q[2];
sx q[2];
rz(-1.4008998) q[2];
sx q[2];
rz(2.1263009) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.037307449) q[1];
sx q[1];
rz(-2.1530495) q[1];
sx q[1];
rz(-2.9557455) q[1];
rz(-2.6456804) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(-1.0752614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.26943794) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(-2.9906452) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(-0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040366216) q[0];
sx q[0];
rz(-1.7352312) q[0];
sx q[0];
rz(0.40584392) q[0];
rz(-pi) q[1];
rz(2.1951139) q[2];
sx q[2];
rz(-0.4779856) q[2];
sx q[2];
rz(1.9566655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90244251) q[1];
sx q[1];
rz(-2.855774) q[1];
sx q[1];
rz(-1.6509389) q[1];
rz(-pi) q[2];
rz(-1.4530573) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(-1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3016004) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(2.278573) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7097968) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(-2.5298932) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69458354) q[0];
sx q[0];
rz(-2.3059079) q[0];
sx q[0];
rz(1.6741333) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8407211) q[2];
sx q[2];
rz(-0.30756018) q[2];
sx q[2];
rz(-0.077066271) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1830814) q[1];
sx q[1];
rz(-1.6861048) q[1];
sx q[1];
rz(-0.99460852) q[1];
rz(-pi) q[2];
rz(1.6226193) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(0.27553976) q[2];
rz(3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(-2.2263288) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1268069) q[0];
sx q[0];
rz(-2.1984587) q[0];
sx q[0];
rz(0.058237596) q[0];
x q[1];
rz(-2.6164242) q[2];
sx q[2];
rz(-1.5319954) q[2];
sx q[2];
rz(-0.9135439) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1716869) q[1];
sx q[1];
rz(-1.7261506) q[1];
sx q[1];
rz(-3.1093662) q[1];
rz(-pi) q[2];
rz(0.74113412) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(2.4853064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(2.6780224) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-0.39658305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81178938) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(2.0407709) q[0];
rz(2.5846892) q[2];
sx q[2];
rz(-2.5957426) q[2];
sx q[2];
rz(-2.5715373) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.035316) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(-1.7229401) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58432213) q[3];
sx q[3];
rz(-1.366368) q[3];
sx q[3];
rz(2.1401329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(-1.5552103) q[2];
rz(-1.4554626) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-0.78480762) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(1.126359) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5267245) q[0];
sx q[0];
rz(-1.4818026) q[0];
sx q[0];
rz(0.29863775) q[0];
rz(-pi) q[1];
rz(-0.14234219) q[2];
sx q[2];
rz(-1.0668313) q[2];
sx q[2];
rz(-0.0056643639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4451051) q[1];
sx q[1];
rz(-1.3828053) q[1];
sx q[1];
rz(-2.6871215) q[1];
rz(-pi) q[2];
rz(-1.2231636) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(-1.7355433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(3.0287108) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927133) q[0];
sx q[0];
rz(-1.8939051) q[0];
sx q[0];
rz(-0.26922853) q[0];
rz(-pi) q[1];
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
rz(2.7777562) q[1];
sx q[1];
rz(-1.6297852) q[1];
sx q[1];
rz(-1.1532409) q[1];
rz(-pi) q[2];
x q[2];
rz(1.950374) q[3];
sx q[3];
rz(-1.4879585) q[3];
sx q[3];
rz(-2.6759202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-2.5194871) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(0.48535767) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94333121) q[0];
sx q[0];
rz(-1.7051464) q[0];
sx q[0];
rz(1.5166548) q[0];
rz(-pi) q[1];
x q[1];
rz(2.543407) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(-2.8796632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24385142) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(0.41390151) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58166196) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(-1.4130842) q[3];
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
rz(0.1594485) q[2];
rz(-1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52206802) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(0.5982582) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3529417) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.7822595) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5355849) q[2];
sx q[2];
rz(-2.8851644) q[2];
sx q[2];
rz(2.3026349) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.88565247) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(0.10140681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4062905) q[3];
sx q[3];
rz(-1.7269772) q[3];
sx q[3];
rz(2.3567049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.2364173) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-2.5219953) q[2];
sx q[2];
rz(-1.6710812) q[2];
sx q[2];
rz(-2.8166213) q[2];
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
