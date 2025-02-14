OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5496991) q[0];
sx q[0];
rz(-1.1130207) q[0];
sx q[0];
rz(9.6564403) q[0];
rz(0.35056937) q[1];
sx q[1];
rz(6.0729519) q[1];
sx q[1];
rz(8.6624866) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0889083) q[0];
sx q[0];
rz(-1.9885953) q[0];
sx q[0];
rz(2.4521928) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0341124) q[2];
sx q[2];
rz(-1.1321403) q[2];
sx q[2];
rz(-0.92213501) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2628535) q[1];
sx q[1];
rz(-1.7567506) q[1];
sx q[1];
rz(1.5145709) q[1];
rz(2.2985994) q[3];
sx q[3];
rz(-2.5268433) q[3];
sx q[3];
rz(0.45045567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1594557) q[2];
sx q[2];
rz(-0.60672131) q[2];
sx q[2];
rz(-0.43178001) q[2];
rz(-0.86709705) q[3];
sx q[3];
rz(-2.1853787) q[3];
sx q[3];
rz(-1.506327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1113488) q[0];
sx q[0];
rz(-2.3115277) q[0];
sx q[0];
rz(-1.1907955) q[0];
rz(1.2442773) q[1];
sx q[1];
rz(-1.1294653) q[1];
sx q[1];
rz(2.5608889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5842019) q[0];
sx q[0];
rz(-1.4317436) q[0];
sx q[0];
rz(-1.1625421) q[0];
rz(-pi) q[1];
rz(2.2207308) q[2];
sx q[2];
rz(-2.5061786) q[2];
sx q[2];
rz(2.970545) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8088219) q[1];
sx q[1];
rz(-2.1968107) q[1];
sx q[1];
rz(-1.9653734) q[1];
rz(-pi) q[2];
x q[2];
rz(2.689183) q[3];
sx q[3];
rz(-0.93660855) q[3];
sx q[3];
rz(0.89434552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1397436) q[2];
sx q[2];
rz(-1.7364343) q[2];
sx q[2];
rz(1.5025567) q[2];
rz(-2.9338845) q[3];
sx q[3];
rz(-1.8307056) q[3];
sx q[3];
rz(1.0072072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.945632) q[0];
sx q[0];
rz(-2.0646844) q[0];
sx q[0];
rz(-0.14863241) q[0];
rz(2.0678988) q[1];
sx q[1];
rz(-0.37207347) q[1];
sx q[1];
rz(2.3587842) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40227213) q[0];
sx q[0];
rz(-1.625652) q[0];
sx q[0];
rz(1.7412454) q[0];
rz(-pi) q[1];
rz(0.53289588) q[2];
sx q[2];
rz(-1.6072465) q[2];
sx q[2];
rz(-2.8849071) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2116993) q[1];
sx q[1];
rz(-0.5553588) q[1];
sx q[1];
rz(0.38513215) q[1];
rz(-pi) q[2];
rz(1.820676) q[3];
sx q[3];
rz(-0.80554038) q[3];
sx q[3];
rz(-0.66242927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17915501) q[2];
sx q[2];
rz(-0.97426263) q[2];
sx q[2];
rz(-1.1243189) q[2];
rz(3.1149241) q[3];
sx q[3];
rz(-0.6921488) q[3];
sx q[3];
rz(-2.8587225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7836595) q[0];
sx q[0];
rz(-1.986035) q[0];
sx q[0];
rz(1.5643157) q[0];
rz(-2.0951001) q[1];
sx q[1];
rz(-1.0341897) q[1];
sx q[1];
rz(1.1276721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2593721) q[0];
sx q[0];
rz(-1.581107) q[0];
sx q[0];
rz(-2.5502648) q[0];
x q[1];
rz(0.079795795) q[2];
sx q[2];
rz(-1.2376806) q[2];
sx q[2];
rz(3.0175957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2232957) q[1];
sx q[1];
rz(-2.363829) q[1];
sx q[1];
rz(0.4593846) q[1];
x q[2];
rz(-2.2644482) q[3];
sx q[3];
rz(-2.9877285) q[3];
sx q[3];
rz(2.2779704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6005738) q[2];
sx q[2];
rz(-1.270673) q[2];
sx q[2];
rz(-0.18826558) q[2];
rz(2.6514734) q[3];
sx q[3];
rz(-1.6907588) q[3];
sx q[3];
rz(1.274444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.984943) q[0];
sx q[0];
rz(-1.5644512) q[0];
sx q[0];
rz(1.6749325) q[0];
rz(-0.4725534) q[1];
sx q[1];
rz(-0.87635374) q[1];
sx q[1];
rz(-2.1579425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0839) q[0];
sx q[0];
rz(-2.6760457) q[0];
sx q[0];
rz(-2.7280612) q[0];
rz(-1.7631268) q[2];
sx q[2];
rz(-2.0266313) q[2];
sx q[2];
rz(-2.6050732) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3762665) q[1];
sx q[1];
rz(-1.5137356) q[1];
sx q[1];
rz(2.0817134) q[1];
x q[2];
rz(-1.8365588) q[3];
sx q[3];
rz(-1.467622) q[3];
sx q[3];
rz(-2.3931695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5094362) q[2];
sx q[2];
rz(-1.7551273) q[2];
sx q[2];
rz(-2.7344088) q[2];
rz(1.2716278) q[3];
sx q[3];
rz(-0.41912246) q[3];
sx q[3];
rz(2.4841888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42095175) q[0];
sx q[0];
rz(-1.9176418) q[0];
sx q[0];
rz(1.2366914) q[0];
rz(-1.1300794) q[1];
sx q[1];
rz(-1.7944444) q[1];
sx q[1];
rz(-2.0225661) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2988855) q[0];
sx q[0];
rz(-1.7241553) q[0];
sx q[0];
rz(1.6394079) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97864464) q[2];
sx q[2];
rz(-0.8675163) q[2];
sx q[2];
rz(2.5539309) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31737008) q[1];
sx q[1];
rz(-0.44840773) q[1];
sx q[1];
rz(0.3387785) q[1];
rz(-pi) q[2];
rz(-1.9672333) q[3];
sx q[3];
rz(-1.160495) q[3];
sx q[3];
rz(1.0956216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59469026) q[2];
sx q[2];
rz(-2.0141527) q[2];
sx q[2];
rz(0.31000578) q[2];
rz(1.7075432) q[3];
sx q[3];
rz(-1.6238345) q[3];
sx q[3];
rz(1.8967418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2731648) q[0];
sx q[0];
rz(-2.7719066) q[0];
sx q[0];
rz(2.0694859) q[0];
rz(1.781069) q[1];
sx q[1];
rz(-0.37630263) q[1];
sx q[1];
rz(1.2881813) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53258713) q[0];
sx q[0];
rz(-0.83206785) q[0];
sx q[0];
rz(-0.84858507) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.200618) q[2];
sx q[2];
rz(-0.67747203) q[2];
sx q[2];
rz(-1.7446454) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.6009734) q[1];
sx q[1];
rz(-0.29635591) q[1];
sx q[1];
rz(-0.61295103) q[1];
rz(-2.6374972) q[3];
sx q[3];
rz(-2.305784) q[3];
sx q[3];
rz(0.1288165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0518904) q[2];
sx q[2];
rz(-1.6334198) q[2];
sx q[2];
rz(-1.4349597) q[2];
rz(0.47576225) q[3];
sx q[3];
rz(-0.69403726) q[3];
sx q[3];
rz(2.7580875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39110228) q[0];
sx q[0];
rz(-0.35677156) q[0];
sx q[0];
rz(-1.9619036) q[0];
rz(-2.7791595) q[1];
sx q[1];
rz(-1.06203) q[1];
sx q[1];
rz(-1.5550044) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060136) q[0];
sx q[0];
rz(-1.4188674) q[0];
sx q[0];
rz(-2.8169027) q[0];
x q[1];
rz(-0.68029515) q[2];
sx q[2];
rz(-1.1914502) q[2];
sx q[2];
rz(2.2915886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5240979) q[1];
sx q[1];
rz(-0.84071181) q[1];
sx q[1];
rz(0.9690271) q[1];
x q[2];
rz(1.8170023) q[3];
sx q[3];
rz(-2.559767) q[3];
sx q[3];
rz(-1.8210321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4156551) q[2];
sx q[2];
rz(-1.1154117) q[2];
sx q[2];
rz(1.1768781) q[2];
rz(-0.40500179) q[3];
sx q[3];
rz(-2.160725) q[3];
sx q[3];
rz(-0.39308959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28615752) q[0];
sx q[0];
rz(-3.0431008) q[0];
sx q[0];
rz(-1.3354907) q[0];
rz(-0.91241854) q[1];
sx q[1];
rz(-1.7203947) q[1];
sx q[1];
rz(-3.1013427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6326014) q[0];
sx q[0];
rz(-0.17290056) q[0];
sx q[0];
rz(1.7937051) q[0];
rz(2.5435872) q[2];
sx q[2];
rz(-1.0710508) q[2];
sx q[2];
rz(2.6464484) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9736824) q[1];
sx q[1];
rz(-1.0492322) q[1];
sx q[1];
rz(0.48814623) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3586799) q[3];
sx q[3];
rz(-0.81394845) q[3];
sx q[3];
rz(2.1569686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6054298) q[2];
sx q[2];
rz(-1.9998113) q[2];
sx q[2];
rz(-2.6119168) q[2];
rz(-1.4179432) q[3];
sx q[3];
rz(-0.46897408) q[3];
sx q[3];
rz(2.6403707) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0025075992) q[0];
sx q[0];
rz(-1.7739828) q[0];
sx q[0];
rz(-0.031877192) q[0];
rz(-1.2079976) q[1];
sx q[1];
rz(-0.54787689) q[1];
sx q[1];
rz(-0.30778232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14292158) q[0];
sx q[0];
rz(-1.0772155) q[0];
sx q[0];
rz(0.30360766) q[0];
x q[1];
rz(1.1319699) q[2];
sx q[2];
rz(-2.8487051) q[2];
sx q[2];
rz(-1.9856118) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.44445723) q[1];
sx q[1];
rz(-0.83205497) q[1];
sx q[1];
rz(0.8529128) q[1];
x q[2];
rz(0.73460893) q[3];
sx q[3];
rz(-0.19140581) q[3];
sx q[3];
rz(-0.020343971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.53177437) q[2];
sx q[2];
rz(-0.69963026) q[2];
sx q[2];
rz(-3.0979274) q[2];
rz(-1.6627436) q[3];
sx q[3];
rz(-2.6378529) q[3];
sx q[3];
rz(1.7536722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26535784) q[0];
sx q[0];
rz(-1.6725412) q[0];
sx q[0];
rz(-1.5681736) q[0];
rz(-0.028989446) q[1];
sx q[1];
rz(-1.1075167) q[1];
sx q[1];
rz(0.97069244) q[1];
rz(-0.12262298) q[2];
sx q[2];
rz(-1.7896252) q[2];
sx q[2];
rz(0.26364506) q[2];
rz(0.51359691) q[3];
sx q[3];
rz(-0.83181341) q[3];
sx q[3];
rz(0.60461525) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
