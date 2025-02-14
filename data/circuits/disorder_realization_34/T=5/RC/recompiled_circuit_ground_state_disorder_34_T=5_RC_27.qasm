OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(4.0255044) q[0];
sx q[0];
rz(8.5691353) q[0];
rz(2.9698676) q[1];
sx q[1];
rz(-3.0260234) q[1];
sx q[1];
rz(-2.5791383) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0732361) q[0];
sx q[0];
rz(-2.773065) q[0];
sx q[0];
rz(-1.7191528) q[0];
rz(-2.0777656) q[2];
sx q[2];
rz(-1.4344782) q[2];
sx q[2];
rz(2.9796114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6537498) q[1];
sx q[1];
rz(-1.1636793) q[1];
sx q[1];
rz(-2.6832677) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7482618) q[3];
sx q[3];
rz(-1.1926023) q[3];
sx q[3];
rz(0.26547576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2354551) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(-2.3423024) q[2];
rz(2.6702787) q[3];
sx q[3];
rz(-2.2282232) q[3];
sx q[3];
rz(-2.2057064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1021295) q[0];
sx q[0];
rz(-1.8164182) q[0];
sx q[0];
rz(-1.7934196) q[0];
rz(-1.7680291) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(-1.1522393) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7822398) q[0];
sx q[0];
rz(-2.3639661) q[0];
sx q[0];
rz(-1.1471143) q[0];
rz(0.94448467) q[2];
sx q[2];
rz(-0.99530333) q[2];
sx q[2];
rz(1.6863509) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5630398) q[1];
sx q[1];
rz(-2.4336947) q[1];
sx q[1];
rz(0.25917525) q[1];
rz(0.82291921) q[3];
sx q[3];
rz(-2.4429818) q[3];
sx q[3];
rz(-0.60958344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79933244) q[2];
sx q[2];
rz(-0.41877425) q[2];
sx q[2];
rz(2.2104134) q[2];
rz(0.10654199) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0837285) q[0];
sx q[0];
rz(-1.7513542) q[0];
sx q[0];
rz(-0.51783836) q[0];
rz(-0.88874108) q[1];
sx q[1];
rz(-0.70777142) q[1];
sx q[1];
rz(-0.4471561) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6155375) q[0];
sx q[0];
rz(-1.8319905) q[0];
sx q[0];
rz(-1.5047856) q[0];
rz(1.0165527) q[2];
sx q[2];
rz(-2.2187761) q[2];
sx q[2];
rz(-0.78314834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67990002) q[1];
sx q[1];
rz(-2.5352806) q[1];
sx q[1];
rz(1.3316505) q[1];
rz(-1.8497883) q[3];
sx q[3];
rz(-2.0800262) q[3];
sx q[3];
rz(1.694547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7670224) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(-2.6089597) q[2];
rz(-2.8940708) q[3];
sx q[3];
rz(-2.4041924) q[3];
sx q[3];
rz(-1.0386764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5175051) q[0];
sx q[0];
rz(-3.0563323) q[0];
sx q[0];
rz(-0.10661539) q[0];
rz(2.7952349) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(-1.9812298) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24460565) q[0];
sx q[0];
rz(-1.7443313) q[0];
sx q[0];
rz(-0.24994295) q[0];
rz(0.82616634) q[2];
sx q[2];
rz(-2.1414504) q[2];
sx q[2];
rz(-2.014239) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2109769) q[1];
sx q[1];
rz(-1.9257571) q[1];
sx q[1];
rz(-0.32101722) q[1];
rz(-2.6909573) q[3];
sx q[3];
rz(-2.051578) q[3];
sx q[3];
rz(-1.6893051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13566636) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(-3.0653817) q[2];
rz(0.58586079) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(1.7990254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0214486) q[0];
sx q[0];
rz(-1.8360538) q[0];
sx q[0];
rz(2.7914877) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.8616385) q[1];
sx q[1];
rz(1.9116481) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1317539) q[0];
sx q[0];
rz(-1.3309877) q[0];
sx q[0];
rz(1.1935704) q[0];
rz(-pi) q[1];
rz(-0.036418865) q[2];
sx q[2];
rz(-1.4055739) q[2];
sx q[2];
rz(-2.9584004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38795162) q[1];
sx q[1];
rz(-1.759583) q[1];
sx q[1];
rz(-1.2354047) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8660981) q[3];
sx q[3];
rz(-2.5341138) q[3];
sx q[3];
rz(-0.5831843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89796394) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.492307) q[2];
rz(-0.40488511) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(-2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355857) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(0.58240044) q[0];
rz(0.68663418) q[1];
sx q[1];
rz(-1.735894) q[1];
sx q[1];
rz(-1.883421) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95520407) q[0];
sx q[0];
rz(-1.135728) q[0];
sx q[0];
rz(-3.0032936) q[0];
rz(2.1297087) q[2];
sx q[2];
rz(-0.61043533) q[2];
sx q[2];
rz(-1.287078) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68738378) q[1];
sx q[1];
rz(-1.6404248) q[1];
sx q[1];
rz(-1.5147665) q[1];
rz(-pi) q[2];
rz(1.2518796) q[3];
sx q[3];
rz(-2.5237759) q[3];
sx q[3];
rz(1.8998673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76576343) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(-1.0180391) q[2];
rz(0.24614075) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(1.6233981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4356284) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(2.9980581) q[1];
sx q[1];
rz(-0.48964557) q[1];
sx q[1];
rz(-0.20763436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.62791) q[0];
sx q[0];
rz(-1.2259036) q[0];
sx q[0];
rz(2.3570127) q[0];
rz(-pi) q[1];
rz(2.1499499) q[2];
sx q[2];
rz(-1.275838) q[2];
sx q[2];
rz(0.3504914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0032867) q[1];
sx q[1];
rz(-2.3992743) q[1];
sx q[1];
rz(1.9600541) q[1];
rz(-pi) q[2];
rz(1.666388) q[3];
sx q[3];
rz(-1.2089653) q[3];
sx q[3];
rz(1.7887448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2988854) q[2];
sx q[2];
rz(-1.8536114) q[2];
sx q[2];
rz(-0.60619727) q[2];
rz(0.68743622) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(0.70639759) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48924482) q[0];
sx q[0];
rz(-3.1135961) q[0];
sx q[0];
rz(-2.0835173) q[0];
rz(3.0319013) q[1];
sx q[1];
rz(-2.0249764) q[1];
sx q[1];
rz(-1.441997) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88272754) q[0];
sx q[0];
rz(-1.8804714) q[0];
sx q[0];
rz(-2.726) q[0];
rz(1.8462366) q[2];
sx q[2];
rz(-2.0552962) q[2];
sx q[2];
rz(0.4415919) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.40934) q[1];
sx q[1];
rz(-0.85459083) q[1];
sx q[1];
rz(-2.2316037) q[1];
x q[2];
rz(2.9061716) q[3];
sx q[3];
rz(-1.306433) q[3];
sx q[3];
rz(-2.2478769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.480392) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(-1.9332168) q[2];
rz(-2.4783573) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(-1.2164345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97063589) q[0];
sx q[0];
rz(-2.1747776) q[0];
sx q[0];
rz(-3.048625) q[0];
rz(1.2804821) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(3.0063937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2743535) q[0];
sx q[0];
rz(-1.6238191) q[0];
sx q[0];
rz(-1.6163325) q[0];
rz(-pi) q[1];
rz(-2.6117295) q[2];
sx q[2];
rz(-0.29680291) q[2];
sx q[2];
rz(2.5713845) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0307903) q[1];
sx q[1];
rz(-0.25486481) q[1];
sx q[1];
rz(-3.0569264) q[1];
x q[2];
rz(1.7235123) q[3];
sx q[3];
rz(-1.4974606) q[3];
sx q[3];
rz(-0.45460816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4325503) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(-0.77152983) q[2];
rz(-2.6500474) q[3];
sx q[3];
rz(-1.8621657) q[3];
sx q[3];
rz(0.5947203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5563357) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(-2.0167895) q[0];
rz(-1.8428165) q[1];
sx q[1];
rz(-2.5262084) q[1];
sx q[1];
rz(-2.3769456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2128396) q[0];
sx q[0];
rz(-3.0055586) q[0];
sx q[0];
rz(-2.4555951) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49955274) q[2];
sx q[2];
rz(-2.1883114) q[2];
sx q[2];
rz(2.2411186) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99801714) q[1];
sx q[1];
rz(-1.7735266) q[1];
sx q[1];
rz(0.12938093) q[1];
rz(1.1922791) q[3];
sx q[3];
rz(-1.0487742) q[3];
sx q[3];
rz(-1.4498364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7654968) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(0.20720227) q[2];
rz(2.1616705) q[3];
sx q[3];
rz(-2.0082974) q[3];
sx q[3];
rz(0.19237147) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021066396) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(-0.5008685) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(2.4293368) q[2];
sx q[2];
rz(-2.9604572) q[2];
sx q[2];
rz(-1.0618718) q[2];
rz(-1.408314) q[3];
sx q[3];
rz(-1.8058895) q[3];
sx q[3];
rz(-1.4598082) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
