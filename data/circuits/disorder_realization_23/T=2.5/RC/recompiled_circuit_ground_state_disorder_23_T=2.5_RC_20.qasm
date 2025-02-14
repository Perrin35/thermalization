OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.58148122) q[0];
sx q[0];
rz(-1.886263) q[0];
sx q[0];
rz(-0.49129593) q[0];
rz(2.8331941) q[1];
sx q[1];
rz(-1.3003132) q[1];
sx q[1];
rz(0.058606776) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97293905) q[0];
sx q[0];
rz(-2.0117451) q[0];
sx q[0];
rz(2.1827374) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.357448) q[2];
sx q[2];
rz(-1.9040739) q[2];
sx q[2];
rz(0.12034697) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9323001) q[1];
sx q[1];
rz(-1.8718616) q[1];
sx q[1];
rz(-1.73897) q[1];
rz(-pi) q[2];
rz(1.8901615) q[3];
sx q[3];
rz(-1.6465558) q[3];
sx q[3];
rz(-3.0204242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0356174) q[2];
sx q[2];
rz(-2.3507037) q[2];
sx q[2];
rz(0.38113511) q[2];
rz(3.1086339) q[3];
sx q[3];
rz(-1.4573174) q[3];
sx q[3];
rz(1.0491252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9829262) q[0];
sx q[0];
rz(-0.52678147) q[0];
sx q[0];
rz(-2.7016933) q[0];
rz(-0.061821763) q[1];
sx q[1];
rz(-1.6114176) q[1];
sx q[1];
rz(2.8672112) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3456164) q[0];
sx q[0];
rz(-1.1152949) q[0];
sx q[0];
rz(2.2296395) q[0];
rz(1.9977278) q[2];
sx q[2];
rz(-1.3620268) q[2];
sx q[2];
rz(-2.9627299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0068757) q[1];
sx q[1];
rz(-1.1235079) q[1];
sx q[1];
rz(-0.083880977) q[1];
rz(0.37526826) q[3];
sx q[3];
rz(-0.57884848) q[3];
sx q[3];
rz(1.8414094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8476094) q[2];
sx q[2];
rz(-3.1276438) q[2];
sx q[2];
rz(-0.070276109) q[2];
rz(2.516732) q[3];
sx q[3];
rz(-1.8127541) q[3];
sx q[3];
rz(0.074987324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4570419) q[0];
sx q[0];
rz(-1.6383189) q[0];
sx q[0];
rz(0.3918089) q[0];
rz(-2.1458972) q[1];
sx q[1];
rz(-2.3052146) q[1];
sx q[1];
rz(-3.0973184) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05737722) q[0];
sx q[0];
rz(-0.65859933) q[0];
sx q[0];
rz(-2.088571) q[0];
rz(-pi) q[1];
rz(2.2703553) q[2];
sx q[2];
rz(-0.92513621) q[2];
sx q[2];
rz(-1.439656) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7024732) q[1];
sx q[1];
rz(-0.28077341) q[1];
sx q[1];
rz(-0.72661119) q[1];
rz(-1.6674394) q[3];
sx q[3];
rz(-2.1657888) q[3];
sx q[3];
rz(-2.4293824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5135045) q[2];
sx q[2];
rz(-2.2254483) q[2];
sx q[2];
rz(-1.8872895) q[2];
rz(0.11145505) q[3];
sx q[3];
rz(-2.1696551) q[3];
sx q[3];
rz(0.85533065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3067538) q[0];
sx q[0];
rz(-0.67973891) q[0];
sx q[0];
rz(0.87598959) q[0];
rz(0.39729473) q[1];
sx q[1];
rz(-1.5700424) q[1];
sx q[1];
rz(-1.809459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5529163) q[0];
sx q[0];
rz(-2.0400926) q[0];
sx q[0];
rz(1.3568715) q[0];
x q[1];
rz(0.45355637) q[2];
sx q[2];
rz(-0.95090129) q[2];
sx q[2];
rz(-0.93728055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5046103) q[1];
sx q[1];
rz(-0.17269293) q[1];
sx q[1];
rz(-0.47422282) q[1];
x q[2];
rz(-1.0334084) q[3];
sx q[3];
rz(-1.3674481) q[3];
sx q[3];
rz(-0.56768473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9925655) q[2];
sx q[2];
rz(-1.488648) q[2];
sx q[2];
rz(-1.7472902) q[2];
rz(-2.1483138) q[3];
sx q[3];
rz(-1.9781878) q[3];
sx q[3];
rz(-1.2629898) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8472327) q[0];
sx q[0];
rz(-2.9672186) q[0];
sx q[0];
rz(0.85884035) q[0];
rz(2.7186588) q[1];
sx q[1];
rz(-0.49086389) q[1];
sx q[1];
rz(1.6805958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.925526) q[0];
sx q[0];
rz(-0.51617814) q[0];
sx q[0];
rz(-0.25052089) q[0];
x q[1];
rz(-1.1255477) q[2];
sx q[2];
rz(-1.5684874) q[2];
sx q[2];
rz(-0.18773676) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8940058) q[1];
sx q[1];
rz(-1.3787706) q[1];
sx q[1];
rz(-2.9321666) q[1];
x q[2];
rz(2.0764396) q[3];
sx q[3];
rz(-1.6958456) q[3];
sx q[3];
rz(1.3920543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6964263) q[2];
sx q[2];
rz(-0.40881613) q[2];
sx q[2];
rz(-1.6430631) q[2];
rz(2.0131352) q[3];
sx q[3];
rz(-1.107629) q[3];
sx q[3];
rz(-2.4988153) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6978825) q[0];
sx q[0];
rz(-1.8118129) q[0];
sx q[0];
rz(-2.4238996) q[0];
rz(0.36092654) q[1];
sx q[1];
rz(-2.2310427) q[1];
sx q[1];
rz(-2.8275729) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0180394) q[0];
sx q[0];
rz(-1.529487) q[0];
sx q[0];
rz(-0.49561758) q[0];
rz(-pi) q[1];
rz(-1.5271565) q[2];
sx q[2];
rz(-0.72481643) q[2];
sx q[2];
rz(-0.61312719) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5962299) q[1];
sx q[1];
rz(-2.0343421) q[1];
sx q[1];
rz(0.35407663) q[1];
rz(-1.3889503) q[3];
sx q[3];
rz(-0.80917984) q[3];
sx q[3];
rz(1.941118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47739425) q[2];
sx q[2];
rz(-1.7919284) q[2];
sx q[2];
rz(-2.9273709) q[2];
rz(-2.2231806) q[3];
sx q[3];
rz(-2.1954506) q[3];
sx q[3];
rz(-1.7414179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4825831) q[0];
sx q[0];
rz(-2.5795586) q[0];
sx q[0];
rz(-0.62579489) q[0];
rz(1.91045) q[1];
sx q[1];
rz(-1.2673255) q[1];
sx q[1];
rz(2.321718) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5441262) q[0];
sx q[0];
rz(-2.1334927) q[0];
sx q[0];
rz(1.1981127) q[0];
rz(-0.71086871) q[2];
sx q[2];
rz(-1.2307302) q[2];
sx q[2];
rz(0.96300551) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7149799) q[1];
sx q[1];
rz(-2.8016114) q[1];
sx q[1];
rz(-1.7991542) q[1];
rz(-0.094535703) q[3];
sx q[3];
rz(-1.5526287) q[3];
sx q[3];
rz(-1.2430735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1890373) q[2];
sx q[2];
rz(-1.2574715) q[2];
sx q[2];
rz(2.7541568) q[2];
rz(-2.6881325) q[3];
sx q[3];
rz(-0.87415868) q[3];
sx q[3];
rz(2.9816154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74359918) q[0];
sx q[0];
rz(-1.1389808) q[0];
sx q[0];
rz(2.0177662) q[0];
rz(-2.3881312) q[1];
sx q[1];
rz(-1.1898142) q[1];
sx q[1];
rz(1.9728647) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982464) q[0];
sx q[0];
rz(-0.69730824) q[0];
sx q[0];
rz(-1.7137609) q[0];
x q[1];
rz(-1.1134603) q[2];
sx q[2];
rz(-1.5654608) q[2];
sx q[2];
rz(2.920546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.67186165) q[1];
sx q[1];
rz(-0.97588343) q[1];
sx q[1];
rz(-2.7487526) q[1];
rz(-2.6791689) q[3];
sx q[3];
rz(-2.1743757) q[3];
sx q[3];
rz(1.0360595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1184825) q[2];
sx q[2];
rz(-1.0175984) q[2];
sx q[2];
rz(-2.6712096) q[2];
rz(-1.4605626) q[3];
sx q[3];
rz(-1.8813671) q[3];
sx q[3];
rz(-2.1605261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7505782) q[0];
sx q[0];
rz(-1.7711201) q[0];
sx q[0];
rz(-1.6455261) q[0];
rz(1.1356614) q[1];
sx q[1];
rz(-1.8541502) q[1];
sx q[1];
rz(-2.6527203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58502349) q[0];
sx q[0];
rz(-2.9980066) q[0];
sx q[0];
rz(-1.344259) q[0];
rz(-pi) q[1];
rz(-0.096285419) q[2];
sx q[2];
rz(-2.0803148) q[2];
sx q[2];
rz(2.4063039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1355537) q[1];
sx q[1];
rz(-3.0873489) q[1];
sx q[1];
rz(0.42363055) q[1];
rz(0.61967697) q[3];
sx q[3];
rz(-1.5159199) q[3];
sx q[3];
rz(2.5375199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6966072) q[2];
sx q[2];
rz(-2.5312436) q[2];
sx q[2];
rz(-1.102591) q[2];
rz(1.6164814) q[3];
sx q[3];
rz(-0.81086719) q[3];
sx q[3];
rz(1.471224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1333604) q[0];
sx q[0];
rz(-2.6120549) q[0];
sx q[0];
rz(1.8089005) q[0];
rz(0.38707271) q[1];
sx q[1];
rz(-1.2553071) q[1];
sx q[1];
rz(1.056384) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0239733) q[0];
sx q[0];
rz(-1.4314382) q[0];
sx q[0];
rz(-1.3246714) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4459287) q[2];
sx q[2];
rz(-1.9859295) q[2];
sx q[2];
rz(-1.268569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3247128) q[1];
sx q[1];
rz(-1.190328) q[1];
sx q[1];
rz(2.3808384) q[1];
x q[2];
rz(2.4293618) q[3];
sx q[3];
rz(-0.81347403) q[3];
sx q[3];
rz(1.939919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.010783823) q[2];
sx q[2];
rz(-2.2937412) q[2];
sx q[2];
rz(-1.3991785) q[2];
rz(1.0731267) q[3];
sx q[3];
rz(-1.9173887) q[3];
sx q[3];
rz(2.6789902) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.241093) q[0];
sx q[0];
rz(-1.6443962) q[0];
sx q[0];
rz(-1.6045438) q[0];
rz(2.4201139) q[1];
sx q[1];
rz(-2.571048) q[1];
sx q[1];
rz(-1.2068988) q[1];
rz(-1.1804562) q[2];
sx q[2];
rz(-1.1710288) q[2];
sx q[2];
rz(-2.2784497) q[2];
rz(0.10151699) q[3];
sx q[3];
rz(-0.77407594) q[3];
sx q[3];
rz(-1.414513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
