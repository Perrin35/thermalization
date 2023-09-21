OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.43006858) q[0];
sx q[0];
rz(-3.0741337) q[0];
sx q[0];
rz(2.467632) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2406143) q[0];
sx q[0];
rz(-0.52871791) q[0];
sx q[0];
rz(2.3762977) q[0];
rz(-pi) q[1];
rz(-1.4235714) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(-1.1046315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4218581) q[1];
sx q[1];
rz(-2.0761479) q[1];
sx q[1];
rz(1.964633) q[1];
x q[2];
rz(0.13866339) q[3];
sx q[3];
rz(-1.1274459) q[3];
sx q[3];
rz(1.5821622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6661466) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(-2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1871724) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(-2.1818838) q[0];
rz(1.6540487) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(2.5141292) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42022959) q[0];
sx q[0];
rz(-0.010275928) q[0];
sx q[0];
rz(-0.72364877) q[0];
rz(0.21164125) q[2];
sx q[2];
rz(-1.5492808) q[2];
sx q[2];
rz(0.053552901) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5261425) q[1];
sx q[1];
rz(-0.28029385) q[1];
sx q[1];
rz(-2.1089094) q[1];
rz(-pi) q[2];
rz(-3.1291817) q[3];
sx q[3];
rz(-1.59613) q[3];
sx q[3];
rz(0.94143553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6015357) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(-2.1616948) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.56101218) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(1.3797492) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.7938991) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4189258) q[0];
sx q[0];
rz(-1.8113266) q[0];
sx q[0];
rz(-2.6618883) q[0];
rz(-1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(2.3929838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8828391) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(3.019481) q[1];
rz(3.0511608) q[3];
sx q[3];
rz(-1.3867497) q[3];
sx q[3];
rz(0.66891608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0195062) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(2.7518318) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6692114) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(0.56458449) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(0.21534236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9504844) q[0];
sx q[0];
rz(-1.9035625) q[0];
sx q[0];
rz(0.82975181) q[0];
rz(-pi) q[1];
rz(1.3232857) q[2];
sx q[2];
rz(-2.5396721) q[2];
sx q[2];
rz(1.1571231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54289651) q[1];
sx q[1];
rz(-0.22194949) q[1];
sx q[1];
rz(-2.0680244) q[1];
rz(0.54939778) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(-2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85201207) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(0.96486282) q[0];
rz(0.016618641) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(-2.4749277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0502888) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(1.9712649) q[0];
x q[1];
rz(1.0275074) q[2];
sx q[2];
rz(-2.5230061) q[2];
sx q[2];
rz(1.0416043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.40068377) q[1];
sx q[1];
rz(-2.2212257) q[1];
sx q[1];
rz(1.2121735) q[1];
x q[2];
rz(-0.93249647) q[3];
sx q[3];
rz(-2.7780048) q[3];
sx q[3];
rz(0.61908412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91606402) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(1.8886245) q[2];
rz(-1.8270252) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1154293) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(1.7793659) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.5302352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513072) q[0];
sx q[0];
rz(-1.4978652) q[0];
sx q[0];
rz(-0.42071995) q[0];
rz(-0.43892626) q[2];
sx q[2];
rz(-1.274684) q[2];
sx q[2];
rz(2.211328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.02009) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(-2.9083088) q[1];
rz(-pi) q[2];
rz(1.3774032) q[3];
sx q[3];
rz(-0.25771991) q[3];
sx q[3];
rz(-1.5148439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(-2.706066) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85833997) q[0];
sx q[0];
rz(-2.8562299) q[0];
sx q[0];
rz(-2.7417389) q[0];
rz(-pi) q[1];
rz(0.73166087) q[2];
sx q[2];
rz(-2.1795142) q[2];
sx q[2];
rz(-0.70639709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4527013) q[1];
sx q[1];
rz(-2.211262) q[1];
sx q[1];
rz(1.0673317) q[1];
x q[2];
rz(-1.0484344) q[3];
sx q[3];
rz(-1.8494542) q[3];
sx q[3];
rz(-0.1230965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(2.506315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104722) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(1.9327823) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8842472) q[2];
sx q[2];
rz(-0.51856504) q[2];
sx q[2];
rz(-2.5093362) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.852408) q[1];
sx q[1];
rz(-1.305456) q[1];
sx q[1];
rz(-2.5177588) q[1];
x q[2];
rz(-2.9397804) q[3];
sx q[3];
rz(-2.1370071) q[3];
sx q[3];
rz(-1.7081529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0044331) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(-0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(-0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(-0.6566748) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(2.231853) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7649987) q[0];
sx q[0];
rz(-1.8928327) q[0];
sx q[0];
rz(-1.2850719) q[0];
rz(1.3526731) q[2];
sx q[2];
rz(-1.911273) q[2];
sx q[2];
rz(-3.043963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4917131) q[1];
sx q[1];
rz(-2.8976106) q[1];
sx q[1];
rz(2.6073543) q[1];
rz(-1.886456) q[3];
sx q[3];
rz(-1.3254032) q[3];
sx q[3];
rz(2.0592225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83071128) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(3.0330372) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-0.51668984) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(-0.89458481) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0827732) q[0];
sx q[0];
rz(-1.528577) q[0];
sx q[0];
rz(-0.93725462) q[0];
rz(-pi) q[1];
rz(1.9616227) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(0.48197907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.48764187) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(-1.3932863) q[1];
rz(-pi) q[2];
rz(-2.22088) q[3];
sx q[3];
rz(-1.4462785) q[3];
sx q[3];
rz(1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(1.0220035) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.548303) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(1.0409566) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(2.2330877) q[2];
sx q[2];
rz(-2.0178595) q[2];
sx q[2];
rz(3.0541228) q[2];
rz(3.0929052) q[3];
sx q[3];
rz(-1.8713453) q[3];
sx q[3];
rz(0.71181675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];