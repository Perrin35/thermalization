OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7115241) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(0.67396069) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(2.34692) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90097839) q[0];
sx q[0];
rz(-0.52871791) q[0];
sx q[0];
rz(-0.76529495) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27172471) q[2];
sx q[2];
rz(-2.6368015) q[2];
sx q[2];
rz(1.4128078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1311878) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(2.5351934) q[1];
rz(-pi) q[2];
rz(-1.2875597) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(-1.8739665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6661466) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(0.90707183) q[3];
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
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1871724) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(0.95970884) q[0];
rz(1.6540487) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(0.62746343) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9976881) q[0];
sx q[0];
rz(-1.5630957) q[0];
sx q[0];
rz(-1.5776004) q[0];
x q[1];
rz(1.5487899) q[2];
sx q[2];
rz(-1.7823879) q[2];
sx q[2];
rz(-1.6197268) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6154502) q[1];
sx q[1];
rz(-0.28029385) q[1];
sx q[1];
rz(2.1089094) q[1];
x q[2];
rz(-0.012410951) q[3];
sx q[3];
rz(-1.59613) q[3];
sx q[3];
rz(2.2001571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54005694) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56101218) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130226) q[0];
sx q[0];
rz(-2.0355814) q[0];
sx q[0];
rz(1.3010498) q[0];
rz(-pi) q[1];
rz(-1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(2.3929838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2587535) q[1];
sx q[1];
rz(-1.6707318) q[1];
sx q[1];
rz(-0.12211166) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1191145) q[3];
sx q[3];
rz(-2.936755) q[3];
sx q[3];
rz(2.0126437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0195062) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(1.1536095) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(0.78891689) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-0.21534236) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0523473) q[0];
sx q[0];
rz(-2.2627292) q[0];
sx q[0];
rz(-0.4381051) q[0];
x q[1];
rz(-2.9748561) q[2];
sx q[2];
rz(-2.1519289) q[2];
sx q[2];
rz(-2.2819448) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.035186471) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(0.10722843) q[1];
x q[2];
rz(2.5921949) q[3];
sx q[3];
rz(-1.1664207) q[3];
sx q[3];
rz(-2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(0.96486282) q[0];
rz(-3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(-2.4749277) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.704741) q[0];
sx q[0];
rz(-1.1724823) q[0];
sx q[0];
rz(-0.10956357) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0235734) q[2];
sx q[2];
rz(-1.2663411) q[2];
sx q[2];
rz(0.98642245) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1946698) q[1];
sx q[1];
rz(-1.2877081) q[1];
sx q[1];
rz(-2.4592295) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2090962) q[3];
sx q[3];
rz(-0.36358788) q[3];
sx q[3];
rz(2.5225085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2255286) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(-1.2529681) q[2];
rz(1.8270252) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(1.7793659) q[0];
rz(-1.7419787) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(-1.6113575) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0936733) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(-1.6506667) q[0];
rz(-pi) q[1];
rz(-2.7026664) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(-0.93026464) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1215026) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(0.23328383) q[1];
x q[2];
rz(1.3774032) q[3];
sx q[3];
rz(-2.8838727) q[3];
sx q[3];
rz(-1.6267488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(-1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2832527) q[0];
sx q[0];
rz(-0.28536277) q[0];
sx q[0];
rz(2.7417389) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3235544) q[2];
sx q[2];
rz(-0.99070264) q[2];
sx q[2];
rz(-1.3387418) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6888914) q[1];
sx q[1];
rz(-0.93033067) q[1];
sx q[1];
rz(-1.0673317) q[1];
rz(-0.31886153) q[3];
sx q[3];
rz(-1.0705035) q[3];
sx q[3];
rz(1.5368411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(0.62057173) q[2];
rz(0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(1.0808806) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-2.506315) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95406336) q[0];
sx q[0];
rz(-0.39217338) q[0];
sx q[0];
rz(-1.9847045) q[0];
rz(-pi) q[1];
rz(-2.0681357) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(-1.2129601) q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
x q[2];
rz(0.99526309) q[3];
sx q[3];
rz(-1.7407773) q[3];
sx q[3];
rz(-0.02804027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2962608) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(-0.6566748) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(-2.231853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0168403) q[0];
sx q[0];
rz(-0.42718712) q[0];
sx q[0];
rz(2.4401526) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7934974) q[2];
sx q[2];
rz(-1.3653792) q[2];
sx q[2];
rz(-1.5470488) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.197387) q[1];
sx q[1];
rz(-1.7802317) q[1];
sx q[1];
rz(-1.4447114) q[1];
x q[2];
rz(2.2496201) q[3];
sx q[3];
rz(-2.744305) q[3];
sx q[3];
rz(2.9904423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(-0.89458481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0827732) q[0];
sx q[0];
rz(-1.6130157) q[0];
sx q[0];
rz(2.204338) q[0];
rz(-pi) q[1];
rz(-1.1799699) q[2];
sx q[2];
rz(-2.639689) q[2];
sx q[2];
rz(2.6596136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6539508) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(-1.7483064) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3668725) q[3];
sx q[3];
rz(-0.66019928) q[3];
sx q[3];
rz(-3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(2.1195892) q[2];
rz(1.3379124) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-0.9085761) q[2];
sx q[2];
rz(-2.3618345) q[2];
sx q[2];
rz(0.97710412) q[2];
rz(1.4150423) q[3];
sx q[3];
rz(-0.30434904) q[3];
sx q[3];
rz(0.54868922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
