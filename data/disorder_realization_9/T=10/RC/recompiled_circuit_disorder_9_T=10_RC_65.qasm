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
rz(-0.67396069) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062382467) q[0];
sx q[0];
rz(-1.9431207) q[0];
sx q[0];
rz(-1.9553493) q[0];
rz(-pi) q[1];
rz(-0.27172471) q[2];
sx q[2];
rz(-2.6368015) q[2];
sx q[2];
rz(-1.7287849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.049584576) q[1];
sx q[1];
rz(-1.2283748) q[1];
sx q[1];
rz(-0.53978668) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(-0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(-1.3249935) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871724) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(0.95970884) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(-0.62746343) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7146485) q[0];
sx q[0];
rz(-1.5639925) q[0];
sx q[0];
rz(-0.0077008458) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0395095) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(1.724147) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6154502) q[1];
sx q[1];
rz(-2.8612988) q[1];
sx q[1];
rz(1.0326833) q[1];
x q[2];
rz(0.012410951) q[3];
sx q[3];
rz(-1.59613) q[3];
sx q[3];
rz(-2.2001571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(-2.1616948) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56101218) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(1.7618435) q[0];
rz(2.9648932) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028570024) q[0];
sx q[0];
rz(-2.0355814) q[0];
sx q[0];
rz(1.3010498) q[0];
rz(0.20027192) q[2];
sx q[2];
rz(-1.592604) q[2];
sx q[2];
rz(0.81776103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8417931) q[1];
sx q[1];
rz(-1.6922957) q[1];
sx q[1];
rz(-1.4701162) q[1];
rz(-3.0511608) q[3];
sx q[3];
rz(-1.7548429) q[3];
sx q[3];
rz(-2.4726766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1220864) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(1.9879831) q[2];
rz(2.7518318) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(2.5770082) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-2.9262503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72258654) q[0];
sx q[0];
rz(-0.79918396) q[0];
sx q[0];
rz(1.0976085) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3232857) q[2];
sx q[2];
rz(-0.60192054) q[2];
sx q[2];
rz(1.9844696) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54289651) q[1];
sx q[1];
rz(-0.22194949) q[1];
sx q[1];
rz(-2.0680244) q[1];
rz(-pi) q[2];
rz(-2.0358884) q[3];
sx q[3];
rz(-1.0700738) q[3];
sx q[3];
rz(0.53946686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(0.016618641) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(0.66666493) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.704741) q[0];
sx q[0];
rz(-1.9691103) q[0];
sx q[0];
rz(-0.10956357) q[0];
x q[1];
rz(-1.0235734) q[2];
sx q[2];
rz(-1.8752516) q[2];
sx q[2];
rz(-2.1551702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95496817) q[1];
sx q[1];
rz(-2.4116227) q[1];
sx q[1];
rz(2.7093922) q[1];
x q[2];
rz(-1.2742217) q[3];
sx q[3];
rz(-1.7843102) q[3];
sx q[3];
rz(-0.34542686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91606402) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.2529681) q[2];
rz(1.3145674) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(-1.7793659) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.5302352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8996457) q[0];
sx q[0];
rz(-0.42662222) q[0];
sx q[0];
rz(0.17701478) q[0];
x q[1];
rz(1.2457232) q[2];
sx q[2];
rz(-1.9893861) q[2];
sx q[2];
rz(-2.3649154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1215026) q[1];
sx q[1];
rz(-1.964718) q[1];
sx q[1];
rz(0.23328383) q[1];
rz(-1.8239162) q[3];
sx q[3];
rz(-1.6198006) q[3];
sx q[3];
rz(-0.24310902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.2221378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-2.9587865) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683838) q[0];
sx q[0];
rz(-1.3084992) q[0];
sx q[0];
rz(1.457085) q[0];
rz(-pi) q[1];
rz(0.80656959) q[2];
sx q[2];
rz(-2.2273846) q[2];
sx q[2];
rz(0.29733959) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4527013) q[1];
sx q[1];
rz(-2.211262) q[1];
sx q[1];
rz(-1.0673317) q[1];
x q[2];
rz(1.0501409) q[3];
sx q[3];
rz(-0.58590349) q[3];
sx q[3];
rz(1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(-1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(2.0525232) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-0.63527766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9104722) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(1.2088103) q[0];
x q[1];
rz(-0.17417553) q[2];
sx q[2];
rz(-2.0617699) q[2];
sx q[2];
rz(2.8665286) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0465225) q[1];
sx q[1];
rz(-2.1696643) q[1];
sx q[1];
rz(-1.8938766) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1463296) q[3];
sx q[3];
rz(-1.4008153) q[3];
sx q[3];
rz(-0.02804027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(-0.00017246406) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(0.85787684) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8453318) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(-2.4849179) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(-2.231853) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37659392) q[0];
sx q[0];
rz(-1.8928327) q[0];
sx q[0];
rz(1.8565208) q[0];
rz(-pi) q[1];
rz(-0.54833834) q[2];
sx q[2];
rz(-0.40204918) q[2];
sx q[2];
rz(2.6532432) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64987953) q[1];
sx q[1];
rz(-2.8976106) q[1];
sx q[1];
rz(0.53423832) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8839888) q[3];
sx q[3];
rz(-1.8766878) q[3];
sx q[3];
rz(-0.56759873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(-0.48661423) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(2.6249028) q[0];
rz(-1.5453045) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(-2.2470078) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6870118) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(-1.6420341) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9616227) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(2.6596136) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48764187) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(-1.7483064) q[1];
rz(2.9856332) q[3];
sx q[3];
rz(-2.2150063) q[3];
sx q[3];
rz(0.19176602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14780012) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(-1.0220035) q[2];
rz(-1.3379124) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(0.90850496) q[2];
sx q[2];
rz(-1.1237332) q[2];
sx q[2];
rz(-0.087469812) q[2];
rz(-1.7265504) q[3];
sx q[3];
rz(-0.30434904) q[3];
sx q[3];
rz(0.54868922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];