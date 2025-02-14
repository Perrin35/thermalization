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
rz(1.9795228) q[0];
sx q[0];
rz(-2.906769) q[0];
sx q[0];
rz(1.809037) q[0];
rz(1.3136343) q[1];
sx q[1];
rz(-1.87275) q[1];
sx q[1];
rz(0.76729265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.965897) q[0];
sx q[0];
rz(-1.5691557) q[0];
sx q[0];
rz(-2.7465435) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1512773) q[2];
sx q[2];
rz(-1.1184208) q[2];
sx q[2];
rz(-2.4828165) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.060081944) q[1];
sx q[1];
rz(-2.7742371) q[1];
sx q[1];
rz(-2.0717595) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6428623) q[3];
sx q[3];
rz(-1.2045897) q[3];
sx q[3];
rz(-1.5281475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8511429) q[2];
sx q[2];
rz(-2.3252611) q[2];
sx q[2];
rz(-1.8905224) q[2];
rz(-0.81392455) q[3];
sx q[3];
rz(-1.1000752) q[3];
sx q[3];
rz(1.4936911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.587065) q[0];
sx q[0];
rz(-1.4545414) q[0];
sx q[0];
rz(-2.5982507) q[0];
rz(-2.4045565) q[1];
sx q[1];
rz(-0.68865028) q[1];
sx q[1];
rz(-3.078281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2800326) q[0];
sx q[0];
rz(-0.82516042) q[0];
sx q[0];
rz(0.78365032) q[0];
rz(-pi) q[1];
x q[1];
rz(2.43119) q[2];
sx q[2];
rz(-2.7766529) q[2];
sx q[2];
rz(2.0055564) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1166414) q[1];
sx q[1];
rz(-2.4519741) q[1];
sx q[1];
rz(-1.0633797) q[1];
rz(1.7821941) q[3];
sx q[3];
rz(-1.4305887) q[3];
sx q[3];
rz(-1.141215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43725499) q[2];
sx q[2];
rz(-2.0180118) q[2];
sx q[2];
rz(2.227318) q[2];
rz(0.28690139) q[3];
sx q[3];
rz(-2.5723781) q[3];
sx q[3];
rz(-0.60681075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8617926) q[0];
sx q[0];
rz(-1.8656116) q[0];
sx q[0];
rz(-1.7472501) q[0];
rz(-0.86976784) q[1];
sx q[1];
rz(-2.6151147) q[1];
sx q[1];
rz(2.8960752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9704392) q[0];
sx q[0];
rz(-1.309179) q[0];
sx q[0];
rz(-0.39830504) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0209544) q[2];
sx q[2];
rz(-1.5392973) q[2];
sx q[2];
rz(-1.3481639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.989262) q[1];
sx q[1];
rz(-0.50571364) q[1];
sx q[1];
rz(-0.590041) q[1];
rz(-3.136335) q[3];
sx q[3];
rz(-0.80734423) q[3];
sx q[3];
rz(-0.30247363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3123582) q[2];
sx q[2];
rz(-1.1623397) q[2];
sx q[2];
rz(-1.3060298) q[2];
rz(2.8570789) q[3];
sx q[3];
rz(-1.5008711) q[3];
sx q[3];
rz(1.3792926) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156085) q[0];
sx q[0];
rz(-2.5475976) q[0];
sx q[0];
rz(-2.3748412) q[0];
rz(1.426567) q[1];
sx q[1];
rz(-1.7568935) q[1];
sx q[1];
rz(-1.0624622) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9869558) q[0];
sx q[0];
rz(-1.5895588) q[0];
sx q[0];
rz(-0.015270998) q[0];
rz(-pi) q[1];
rz(-1.3704186) q[2];
sx q[2];
rz(-0.60094072) q[2];
sx q[2];
rz(2.9487425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4442788) q[1];
sx q[1];
rz(-2.0525888) q[1];
sx q[1];
rz(2.8413111) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7492529) q[3];
sx q[3];
rz(-2.9846387) q[3];
sx q[3];
rz(0.63877393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1876424) q[2];
sx q[2];
rz(-3.1270449) q[2];
sx q[2];
rz(0.1178096) q[2];
rz(0.5438965) q[3];
sx q[3];
rz(-2.1102648) q[3];
sx q[3];
rz(-2.4081374) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8796922) q[0];
sx q[0];
rz(-1.4786792) q[0];
sx q[0];
rz(-0.53801584) q[0];
rz(-3.0777439) q[1];
sx q[1];
rz(-1.0212746) q[1];
sx q[1];
rz(1.341238) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1088561) q[0];
sx q[0];
rz(-1.0633611) q[0];
sx q[0];
rz(-0.638559) q[0];
x q[1];
rz(-1.2476498) q[2];
sx q[2];
rz(-2.1003336) q[2];
sx q[2];
rz(-2.8036237) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3772014) q[1];
sx q[1];
rz(-2.1294597) q[1];
sx q[1];
rz(-1.4989946) q[1];
rz(-pi) q[2];
rz(2.8654534) q[3];
sx q[3];
rz(-0.44790972) q[3];
sx q[3];
rz(-0.78820273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4591879) q[2];
sx q[2];
rz(-0.29067278) q[2];
sx q[2];
rz(2.235137) q[2];
rz(-1.2022379) q[3];
sx q[3];
rz(-1.5463444) q[3];
sx q[3];
rz(1.7513587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15659139) q[0];
sx q[0];
rz(-2.5522794) q[0];
sx q[0];
rz(0.44664788) q[0];
rz(-2.7062972) q[1];
sx q[1];
rz(-2.5050102) q[1];
sx q[1];
rz(0.84904233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9049677) q[0];
sx q[0];
rz(-1.6456438) q[0];
sx q[0];
rz(1.4017795) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4150494) q[2];
sx q[2];
rz(-1.0439936) q[2];
sx q[2];
rz(-2.0447363) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0739797) q[1];
sx q[1];
rz(-0.22691209) q[1];
sx q[1];
rz(-0.048048419) q[1];
rz(-pi) q[2];
rz(-1.3983512) q[3];
sx q[3];
rz(-1.9151805) q[3];
sx q[3];
rz(-1.9454959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1242096) q[2];
sx q[2];
rz(-2.1653192) q[2];
sx q[2];
rz(-1.3503831) q[2];
rz(-0.24556686) q[3];
sx q[3];
rz(-0.96846247) q[3];
sx q[3];
rz(0.33026162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46828142) q[0];
sx q[0];
rz(-0.89235965) q[0];
sx q[0];
rz(0.44912502) q[0];
rz(-1.4324073) q[1];
sx q[1];
rz(-0.9811554) q[1];
sx q[1];
rz(1.2264576) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5275636) q[0];
sx q[0];
rz(-1.6467268) q[0];
sx q[0];
rz(2.1873104) q[0];
x q[1];
rz(-0.50275393) q[2];
sx q[2];
rz(-1.9079676) q[2];
sx q[2];
rz(-1.7751116) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.12237182) q[1];
sx q[1];
rz(-1.6748412) q[1];
sx q[1];
rz(1.6582483) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.874648) q[3];
sx q[3];
rz(-1.2011853) q[3];
sx q[3];
rz(3.1351442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2655098) q[2];
sx q[2];
rz(-2.2777057) q[2];
sx q[2];
rz(2.6465805) q[2];
rz(2.072295) q[3];
sx q[3];
rz(-0.73136955) q[3];
sx q[3];
rz(0.67720145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4070364) q[0];
sx q[0];
rz(-0.41663909) q[0];
sx q[0];
rz(-0.46105841) q[0];
rz(3.0126493) q[1];
sx q[1];
rz(-0.72595969) q[1];
sx q[1];
rz(2.2393548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26048294) q[0];
sx q[0];
rz(-0.04180464) q[0];
sx q[0];
rz(-2.9737472) q[0];
x q[1];
rz(3.0831218) q[2];
sx q[2];
rz(-1.2035511) q[2];
sx q[2];
rz(-2.1580687) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5758491) q[1];
sx q[1];
rz(-1.1593124) q[1];
sx q[1];
rz(-2.7516512) q[1];
x q[2];
rz(-1.447313) q[3];
sx q[3];
rz(-0.65992814) q[3];
sx q[3];
rz(2.8483158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43877131) q[2];
sx q[2];
rz(-2.2192025) q[2];
sx q[2];
rz(-1.2141466) q[2];
rz(-2.7625648) q[3];
sx q[3];
rz(-1.6774079) q[3];
sx q[3];
rz(1.3488784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9146861) q[0];
sx q[0];
rz(-1.560697) q[0];
sx q[0];
rz(3.1353986) q[0];
rz(2.6365623) q[1];
sx q[1];
rz(-2.6513702) q[1];
sx q[1];
rz(2.6944611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762334) q[0];
sx q[0];
rz(-2.0885944) q[0];
sx q[0];
rz(-0.8014265) q[0];
rz(-pi) q[1];
rz(0.51905175) q[2];
sx q[2];
rz(-1.468892) q[2];
sx q[2];
rz(-0.48619147) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9046482) q[1];
sx q[1];
rz(-1.5997026) q[1];
sx q[1];
rz(-1.7224599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4486964) q[3];
sx q[3];
rz(-1.7432392) q[3];
sx q[3];
rz(1.8320097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7827683) q[2];
sx q[2];
rz(-2.0147822) q[2];
sx q[2];
rz(2.6940441) q[2];
rz(1.1131845) q[3];
sx q[3];
rz(-1.0380849) q[3];
sx q[3];
rz(-0.32570496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96414763) q[0];
sx q[0];
rz(-0.064346813) q[0];
sx q[0];
rz(0.84661761) q[0];
rz(2.8586491) q[1];
sx q[1];
rz(-1.6694262) q[1];
sx q[1];
rz(2.4683594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2344674) q[0];
sx q[0];
rz(-1.5801718) q[0];
sx q[0];
rz(-3.1251291) q[0];
rz(-pi) q[1];
x q[1];
rz(2.368957) q[2];
sx q[2];
rz(-1.7353206) q[2];
sx q[2];
rz(-2.5009837) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26499012) q[1];
sx q[1];
rz(-0.68176378) q[1];
sx q[1];
rz(1.216128) q[1];
rz(-0.29127647) q[3];
sx q[3];
rz(-1.2493973) q[3];
sx q[3];
rz(2.0673366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29413572) q[2];
sx q[2];
rz(-0.86926502) q[2];
sx q[2];
rz(2.4693303) q[2];
rz(1.8654035) q[3];
sx q[3];
rz(-1.0381235) q[3];
sx q[3];
rz(2.7289895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1558253) q[0];
sx q[0];
rz(-1.1341996) q[0];
sx q[0];
rz(-1.2962935) q[0];
rz(-3.1271707) q[1];
sx q[1];
rz(-1.8564693) q[1];
sx q[1];
rz(1.2546702) q[1];
rz(-1.1070743) q[2];
sx q[2];
rz(-2.2131789) q[2];
sx q[2];
rz(0.65050355) q[2];
rz(1.779535) q[3];
sx q[3];
rz(-1.345428) q[3];
sx q[3];
rz(-0.71161436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
