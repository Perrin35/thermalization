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
rz(3.3764163) q[0];
sx q[0];
rz(11.233815) q[0];
rz(1.3136343) q[1];
sx q[1];
rz(-1.87275) q[1];
sx q[1];
rz(0.76729265) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7504265) q[0];
sx q[0];
rz(-0.39505233) q[0];
sx q[0];
rz(0.0042628757) q[0];
rz(-1.9903153) q[2];
sx q[2];
rz(-1.1184208) q[2];
sx q[2];
rz(0.6587761) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1033039) q[1];
sx q[1];
rz(-1.3974408) q[1];
sx q[1];
rz(1.2452673) q[1];
rz(-pi) q[2];
rz(2.4657605) q[3];
sx q[3];
rz(-2.5321372) q[3];
sx q[3];
rz(-2.6026562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2904498) q[2];
sx q[2];
rz(-2.3252611) q[2];
sx q[2];
rz(-1.8905224) q[2];
rz(-2.3276681) q[3];
sx q[3];
rz(-2.0415174) q[3];
sx q[3];
rz(-1.6479015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.587065) q[0];
sx q[0];
rz(-1.4545414) q[0];
sx q[0];
rz(-2.5982507) q[0];
rz(0.73703611) q[1];
sx q[1];
rz(-0.68865028) q[1];
sx q[1];
rz(0.063311689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2800326) q[0];
sx q[0];
rz(-0.82516042) q[0];
sx q[0];
rz(-2.3579423) q[0];
rz(-0.28191992) q[2];
sx q[2];
rz(-1.3358982) q[2];
sx q[2];
rz(0.24215936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0009389) q[1];
sx q[1];
rz(-1.2564827) q[1];
sx q[1];
rz(0.94625603) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7821941) q[3];
sx q[3];
rz(-1.4305887) q[3];
sx q[3];
rz(2.0003776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43725499) q[2];
sx q[2];
rz(-1.1235808) q[2];
sx q[2];
rz(0.91427461) q[2];
rz(-2.8546913) q[3];
sx q[3];
rz(-0.56921452) q[3];
sx q[3];
rz(-2.5347819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(0.86976784) q[1];
sx q[1];
rz(-0.52647796) q[1];
sx q[1];
rz(2.8960752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29123339) q[0];
sx q[0];
rz(-1.1867674) q[0];
sx q[0];
rz(-1.2880833) q[0];
rz(-3.1046529) q[2];
sx q[2];
rz(-2.1203342) q[2];
sx q[2];
rz(-2.93826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88855381) q[1];
sx q[1];
rz(-1.843707) q[1];
sx q[1];
rz(-0.43123475) q[1];
rz(-1.5762899) q[3];
sx q[3];
rz(-0.76346654) q[3];
sx q[3];
rz(-0.31007773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3123582) q[2];
sx q[2];
rz(-1.979253) q[2];
sx q[2];
rz(-1.8355628) q[2];
rz(-0.28451377) q[3];
sx q[3];
rz(-1.6407216) q[3];
sx q[3];
rz(-1.3792926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5259842) q[0];
sx q[0];
rz(-2.5475976) q[0];
sx q[0];
rz(2.3748412) q[0];
rz(1.426567) q[1];
sx q[1];
rz(-1.7568935) q[1];
sx q[1];
rz(-1.0624622) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7257197) q[0];
sx q[0];
rz(-1.555528) q[0];
sx q[0];
rz(1.589561) q[0];
x q[1];
rz(1.771174) q[2];
sx q[2];
rz(-0.60094072) q[2];
sx q[2];
rz(-0.19285017) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4442788) q[1];
sx q[1];
rz(-1.0890038) q[1];
sx q[1];
rz(-0.30028157) q[1];
rz(-0.39233971) q[3];
sx q[3];
rz(-2.9846387) q[3];
sx q[3];
rz(-0.63877393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95395025) q[2];
sx q[2];
rz(-0.01454777) q[2];
sx q[2];
rz(-0.1178096) q[2];
rz(-0.5438965) q[3];
sx q[3];
rz(-2.1102648) q[3];
sx q[3];
rz(-0.73345524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2619005) q[0];
sx q[0];
rz(-1.4786792) q[0];
sx q[0];
rz(2.6035768) q[0];
rz(-3.0777439) q[1];
sx q[1];
rz(-2.1203181) q[1];
sx q[1];
rz(1.8003546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1088561) q[0];
sx q[0];
rz(-1.0633611) q[0];
sx q[0];
rz(-2.5030337) q[0];
rz(-pi) q[1];
rz(2.5885905) q[2];
sx q[2];
rz(-1.2931839) q[2];
sx q[2];
rz(1.0652519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3732934) q[1];
sx q[1];
rz(-1.5099257) q[1];
sx q[1];
rz(2.5817691) q[1];
x q[2];
rz(-2.7085764) q[3];
sx q[3];
rz(-1.6891494) q[3];
sx q[3];
rz(-2.1089212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6824048) q[2];
sx q[2];
rz(-2.8509199) q[2];
sx q[2];
rz(0.9064557) q[2];
rz(-1.2022379) q[3];
sx q[3];
rz(-1.5463444) q[3];
sx q[3];
rz(1.7513587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15659139) q[0];
sx q[0];
rz(-2.5522794) q[0];
sx q[0];
rz(-0.44664788) q[0];
rz(-0.43529549) q[1];
sx q[1];
rz(-2.5050102) q[1];
sx q[1];
rz(2.2925503) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0627611) q[0];
sx q[0];
rz(-2.9568892) q[0];
sx q[0];
rz(1.1514501) q[0];
x q[1];
rz(-0.53211114) q[2];
sx q[2];
rz(-1.4363043) q[2];
sx q[2];
rz(-2.7464339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.067612996) q[1];
sx q[1];
rz(-2.9146806) q[1];
sx q[1];
rz(-3.0935442) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44621991) q[3];
sx q[3];
rz(-0.38360197) q[3];
sx q[3];
rz(-1.6723796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0173831) q[2];
sx q[2];
rz(-2.1653192) q[2];
sx q[2];
rz(-1.7912095) q[2];
rz(0.24556686) q[3];
sx q[3];
rz(-2.1731302) q[3];
sx q[3];
rz(0.33026162) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46828142) q[0];
sx q[0];
rz(-0.89235965) q[0];
sx q[0];
rz(2.6924676) q[0];
rz(-1.7091854) q[1];
sx q[1];
rz(-2.1604373) q[1];
sx q[1];
rz(1.2264576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5275636) q[0];
sx q[0];
rz(-1.4948658) q[0];
sx q[0];
rz(2.1873104) q[0];
rz(-pi) q[1];
rz(-0.62896987) q[2];
sx q[2];
rz(-2.5444053) q[2];
sx q[2];
rz(0.33729182) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0192208) q[1];
sx q[1];
rz(-1.6748412) q[1];
sx q[1];
rz(-1.6582483) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.874648) q[3];
sx q[3];
rz(-1.9404074) q[3];
sx q[3];
rz(0.0064484869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8760828) q[2];
sx q[2];
rz(-2.2777057) q[2];
sx q[2];
rz(-2.6465805) q[2];
rz(1.0692976) q[3];
sx q[3];
rz(-2.4102231) q[3];
sx q[3];
rz(0.67720145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7345562) q[0];
sx q[0];
rz(-0.41663909) q[0];
sx q[0];
rz(-0.46105841) q[0];
rz(-3.0126493) q[1];
sx q[1];
rz(-2.415633) q[1];
sx q[1];
rz(-0.9022378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092493499) q[0];
sx q[0];
rz(-1.6120132) q[0];
sx q[0];
rz(-1.5638086) q[0];
x q[1];
rz(-1.721549) q[2];
sx q[2];
rz(-0.37166212) q[2];
sx q[2];
rz(-1.1451385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.56574351) q[1];
sx q[1];
rz(-1.1593124) q[1];
sx q[1];
rz(-0.38994148) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.227026) q[3];
sx q[3];
rz(-1.4952139) q[3];
sx q[3];
rz(1.9618159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43877131) q[2];
sx q[2];
rz(-2.2192025) q[2];
sx q[2];
rz(-1.2141466) q[2];
rz(2.7625648) q[3];
sx q[3];
rz(-1.4641848) q[3];
sx q[3];
rz(1.3488784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9146861) q[0];
sx q[0];
rz(-1.5808957) q[0];
sx q[0];
rz(-3.1353986) q[0];
rz(-2.6365623) q[1];
sx q[1];
rz(-0.49022245) q[1];
sx q[1];
rz(-0.44713155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4775708) q[0];
sx q[0];
rz(-0.89671269) q[0];
sx q[0];
rz(-0.88468219) q[0];
rz(-pi) q[1];
rz(1.4535663) q[2];
sx q[2];
rz(-2.0868868) q[2];
sx q[2];
rz(-1.0265526) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3382692) q[1];
sx q[1];
rz(-1.7223961) q[1];
sx q[1];
rz(3.1123509) q[1];
rz(-pi) q[2];
rz(0.38180967) q[3];
sx q[3];
rz(-2.6630132) q[3];
sx q[3];
rz(-0.081153966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7827683) q[2];
sx q[2];
rz(-1.1268104) q[2];
sx q[2];
rz(-2.6940441) q[2];
rz(-2.0284082) q[3];
sx q[3];
rz(-1.0380849) q[3];
sx q[3];
rz(2.8158877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96414763) q[0];
sx q[0];
rz(-0.064346813) q[0];
sx q[0];
rz(-0.84661761) q[0];
rz(-2.8586491) q[1];
sx q[1];
rz(-1.6694262) q[1];
sx q[1];
rz(0.67323321) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1812912) q[0];
sx q[0];
rz(-3.1226469) q[0];
sx q[0];
rz(-2.6238954) q[0];
rz(-pi) q[1];
rz(-1.3429672) q[2];
sx q[2];
rz(-0.81124094) q[2];
sx q[2];
rz(1.0885061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1800419) q[1];
sx q[1];
rz(-2.2030239) q[1];
sx q[1];
rz(-2.8668731) q[1];
rz(-pi) q[2];
rz(0.29127647) q[3];
sx q[3];
rz(-1.8921953) q[3];
sx q[3];
rz(2.0673366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8474569) q[2];
sx q[2];
rz(-2.2723276) q[2];
sx q[2];
rz(2.4693303) q[2];
rz(1.8654035) q[3];
sx q[3];
rz(-1.0381235) q[3];
sx q[3];
rz(-0.41260317) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1558253) q[0];
sx q[0];
rz(-1.1341996) q[0];
sx q[0];
rz(-1.2962935) q[0];
rz(3.1271707) q[1];
sx q[1];
rz(-1.2851234) q[1];
sx q[1];
rz(-1.8869225) q[1];
rz(0.6966656) q[2];
sx q[2];
rz(-1.2045384) q[2];
sx q[2];
rz(-0.62919553) q[2];
rz(-1.779535) q[3];
sx q[3];
rz(-1.7961647) q[3];
sx q[3];
rz(2.4299783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
