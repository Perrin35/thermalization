OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8310175) q[0];
sx q[0];
rz(-2.0895045) q[0];
sx q[0];
rz(-1.6488099) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(2.066943) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5045835) q[0];
sx q[0];
rz(-1.5890117) q[0];
sx q[0];
rz(-1.5760742) q[0];
rz(1.4104112) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(-0.10057848) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1638012) q[1];
sx q[1];
rz(-1.4492387) q[1];
sx q[1];
rz(1.6519283) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1172469) q[3];
sx q[3];
rz(-2.0587066) q[3];
sx q[3];
rz(0.32353668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.6281698) q[2];
rz(-1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(-1.9182385) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.6404023) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6340856) q[0];
sx q[0];
rz(-1.5459783) q[0];
sx q[0];
rz(3.0188942) q[0];
rz(-pi) q[1];
rz(0.31200774) q[2];
sx q[2];
rz(-2.0199676) q[2];
sx q[2];
rz(2.1829407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0928287) q[1];
sx q[1];
rz(-1.7424912) q[1];
sx q[1];
rz(-2.462639) q[1];
rz(-1.5562015) q[3];
sx q[3];
rz(-0.39108927) q[3];
sx q[3];
rz(2.1434243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(1.6274874) q[2];
rz(-1.1668011) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353772) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.6378145) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5119748) q[0];
sx q[0];
rz(-3.0324728) q[0];
sx q[0];
rz(0.93018053) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2321646) q[2];
sx q[2];
rz(-2.4436908) q[2];
sx q[2];
rz(0.94240377) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1551731) q[1];
sx q[1];
rz(-1.5281614) q[1];
sx q[1];
rz(-0.30100545) q[1];
x q[2];
rz(-3.010473) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(0.91344792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(2.3977996) q[2];
rz(0.25742325) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-2.0898576) q[0];
rz(-0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(-1.1490885) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968397) q[0];
sx q[0];
rz(-1.4595932) q[0];
sx q[0];
rz(2.8993594) q[0];
rz(-2.1707702) q[2];
sx q[2];
rz(-2.4568825) q[2];
sx q[2];
rz(-0.1240571) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.035605343) q[1];
sx q[1];
rz(-0.78250256) q[1];
sx q[1];
rz(-1.1483907) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2265008) q[3];
sx q[3];
rz(-0.85840423) q[3];
sx q[3];
rz(0.32574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(-0.66876283) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-0.016383735) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(-2.2573684) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852172) q[0];
sx q[0];
rz(-1.7612805) q[0];
sx q[0];
rz(0.18116829) q[0];
rz(-0.82489478) q[2];
sx q[2];
rz(-0.56606228) q[2];
sx q[2];
rz(-0.60264665) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9235232) q[1];
sx q[1];
rz(-1.7150208) q[1];
sx q[1];
rz(-1.5037392) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0377117) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(-2.0562293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23652442) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(2.0489676) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(1.3177692) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(0.96907369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9861525) q[0];
sx q[0];
rz(-1.4618205) q[0];
sx q[0];
rz(-0.040779671) q[0];
x q[1];
rz(3.063383) q[2];
sx q[2];
rz(-0.71560301) q[2];
sx q[2];
rz(-0.90925928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5507257) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(-0.73901891) q[1];
rz(-pi) q[2];
x q[2];
rz(3*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(-1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(-1.3597885) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.70599) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.766151) q[0];
sx q[0];
rz(-2.0631844) q[0];
sx q[0];
rz(-2.173461) q[0];
rz(-pi) q[1];
rz(-3.1232386) q[2];
sx q[2];
rz(-1.5905426) q[2];
sx q[2];
rz(2.8543618) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.095852921) q[1];
sx q[1];
rz(-0.47912712) q[1];
sx q[1];
rz(-0.87289255) q[1];
x q[2];
rz(2.0766047) q[3];
sx q[3];
rz(-1.9516264) q[3];
sx q[3];
rz(1.8246458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(0.15667668) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.3652323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5196461) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(-1.4343065) q[0];
rz(-3.119602) q[2];
sx q[2];
rz(-1.6985661) q[2];
sx q[2];
rz(2.2760504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8147162) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(-0.53175064) q[1];
rz(2.7435999) q[3];
sx q[3];
rz(-0.45106217) q[3];
sx q[3];
rz(0.034702452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.075295538) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(2.1832809) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(2.9039834) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(2.8839135) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(2.2176567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7565219) q[0];
sx q[0];
rz(-2.0002881) q[0];
sx q[0];
rz(-0.71682741) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.02248259) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(-1.2475916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82211557) q[1];
sx q[1];
rz(-1.176782) q[1];
sx q[1];
rz(1.0072717) q[1];
x q[2];
rz(-0.68848916) q[3];
sx q[3];
rz(-2.2985035) q[3];
sx q[3];
rz(1.1779932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6442287) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(-1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.37344638) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-3.1034234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.569012) q[0];
sx q[0];
rz(-0.60927143) q[0];
sx q[0];
rz(1.5241745) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9534692) q[2];
sx q[2];
rz(-1.4854991) q[2];
sx q[2];
rz(1.485699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8074179) q[1];
sx q[1];
rz(-1.6125229) q[1];
sx q[1];
rz(-1.9445022) q[1];
rz(-0.84382236) q[3];
sx q[3];
rz(-0.80599621) q[3];
sx q[3];
rz(0.7741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.96597153) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(-1.9434628) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(-1.6434796) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-2.6565068) q[2];
sx q[2];
rz(-3.017364) q[2];
sx q[2];
rz(-0.45807522) q[2];
rz(2.8244143) q[3];
sx q[3];
rz(-0.89430292) q[3];
sx q[3];
rz(3.111307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
