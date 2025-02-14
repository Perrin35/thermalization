OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7118536) q[0];
sx q[0];
rz(5.138152) q[0];
sx q[0];
rz(10.44147) q[0];
rz(-0.84438762) q[1];
sx q[1];
rz(-0.8165741) q[1];
sx q[1];
rz(0.56632298) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1725572) q[0];
sx q[0];
rz(-2.175967) q[0];
sx q[0];
rz(-3.023554) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3728605) q[2];
sx q[2];
rz(-1.6467484) q[2];
sx q[2];
rz(1.6449442) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.50496347) q[1];
sx q[1];
rz(-2.7111021) q[1];
sx q[1];
rz(2.9477547) q[1];
x q[2];
rz(1.3641961) q[3];
sx q[3];
rz(-2.1461764) q[3];
sx q[3];
rz(2.8386338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6301959) q[2];
sx q[2];
rz(-0.87367311) q[2];
sx q[2];
rz(2.013618) q[2];
rz(1.6389716) q[3];
sx q[3];
rz(-0.84533826) q[3];
sx q[3];
rz(0.42475167) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2750435) q[0];
sx q[0];
rz(-0.4568704) q[0];
sx q[0];
rz(3.1006815) q[0];
rz(2.5919137) q[1];
sx q[1];
rz(-2.7626541) q[1];
sx q[1];
rz(-0.080370195) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98438063) q[0];
sx q[0];
rz(-1.5750021) q[0];
sx q[0];
rz(-2.4192817) q[0];
rz(-pi) q[1];
rz(-0.92795579) q[2];
sx q[2];
rz(-2.0242226) q[2];
sx q[2];
rz(0.7446028) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7521022) q[1];
sx q[1];
rz(-0.48797777) q[1];
sx q[1];
rz(1.6882903) q[1];
x q[2];
rz(-2.8714058) q[3];
sx q[3];
rz(-0.98434292) q[3];
sx q[3];
rz(-1.0932066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72767672) q[2];
sx q[2];
rz(-2.2985986) q[2];
sx q[2];
rz(2.8576039) q[2];
rz(1.6207638) q[3];
sx q[3];
rz(-1.5903383) q[3];
sx q[3];
rz(-0.76333299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82336998) q[0];
sx q[0];
rz(-2.9871873) q[0];
sx q[0];
rz(0.41686091) q[0];
rz(0.93337026) q[1];
sx q[1];
rz(-2.2552172) q[1];
sx q[1];
rz(-2.0534168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8899208) q[0];
sx q[0];
rz(-1.8150738) q[0];
sx q[0];
rz(-2.1273462) q[0];
x q[1];
rz(0.56060411) q[2];
sx q[2];
rz(-1.6050771) q[2];
sx q[2];
rz(-0.83002582) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58055604) q[1];
sx q[1];
rz(-1.2966178) q[1];
sx q[1];
rz(-2.8304173) q[1];
rz(-pi) q[2];
rz(-1.5059286) q[3];
sx q[3];
rz(-1.2838351) q[3];
sx q[3];
rz(2.502041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1209968) q[2];
sx q[2];
rz(-0.72046295) q[2];
sx q[2];
rz(-2.8159115) q[2];
rz(-2.7459512) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(2.4237848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91738236) q[0];
sx q[0];
rz(-0.93467394) q[0];
sx q[0];
rz(3.000946) q[0];
rz(-0.39528254) q[1];
sx q[1];
rz(-1.544516) q[1];
sx q[1];
rz(-2.7093754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2435139) q[0];
sx q[0];
rz(-1.167341) q[0];
sx q[0];
rz(1.0591169) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.08608685) q[2];
sx q[2];
rz(-2.0755092) q[2];
sx q[2];
rz(1.9534257) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85104698) q[1];
sx q[1];
rz(-1.49069) q[1];
sx q[1];
rz(0.27388957) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3234576) q[3];
sx q[3];
rz(-0.77518089) q[3];
sx q[3];
rz(1.1213746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29493368) q[2];
sx q[2];
rz(-2.0294971) q[2];
sx q[2];
rz(-0.46889949) q[2];
rz(-0.15639671) q[3];
sx q[3];
rz(-0.96894914) q[3];
sx q[3];
rz(1.8739353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0676607) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(0.78039783) q[0];
rz(0.91267768) q[1];
sx q[1];
rz(-0.78881216) q[1];
sx q[1];
rz(1.7524293) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71043832) q[0];
sx q[0];
rz(-1.2200933) q[0];
sx q[0];
rz(1.8560657) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58163403) q[2];
sx q[2];
rz(-1.6791145) q[2];
sx q[2];
rz(-1.9519639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.011574419) q[1];
sx q[1];
rz(-0.60638529) q[1];
sx q[1];
rz(2.9467086) q[1];
rz(-pi) q[2];
rz(1.0205998) q[3];
sx q[3];
rz(-0.76348272) q[3];
sx q[3];
rz(-2.965345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7954365) q[2];
sx q[2];
rz(-0.6554335) q[2];
sx q[2];
rz(2.9620841) q[2];
rz(1.6620212) q[3];
sx q[3];
rz(-2.447465) q[3];
sx q[3];
rz(3.1365385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805304) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(-1.4638715) q[0];
rz(0.013710984) q[1];
sx q[1];
rz(-1.6746215) q[1];
sx q[1];
rz(-2.1965068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6701407) q[0];
sx q[0];
rz(-1.2218804) q[0];
sx q[0];
rz(0.80100153) q[0];
rz(-pi) q[1];
rz(1.361998) q[2];
sx q[2];
rz(-1.2660625) q[2];
sx q[2];
rz(1.7727675) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7146804) q[1];
sx q[1];
rz(-1.1537997) q[1];
sx q[1];
rz(0.51583536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1626445) q[3];
sx q[3];
rz(-1.1148387) q[3];
sx q[3];
rz(0.78612529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1549687) q[2];
sx q[2];
rz(-2.0782317) q[2];
sx q[2];
rz(-1.0490136) q[2];
rz(-1.6074041) q[3];
sx q[3];
rz(-1.0020071) q[3];
sx q[3];
rz(1.3656507) q[3];
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
rz(-0.55117115) q[0];
sx q[0];
rz(-2.6541002) q[0];
sx q[0];
rz(0.76229873) q[0];
rz(2.6679692) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(0.90829888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1414836) q[0];
sx q[0];
rz(-2.3628652) q[0];
sx q[0];
rz(-1.578192) q[0];
x q[1];
rz(-1.3495096) q[2];
sx q[2];
rz(-1.7402116) q[2];
sx q[2];
rz(-2.7986023) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1181284) q[1];
sx q[1];
rz(-0.23734094) q[1];
sx q[1];
rz(-1.6327536) q[1];
x q[2];
rz(-1.5987464) q[3];
sx q[3];
rz(-1.5833491) q[3];
sx q[3];
rz(-2.8421228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0198387) q[2];
sx q[2];
rz(-0.25889954) q[2];
sx q[2];
rz(0.961595) q[2];
rz(2.614295) q[3];
sx q[3];
rz(-1.7617825) q[3];
sx q[3];
rz(-0.56078792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90802646) q[0];
sx q[0];
rz(-2.5416424) q[0];
sx q[0];
rz(2.7954234) q[0];
rz(1.2973805) q[1];
sx q[1];
rz(-2.4751414) q[1];
sx q[1];
rz(0.60874879) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.611871) q[0];
sx q[0];
rz(-0.68096113) q[0];
sx q[0];
rz(1.7044071) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26774755) q[2];
sx q[2];
rz(-1.367298) q[2];
sx q[2];
rz(0.45058077) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.765878) q[1];
sx q[1];
rz(-2.8676413) q[1];
sx q[1];
rz(2.8841444) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5777588) q[3];
sx q[3];
rz(-2.2632416) q[3];
sx q[3];
rz(-3.0094224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.346647) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(-2.4897599) q[2];
rz(2.4018304) q[3];
sx q[3];
rz(-2.0756105) q[3];
sx q[3];
rz(2.2447926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5480492) q[0];
sx q[0];
rz(-2.6683922) q[0];
sx q[0];
rz(-2.5439673) q[0];
rz(1.301544) q[1];
sx q[1];
rz(-1.5751244) q[1];
sx q[1];
rz(-2.8155933) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2151011) q[0];
sx q[0];
rz(-2.282906) q[0];
sx q[0];
rz(2.7551921) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8284441) q[2];
sx q[2];
rz(-1.2533873) q[2];
sx q[2];
rz(-0.35955113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44931901) q[1];
sx q[1];
rz(-2.2137768) q[1];
sx q[1];
rz(2.8161956) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11362152) q[3];
sx q[3];
rz(-1.4534165) q[3];
sx q[3];
rz(-2.6746933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43805435) q[2];
sx q[2];
rz(-1.9376829) q[2];
sx q[2];
rz(2.8013012) q[2];
rz(1.5002286) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(0.59804183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21094766) q[0];
sx q[0];
rz(-1.1177381) q[0];
sx q[0];
rz(3.0979544) q[0];
rz(-2.4234407) q[1];
sx q[1];
rz(-1.3190045) q[1];
sx q[1];
rz(-1.7125548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14075101) q[0];
sx q[0];
rz(-1.4006097) q[0];
sx q[0];
rz(-0.19077459) q[0];
rz(-pi) q[1];
rz(-2.9731644) q[2];
sx q[2];
rz(-0.68466869) q[2];
sx q[2];
rz(1.9569966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.43914) q[1];
sx q[1];
rz(-2.0589152) q[1];
sx q[1];
rz(-1.4737275) q[1];
rz(-2.7765482) q[3];
sx q[3];
rz(-1.2170346) q[3];
sx q[3];
rz(0.34788528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67937294) q[2];
sx q[2];
rz(-0.93928176) q[2];
sx q[2];
rz(-0.77793724) q[2];
rz(-2.098162) q[3];
sx q[3];
rz(-1.5305887) q[3];
sx q[3];
rz(1.9061609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048024561) q[0];
sx q[0];
rz(-2.1401736) q[0];
sx q[0];
rz(-2.370136) q[0];
rz(1.363516) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(-0.58495782) q[2];
sx q[2];
rz(-0.81333209) q[2];
sx q[2];
rz(1.3438136) q[2];
rz(-2.1322973) q[3];
sx q[3];
rz(-0.75405585) q[3];
sx q[3];
rz(2.8094359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
