OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(1.3774435) q[0];
sx q[0];
rz(5.8537771) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(-1.5996999) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854092) q[0];
sx q[0];
rz(-1.2852291) q[0];
sx q[0];
rz(-2.9182316) q[0];
rz(2.9028106) q[2];
sx q[2];
rz(-2.9936643) q[2];
sx q[2];
rz(1.5419568) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0059389) q[1];
sx q[1];
rz(-0.42891132) q[1];
sx q[1];
rz(0.072838293) q[1];
rz(-pi) q[2];
rz(-1.2774994) q[3];
sx q[3];
rz(-1.9175005) q[3];
sx q[3];
rz(1.6270005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30449197) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(-0.692918) q[2];
rz(-0.20017008) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(-2.0076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8189341) q[0];
sx q[0];
rz(-0.19911961) q[0];
sx q[0];
rz(2.0767427) q[0];
rz(0.62243593) q[1];
sx q[1];
rz(-1.7935926) q[1];
sx q[1];
rz(-2.650824) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7642794) q[0];
sx q[0];
rz(-1.5945596) q[0];
sx q[0];
rz(1.5800493) q[0];
x q[1];
rz(2.7193473) q[2];
sx q[2];
rz(-1.7755055) q[2];
sx q[2];
rz(1.7038356) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1176396) q[1];
sx q[1];
rz(-2.3815386) q[1];
sx q[1];
rz(3.0031086) q[1];
rz(-pi) q[2];
x q[2];
rz(0.048100483) q[3];
sx q[3];
rz(-1.7399746) q[3];
sx q[3];
rz(-0.98911232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6530767) q[2];
sx q[2];
rz(-1.642903) q[2];
sx q[2];
rz(0.24031362) q[2];
rz(-0.49992418) q[3];
sx q[3];
rz(-2.5559055) q[3];
sx q[3];
rz(1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.2568473) q[0];
sx q[0];
rz(-2.186543) q[0];
sx q[0];
rz(-0.98168674) q[0];
rz(2.6495972) q[1];
sx q[1];
rz(-1.0849846) q[1];
sx q[1];
rz(1.6384151) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7273151) q[0];
sx q[0];
rz(-2.2637667) q[0];
sx q[0];
rz(-2.5056865) q[0];
x q[1];
rz(-1.515834) q[2];
sx q[2];
rz(-1.171249) q[2];
sx q[2];
rz(-1.9683226) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43736378) q[1];
sx q[1];
rz(-1.2384602) q[1];
sx q[1];
rz(1.1613013) q[1];
x q[2];
rz(1.668456) q[3];
sx q[3];
rz(-1.5499695) q[3];
sx q[3];
rz(1.5725003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4846399) q[2];
sx q[2];
rz(-2.6159365) q[2];
sx q[2];
rz(-0.12118113) q[2];
rz(1.0080522) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(-1.0813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0722395) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(-0.51373154) q[0];
rz(-0.95798245) q[1];
sx q[1];
rz(-1.999141) q[1];
sx q[1];
rz(-1.6987919) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0009036) q[0];
sx q[0];
rz(-1.1620518) q[0];
sx q[0];
rz(-1.391414) q[0];
x q[1];
rz(-0.36607217) q[2];
sx q[2];
rz(-2.5806081) q[2];
sx q[2];
rz(1.6268886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.629238) q[1];
sx q[1];
rz(-1.6529473) q[1];
sx q[1];
rz(2.1869529) q[1];
rz(-pi) q[2];
rz(1.6677271) q[3];
sx q[3];
rz(-2.2301794) q[3];
sx q[3];
rz(-1.7771378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5865667) q[2];
sx q[2];
rz(-0.96379605) q[2];
sx q[2];
rz(-1.8761934) q[2];
rz(-1.1749367) q[3];
sx q[3];
rz(-0.73052162) q[3];
sx q[3];
rz(-1.1635715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88605276) q[0];
sx q[0];
rz(-2.9386254) q[0];
sx q[0];
rz(-1.3245921) q[0];
rz(-0.96877226) q[1];
sx q[1];
rz(-1.4854393) q[1];
sx q[1];
rz(-2.103215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85920716) q[0];
sx q[0];
rz(-0.41658724) q[0];
sx q[0];
rz(-0.4075012) q[0];
rz(-pi) q[1];
rz(-2.3164301) q[2];
sx q[2];
rz(-1.0010011) q[2];
sx q[2];
rz(1.2078326) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9818791) q[1];
sx q[1];
rz(-2.7687589) q[1];
sx q[1];
rz(-3.0283958) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0821003) q[3];
sx q[3];
rz(-0.97701525) q[3];
sx q[3];
rz(0.64835129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6394627) q[2];
sx q[2];
rz(-0.30439964) q[2];
sx q[2];
rz(2.2501865) q[2];
rz(1.8799479) q[3];
sx q[3];
rz(-1.5104048) q[3];
sx q[3];
rz(0.6663028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0618184) q[0];
sx q[0];
rz(-2.6462055) q[0];
sx q[0];
rz(-0.99639446) q[0];
rz(0.90006104) q[1];
sx q[1];
rz(-0.78949094) q[1];
sx q[1];
rz(-2.9642504) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19939199) q[0];
sx q[0];
rz(-2.9900402) q[0];
sx q[0];
rz(1.5723537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93219535) q[2];
sx q[2];
rz(-1.5538486) q[2];
sx q[2];
rz(1.0826966) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85834938) q[1];
sx q[1];
rz(-2.5734786) q[1];
sx q[1];
rz(2.1217705) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3979891) q[3];
sx q[3];
rz(-0.51522845) q[3];
sx q[3];
rz(-0.15095161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29048723) q[2];
sx q[2];
rz(-1.1815716) q[2];
sx q[2];
rz(-0.30430749) q[2];
rz(-2.815222) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(2.2296026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.07311634) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(-2.2090744) q[0];
rz(3.0724691) q[1];
sx q[1];
rz(-2.9239475) q[1];
sx q[1];
rz(1.0401915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067834082) q[0];
sx q[0];
rz(-0.76601765) q[0];
sx q[0];
rz(-0.037184663) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1625453) q[2];
sx q[2];
rz(-0.30610105) q[2];
sx q[2];
rz(-0.09397587) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55807796) q[1];
sx q[1];
rz(-2.3574074) q[1];
sx q[1];
rz(0.23434831) q[1];
rz(-1.0024985) q[3];
sx q[3];
rz(-0.95469427) q[3];
sx q[3];
rz(1.440986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0065464) q[2];
sx q[2];
rz(-2.4032205) q[2];
sx q[2];
rz(-2.4086003) q[2];
rz(-1.0829571) q[3];
sx q[3];
rz(-2.0854009) q[3];
sx q[3];
rz(2.1347031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55423823) q[0];
sx q[0];
rz(-3.0576958) q[0];
sx q[0];
rz(-2.6485637) q[0];
rz(-0.21656491) q[1];
sx q[1];
rz(-2.000587) q[1];
sx q[1];
rz(-1.0239748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42070779) q[0];
sx q[0];
rz(-2.5455089) q[0];
sx q[0];
rz(0.021230265) q[0];
x q[1];
rz(-2.8093084) q[2];
sx q[2];
rz(-1.6783251) q[2];
sx q[2];
rz(2.8018453) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.4272642) q[1];
sx q[1];
rz(-1.9275394) q[1];
sx q[1];
rz(-1.7206232) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1795565) q[3];
sx q[3];
rz(-2.4613698) q[3];
sx q[3];
rz(0.30339751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0563125) q[2];
sx q[2];
rz(-1.6322501) q[2];
sx q[2];
rz(-2.4658266) q[2];
rz(2.7285649) q[3];
sx q[3];
rz(-1.4512117) q[3];
sx q[3];
rz(-0.54615027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6290879) q[0];
sx q[0];
rz(-0.65852037) q[0];
sx q[0];
rz(-3.107048) q[0];
rz(3.0870364) q[1];
sx q[1];
rz(-1.1628393) q[1];
sx q[1];
rz(-0.82130718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9050582) q[0];
sx q[0];
rz(-2.5590959) q[0];
sx q[0];
rz(2.9597917) q[0];
rz(-pi) q[1];
rz(-3.0376833) q[2];
sx q[2];
rz(-3.0523411) q[2];
sx q[2];
rz(2.1572942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7682338) q[1];
sx q[1];
rz(-1.4211417) q[1];
sx q[1];
rz(1.6878769) q[1];
rz(-pi) q[2];
rz(0.91854878) q[3];
sx q[3];
rz(-1.1763185) q[3];
sx q[3];
rz(-0.6834417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36046946) q[2];
sx q[2];
rz(-2.1253773) q[2];
sx q[2];
rz(3.0086369) q[2];
rz(-0.41108701) q[3];
sx q[3];
rz(-1.6229595) q[3];
sx q[3];
rz(-1.0857922) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046831176) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(1.2588311) q[0];
rz(-1.5065441) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(-0.81319317) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55594873) q[0];
sx q[0];
rz(-3.0559982) q[0];
sx q[0];
rz(0.42285796) q[0];
rz(-2.2426064) q[2];
sx q[2];
rz(-1.4788663) q[2];
sx q[2];
rz(-1.7882333) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.836536) q[1];
sx q[1];
rz(-0.99545762) q[1];
sx q[1];
rz(1.025455) q[1];
rz(-0.45222262) q[3];
sx q[3];
rz(-1.4503696) q[3];
sx q[3];
rz(-2.6129006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64858156) q[2];
sx q[2];
rz(-2.0040671) q[2];
sx q[2];
rz(1.0515593) q[2];
rz(2.1227396) q[3];
sx q[3];
rz(-0.84210432) q[3];
sx q[3];
rz(2.4837608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9857585) q[0];
sx q[0];
rz(-0.65407615) q[0];
sx q[0];
rz(0.88055897) q[0];
rz(1.3195994) q[1];
sx q[1];
rz(-1.1450014) q[1];
sx q[1];
rz(-3.0897279) q[1];
rz(-1.6435087) q[2];
sx q[2];
rz(-1.3617392) q[2];
sx q[2];
rz(-2.9206252) q[2];
rz(-1.6307835) q[3];
sx q[3];
rz(-1.7675478) q[3];
sx q[3];
rz(0.02789733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
