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
rz(-3.0847235) q[0];
sx q[0];
rz(-2.9480204) q[0];
sx q[0];
rz(-2.733736) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(4.3932228) q[1];
sx q[1];
rz(7.8235758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33690573) q[0];
sx q[0];
rz(-1.1338455) q[0];
sx q[0];
rz(-3.0880465) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7061739) q[2];
sx q[2];
rz(-1.1203566) q[2];
sx q[2];
rz(-0.13794402) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.08886628) q[1];
sx q[1];
rz(-1.7824934) q[1];
sx q[1];
rz(-0.39525169) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61397378) q[3];
sx q[3];
rz(-2.2212265) q[3];
sx q[3];
rz(-0.87343317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5876329) q[2];
sx q[2];
rz(-2.0822058) q[2];
sx q[2];
rz(0.21767347) q[2];
rz(-0.31146464) q[3];
sx q[3];
rz(-0.61190999) q[3];
sx q[3];
rz(-1.9757087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.72164732) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(-2.4496147) q[0];
rz(-0.75633374) q[1];
sx q[1];
rz(-2.6192009) q[1];
sx q[1];
rz(-1.5914894) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1082008) q[0];
sx q[0];
rz(-2.4188359) q[0];
sx q[0];
rz(0.18527822) q[0];
x q[1];
rz(1.9687551) q[2];
sx q[2];
rz(-1.474547) q[2];
sx q[2];
rz(0.41589662) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.31526947) q[1];
sx q[1];
rz(-2.9401952) q[1];
sx q[1];
rz(0.65314318) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1308998) q[3];
sx q[3];
rz(-2.0030795) q[3];
sx q[3];
rz(-0.11350693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1456566) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(-0.23925979) q[2];
rz(1.650882) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(1.2283121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8656798) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(3.1269585) q[0];
rz(-0.89302653) q[1];
sx q[1];
rz(-0.97517401) q[1];
sx q[1];
rz(2.9258974) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94616643) q[0];
sx q[0];
rz(-3.0901244) q[0];
sx q[0];
rz(-1.6610751) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9142308) q[2];
sx q[2];
rz(-0.94007713) q[2];
sx q[2];
rz(-1.7523927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0862992) q[1];
sx q[1];
rz(-2.8759416) q[1];
sx q[1];
rz(-0.50528996) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78554947) q[3];
sx q[3];
rz(-0.86570287) q[3];
sx q[3];
rz(2.7712703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0594242) q[2];
sx q[2];
rz(-1.9614204) q[2];
sx q[2];
rz(0.99349418) q[2];
rz(0.70837402) q[3];
sx q[3];
rz(-1.1185442) q[3];
sx q[3];
rz(0.12349252) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941536) q[0];
sx q[0];
rz(-0.49462947) q[0];
sx q[0];
rz(0.69394773) q[0];
rz(2.8548062) q[1];
sx q[1];
rz(-1.5319805) q[1];
sx q[1];
rz(-2.0538816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0407242) q[0];
sx q[0];
rz(-0.52170149) q[0];
sx q[0];
rz(-1.3203599) q[0];
rz(0.7544341) q[2];
sx q[2];
rz(-0.70068073) q[2];
sx q[2];
rz(0.59631077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4162594) q[1];
sx q[1];
rz(-0.97914588) q[1];
sx q[1];
rz(-2.7736431) q[1];
x q[2];
rz(-3.089043) q[3];
sx q[3];
rz(-0.71501117) q[3];
sx q[3];
rz(-0.258436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9168758) q[2];
sx q[2];
rz(-1.9852873) q[2];
sx q[2];
rz(1.7737596) q[2];
rz(1.9035089) q[3];
sx q[3];
rz(-1.5690469) q[3];
sx q[3];
rz(0.21627538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5991768) q[0];
sx q[0];
rz(-1.2681862) q[0];
sx q[0];
rz(2.5237778) q[0];
rz(1.1082331) q[1];
sx q[1];
rz(-0.30312678) q[1];
sx q[1];
rz(-0.88315001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41048356) q[0];
sx q[0];
rz(-1.5781814) q[0];
sx q[0];
rz(-2.3855668) q[0];
rz(-2.8034535) q[2];
sx q[2];
rz(-0.45881841) q[2];
sx q[2];
rz(-0.79566075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7204764) q[1];
sx q[1];
rz(-1.5052374) q[1];
sx q[1];
rz(2.0620729) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95728504) q[3];
sx q[3];
rz(-2.4924879) q[3];
sx q[3];
rz(-1.7136991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21657476) q[2];
sx q[2];
rz(-0.25187945) q[2];
sx q[2];
rz(-1.3612755) q[2];
rz(1.2733634) q[3];
sx q[3];
rz(-1.7073771) q[3];
sx q[3];
rz(2.1586965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(3.0971766) q[0];
sx q[0];
rz(-1.2460848) q[0];
sx q[0];
rz(-2.521305) q[0];
rz(2.7445131) q[1];
sx q[1];
rz(-0.47080165) q[1];
sx q[1];
rz(-1.1995859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9387377) q[0];
sx q[0];
rz(-1.5498383) q[0];
sx q[0];
rz(-1.6771132) q[0];
x q[1];
rz(1.7825837) q[2];
sx q[2];
rz(-1.4265991) q[2];
sx q[2];
rz(1.5590925) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37383693) q[1];
sx q[1];
rz(-1.7636313) q[1];
sx q[1];
rz(0.53044935) q[1];
rz(-pi) q[2];
rz(1.4523466) q[3];
sx q[3];
rz(-1.9129686) q[3];
sx q[3];
rz(-2.7717048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8951796) q[2];
sx q[2];
rz(-1.8975039) q[2];
sx q[2];
rz(-0.67071521) q[2];
rz(2.8504168) q[3];
sx q[3];
rz(-0.61101919) q[3];
sx q[3];
rz(-2.5197855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42596844) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(2.9631462) q[0];
rz(-0.87002358) q[1];
sx q[1];
rz(-0.88302892) q[1];
sx q[1];
rz(2.3387486) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69235301) q[0];
sx q[0];
rz(-0.62004161) q[0];
sx q[0];
rz(-0.043278261) q[0];
rz(0.58391352) q[2];
sx q[2];
rz(-1.7758992) q[2];
sx q[2];
rz(2.5896304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4124914) q[1];
sx q[1];
rz(-1.8961398) q[1];
sx q[1];
rz(-2.4304683) q[1];
rz(-pi) q[2];
rz(-1.9196627) q[3];
sx q[3];
rz(-0.65769821) q[3];
sx q[3];
rz(1.5546834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37658438) q[2];
sx q[2];
rz(-1.2976126) q[2];
sx q[2];
rz(0.89361781) q[2];
rz(1.5997959) q[3];
sx q[3];
rz(-1.4897646) q[3];
sx q[3];
rz(-2.5347575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97745085) q[0];
sx q[0];
rz(-2.7315388) q[0];
sx q[0];
rz(-0.2151016) q[0];
rz(-2.2970301) q[1];
sx q[1];
rz(-1.8212049) q[1];
sx q[1];
rz(2.3473158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1298908) q[0];
sx q[0];
rz(-3.112535) q[0];
sx q[0];
rz(0.41883166) q[0];
x q[1];
rz(-1.2634981) q[2];
sx q[2];
rz(-1.1329044) q[2];
sx q[2];
rz(1.7192993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25674838) q[1];
sx q[1];
rz(-2.1089541) q[1];
sx q[1];
rz(3.0358008) q[1];
rz(-2.3598757) q[3];
sx q[3];
rz(-1.8060883) q[3];
sx q[3];
rz(-0.58636753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27791417) q[2];
sx q[2];
rz(-1.9387551) q[2];
sx q[2];
rz(-1.4097144) q[2];
rz(-2.388741) q[3];
sx q[3];
rz(-1.9949621) q[3];
sx q[3];
rz(2.1394155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79733217) q[0];
sx q[0];
rz(-0.40818885) q[0];
sx q[0];
rz(-0.97287384) q[0];
rz(1.5219888) q[1];
sx q[1];
rz(-1.3643967) q[1];
sx q[1];
rz(0.63046986) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3469543) q[0];
sx q[0];
rz(-1.7525867) q[0];
sx q[0];
rz(-1.4304881) q[0];
rz(-1.7277355) q[2];
sx q[2];
rz(-1.5782159) q[2];
sx q[2];
rz(-2.3049624) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18955293) q[1];
sx q[1];
rz(-1.5209043) q[1];
sx q[1];
rz(-0.76811172) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6338324) q[3];
sx q[3];
rz(-1.3975289) q[3];
sx q[3];
rz(-1.8484556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.015532739) q[2];
sx q[2];
rz(-1.1522747) q[2];
sx q[2];
rz(1.8401592) q[2];
rz(1.556501) q[3];
sx q[3];
rz(-0.82587487) q[3];
sx q[3];
rz(-0.41527709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5179317) q[0];
sx q[0];
rz(-0.15075891) q[0];
sx q[0];
rz(-0.49322042) q[0];
rz(1.2552235) q[1];
sx q[1];
rz(-1.247765) q[1];
sx q[1];
rz(0.36466041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5887774) q[0];
sx q[0];
rz(-1.7087666) q[0];
sx q[0];
rz(1.9068043) q[0];
x q[1];
rz(0.45367806) q[2];
sx q[2];
rz(-0.96900207) q[2];
sx q[2];
rz(1.8555249) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4019822) q[1];
sx q[1];
rz(-1.6270042) q[1];
sx q[1];
rz(-0.79938857) q[1];
rz(-pi) q[2];
x q[2];
rz(0.040382645) q[3];
sx q[3];
rz(-1.2546439) q[3];
sx q[3];
rz(0.97163661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4517335) q[2];
sx q[2];
rz(-2.1351337) q[2];
sx q[2];
rz(0.9922007) q[2];
rz(-2.774488) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(0.67830694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51919666) q[0];
sx q[0];
rz(-1.4764897) q[0];
sx q[0];
rz(1.3523703) q[0];
rz(-0.92338152) q[1];
sx q[1];
rz(-2.2877749) q[1];
sx q[1];
rz(1.9705082) q[1];
rz(-0.95465701) q[2];
sx q[2];
rz(-2.8980394) q[2];
sx q[2];
rz(0.63799636) q[2];
rz(-0.67305858) q[3];
sx q[3];
rz(-1.0792427) q[3];
sx q[3];
rz(-1.4918809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
