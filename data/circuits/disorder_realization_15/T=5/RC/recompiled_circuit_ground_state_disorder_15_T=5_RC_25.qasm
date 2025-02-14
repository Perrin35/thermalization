OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.48159596) q[0];
sx q[0];
rz(3.5591535) q[0];
sx q[0];
rz(9.4661718) q[0];
rz(0.59250915) q[1];
sx q[1];
rz(-0.23696391) q[1];
sx q[1];
rz(-0.89495975) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76413078) q[0];
sx q[0];
rz(-1.9870166) q[0];
sx q[0];
rz(-0.14472117) q[0];
x q[1];
rz(-2.1092806) q[2];
sx q[2];
rz(-2.5729524) q[2];
sx q[2];
rz(1.2775354) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30679157) q[1];
sx q[1];
rz(-0.39265206) q[1];
sx q[1];
rz(0.77775915) q[1];
rz(2.7616463) q[3];
sx q[3];
rz(-0.11726876) q[3];
sx q[3];
rz(-0.46739551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71901739) q[2];
sx q[2];
rz(-2.1696679) q[2];
sx q[2];
rz(1.0389339) q[2];
rz(-0.59764189) q[3];
sx q[3];
rz(-0.20331764) q[3];
sx q[3];
rz(1.982127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1379625) q[0];
sx q[0];
rz(-2.2901386) q[0];
sx q[0];
rz(2.8077069) q[0];
rz(0.69260287) q[1];
sx q[1];
rz(-1.875016) q[1];
sx q[1];
rz(2.0803221) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675628) q[0];
sx q[0];
rz(-2.2252796) q[0];
sx q[0];
rz(-0.45463134) q[0];
rz(-pi) q[1];
rz(-1.3150923) q[2];
sx q[2];
rz(-1.1785011) q[2];
sx q[2];
rz(-0.43161156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0063707) q[1];
sx q[1];
rz(-1.3961482) q[1];
sx q[1];
rz(2.6937301) q[1];
x q[2];
rz(1.3676104) q[3];
sx q[3];
rz(-1.7082761) q[3];
sx q[3];
rz(0.41609496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4873203) q[2];
sx q[2];
rz(-1.4620179) q[2];
sx q[2];
rz(-0.1839323) q[2];
rz(0.0025302689) q[3];
sx q[3];
rz(-2.3380184) q[3];
sx q[3];
rz(2.7185503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84394395) q[0];
sx q[0];
rz(-0.14744814) q[0];
sx q[0];
rz(2.3609128) q[0];
rz(-2.2536904) q[1];
sx q[1];
rz(-2.0474515) q[1];
sx q[1];
rz(0.81781864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6827777) q[0];
sx q[0];
rz(-1.11894) q[0];
sx q[0];
rz(2.8917612) q[0];
rz(-pi) q[1];
rz(2.3505402) q[2];
sx q[2];
rz(-1.0246236) q[2];
sx q[2];
rz(0.56529048) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1314268) q[1];
sx q[1];
rz(-1.5714116) q[1];
sx q[1];
rz(-2.5593131) q[1];
x q[2];
rz(0.85237827) q[3];
sx q[3];
rz(-0.63567978) q[3];
sx q[3];
rz(-1.6851684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4320977) q[2];
sx q[2];
rz(-3.0165387) q[2];
sx q[2];
rz(-1.032426) q[2];
rz(0.89140511) q[3];
sx q[3];
rz(-0.67474198) q[3];
sx q[3];
rz(0.83511043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8472403) q[0];
sx q[0];
rz(-2.8645741) q[0];
sx q[0];
rz(2.5936122) q[0];
rz(-0.056593865) q[1];
sx q[1];
rz(-1.2013925) q[1];
sx q[1];
rz(-3.1170381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5540172) q[0];
sx q[0];
rz(-1.8880821) q[0];
sx q[0];
rz(-1.1486828) q[0];
rz(-pi) q[1];
rz(-1.9664832) q[2];
sx q[2];
rz(-1.8818237) q[2];
sx q[2];
rz(-2.0360586) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.430342) q[1];
sx q[1];
rz(-0.81626662) q[1];
sx q[1];
rz(-1.0984983) q[1];
x q[2];
rz(-1.914045) q[3];
sx q[3];
rz(-1.5446179) q[3];
sx q[3];
rz(0.13939339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6294488) q[2];
sx q[2];
rz(-2.5089743) q[2];
sx q[2];
rz(0.72652417) q[2];
rz(-1.7065382) q[3];
sx q[3];
rz(-1.8972998) q[3];
sx q[3];
rz(-1.6632891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59359819) q[0];
sx q[0];
rz(-1.7578121) q[0];
sx q[0];
rz(1.3872248) q[0];
rz(2.6902426) q[1];
sx q[1];
rz(-2.4001922) q[1];
sx q[1];
rz(-2.6373934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6071753) q[0];
sx q[0];
rz(-1.4481521) q[0];
sx q[0];
rz(-2.9706435) q[0];
rz(-pi) q[1];
rz(-2.4455657) q[2];
sx q[2];
rz(-2.2559204) q[2];
sx q[2];
rz(-2.7643577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5686149) q[1];
sx q[1];
rz(-2.0842127) q[1];
sx q[1];
rz(-0.77691369) q[1];
rz(2.3285015) q[3];
sx q[3];
rz(-0.92223863) q[3];
sx q[3];
rz(0.5483707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35895687) q[2];
sx q[2];
rz(-2.2613596) q[2];
sx q[2];
rz(3.1202988) q[2];
rz(3.0535789) q[3];
sx q[3];
rz(-2.8787677) q[3];
sx q[3];
rz(-2.9737441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39040318) q[0];
sx q[0];
rz(-2.5683371) q[0];
sx q[0];
rz(2.7943352) q[0];
rz(-1.685453) q[1];
sx q[1];
rz(-2.8442454) q[1];
sx q[1];
rz(-2.8567294) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96831095) q[0];
sx q[0];
rz(-1.8704458) q[0];
sx q[0];
rz(-2.8408627) q[0];
x q[1];
rz(0.059528298) q[2];
sx q[2];
rz(-0.7452508) q[2];
sx q[2];
rz(1.1203114) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9544475) q[1];
sx q[1];
rz(-1.3055077) q[1];
sx q[1];
rz(-2.250606) q[1];
x q[2];
rz(1.2142193) q[3];
sx q[3];
rz(-1.8468401) q[3];
sx q[3];
rz(2.2726633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0028637) q[2];
sx q[2];
rz(-1.8314654) q[2];
sx q[2];
rz(2.1076473) q[2];
rz(2.406534) q[3];
sx q[3];
rz(-1.4933973) q[3];
sx q[3];
rz(-0.66966164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1351778) q[0];
sx q[0];
rz(-2.4786351) q[0];
sx q[0];
rz(2.2354777) q[0];
rz(-0.16618973) q[1];
sx q[1];
rz(-2.6527185) q[1];
sx q[1];
rz(2.8320872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75342732) q[0];
sx q[0];
rz(-1.7953402) q[0];
sx q[0];
rz(-2.9675447) q[0];
x q[1];
rz(0.1997582) q[2];
sx q[2];
rz(-0.80325489) q[2];
sx q[2];
rz(-1.8920395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5823203) q[1];
sx q[1];
rz(-2.4666335) q[1];
sx q[1];
rz(-1.9155424) q[1];
rz(-pi) q[2];
rz(1.0528492) q[3];
sx q[3];
rz(-3.042683) q[3];
sx q[3];
rz(-0.78341752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4497778) q[2];
sx q[2];
rz(-2.4725627) q[2];
sx q[2];
rz(1.6548033) q[2];
rz(2.6438223) q[3];
sx q[3];
rz(-2.3064752) q[3];
sx q[3];
rz(0.2275137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787702) q[0];
sx q[0];
rz(-2.4202122) q[0];
sx q[0];
rz(-0.26807868) q[0];
rz(-0.90191853) q[1];
sx q[1];
rz(-2.0733158) q[1];
sx q[1];
rz(1.8394151) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64619795) q[0];
sx q[0];
rz(-1.7359635) q[0];
sx q[0];
rz(2.5423869) q[0];
rz(1.5728358) q[2];
sx q[2];
rz(-2.0173948) q[2];
sx q[2];
rz(1.0616605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0331653) q[1];
sx q[1];
rz(-1.4301586) q[1];
sx q[1];
rz(2.161639) q[1];
rz(-pi) q[2];
rz(1.4623649) q[3];
sx q[3];
rz(-2.83395) q[3];
sx q[3];
rz(0.82266146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34182081) q[2];
sx q[2];
rz(-1.2398961) q[2];
sx q[2];
rz(-0.58490252) q[2];
rz(-1.9386442) q[3];
sx q[3];
rz(-1.3936309) q[3];
sx q[3];
rz(1.5929619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0840266) q[0];
sx q[0];
rz(-2.6507222) q[0];
sx q[0];
rz(0.67114818) q[0];
rz(0.28823832) q[1];
sx q[1];
rz(-2.3858374) q[1];
sx q[1];
rz(0.63405687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8347522) q[0];
sx q[0];
rz(-3.1001087) q[0];
sx q[0];
rz(1.1783429) q[0];
rz(-0.05001683) q[2];
sx q[2];
rz(-1.4054308) q[2];
sx q[2];
rz(-0.53420137) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1131192) q[1];
sx q[1];
rz(-0.32580842) q[1];
sx q[1];
rz(0.079976639) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4692772) q[3];
sx q[3];
rz(-1.3642715) q[3];
sx q[3];
rz(0.12609161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5481698) q[2];
sx q[2];
rz(-1.5680743) q[2];
sx q[2];
rz(-2.8578952) q[2];
rz(0.13188322) q[3];
sx q[3];
rz(-2.7851084) q[3];
sx q[3];
rz(-0.62294817) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.236096) q[0];
sx q[0];
rz(-2.997969) q[0];
sx q[0];
rz(-0.96963257) q[0];
rz(-0.49599221) q[1];
sx q[1];
rz(-2.453936) q[1];
sx q[1];
rz(0.86404854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43661731) q[0];
sx q[0];
rz(-0.71067315) q[0];
sx q[0];
rz(-1.4628167) q[0];
x q[1];
rz(-2.0798565) q[2];
sx q[2];
rz(-2.6531918) q[2];
sx q[2];
rz(-2.050972) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38258994) q[1];
sx q[1];
rz(-2.14258) q[1];
sx q[1];
rz(2.5409805) q[1];
rz(0.91409036) q[3];
sx q[3];
rz(-1.9352846) q[3];
sx q[3];
rz(1.2165704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25829092) q[2];
sx q[2];
rz(-2.624888) q[2];
sx q[2];
rz(-1.0374163) q[2];
rz(2.9344905) q[3];
sx q[3];
rz(-2.243302) q[3];
sx q[3];
rz(-0.5894388) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265738) q[0];
sx q[0];
rz(-1.6828231) q[0];
sx q[0];
rz(-1.9376301) q[0];
rz(-0.013982458) q[1];
sx q[1];
rz(-1.7054066) q[1];
sx q[1];
rz(-1.1048497) q[1];
rz(0.034910708) q[2];
sx q[2];
rz(-1.3401122) q[2];
sx q[2];
rz(2.050436) q[2];
rz(2.8057475) q[3];
sx q[3];
rz(-0.8728948) q[3];
sx q[3];
rz(-0.030073312) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
