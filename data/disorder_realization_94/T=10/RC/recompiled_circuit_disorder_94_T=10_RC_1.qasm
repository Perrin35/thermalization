OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(-1.5885408) q[0];
sx q[0];
rz(1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(0.61520666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9278487) q[0];
sx q[0];
rz(-0.97457492) q[0];
sx q[0];
rz(-0.56791373) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46703672) q[2];
sx q[2];
rz(-0.39614284) q[2];
sx q[2];
rz(-0.064183891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93278904) q[1];
sx q[1];
rz(-0.99901576) q[1];
sx q[1];
rz(0.51704452) q[1];
rz(-2.8350713) q[3];
sx q[3];
rz(-1.74311) q[3];
sx q[3];
rz(0.70787187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(2.8033076) q[2];
rz(1.4398549) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.171339) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(-1.617584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.679927) q[0];
sx q[0];
rz(-1.7704417) q[0];
sx q[0];
rz(3.1398849) q[0];
rz(1.6127869) q[2];
sx q[2];
rz(-0.4547555) q[2];
sx q[2];
rz(3.0595879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6807032) q[1];
sx q[1];
rz(-1.6398805) q[1];
sx q[1];
rz(0.99980385) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1094692) q[3];
sx q[3];
rz(-0.20878775) q[3];
sx q[3];
rz(1.8148282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.1478109) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-2.8258064) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(-0.25207239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.789327) q[0];
sx q[0];
rz(-2.0354712) q[0];
sx q[0];
rz(-2.4436823) q[0];
rz(-pi) q[1];
rz(2.2601068) q[2];
sx q[2];
rz(-0.93886095) q[2];
sx q[2];
rz(1.2607247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6985059) q[1];
sx q[1];
rz(-1.7585187) q[1];
sx q[1];
rz(2.1952573) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39450816) q[3];
sx q[3];
rz(-1.2393701) q[3];
sx q[3];
rz(3.1098207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1217653) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(0.034051731) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5144192) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(1.5267641) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(0.70708752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2061545) q[0];
sx q[0];
rz(-1.1797138) q[0];
sx q[0];
rz(0.48396707) q[0];
rz(0.25755067) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(-0.79007733) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4428963) q[1];
sx q[1];
rz(-0.92123182) q[1];
sx q[1];
rz(-0.15028468) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96890038) q[3];
sx q[3];
rz(-2.7728191) q[3];
sx q[3];
rz(0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2924071) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(2.6814931) q[2];
rz(1.397331) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8206772) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(1.0955411) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(-2.8869693) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8481962) q[0];
sx q[0];
rz(-3.0612429) q[0];
sx q[0];
rz(-1.3991762) q[0];
rz(-pi) q[1];
x q[1];
rz(0.014572797) q[2];
sx q[2];
rz(-0.99327786) q[2];
sx q[2];
rz(0.84469675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9533206) q[1];
sx q[1];
rz(-2.4362262) q[1];
sx q[1];
rz(-0.377368) q[1];
x q[2];
rz(2.280064) q[3];
sx q[3];
rz(-1.9083175) q[3];
sx q[3];
rz(2.5634114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.022481) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6901533) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(0.29944637) q[0];
rz(2.1014138) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(0.20656955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5305938) q[0];
sx q[0];
rz(-1.3337413) q[0];
sx q[0];
rz(-0.80942746) q[0];
x q[1];
rz(2.0986404) q[2];
sx q[2];
rz(-1.9163418) q[2];
sx q[2];
rz(-2.8258459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8787074) q[1];
sx q[1];
rz(-0.72808121) q[1];
sx q[1];
rz(1.3036149) q[1];
rz(-pi) q[2];
rz(-2.2264678) q[3];
sx q[3];
rz(-1.3279928) q[3];
sx q[3];
rz(1.9608998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(-2.3287866) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92641002) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(0.38152951) q[0];
rz(0.58386699) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.8136224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8570003) q[0];
sx q[0];
rz(-0.45495957) q[0];
sx q[0];
rz(-1.8332464) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1778657) q[2];
sx q[2];
rz(-2.4237195) q[2];
sx q[2];
rz(1.8075862) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2139637) q[1];
sx q[1];
rz(-0.8968401) q[1];
sx q[1];
rz(-0.2143292) q[1];
rz(-pi) q[2];
rz(2.6286969) q[3];
sx q[3];
rz(-1.1937871) q[3];
sx q[3];
rz(2.5412154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(0.84806228) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(0.18187901) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(0.95091933) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88175854) q[0];
sx q[0];
rz(-1.8108597) q[0];
sx q[0];
rz(1.865922) q[0];
x q[1];
rz(1.3457001) q[2];
sx q[2];
rz(-0.27268073) q[2];
sx q[2];
rz(2.8358012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7656895) q[1];
sx q[1];
rz(-1.5058869) q[1];
sx q[1];
rz(2.9325571) q[1];
x q[2];
rz(0.82138737) q[3];
sx q[3];
rz(-2.6139724) q[3];
sx q[3];
rz(-1.9963095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(0.075909464) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(0.63265911) q[3];
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
rz(-2.6178745) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(0.41704047) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(-0.65972796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035559) q[0];
sx q[0];
rz(-0.8350026) q[0];
sx q[0];
rz(0.8291709) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21978901) q[2];
sx q[2];
rz(-0.79384365) q[2];
sx q[2];
rz(1.1500051) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8068741) q[1];
sx q[1];
rz(-1.2804619) q[1];
sx q[1];
rz(-2.2335386) q[1];
rz(3.0845853) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(-1.9407335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(-1.0260322) q[2];
rz(-0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.898107) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.8359258) q[0];
rz(-1.1765515) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(-2.1059039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0337692) q[0];
sx q[0];
rz(-2.6080837) q[0];
sx q[0];
rz(-0.67722042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91949384) q[2];
sx q[2];
rz(-0.80923015) q[2];
sx q[2];
rz(-2.2724255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4020821) q[1];
sx q[1];
rz(-1.2716736) q[1];
sx q[1];
rz(-1.3189391) q[1];
rz(-pi) q[2];
rz(1.3003179) q[3];
sx q[3];
rz(-0.83993739) q[3];
sx q[3];
rz(-0.65281103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(-1.1516085) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37968996) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(-1.5079386) q[1];
sx q[1];
rz(-2.5588551) q[1];
sx q[1];
rz(-0.48164639) q[1];
rz(-1.528231) q[2];
sx q[2];
rz(-1.1602989) q[2];
sx q[2];
rz(1.5699671) q[2];
rz(-0.012398331) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
