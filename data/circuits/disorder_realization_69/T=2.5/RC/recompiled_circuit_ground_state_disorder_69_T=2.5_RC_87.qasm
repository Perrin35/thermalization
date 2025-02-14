OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2332239) q[0];
sx q[0];
rz(-1.5492735) q[0];
sx q[0];
rz(-0.17750658) q[0];
rz(-1.5818051) q[1];
sx q[1];
rz(-1.3270562) q[1];
sx q[1];
rz(1.1787193) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27025199) q[0];
sx q[0];
rz(-0.9862904) q[0];
sx q[0];
rz(0.22342213) q[0];
rz(2.841973) q[2];
sx q[2];
rz(-2.1914406) q[2];
sx q[2];
rz(0.52052467) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70397729) q[1];
sx q[1];
rz(-2.0052138) q[1];
sx q[1];
rz(-2.2056286) q[1];
rz(-3.0079675) q[3];
sx q[3];
rz(-2.3006064) q[3];
sx q[3];
rz(1.3624024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.061964758) q[2];
sx q[2];
rz(-2.2140333) q[2];
sx q[2];
rz(0.24688841) q[2];
rz(-0.43761349) q[3];
sx q[3];
rz(-1.0381617) q[3];
sx q[3];
rz(-0.26052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058411773) q[0];
sx q[0];
rz(-2.658168) q[0];
sx q[0];
rz(-1.1046945) q[0];
rz(2.3339234) q[1];
sx q[1];
rz(-0.75391114) q[1];
sx q[1];
rz(0.63899904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09428962) q[0];
sx q[0];
rz(-1.847379) q[0];
sx q[0];
rz(-1.7533789) q[0];
rz(-pi) q[1];
rz(2.866687) q[2];
sx q[2];
rz(-1.6187117) q[2];
sx q[2];
rz(2.3570182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2874291) q[1];
sx q[1];
rz(-0.87814444) q[1];
sx q[1];
rz(1.0911029) q[1];
rz(-pi) q[2];
rz(-2.2792729) q[3];
sx q[3];
rz(-0.42565027) q[3];
sx q[3];
rz(-1.6274393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69940007) q[2];
sx q[2];
rz(-1.3905808) q[2];
sx q[2];
rz(-0.81713444) q[2];
rz(-2.6160431) q[3];
sx q[3];
rz(-2.3216129) q[3];
sx q[3];
rz(-1.9513963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84080559) q[0];
sx q[0];
rz(-1.5910933) q[0];
sx q[0];
rz(2.9817885) q[0];
rz(2.6218759) q[1];
sx q[1];
rz(-2.0964041) q[1];
sx q[1];
rz(0.91667485) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5477429) q[0];
sx q[0];
rz(-1.5968939) q[0];
sx q[0];
rz(-2.5711077) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4009928) q[2];
sx q[2];
rz(-1.3889905) q[2];
sx q[2];
rz(-1.2642853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8868272) q[1];
sx q[1];
rz(-1.6187833) q[1];
sx q[1];
rz(1.7601723) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9514922) q[3];
sx q[3];
rz(-0.62119609) q[3];
sx q[3];
rz(2.7060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2436169) q[2];
sx q[2];
rz(-2.6349082) q[2];
sx q[2];
rz(-1.766073) q[2];
rz(-0.033179387) q[3];
sx q[3];
rz(-1.5876007) q[3];
sx q[3];
rz(2.3110265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0783949) q[0];
sx q[0];
rz(-2.4802408) q[0];
sx q[0];
rz(-2.1589494) q[0];
rz(0.64535087) q[1];
sx q[1];
rz(-1.9159562) q[1];
sx q[1];
rz(0.4172079) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44051927) q[0];
sx q[0];
rz(-2.5357781) q[0];
sx q[0];
rz(0.54744253) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9511767) q[2];
sx q[2];
rz(-1.6350758) q[2];
sx q[2];
rz(0.96997875) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.55298333) q[1];
sx q[1];
rz(-0.89834301) q[1];
sx q[1];
rz(-1.49468) q[1];
x q[2];
rz(2.5933635) q[3];
sx q[3];
rz(-1.3286202) q[3];
sx q[3];
rz(0.13270031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30503201) q[2];
sx q[2];
rz(-0.3362259) q[2];
sx q[2];
rz(-1.9545371) q[2];
rz(-0.120397) q[3];
sx q[3];
rz(-1.7359366) q[3];
sx q[3];
rz(-0.19336893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8315941) q[0];
sx q[0];
rz(-0.11045063) q[0];
sx q[0];
rz(2.4647392) q[0];
rz(-1.5325158) q[1];
sx q[1];
rz(-2.0307505) q[1];
sx q[1];
rz(-0.19930509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53056419) q[0];
sx q[0];
rz(-0.26913759) q[0];
sx q[0];
rz(-1.2236973) q[0];
rz(-pi) q[1];
rz(-0.32782741) q[2];
sx q[2];
rz(-2.77711) q[2];
sx q[2];
rz(2.5555697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.052303132) q[1];
sx q[1];
rz(-0.48434231) q[1];
sx q[1];
rz(-0.67708309) q[1];
rz(-pi) q[2];
rz(1.293888) q[3];
sx q[3];
rz(-1.5770638) q[3];
sx q[3];
rz(-0.66440493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4946332) q[2];
sx q[2];
rz(-0.55611742) q[2];
sx q[2];
rz(2.6456918) q[2];
rz(2.2206709) q[3];
sx q[3];
rz(-2.0995188) q[3];
sx q[3];
rz(-2.8310149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344236) q[0];
sx q[0];
rz(-1.9155707) q[0];
sx q[0];
rz(2.7219211) q[0];
rz(-2.2113129) q[1];
sx q[1];
rz(-2.1280961) q[1];
sx q[1];
rz(-2.1297661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7714846) q[0];
sx q[0];
rz(-0.91910942) q[0];
sx q[0];
rz(-1.6328567) q[0];
x q[1];
rz(2.490245) q[2];
sx q[2];
rz(-2.911951) q[2];
sx q[2];
rz(-0.95637134) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8817084) q[1];
sx q[1];
rz(-1.0233423) q[1];
sx q[1];
rz(0.24357067) q[1];
x q[2];
rz(2.1675918) q[3];
sx q[3];
rz(-1.7241641) q[3];
sx q[3];
rz(-2.0917348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4100627) q[2];
sx q[2];
rz(-2.5750934) q[2];
sx q[2];
rz(-2.0274963) q[2];
rz(-1.7690432) q[3];
sx q[3];
rz(-1.0053582) q[3];
sx q[3];
rz(-0.70926386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0396742) q[0];
sx q[0];
rz(-1.2693951) q[0];
sx q[0];
rz(-1.6495548) q[0];
rz(2.6986625) q[1];
sx q[1];
rz(-1.8699173) q[1];
sx q[1];
rz(-0.84239229) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7105167) q[0];
sx q[0];
rz(-2.1177409) q[0];
sx q[0];
rz(0.7677107) q[0];
rz(-2.2350253) q[2];
sx q[2];
rz(-1.3593055) q[2];
sx q[2];
rz(1.4480526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4381899) q[1];
sx q[1];
rz(-2.3329607) q[1];
sx q[1];
rz(2.1548197) q[1];
x q[2];
rz(3.0987175) q[3];
sx q[3];
rz(-1.4673641) q[3];
sx q[3];
rz(-0.34260157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81279057) q[2];
sx q[2];
rz(-0.91965681) q[2];
sx q[2];
rz(-3.1205422) q[2];
rz(2.1444495) q[3];
sx q[3];
rz(-0.54233426) q[3];
sx q[3];
rz(-1.5525345) q[3];
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
rz(-1.7954623) q[0];
sx q[0];
rz(-1.2393476) q[0];
sx q[0];
rz(1.963266) q[0];
rz(2.1306253) q[1];
sx q[1];
rz(-1.6594454) q[1];
sx q[1];
rz(2.7706026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3431992) q[0];
sx q[0];
rz(-1.6863375) q[0];
sx q[0];
rz(-0.2327524) q[0];
x q[1];
rz(0.095711555) q[2];
sx q[2];
rz(-1.40279) q[2];
sx q[2];
rz(-1.8883349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1474516) q[1];
sx q[1];
rz(-0.32869654) q[1];
sx q[1];
rz(2.0533153) q[1];
x q[2];
rz(2.1493672) q[3];
sx q[3];
rz(-0.5360837) q[3];
sx q[3];
rz(1.9158196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79169881) q[2];
sx q[2];
rz(-1.0827622) q[2];
sx q[2];
rz(-2.2278111) q[2];
rz(-2.6036116) q[3];
sx q[3];
rz(-0.83298433) q[3];
sx q[3];
rz(2.7624281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051788483) q[0];
sx q[0];
rz(-1.7772978) q[0];
sx q[0];
rz(2.9096933) q[0];
rz(-1.505835) q[1];
sx q[1];
rz(-1.4488522) q[1];
sx q[1];
rz(-2.4969782) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58471352) q[0];
sx q[0];
rz(-1.2290092) q[0];
sx q[0];
rz(0.36447592) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5368136) q[2];
sx q[2];
rz(-1.591601) q[2];
sx q[2];
rz(0.60415798) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.3642824) q[1];
sx q[1];
rz(-2.9743618) q[1];
sx q[1];
rz(-1.689453) q[1];
x q[2];
rz(-0.50231825) q[3];
sx q[3];
rz(-0.76767477) q[3];
sx q[3];
rz(-1.8319825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5171648) q[2];
sx q[2];
rz(-1.5394053) q[2];
sx q[2];
rz(1.8713162) q[2];
rz(-2.7020448) q[3];
sx q[3];
rz(-2.5054273) q[3];
sx q[3];
rz(-0.51072454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9257833) q[0];
sx q[0];
rz(-1.9010811) q[0];
sx q[0];
rz(2.1155758) q[0];
rz(2.6846474) q[1];
sx q[1];
rz(-0.89070717) q[1];
sx q[1];
rz(0.93422008) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1608604) q[0];
sx q[0];
rz(-2.6965158) q[0];
sx q[0];
rz(-2.829192) q[0];
x q[1];
rz(1.9719741) q[2];
sx q[2];
rz(-2.8658762) q[2];
sx q[2];
rz(-1.3324312) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6187541) q[1];
sx q[1];
rz(-1.4013023) q[1];
sx q[1];
rz(-2.3970669) q[1];
x q[2];
rz(2.281743) q[3];
sx q[3];
rz(-0.78892869) q[3];
sx q[3];
rz(-3.0956059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43941867) q[2];
sx q[2];
rz(-2.4472523) q[2];
sx q[2];
rz(2.0806506) q[2];
rz(-2.151978) q[3];
sx q[3];
rz(-1.7695534) q[3];
sx q[3];
rz(2.8770679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92985741) q[0];
sx q[0];
rz(-2.0462357) q[0];
sx q[0];
rz(-0.76475058) q[0];
rz(-1.1872956) q[1];
sx q[1];
rz(-1.4359513) q[1];
sx q[1];
rz(-0.5618701) q[1];
rz(-0.019657481) q[2];
sx q[2];
rz(-1.622762) q[2];
sx q[2];
rz(0.30773559) q[2];
rz(-2.017631) q[3];
sx q[3];
rz(-1.1458967) q[3];
sx q[3];
rz(0.55895373) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
