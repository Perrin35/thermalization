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
rz(1.0072768) q[0];
sx q[0];
rz(2.4957823) q[0];
sx q[0];
rz(9.0584005) q[0];
rz(-2.3925048) q[1];
sx q[1];
rz(-1.6269416) q[1];
sx q[1];
rz(-2.9989624) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8170022) q[0];
sx q[0];
rz(-2.3329371) q[0];
sx q[0];
rz(2.3456744) q[0];
rz(-pi) q[1];
rz(-2.5900625) q[2];
sx q[2];
rz(-0.26187927) q[2];
sx q[2];
rz(2.8079845) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2995404) q[1];
sx q[1];
rz(-1.6900542) q[1];
sx q[1];
rz(-2.0485318) q[1];
x q[2];
rz(-0.47358114) q[3];
sx q[3];
rz(-1.2788805) q[3];
sx q[3];
rz(1.1137258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7819405) q[2];
sx q[2];
rz(-2.6617229) q[2];
sx q[2];
rz(-1.5894319) q[2];
rz(1.747067) q[3];
sx q[3];
rz(-2.7701869) q[3];
sx q[3];
rz(-2.4594405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29779103) q[0];
sx q[0];
rz(-0.60972917) q[0];
sx q[0];
rz(0.30526701) q[0];
rz(1.3097395) q[1];
sx q[1];
rz(-2.7204308) q[1];
sx q[1];
rz(0.052068204) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16927707) q[0];
sx q[0];
rz(-0.18405633) q[0];
sx q[0];
rz(2.4993624) q[0];
rz(-pi) q[1];
rz(-1.270625) q[2];
sx q[2];
rz(-1.2121634) q[2];
sx q[2];
rz(1.7110362) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1250097) q[1];
sx q[1];
rz(-2.8270611) q[1];
sx q[1];
rz(2.1279863) q[1];
rz(-0.90272119) q[3];
sx q[3];
rz(-1.4666712) q[3];
sx q[3];
rz(-1.6834761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0569968) q[2];
sx q[2];
rz(-1.3330385) q[2];
sx q[2];
rz(1.4196654) q[2];
rz(2.2325884) q[3];
sx q[3];
rz(-0.82908583) q[3];
sx q[3];
rz(-1.4066633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97515714) q[0];
sx q[0];
rz(-1.1844013) q[0];
sx q[0];
rz(-1.3982406) q[0];
rz(2.1947529) q[1];
sx q[1];
rz(-1.0942752) q[1];
sx q[1];
rz(1.2907226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4104267) q[0];
sx q[0];
rz(-1.435758) q[0];
sx q[0];
rz(-3.0655087) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3554736) q[2];
sx q[2];
rz(-1.9232456) q[2];
sx q[2];
rz(-2.4967683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4358999) q[1];
sx q[1];
rz(-1.5997643) q[1];
sx q[1];
rz(1.5040008) q[1];
x q[2];
rz(1.3685143) q[3];
sx q[3];
rz(-2.0401388) q[3];
sx q[3];
rz(0.795006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0134086) q[2];
sx q[2];
rz(-0.42542294) q[2];
sx q[2];
rz(-2.181633) q[2];
rz(0.32402447) q[3];
sx q[3];
rz(-0.5158546) q[3];
sx q[3];
rz(0.44017756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1337166) q[0];
sx q[0];
rz(-0.77819264) q[0];
sx q[0];
rz(2.9414951) q[0];
rz(2.5307185) q[1];
sx q[1];
rz(-0.69823825) q[1];
sx q[1];
rz(-2.3470338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2539729) q[0];
sx q[0];
rz(-1.4591503) q[0];
sx q[0];
rz(2.991204) q[0];
rz(2.9616293) q[2];
sx q[2];
rz(-2.2359701) q[2];
sx q[2];
rz(2.6106204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.17758372) q[1];
sx q[1];
rz(-1.4191886) q[1];
sx q[1];
rz(-1.2535918) q[1];
rz(-2.4684687) q[3];
sx q[3];
rz(-0.84852058) q[3];
sx q[3];
rz(-0.22658843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2574629) q[2];
sx q[2];
rz(-0.56371671) q[2];
sx q[2];
rz(1.5719315) q[2];
rz(-1.7852768) q[3];
sx q[3];
rz(-1.9808199) q[3];
sx q[3];
rz(0.13264382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86303478) q[0];
sx q[0];
rz(-0.25535169) q[0];
sx q[0];
rz(-0.72364664) q[0];
rz(3.0408995) q[1];
sx q[1];
rz(-1.0483402) q[1];
sx q[1];
rz(3.0752799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6322256) q[0];
sx q[0];
rz(-0.57183391) q[0];
sx q[0];
rz(-1.454005) q[0];
x q[1];
rz(-2.5939717) q[2];
sx q[2];
rz(-0.98973083) q[2];
sx q[2];
rz(0.18085322) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3795586) q[1];
sx q[1];
rz(-0.67081106) q[1];
sx q[1];
rz(-1.2993774) q[1];
x q[2];
rz(-0.39855115) q[3];
sx q[3];
rz(-1.8722187) q[3];
sx q[3];
rz(0.66661807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2012607) q[2];
sx q[2];
rz(-1.2137493) q[2];
sx q[2];
rz(-2.9230389) q[2];
rz(-2.7471733) q[3];
sx q[3];
rz(-0.68587488) q[3];
sx q[3];
rz(2.7421537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0938996) q[0];
sx q[0];
rz(-0.91966367) q[0];
sx q[0];
rz(-0.0041051824) q[0];
rz(2.507569) q[1];
sx q[1];
rz(-1.5019633) q[1];
sx q[1];
rz(2.1852469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0633558) q[0];
sx q[0];
rz(-0.18785297) q[0];
sx q[0];
rz(2.9211723) q[0];
rz(1.9955098) q[2];
sx q[2];
rz(-2.0419952) q[2];
sx q[2];
rz(1.6933407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33935336) q[1];
sx q[1];
rz(-1.5652854) q[1];
sx q[1];
rz(0.27779419) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31305571) q[3];
sx q[3];
rz(-0.36479539) q[3];
sx q[3];
rz(-0.49653253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6213106) q[2];
sx q[2];
rz(-1.4352398) q[2];
sx q[2];
rz(0.83462805) q[2];
rz(-0.37072119) q[3];
sx q[3];
rz(-2.4407237) q[3];
sx q[3];
rz(2.4247775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653559) q[0];
sx q[0];
rz(-3.1270202) q[0];
sx q[0];
rz(-0.45005774) q[0];
rz(2.693148) q[1];
sx q[1];
rz(-0.83201718) q[1];
sx q[1];
rz(0.73743302) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0919439) q[0];
sx q[0];
rz(-1.7554325) q[0];
sx q[0];
rz(-2.0020657) q[0];
rz(-pi) q[1];
rz(-1.8476494) q[2];
sx q[2];
rz(-2.9947318) q[2];
sx q[2];
rz(3.0727149) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1597643) q[1];
sx q[1];
rz(-0.48063403) q[1];
sx q[1];
rz(-0.7483866) q[1];
x q[2];
rz(1.3025018) q[3];
sx q[3];
rz(-0.80888575) q[3];
sx q[3];
rz(2.9380662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2927148) q[2];
sx q[2];
rz(-0.60056168) q[2];
sx q[2];
rz(0.71504492) q[2];
rz(2.6333366) q[3];
sx q[3];
rz(-1.4628937) q[3];
sx q[3];
rz(1.7432632) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4051064) q[0];
sx q[0];
rz(-2.5953601) q[0];
sx q[0];
rz(-0.35838321) q[0];
rz(-0.37250039) q[1];
sx q[1];
rz(-1.7820216) q[1];
sx q[1];
rz(0.43982664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0956235) q[0];
sx q[0];
rz(-1.5295424) q[0];
sx q[0];
rz(0.79605632) q[0];
rz(2.0310287) q[2];
sx q[2];
rz(-1.2046736) q[2];
sx q[2];
rz(1.5879968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.81573) q[1];
sx q[1];
rz(-1.3943895) q[1];
sx q[1];
rz(-2.8730064) q[1];
x q[2];
rz(-2.8183455) q[3];
sx q[3];
rz(-1.519299) q[3];
sx q[3];
rz(0.53158376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0899352) q[2];
sx q[2];
rz(-3.0445485) q[2];
sx q[2];
rz(-0.70866054) q[2];
rz(-2.8900201) q[3];
sx q[3];
rz(-2.1422062) q[3];
sx q[3];
rz(3.1075509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8468903) q[0];
sx q[0];
rz(-2.1069694) q[0];
sx q[0];
rz(2.5501472) q[0];
rz(-3.1239608) q[1];
sx q[1];
rz(-2.089274) q[1];
sx q[1];
rz(-0.10094053) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2786733) q[0];
sx q[0];
rz(-1.8846719) q[0];
sx q[0];
rz(0.14213965) q[0];
rz(-0.63669651) q[2];
sx q[2];
rz(-1.6023738) q[2];
sx q[2];
rz(0.65177887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7314607) q[1];
sx q[1];
rz(-1.6900151) q[1];
sx q[1];
rz(-1.340926) q[1];
rz(-1.798257) q[3];
sx q[3];
rz(-2.449027) q[3];
sx q[3];
rz(0.3161968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4883604) q[2];
sx q[2];
rz(-2.4182352) q[2];
sx q[2];
rz(1.4867268) q[2];
rz(0.15445736) q[3];
sx q[3];
rz(-2.0363225) q[3];
sx q[3];
rz(-0.36623335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79009295) q[0];
sx q[0];
rz(-0.14227754) q[0];
sx q[0];
rz(1.0737786) q[0];
rz(-1.0723266) q[1];
sx q[1];
rz(-0.25389478) q[1];
sx q[1];
rz(3.0825739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.260215) q[0];
sx q[0];
rz(-2.9395967) q[0];
sx q[0];
rz(-2.9750573) q[0];
rz(-2.0201276) q[2];
sx q[2];
rz(-2.5856433) q[2];
sx q[2];
rz(-0.96937552) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5328588) q[1];
sx q[1];
rz(-0.97670943) q[1];
sx q[1];
rz(2.2284177) q[1];
x q[2];
rz(2.2667472) q[3];
sx q[3];
rz(-2.2922462) q[3];
sx q[3];
rz(1.9375999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.8129639) q[2];
sx q[2];
rz(-1.8495411) q[2];
sx q[2];
rz(-0.22895075) q[2];
rz(-1.9479343) q[3];
sx q[3];
rz(-2.8173859) q[3];
sx q[3];
rz(1.6875904) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0070294587) q[0];
sx q[0];
rz(-2.3352191) q[0];
sx q[0];
rz(2.4051608) q[0];
rz(-1.3855343) q[1];
sx q[1];
rz(-1.3722739) q[1];
sx q[1];
rz(-1.3442232) q[1];
rz(-2.110021) q[2];
sx q[2];
rz(-1.5466718) q[2];
sx q[2];
rz(-0.9654733) q[2];
rz(-0.28382873) q[3];
sx q[3];
rz(-2.5093439) q[3];
sx q[3];
rz(0.96998246) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
