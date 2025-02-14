OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9417579) q[0];
sx q[0];
rz(-2.4500442) q[0];
sx q[0];
rz(-2.3681695) q[0];
rz(3.050488) q[1];
sx q[1];
rz(-0.79371047) q[1];
sx q[1];
rz(-1.33338) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42595902) q[0];
sx q[0];
rz(-1.5797401) q[0];
sx q[0];
rz(2.483213) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0716419) q[2];
sx q[2];
rz(-1.3678275) q[2];
sx q[2];
rz(0.90107337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65705196) q[1];
sx q[1];
rz(-0.96772497) q[1];
sx q[1];
rz(0.48274996) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5035787) q[3];
sx q[3];
rz(-1.5878233) q[3];
sx q[3];
rz(2.5650781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29204631) q[2];
sx q[2];
rz(-1.8201733) q[2];
sx q[2];
rz(0.41479659) q[2];
rz(1.8997806) q[3];
sx q[3];
rz(-2.0262227) q[3];
sx q[3];
rz(-0.19150664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5135797) q[0];
sx q[0];
rz(-0.11592557) q[0];
sx q[0];
rz(2.7017748) q[0];
rz(0.10305931) q[1];
sx q[1];
rz(-0.28918806) q[1];
sx q[1];
rz(1.9724847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3713908) q[0];
sx q[0];
rz(-1.4318716) q[0];
sx q[0];
rz(0.94111218) q[0];
rz(-2.1610519) q[2];
sx q[2];
rz(-1.613918) q[2];
sx q[2];
rz(-2.0415963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.00094571908) q[1];
sx q[1];
rz(-2.2123811) q[1];
sx q[1];
rz(-2.916275) q[1];
x q[2];
rz(-1.3430994) q[3];
sx q[3];
rz(-2.4191354) q[3];
sx q[3];
rz(1.3355153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.64572) q[2];
sx q[2];
rz(-1.7503259) q[2];
sx q[2];
rz(-0.24620852) q[2];
rz(1.5496893) q[3];
sx q[3];
rz(-2.527745) q[3];
sx q[3];
rz(-1.7483819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.15568) q[0];
sx q[0];
rz(-1.14013) q[0];
sx q[0];
rz(2.929856) q[0];
rz(3.1302997) q[1];
sx q[1];
rz(-2.3432422) q[1];
sx q[1];
rz(-1.6976154) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31682184) q[0];
sx q[0];
rz(-1.5707234) q[0];
sx q[0];
rz(3.1064338) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.774774) q[2];
sx q[2];
rz(-2.8958671) q[2];
sx q[2];
rz(1.7468921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1423751) q[1];
sx q[1];
rz(-0.96227598) q[1];
sx q[1];
rz(2.0318248) q[1];
x q[2];
rz(0.45076799) q[3];
sx q[3];
rz(-1.6524466) q[3];
sx q[3];
rz(1.1686366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4766562) q[2];
sx q[2];
rz(-2.0659451) q[2];
sx q[2];
rz(2.5962489) q[2];
rz(2.2319345) q[3];
sx q[3];
rz(-2.6792512) q[3];
sx q[3];
rz(1.2130515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.42105168) q[0];
sx q[0];
rz(-0.67890972) q[0];
sx q[0];
rz(2.3140267) q[0];
rz(1.1830117) q[1];
sx q[1];
rz(-1.2748101) q[1];
sx q[1];
rz(1.8756728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7050305) q[0];
sx q[0];
rz(-1.6616735) q[0];
sx q[0];
rz(1.4687125) q[0];
rz(1.363344) q[2];
sx q[2];
rz(-1.5816763) q[2];
sx q[2];
rz(1.2258197) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75531653) q[1];
sx q[1];
rz(-2.1960253) q[1];
sx q[1];
rz(-3.1023854) q[1];
x q[2];
rz(-2.9814441) q[3];
sx q[3];
rz(-1.7584828) q[3];
sx q[3];
rz(1.2040862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2563235) q[2];
sx q[2];
rz(-1.4726535) q[2];
sx q[2];
rz(0.081776865) q[2];
rz(-2.8335588) q[3];
sx q[3];
rz(-0.48397288) q[3];
sx q[3];
rz(1.6035732) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643739) q[0];
sx q[0];
rz(-2.866221) q[0];
sx q[0];
rz(-0.070505738) q[0];
rz(0.33440691) q[1];
sx q[1];
rz(-0.53135482) q[1];
sx q[1];
rz(-2.6883584) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2775947) q[0];
sx q[0];
rz(-2.3454614) q[0];
sx q[0];
rz(0.73426883) q[0];
rz(1.1788751) q[2];
sx q[2];
rz(-1.7818613) q[2];
sx q[2];
rz(-3.0818194) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47667015) q[1];
sx q[1];
rz(-1.3241265) q[1];
sx q[1];
rz(1.6849405) q[1];
rz(-pi) q[2];
rz(-2.2599561) q[3];
sx q[3];
rz(-0.83587468) q[3];
sx q[3];
rz(2.3595032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8411023) q[2];
sx q[2];
rz(-1.6359685) q[2];
sx q[2];
rz(-2.9312768) q[2];
rz(0.39120832) q[3];
sx q[3];
rz(-0.20520964) q[3];
sx q[3];
rz(-0.44262639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5000358) q[0];
sx q[0];
rz(-2.0099202) q[0];
sx q[0];
rz(0.17182194) q[0];
rz(1.3459282) q[1];
sx q[1];
rz(-2.185952) q[1];
sx q[1];
rz(1.1489493) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1412107) q[0];
sx q[0];
rz(-0.11185574) q[0];
sx q[0];
rz(3.0684708) q[0];
rz(-pi) q[1];
rz(2.9093018) q[2];
sx q[2];
rz(-1.8477167) q[2];
sx q[2];
rz(-2.0611219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2408381) q[1];
sx q[1];
rz(-0.60610702) q[1];
sx q[1];
rz(1.89487) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9385185) q[3];
sx q[3];
rz(-1.0725029) q[3];
sx q[3];
rz(-0.6607252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12080869) q[2];
sx q[2];
rz(-0.9641996) q[2];
sx q[2];
rz(2.7210893) q[2];
rz(1.6546107) q[3];
sx q[3];
rz(-1.9872811) q[3];
sx q[3];
rz(-1.1544863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625951) q[0];
sx q[0];
rz(-1.2263466) q[0];
sx q[0];
rz(2.763789) q[0];
rz(1.7516349) q[1];
sx q[1];
rz(-1.708834) q[1];
sx q[1];
rz(-2.7094944) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8209131) q[0];
sx q[0];
rz(-2.8519676) q[0];
sx q[0];
rz(0.65307133) q[0];
x q[1];
rz(1.6040726) q[2];
sx q[2];
rz(-1.6634523) q[2];
sx q[2];
rz(-1.7750545) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3258265) q[1];
sx q[1];
rz(-1.5349421) q[1];
sx q[1];
rz(-0.46200606) q[1];
x q[2];
rz(2.7338846) q[3];
sx q[3];
rz(-0.93824905) q[3];
sx q[3];
rz(-0.11128259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7695339) q[2];
sx q[2];
rz(-0.62166119) q[2];
sx q[2];
rz(-0.2335693) q[2];
rz(-1.4522878) q[3];
sx q[3];
rz(-1.7100311) q[3];
sx q[3];
rz(-1.5805894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772783) q[0];
sx q[0];
rz(-1.0294585) q[0];
sx q[0];
rz(-2.8218063) q[0];
rz(-2.2794967) q[1];
sx q[1];
rz(-1.2513221) q[1];
sx q[1];
rz(-1.7822942) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0727973) q[0];
sx q[0];
rz(-2.2900343) q[0];
sx q[0];
rz(0.79049514) q[0];
x q[1];
rz(0.14131693) q[2];
sx q[2];
rz(-1.3719146) q[2];
sx q[2];
rz(-2.3945253) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8492125) q[1];
sx q[1];
rz(-1.3977726) q[1];
sx q[1];
rz(-0.72489691) q[1];
rz(-1.2008181) q[3];
sx q[3];
rz(-0.66744679) q[3];
sx q[3];
rz(-2.1521309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84997815) q[2];
sx q[2];
rz(-2.7225967) q[2];
sx q[2];
rz(-0.46763793) q[2];
rz(0.48401287) q[3];
sx q[3];
rz(-1.368112) q[3];
sx q[3];
rz(1.2776933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676232) q[0];
sx q[0];
rz(-1.6521709) q[0];
sx q[0];
rz(1.1349348) q[0];
rz(-3.0813772) q[1];
sx q[1];
rz(-0.73679149) q[1];
sx q[1];
rz(1.5312451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1015243) q[0];
sx q[0];
rz(-2.1520242) q[0];
sx q[0];
rz(2.5826497) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0410454) q[2];
sx q[2];
rz(-0.95472958) q[2];
sx q[2];
rz(0.043827961) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.169226) q[1];
sx q[1];
rz(-1.833583) q[1];
sx q[1];
rz(1.6714833) q[1];
rz(2.3594112) q[3];
sx q[3];
rz(-1.5986658) q[3];
sx q[3];
rz(0.079558209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81426364) q[2];
sx q[2];
rz(-1.9244104) q[2];
sx q[2];
rz(2.4129756) q[2];
rz(-0.21555756) q[3];
sx q[3];
rz(-1.39648) q[3];
sx q[3];
rz(1.0936776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88290596) q[0];
sx q[0];
rz(-0.031143324) q[0];
sx q[0];
rz(0.31797847) q[0];
rz(-1.7440354) q[1];
sx q[1];
rz(-1.8122858) q[1];
sx q[1];
rz(-2.8584282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69071373) q[0];
sx q[0];
rz(-2.9895999) q[0];
sx q[0];
rz(-0.99894036) q[0];
rz(-pi) q[1];
rz(3.0599784) q[2];
sx q[2];
rz(-2.3858262) q[2];
sx q[2];
rz(0.72521842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8258543) q[1];
sx q[1];
rz(-2.2675677) q[1];
sx q[1];
rz(0.61751868) q[1];
rz(-3.0491676) q[3];
sx q[3];
rz(-1.3636949) q[3];
sx q[3];
rz(0.3849808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91750034) q[2];
sx q[2];
rz(-1.516569) q[2];
sx q[2];
rz(0.089281233) q[2];
rz(1.3034405) q[3];
sx q[3];
rz(-2.3035514) q[3];
sx q[3];
rz(-0.97344056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3598809) q[0];
sx q[0];
rz(-2.4926873) q[0];
sx q[0];
rz(0.98095184) q[0];
rz(-0.5858865) q[1];
sx q[1];
rz(-1.8460907) q[1];
sx q[1];
rz(1.5150217) q[1];
rz(1.7513314) q[2];
sx q[2];
rz(-1.3611887) q[2];
sx q[2];
rz(1.4088189) q[2];
rz(0.2570493) q[3];
sx q[3];
rz(-1.201073) q[3];
sx q[3];
rz(1.4504688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
