OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0527394) q[0];
sx q[0];
rz(-2.7272447) q[0];
sx q[0];
rz(0.9847087) q[0];
rz(-0.86496487) q[1];
sx q[1];
rz(-0.84671658) q[1];
sx q[1];
rz(-2.4196978) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979209) q[0];
sx q[0];
rz(-1.1211485) q[0];
sx q[0];
rz(1.1051154) q[0];
rz(-pi) q[1];
rz(-0.73562311) q[2];
sx q[2];
rz(-1.8343226) q[2];
sx q[2];
rz(1.5401538) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5923323) q[1];
sx q[1];
rz(-1.6417786) q[1];
sx q[1];
rz(0.029220079) q[1];
rz(-pi) q[2];
rz(-1.9328362) q[3];
sx q[3];
rz(-2.7305121) q[3];
sx q[3];
rz(0.53888881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7559173) q[2];
sx q[2];
rz(-2.6188681) q[2];
sx q[2];
rz(-2.6965466) q[2];
rz(1.880315) q[3];
sx q[3];
rz(-1.4538366) q[3];
sx q[3];
rz(-1.6311579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181013) q[0];
sx q[0];
rz(-2.0913251) q[0];
sx q[0];
rz(2.4871248) q[0];
rz(-2.0596313) q[1];
sx q[1];
rz(-1.8257717) q[1];
sx q[1];
rz(-1.6771603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1337812) q[0];
sx q[0];
rz(-2.0424583) q[0];
sx q[0];
rz(-3.0852484) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0686321) q[2];
sx q[2];
rz(-2.1734997) q[2];
sx q[2];
rz(-0.51034865) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9935181) q[1];
sx q[1];
rz(-1.6119527) q[1];
sx q[1];
rz(2.4012474) q[1];
rz(-pi) q[2];
rz(-1.220318) q[3];
sx q[3];
rz(-1.1247903) q[3];
sx q[3];
rz(-1.9917909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64708853) q[2];
sx q[2];
rz(-0.096179811) q[2];
sx q[2];
rz(2.3829226) q[2];
rz(1.3282954) q[3];
sx q[3];
rz(-1.0814861) q[3];
sx q[3];
rz(-1.3767415) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3674304) q[0];
sx q[0];
rz(-1.6663015) q[0];
sx q[0];
rz(-0.037516315) q[0];
rz(-1.0388177) q[1];
sx q[1];
rz(-0.33375868) q[1];
sx q[1];
rz(-0.85938251) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5512269) q[0];
sx q[0];
rz(-2.1810576) q[0];
sx q[0];
rz(1.3406189) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0065303) q[2];
sx q[2];
rz(-2.3762868) q[2];
sx q[2];
rz(-2.0489592) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.22166477) q[1];
sx q[1];
rz(-1.2499735) q[1];
sx q[1];
rz(-2.0179835) q[1];
rz(2.476781) q[3];
sx q[3];
rz(-0.71352173) q[3];
sx q[3];
rz(2.2647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77217707) q[2];
sx q[2];
rz(-0.93537664) q[2];
sx q[2];
rz(-0.76425648) q[2];
rz(-0.94591004) q[3];
sx q[3];
rz(-1.2063682) q[3];
sx q[3];
rz(-0.8980155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1979444) q[0];
sx q[0];
rz(-2.2479489) q[0];
sx q[0];
rz(-1.3125032) q[0];
rz(2.8370044) q[1];
sx q[1];
rz(-2.1423788) q[1];
sx q[1];
rz(-0.33590683) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2755796) q[0];
sx q[0];
rz(-3.0678425) q[0];
sx q[0];
rz(-1.3763675) q[0];
rz(0.44563771) q[2];
sx q[2];
rz(-1.012371) q[2];
sx q[2];
rz(-1.7281009) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.057317352) q[1];
sx q[1];
rz(-1.0773939) q[1];
sx q[1];
rz(1.1019023) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24728878) q[3];
sx q[3];
rz(-0.22824057) q[3];
sx q[3];
rz(0.3397371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.765982) q[2];
sx q[2];
rz(-0.42372647) q[2];
sx q[2];
rz(-1.83164) q[2];
rz(1.8464108) q[3];
sx q[3];
rz(-0.83596197) q[3];
sx q[3];
rz(1.0999934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64422166) q[0];
sx q[0];
rz(-0.27232429) q[0];
sx q[0];
rz(-2.3261133) q[0];
rz(-0.71715912) q[1];
sx q[1];
rz(-1.6970789) q[1];
sx q[1];
rz(-1.9898293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718147) q[0];
sx q[0];
rz(-1.0564532) q[0];
sx q[0];
rz(1.8215979) q[0];
x q[1];
rz(1.8134591) q[2];
sx q[2];
rz(-1.8662226) q[2];
sx q[2];
rz(2.5385419) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1895441) q[1];
sx q[1];
rz(-0.77922339) q[1];
sx q[1];
rz(2.0019931) q[1];
x q[2];
rz(1.607343) q[3];
sx q[3];
rz(-1.2690971) q[3];
sx q[3];
rz(-2.4100565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.069245) q[2];
sx q[2];
rz(-2.3657511) q[2];
sx q[2];
rz(-1.5838712) q[2];
rz(-2.1690058) q[3];
sx q[3];
rz(-1.9105304) q[3];
sx q[3];
rz(1.7044273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4438181) q[0];
sx q[0];
rz(-2.8688718) q[0];
sx q[0];
rz(0.66993237) q[0];
rz(-2.2386235) q[1];
sx q[1];
rz(-0.65330708) q[1];
sx q[1];
rz(3.0526551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085468) q[0];
sx q[0];
rz(-0.65292299) q[0];
sx q[0];
rz(1.2364975) q[0];
rz(-pi) q[1];
rz(-0.64623364) q[2];
sx q[2];
rz(-1.7719194) q[2];
sx q[2];
rz(1.0758019) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.054259334) q[1];
sx q[1];
rz(-1.8197734) q[1];
sx q[1];
rz(3.0691248) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1286147) q[3];
sx q[3];
rz(-2.0484784) q[3];
sx q[3];
rz(1.1123808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6637806) q[2];
sx q[2];
rz(-2.0974443) q[2];
sx q[2];
rz(-0.021154724) q[2];
rz(2.2145005) q[3];
sx q[3];
rz(-1.5320211) q[3];
sx q[3];
rz(0.65897861) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053719036) q[0];
sx q[0];
rz(-1.3608195) q[0];
sx q[0];
rz(-1.1370283) q[0];
rz(2.673705) q[1];
sx q[1];
rz(-1.7158022) q[1];
sx q[1];
rz(0.55380026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0040941) q[0];
sx q[0];
rz(-1.8589673) q[0];
sx q[0];
rz(-0.96364809) q[0];
x q[1];
rz(-2.2442859) q[2];
sx q[2];
rz(-2.6061432) q[2];
sx q[2];
rz(2.698512) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77417654) q[1];
sx q[1];
rz(-1.2700081) q[1];
sx q[1];
rz(2.1117626) q[1];
rz(-pi) q[2];
rz(-0.99618995) q[3];
sx q[3];
rz(-1.6227727) q[3];
sx q[3];
rz(1.0135723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0849453) q[2];
sx q[2];
rz(-2.2979996) q[2];
sx q[2];
rz(2.4455369) q[2];
rz(-0.36786914) q[3];
sx q[3];
rz(-2.0335679) q[3];
sx q[3];
rz(2.8554816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.87357658) q[0];
sx q[0];
rz(-1.753267) q[0];
sx q[0];
rz(0.8152813) q[0];
rz(2.6878327) q[1];
sx q[1];
rz(-0.90180698) q[1];
sx q[1];
rz(2.544983) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0469358) q[0];
sx q[0];
rz(-0.8917745) q[0];
sx q[0];
rz(-0.89767098) q[0];
x q[1];
rz(-1.6758133) q[2];
sx q[2];
rz(-2.3987282) q[2];
sx q[2];
rz(-0.32203963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4441159) q[1];
sx q[1];
rz(-1.8702862) q[1];
sx q[1];
rz(-1.0018574) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1659307) q[3];
sx q[3];
rz(-1.941276) q[3];
sx q[3];
rz(2.3080829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.10670184) q[2];
sx q[2];
rz(-1.7009578) q[2];
sx q[2];
rz(-1.9218669) q[2];
rz(-1.8969511) q[3];
sx q[3];
rz(-1.605426) q[3];
sx q[3];
rz(2.9419148) q[3];
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
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.225746) q[0];
sx q[0];
rz(-2.2016278) q[0];
sx q[0];
rz(1.655727) q[0];
rz(-1.111221) q[1];
sx q[1];
rz(-1.8467555) q[1];
sx q[1];
rz(1.750754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4877473) q[0];
sx q[0];
rz(-1.5731205) q[0];
sx q[0];
rz(-0.33727686) q[0];
rz(2.7971091) q[2];
sx q[2];
rz(-1.4622258) q[2];
sx q[2];
rz(2.6950633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2272033) q[1];
sx q[1];
rz(-1.3628642) q[1];
sx q[1];
rz(-3.1062011) q[1];
rz(0.81318716) q[3];
sx q[3];
rz(-0.48661806) q[3];
sx q[3];
rz(0.77017654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.86604649) q[2];
sx q[2];
rz(-1.6547497) q[2];
sx q[2];
rz(0.1869959) q[2];
rz(-0.20728076) q[3];
sx q[3];
rz(-0.63396251) q[3];
sx q[3];
rz(-2.0980339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279376) q[0];
sx q[0];
rz(-2.3917103) q[0];
sx q[0];
rz(-0.0023181152) q[0];
rz(1.9193316) q[1];
sx q[1];
rz(-1.5536676) q[1];
sx q[1];
rz(1.7244171) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42444705) q[0];
sx q[0];
rz(-1.3657278) q[0];
sx q[0];
rz(1.4385094) q[0];
rz(-2.6409686) q[2];
sx q[2];
rz(-2.5412895) q[2];
sx q[2];
rz(2.9666117) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.05310381) q[1];
sx q[1];
rz(-1.3657059) q[1];
sx q[1];
rz(0.049909485) q[1];
x q[2];
rz(0.98448344) q[3];
sx q[3];
rz(-1.7147736) q[3];
sx q[3];
rz(-1.7881623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2554539) q[2];
sx q[2];
rz(-1.5172493) q[2];
sx q[2];
rz(2.9296866) q[2];
rz(2.741277) q[3];
sx q[3];
rz(-2.2369657) q[3];
sx q[3];
rz(1.3306085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0935852) q[0];
sx q[0];
rz(-1.2496017) q[0];
sx q[0];
rz(-1.6954419) q[0];
rz(-2.6493337) q[1];
sx q[1];
rz(-1.6620363) q[1];
sx q[1];
rz(-1.5580039) q[1];
rz(-1.3561495) q[2];
sx q[2];
rz(-1.7742021) q[2];
sx q[2];
rz(0.85819744) q[2];
rz(-2.0782804) q[3];
sx q[3];
rz(-2.0290658) q[3];
sx q[3];
rz(3.1334044) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
