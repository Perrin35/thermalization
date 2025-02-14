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
rz(-3.0523678) q[0];
sx q[0];
rz(-2.3951946) q[0];
sx q[0];
rz(2.449692) q[0];
rz(-2.6785985) q[1];
sx q[1];
rz(-1.5249335) q[1];
sx q[1];
rz(-2.1114299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2041516) q[0];
sx q[0];
rz(-2.5290956) q[0];
sx q[0];
rz(2.7322311) q[0];
rz(0.18333034) q[2];
sx q[2];
rz(-1.1784484) q[2];
sx q[2];
rz(-3.0073056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9889415) q[1];
sx q[1];
rz(-1.2458774) q[1];
sx q[1];
rz(-0.12293651) q[1];
x q[2];
rz(-1.6923156) q[3];
sx q[3];
rz(-2.2560824) q[3];
sx q[3];
rz(0.44997643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87972155) q[2];
sx q[2];
rz(-2.3679831) q[2];
sx q[2];
rz(2.775178) q[2];
rz(2.053818) q[3];
sx q[3];
rz(-2.352114) q[3];
sx q[3];
rz(-0.68939775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.3159897) q[0];
sx q[0];
rz(-3.0593384) q[0];
sx q[0];
rz(-2.2218623) q[0];
rz(0.74291825) q[1];
sx q[1];
rz(-1.2580322) q[1];
sx q[1];
rz(1.119335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3920721) q[0];
sx q[0];
rz(-1.5095834) q[0];
sx q[0];
rz(-1.7158414) q[0];
rz(-2.934147) q[2];
sx q[2];
rz(-0.72978696) q[2];
sx q[2];
rz(1.6285673) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7537958) q[1];
sx q[1];
rz(-2.2860074) q[1];
sx q[1];
rz(1.8868576) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1503937) q[3];
sx q[3];
rz(-2.0127735) q[3];
sx q[3];
rz(2.9143765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6433581) q[2];
sx q[2];
rz(-1.9469807) q[2];
sx q[2];
rz(-0.052113459) q[2];
rz(-1.107996) q[3];
sx q[3];
rz(-2.2823157) q[3];
sx q[3];
rz(-3.0767379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68010083) q[0];
sx q[0];
rz(-1.9100186) q[0];
sx q[0];
rz(1.8517866) q[0];
rz(0.75311226) q[1];
sx q[1];
rz(-1.4537922) q[1];
sx q[1];
rz(0.49711102) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.755794) q[0];
sx q[0];
rz(-0.69911861) q[0];
sx q[0];
rz(2.505245) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92612626) q[2];
sx q[2];
rz(-2.1417649) q[2];
sx q[2];
rz(1.0356916) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0785574) q[1];
sx q[1];
rz(-0.37606701) q[1];
sx q[1];
rz(3.0381812) q[1];
rz(-pi) q[2];
rz(-3.0158832) q[3];
sx q[3];
rz(-0.65646986) q[3];
sx q[3];
rz(1.8905115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6734753) q[2];
sx q[2];
rz(-2.868729) q[2];
sx q[2];
rz(2.4172778) q[2];
rz(-2.916548) q[3];
sx q[3];
rz(-1.8658172) q[3];
sx q[3];
rz(2.3465274) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1555136) q[0];
sx q[0];
rz(-0.88548311) q[0];
sx q[0];
rz(-0.28050637) q[0];
rz(1.4836503) q[1];
sx q[1];
rz(-1.2018485) q[1];
sx q[1];
rz(2.6703506) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0126033) q[0];
sx q[0];
rz(-1.5480729) q[0];
sx q[0];
rz(-1.9381592) q[0];
rz(-pi) q[1];
rz(0.066023237) q[2];
sx q[2];
rz(-1.6834661) q[2];
sx q[2];
rz(2.6875935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.9556826) q[1];
sx q[1];
rz(-1.5430527) q[1];
sx q[1];
rz(-1.300059) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97906487) q[3];
sx q[3];
rz(-0.85460288) q[3];
sx q[3];
rz(1.5997052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9435141) q[2];
sx q[2];
rz(-0.55800617) q[2];
sx q[2];
rz(-1.0151218) q[2];
rz(2.4070168) q[3];
sx q[3];
rz(-1.3646804) q[3];
sx q[3];
rz(2.6034897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15767844) q[0];
sx q[0];
rz(-1.2175125) q[0];
sx q[0];
rz(2.1814046) q[0];
rz(-0.045711191) q[1];
sx q[1];
rz(-2.263133) q[1];
sx q[1];
rz(3.1083621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7901781) q[0];
sx q[0];
rz(-2.040967) q[0];
sx q[0];
rz(-2.1432502) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.015093283) q[2];
sx q[2];
rz(-1.8383257) q[2];
sx q[2];
rz(1.0969197) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3840365) q[1];
sx q[1];
rz(-1.9030326) q[1];
sx q[1];
rz(0.13793034) q[1];
rz(-pi) q[2];
rz(-1.1636249) q[3];
sx q[3];
rz(-1.1615175) q[3];
sx q[3];
rz(0.57151505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9783322) q[2];
sx q[2];
rz(-2.1425118) q[2];
sx q[2];
rz(0.31326374) q[2];
rz(2.7836109) q[3];
sx q[3];
rz(-1.0397747) q[3];
sx q[3];
rz(3.0430005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98704308) q[0];
sx q[0];
rz(-1.6313666) q[0];
sx q[0];
rz(-2.6556067) q[0];
rz(1.9588574) q[1];
sx q[1];
rz(-2.4458838) q[1];
sx q[1];
rz(0.33014306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81583038) q[0];
sx q[0];
rz(-1.6005524) q[0];
sx q[0];
rz(1.8682006) q[0];
rz(-pi) q[1];
rz(0.92440549) q[2];
sx q[2];
rz(-0.87714783) q[2];
sx q[2];
rz(2.2104757) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6961665) q[1];
sx q[1];
rz(-1.069134) q[1];
sx q[1];
rz(-1.0234246) q[1];
rz(-pi) q[2];
rz(-2.7948912) q[3];
sx q[3];
rz(-1.7809521) q[3];
sx q[3];
rz(-2.8488248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5064454) q[2];
sx q[2];
rz(-2.5250285) q[2];
sx q[2];
rz(0.35663024) q[2];
rz(0.59456524) q[3];
sx q[3];
rz(-1.6584572) q[3];
sx q[3];
rz(2.4026292) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4643788) q[0];
sx q[0];
rz(-2.3742299) q[0];
sx q[0];
rz(-2.5546524) q[0];
rz(0.37239536) q[1];
sx q[1];
rz(-1.3601235) q[1];
sx q[1];
rz(-0.24758235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42668396) q[0];
sx q[0];
rz(-1.6021906) q[0];
sx q[0];
rz(2.035729) q[0];
rz(-pi) q[1];
rz(-1.3273986) q[2];
sx q[2];
rz(-0.95275022) q[2];
sx q[2];
rz(-2.7400561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0928318) q[1];
sx q[1];
rz(-0.93834651) q[1];
sx q[1];
rz(-2.7238242) q[1];
rz(-pi) q[2];
rz(-2.7917737) q[3];
sx q[3];
rz(-1.8050151) q[3];
sx q[3];
rz(1.2522335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26329142) q[2];
sx q[2];
rz(-1.9204488) q[2];
sx q[2];
rz(0.2090052) q[2];
rz(-2.5448997) q[3];
sx q[3];
rz(-0.34691072) q[3];
sx q[3];
rz(-1.5359115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.7335032) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(-2.7485513) q[0];
rz(0.6706925) q[1];
sx q[1];
rz(-1.9978943) q[1];
sx q[1];
rz(-1.6938946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2358952) q[0];
sx q[0];
rz(-1.6950771) q[0];
sx q[0];
rz(1.0727364) q[0];
rz(-pi) q[1];
rz(-1.1426439) q[2];
sx q[2];
rz(-2.482126) q[2];
sx q[2];
rz(0.028837774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32998521) q[1];
sx q[1];
rz(-1.0701655) q[1];
sx q[1];
rz(2.5536584) q[1];
rz(2.3438934) q[3];
sx q[3];
rz(-1.9142861) q[3];
sx q[3];
rz(2.2960312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4286917) q[2];
sx q[2];
rz(-0.8608326) q[2];
sx q[2];
rz(-0.11670308) q[2];
rz(-1.0137001) q[3];
sx q[3];
rz(-1.2319535) q[3];
sx q[3];
rz(1.3056171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75334615) q[0];
sx q[0];
rz(-1.4660864) q[0];
sx q[0];
rz(-1.4353132) q[0];
rz(-0.995579) q[1];
sx q[1];
rz(-1.3187871) q[1];
sx q[1];
rz(-2.6011655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1216967) q[0];
sx q[0];
rz(-1.833124) q[0];
sx q[0];
rz(-2.645825) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2378906) q[2];
sx q[2];
rz(-1.0882749) q[2];
sx q[2];
rz(-2.7806892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7240043) q[1];
sx q[1];
rz(-0.88694983) q[1];
sx q[1];
rz(1.0761912) q[1];
rz(-pi) q[2];
rz(-2.5164817) q[3];
sx q[3];
rz(-1.907699) q[3];
sx q[3];
rz(-0.57527104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.871375) q[2];
sx q[2];
rz(-1.8367218) q[2];
sx q[2];
rz(0.77667856) q[2];
rz(3.1089879) q[3];
sx q[3];
rz(-1.5332581) q[3];
sx q[3];
rz(1.5743272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.8107574) q[0];
sx q[0];
rz(-0.36892712) q[0];
sx q[0];
rz(-2.708013) q[0];
rz(-2.4724204) q[1];
sx q[1];
rz(-1.9586261) q[1];
sx q[1];
rz(2.7152667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2992914) q[0];
sx q[0];
rz(-0.29336818) q[0];
sx q[0];
rz(2.3175453) q[0];
rz(-1.2157006) q[2];
sx q[2];
rz(-0.85596687) q[2];
sx q[2];
rz(-1.0155884) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5475216) q[1];
sx q[1];
rz(-1.2202377) q[1];
sx q[1];
rz(-0.5188397) q[1];
rz(2.5567708) q[3];
sx q[3];
rz(-1.9859145) q[3];
sx q[3];
rz(-0.41497542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7508037) q[2];
sx q[2];
rz(-0.28826737) q[2];
sx q[2];
rz(-1.1061579) q[2];
rz(3.028051) q[3];
sx q[3];
rz(-2.2083486) q[3];
sx q[3];
rz(0.043702628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.9857367) q[0];
sx q[0];
rz(-0.899213) q[0];
sx q[0];
rz(-0.56094299) q[0];
rz(0.38324311) q[1];
sx q[1];
rz(-1.1226729) q[1];
sx q[1];
rz(2.9164006) q[1];
rz(-0.46950317) q[2];
sx q[2];
rz(-2.532685) q[2];
sx q[2];
rz(-1.572123) q[2];
rz(1.0998955) q[3];
sx q[3];
rz(-0.80977898) q[3];
sx q[3];
rz(1.5017189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
