OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(-2.196329) q[0];
sx q[0];
rz(0.52559108) q[0];
rz(0.2431915) q[1];
sx q[1];
rz(4.3742124) q[1];
sx q[1];
rz(10.32962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63431595) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(2.1202203) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9154539) q[2];
sx q[2];
rz(-1.7625426) q[2];
sx q[2];
rz(-2.6297671) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9368254) q[1];
sx q[1];
rz(-1.1751886) q[1];
sx q[1];
rz(-2.8863465) q[1];
rz(-pi) q[2];
rz(2.1078029) q[3];
sx q[3];
rz(-2.783943) q[3];
sx q[3];
rz(1.8358177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-3.0453483) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(-2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(2.3764215) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239319) q[0];
sx q[0];
rz(-1.5045325) q[0];
sx q[0];
rz(-3.0793385) q[0];
x q[1];
rz(0.40132482) q[2];
sx q[2];
rz(-0.95762816) q[2];
sx q[2];
rz(-0.80712986) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4274802) q[1];
sx q[1];
rz(-0.44891) q[1];
sx q[1];
rz(2.950579) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29942056) q[3];
sx q[3];
rz(-2.6357108) q[3];
sx q[3];
rz(-2.8477856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(2.4760903) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5488141) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(-0.51868784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9496574) q[0];
sx q[0];
rz(-1.2149095) q[0];
sx q[0];
rz(2.3200672) q[0];
rz(0.2661744) q[2];
sx q[2];
rz(-1.6632348) q[2];
sx q[2];
rz(-0.12649525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4244713) q[1];
sx q[1];
rz(-0.79454225) q[1];
sx q[1];
rz(-2.8984927) q[1];
rz(-pi) q[2];
rz(0.054611562) q[3];
sx q[3];
rz(-1.5699937) q[3];
sx q[3];
rz(0.69566788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(-0.60244256) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(2.6113367) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78631567) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(-0.13537188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5330403) q[2];
sx q[2];
rz(-2.2366479) q[2];
sx q[2];
rz(2.7666639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2345703) q[1];
sx q[1];
rz(-1.7090194) q[1];
sx q[1];
rz(-0.45781086) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.030839132) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(-0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693562) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(-2.4647734) q[0];
rz(2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-0.61606032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84232932) q[0];
sx q[0];
rz(-0.063532524) q[0];
sx q[0];
rz(2.1701943) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6904171) q[2];
sx q[2];
rz(-1.2079442) q[2];
sx q[2];
rz(2.5809443) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9628145) q[1];
sx q[1];
rz(-1.5652834) q[1];
sx q[1];
rz(2.5639736) q[1];
x q[2];
rz(0.30619196) q[3];
sx q[3];
rz(-0.91727835) q[3];
sx q[3];
rz(-2.9432076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(-1.6051965) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42246321) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(2.6155112) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(-2.2568259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621983) q[0];
sx q[0];
rz(-1.4967376) q[0];
sx q[0];
rz(1.9522569) q[0];
rz(-pi) q[1];
rz(3.0753291) q[2];
sx q[2];
rz(-0.21050669) q[2];
sx q[2];
rz(1.9312242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5342418) q[1];
sx q[1];
rz(-0.51041767) q[1];
sx q[1];
rz(0.95153248) q[1];
x q[2];
rz(1.9725419) q[3];
sx q[3];
rz(-0.570795) q[3];
sx q[3];
rz(1.4626383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-2.8302144) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7235274) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-0.061070651) q[0];
rz(0.04018499) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(-2.4087002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6717364) q[0];
sx q[0];
rz(-2.0977019) q[0];
sx q[0];
rz(1.1080452) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.272846) q[2];
sx q[2];
rz(-1.1940496) q[2];
sx q[2];
rz(-1.6422611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1454754) q[1];
sx q[1];
rz(-2.5158094) q[1];
sx q[1];
rz(-1.5819509) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1002058) q[3];
sx q[3];
rz(-1.4635411) q[3];
sx q[3];
rz(1.0556575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(-0.51458365) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(2.432166) q[0];
rz(-1.5052694) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-0.27871305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.404971) q[0];
sx q[0];
rz(-1.706089) q[0];
sx q[0];
rz(2.9297329) q[0];
rz(-pi) q[1];
rz(-2.5829685) q[2];
sx q[2];
rz(-1.7785636) q[2];
sx q[2];
rz(-0.34780207) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0875138) q[1];
sx q[1];
rz(-1.2607818) q[1];
sx q[1];
rz(0.22884303) q[1];
x q[2];
rz(0.078019402) q[3];
sx q[3];
rz(-1.7934414) q[3];
sx q[3];
rz(-1.5732461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.30148208) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(2.2864443) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(2.7364031) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(0.96889281) q[0];
rz(-0.47337198) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(-1.1425346) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965217) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(1.2678498) q[0];
x q[1];
rz(-0.70334401) q[2];
sx q[2];
rz(-1.5726657) q[2];
sx q[2];
rz(-1.4700996) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7361703) q[1];
sx q[1];
rz(-2.1076964) q[1];
sx q[1];
rz(1.7169397) q[1];
x q[2];
rz(-0.90060602) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(-2.1538494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(-2.1208105) q[2];
rz(-0.3237237) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97994119) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(-1.9030301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4558444) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(-0.2483764) q[0];
x q[1];
rz(-1.5194703) q[2];
sx q[2];
rz(-2.3231299) q[2];
sx q[2];
rz(2.3527956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0760127) q[1];
sx q[1];
rz(-0.91846839) q[1];
sx q[1];
rz(-2.9933762) q[1];
x q[2];
rz(0.11390399) q[3];
sx q[3];
rz(-0.98131991) q[3];
sx q[3];
rz(1.9437499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4548268) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-0.043838538) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713521) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(2.8364137) q[2];
sx q[2];
rz(-1.1883493) q[2];
sx q[2];
rz(-2.3408008) q[2];
rz(-1.745789) q[3];
sx q[3];
rz(-2.3621515) q[3];
sx q[3];
rz(3.0099517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
