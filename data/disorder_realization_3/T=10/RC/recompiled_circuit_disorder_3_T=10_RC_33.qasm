OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(-3.1103818) q[0];
sx q[0];
rz(-2.6565235) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1026693) q[0];
sx q[0];
rz(-1.2352714) q[0];
sx q[0];
rz(1.194792) q[0];
rz(-pi) q[1];
rz(-0.094751058) q[2];
sx q[2];
rz(-2.13846) q[2];
sx q[2];
rz(-1.3324347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2990883) q[1];
sx q[1];
rz(-2.7111004) q[1];
sx q[1];
rz(-1.4166142) q[1];
x q[2];
rz(0.26478404) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(0.36987723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-0.031575354) q[2];
rz(1.8850373) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(-0.53952113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927412) q[0];
sx q[0];
rz(-0.451938) q[0];
sx q[0];
rz(1.4638204) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72210724) q[2];
sx q[2];
rz(-1.8849843) q[2];
sx q[2];
rz(-2.9177641) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62944618) q[1];
sx q[1];
rz(-1.0735895) q[1];
sx q[1];
rz(2.9260103) q[1];
rz(-pi) q[2];
rz(0.94445618) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(2.8163547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77984017) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(0.2562491) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019779531) q[0];
sx q[0];
rz(-2.461713) q[0];
sx q[0];
rz(0.42827423) q[0];
x q[1];
rz(2.696064) q[2];
sx q[2];
rz(-1.1330714) q[2];
sx q[2];
rz(2.90403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5445404) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(0.87045963) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20817169) q[3];
sx q[3];
rz(-0.18067193) q[3];
sx q[3];
rz(-0.65073035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(2.5668872) q[2];
rz(-1.3556708) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(0.25948778) q[0];
rz(-1.9909987) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-0.73192275) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67128348) q[0];
sx q[0];
rz(-2.0844315) q[0];
sx q[0];
rz(2.1156103) q[0];
rz(1.4540133) q[2];
sx q[2];
rz(-2.2846966) q[2];
sx q[2];
rz(-1.8494948) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54938984) q[1];
sx q[1];
rz(-1.7648186) q[1];
sx q[1];
rz(-2.8787896) q[1];
rz(0.51283522) q[3];
sx q[3];
rz(-0.25179112) q[3];
sx q[3];
rz(-2.4076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(2.990492) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(1.8792101) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-2.343822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3691468) q[0];
sx q[0];
rz(-1.6909084) q[0];
sx q[0];
rz(3.1390879) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7902137) q[2];
sx q[2];
rz(-1.8674388) q[2];
sx q[2];
rz(2.9079633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3623912) q[1];
sx q[1];
rz(-2.5070842) q[1];
sx q[1];
rz(-1.7290551) q[1];
rz(-1.2410774) q[3];
sx q[3];
rz(-1.9922678) q[3];
sx q[3];
rz(0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34565869) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-2.8395555) q[2];
rz(1.1473514) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.1557895) q[0];
rz(-1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(3.0715122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72567155) q[0];
sx q[0];
rz(-1.7416746) q[0];
sx q[0];
rz(2.9861949) q[0];
x q[1];
rz(1.0191392) q[2];
sx q[2];
rz(-0.86691228) q[2];
sx q[2];
rz(-2.4920419) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8225704) q[1];
sx q[1];
rz(-0.59148568) q[1];
sx q[1];
rz(0.2925847) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.035404215) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(2.281718) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(3.133657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82848362) q[0];
sx q[0];
rz(-0.61921739) q[0];
sx q[0];
rz(2.6909268) q[0];
rz(-pi) q[1];
rz(-1.0646348) q[2];
sx q[2];
rz(-2.3828265) q[2];
sx q[2];
rz(2.2366692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7759526) q[1];
sx q[1];
rz(-1.0675634) q[1];
sx q[1];
rz(-2.3535437) q[1];
rz(-0.90585917) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(2.0743899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(0.41637862) q[2];
rz(1.3683866) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(2.2015613) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.6392802) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72909268) q[0];
sx q[0];
rz(-1.3206498) q[0];
sx q[0];
rz(-1.5604707) q[0];
rz(-2.084923) q[2];
sx q[2];
rz(-2.3377315) q[2];
sx q[2];
rz(0.67509292) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18174905) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(0.77397857) q[1];
rz(-pi) q[2];
rz(-1.3650465) q[3];
sx q[3];
rz(-1.660941) q[3];
sx q[3];
rz(-0.84116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6442948) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(-0.30125695) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97312462) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(-0.43564963) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(2.7246144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1736261) q[0];
sx q[0];
rz(-0.91714232) q[0];
sx q[0];
rz(-3.0941512) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66395335) q[2];
sx q[2];
rz(-1.10154) q[2];
sx q[2];
rz(-0.19367733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.13874395) q[1];
sx q[1];
rz(-1.7094304) q[1];
sx q[1];
rz(1.0671875) q[1];
rz(-2.3973893) q[3];
sx q[3];
rz(-2.9508698) q[3];
sx q[3];
rz(-0.71803367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(-2.9337692) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41480961) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(1.5243994) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(1.3226002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1117301) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(1.8882621) q[0];
x q[1];
rz(-2.0199213) q[2];
sx q[2];
rz(-1.3438864) q[2];
sx q[2];
rz(-0.34044468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1982556) q[1];
sx q[1];
rz(-0.88611929) q[1];
sx q[1];
rz(-0.68590045) q[1];
rz(-pi) q[2];
rz(-1.3721458) q[3];
sx q[3];
rz(-2.47654) q[3];
sx q[3];
rz(0.24634758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(-1.2223876) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(1.5608816) q[2];
sx q[2];
rz(-2.0799939) q[2];
sx q[2];
rz(-1.6891198) q[2];
rz(1.3397459) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
