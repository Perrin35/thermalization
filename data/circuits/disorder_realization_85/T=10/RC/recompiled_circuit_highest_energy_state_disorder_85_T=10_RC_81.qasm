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
rz(-2.149481) q[0];
sx q[0];
rz(-0.7465201) q[0];
sx q[0];
rz(2.574805) q[0];
rz(4.3920565) q[1];
sx q[1];
rz(5.1345706) q[1];
sx q[1];
rz(13.486025) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84595478) q[0];
sx q[0];
rz(-1.4104809) q[0];
sx q[0];
rz(2.3491067) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72534277) q[2];
sx q[2];
rz(-1.9646461) q[2];
sx q[2];
rz(-2.3462636) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7417042) q[1];
sx q[1];
rz(-2.2117858) q[1];
sx q[1];
rz(-1.6287032) q[1];
x q[2];
rz(-2.9018974) q[3];
sx q[3];
rz(-0.54876941) q[3];
sx q[3];
rz(-1.5769928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2016242) q[2];
sx q[2];
rz(-0.93279606) q[2];
sx q[2];
rz(-2.0882108) q[2];
rz(-0.71422226) q[3];
sx q[3];
rz(-0.37319365) q[3];
sx q[3];
rz(1.8478954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.5193704) q[0];
sx q[0];
rz(-0.70835963) q[0];
sx q[0];
rz(-2.569662) q[0];
rz(1.2988623) q[1];
sx q[1];
rz(-0.79890257) q[1];
sx q[1];
rz(0.78972185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2970726) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(-1.8734832) q[0];
rz(-0.65424325) q[2];
sx q[2];
rz(-1.6983319) q[2];
sx q[2];
rz(2.5728284) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4724628) q[1];
sx q[1];
rz(-1.2910559) q[1];
sx q[1];
rz(0.6196687) q[1];
x q[2];
rz(-1.0785473) q[3];
sx q[3];
rz(-2.1849602) q[3];
sx q[3];
rz(-1.9191161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85084891) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(-1.3624462) q[2];
rz(0.29081523) q[3];
sx q[3];
rz(-1.7751866) q[3];
sx q[3];
rz(-3.0349558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13880759) q[0];
sx q[0];
rz(-1.0951575) q[0];
sx q[0];
rz(-2.441067) q[0];
rz(0.48577148) q[1];
sx q[1];
rz(-0.829521) q[1];
sx q[1];
rz(-0.22451678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1248847) q[0];
sx q[0];
rz(-1.9210235) q[0];
sx q[0];
rz(-1.4136397) q[0];
rz(-pi) q[1];
rz(2.1491884) q[2];
sx q[2];
rz(-1.6534205) q[2];
sx q[2];
rz(2.3985661) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3697046) q[1];
sx q[1];
rz(-1.4955031) q[1];
sx q[1];
rz(-0.98702191) q[1];
x q[2];
rz(0.84214476) q[3];
sx q[3];
rz(-1.6557459) q[3];
sx q[3];
rz(2.2865975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.756788) q[2];
sx q[2];
rz(-0.89185682) q[2];
sx q[2];
rz(3.0925114) q[2];
rz(3.0745506) q[3];
sx q[3];
rz(-1.1578355) q[3];
sx q[3];
rz(-2.4750347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1753569) q[0];
sx q[0];
rz(-2.5844564) q[0];
sx q[0];
rz(0.019158451) q[0];
rz(1.3592023) q[1];
sx q[1];
rz(-1.7117585) q[1];
sx q[1];
rz(-3.137099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0877197) q[0];
sx q[0];
rz(-1.036265) q[0];
sx q[0];
rz(-2.1561532) q[0];
x q[1];
rz(-0.49764823) q[2];
sx q[2];
rz(-1.9714173) q[2];
sx q[2];
rz(1.7965574) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7178065) q[1];
sx q[1];
rz(-0.39925925) q[1];
sx q[1];
rz(1.1244133) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2035844) q[3];
sx q[3];
rz(-1.4798375) q[3];
sx q[3];
rz(1.1435777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4517453) q[2];
sx q[2];
rz(-2.0514026) q[2];
sx q[2];
rz(-1.6881662) q[2];
rz(-1.2545741) q[3];
sx q[3];
rz(-1.5213608) q[3];
sx q[3];
rz(-0.5622676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8101863) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(1.1619262) q[0];
rz(2.3792073) q[1];
sx q[1];
rz(-2.1547909) q[1];
sx q[1];
rz(-0.42957482) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7976101) q[0];
sx q[0];
rz(-1.6071885) q[0];
sx q[0];
rz(-1.162536) q[0];
x q[1];
rz(2.863071) q[2];
sx q[2];
rz(-2.3633011) q[2];
sx q[2];
rz(-2.8964004) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.030368806) q[1];
sx q[1];
rz(-0.75654034) q[1];
sx q[1];
rz(2.0732422) q[1];
rz(1.1384835) q[3];
sx q[3];
rz(-2.4412182) q[3];
sx q[3];
rz(-2.5770503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21542159) q[2];
sx q[2];
rz(-1.2913707) q[2];
sx q[2];
rz(-1.311709) q[2];
rz(-0.46071509) q[3];
sx q[3];
rz(-2.1381133) q[3];
sx q[3];
rz(-0.13319143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120414) q[0];
sx q[0];
rz(-2.9819745) q[0];
sx q[0];
rz(-1.106369) q[0];
rz(-3.0270992) q[1];
sx q[1];
rz(-1.0527) q[1];
sx q[1];
rz(3.0879424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.69218) q[0];
sx q[0];
rz(-0.25568889) q[0];
sx q[0];
rz(1.5480255) q[0];
rz(-2.3340763) q[2];
sx q[2];
rz(-1.4984594) q[2];
sx q[2];
rz(1.7010207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1785894) q[1];
sx q[1];
rz(-0.51437639) q[1];
sx q[1];
rz(1.6806755) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5074282) q[3];
sx q[3];
rz(-1.0552707) q[3];
sx q[3];
rz(0.17532119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2758808) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(-2.6344521) q[2];
rz(1.3745314) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(1.1486294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.62992612) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(3.1003057) q[0];
rz(0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(0.98947492) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51136298) q[0];
sx q[0];
rz(-1.0697027) q[0];
sx q[0];
rz(0.54176919) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6648952) q[2];
sx q[2];
rz(-0.39604353) q[2];
sx q[2];
rz(-0.11644289) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7478024) q[1];
sx q[1];
rz(-2.248924) q[1];
sx q[1];
rz(-2.3403779) q[1];
x q[2];
rz(0.37240828) q[3];
sx q[3];
rz(-0.95015991) q[3];
sx q[3];
rz(-2.2351976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5993293) q[2];
sx q[2];
rz(-2.0928536) q[2];
sx q[2];
rz(0.85912022) q[2];
rz(3.1298992) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(-0.11277994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3149253) q[0];
sx q[0];
rz(-0.59190094) q[0];
sx q[0];
rz(0.23183204) q[0];
rz(0.58206093) q[1];
sx q[1];
rz(-1.8128017) q[1];
sx q[1];
rz(1.2190762) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9909331) q[0];
sx q[0];
rz(-0.57123371) q[0];
sx q[0];
rz(-0.21749638) q[0];
rz(0.35301669) q[2];
sx q[2];
rz(-1.5215877) q[2];
sx q[2];
rz(-2.4491058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6775359) q[1];
sx q[1];
rz(-1.7221863) q[1];
sx q[1];
rz(-1.38358) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5737765) q[3];
sx q[3];
rz(-2.3006342) q[3];
sx q[3];
rz(-1.0882167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9570738) q[2];
sx q[2];
rz(-1.7902713) q[2];
sx q[2];
rz(2.3547724) q[2];
rz(-1.6400853) q[3];
sx q[3];
rz(-0.4314751) q[3];
sx q[3];
rz(-1.4260346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805434) q[0];
sx q[0];
rz(-1.3732055) q[0];
sx q[0];
rz(2.2013262) q[0];
rz(0.44479784) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(-0.14642265) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0068374182) q[0];
sx q[0];
rz(-1.869258) q[0];
sx q[0];
rz(1.7025856) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74898656) q[2];
sx q[2];
rz(-0.36452132) q[2];
sx q[2];
rz(-2.3836992) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9616016) q[1];
sx q[1];
rz(-0.72979425) q[1];
sx q[1];
rz(-0.37288937) q[1];
x q[2];
rz(1.6723456) q[3];
sx q[3];
rz(-0.78208215) q[3];
sx q[3];
rz(-2.8080432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0545097) q[2];
sx q[2];
rz(-1.6528218) q[2];
sx q[2];
rz(2.3040237) q[2];
rz(-2.2086823) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(2.3605409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8785716) q[0];
sx q[0];
rz(-1.913338) q[0];
sx q[0];
rz(1.7735057) q[0];
rz(-0.076016501) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(2.0215633) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.66066) q[0];
sx q[0];
rz(-1.2613861) q[0];
sx q[0];
rz(-0.27746986) q[0];
rz(-pi) q[1];
rz(2.7897808) q[2];
sx q[2];
rz(-2.4016671) q[2];
sx q[2];
rz(0.49525317) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.489913) q[1];
sx q[1];
rz(-1.2800242) q[1];
sx q[1];
rz(0.099822961) q[1];
x q[2];
rz(1.0088167) q[3];
sx q[3];
rz(-2.9491317) q[3];
sx q[3];
rz(-0.61103067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40176216) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(-2.8945727) q[2];
rz(-0.57010993) q[3];
sx q[3];
rz(-0.48738042) q[3];
sx q[3];
rz(-1.1183687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0433255) q[0];
sx q[0];
rz(-1.387431) q[0];
sx q[0];
rz(-0.36547216) q[0];
rz(3.0430766) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(1.6016207) q[2];
sx q[2];
rz(-2.7285103) q[2];
sx q[2];
rz(2.8693347) q[2];
rz(-0.3803654) q[3];
sx q[3];
rz(-1.3283397) q[3];
sx q[3];
rz(2.0793177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
