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
rz(3.6291549) q[0];
sx q[0];
rz(3.5222375) q[0];
sx q[0];
rz(9.335523) q[0];
rz(-0.46269497) q[1];
sx q[1];
rz(-2.8609639) q[1];
sx q[1];
rz(1.7791003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69169626) q[0];
sx q[0];
rz(-1.5917771) q[0];
sx q[0];
rz(-2.7684661) q[0];
rz(-2.9930395) q[2];
sx q[2];
rz(-1.534621) q[2];
sx q[2];
rz(0.52001563) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0995347) q[1];
sx q[1];
rz(-2.1452234) q[1];
sx q[1];
rz(-1.3045207) q[1];
x q[2];
rz(2.6854679) q[3];
sx q[3];
rz(-1.6509735) q[3];
sx q[3];
rz(-2.2661569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6224299) q[2];
sx q[2];
rz(-0.76202718) q[2];
sx q[2];
rz(-2.2117174) q[2];
rz(0.30396384) q[3];
sx q[3];
rz(-1.3320987) q[3];
sx q[3];
rz(-0.16873321) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794432) q[0];
sx q[0];
rz(-1.5771663) q[0];
sx q[0];
rz(-1.8899348) q[0];
rz(-0.21580639) q[1];
sx q[1];
rz(-1.0390751) q[1];
sx q[1];
rz(3.0024517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9342095) q[0];
sx q[0];
rz(-1.73044) q[0];
sx q[0];
rz(2.2347694) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47454796) q[2];
sx q[2];
rz(-2.5926771) q[2];
sx q[2];
rz(-2.7868556) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3635837) q[1];
sx q[1];
rz(-0.87555779) q[1];
sx q[1];
rz(1.259205) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71897935) q[3];
sx q[3];
rz(-2.1785979) q[3];
sx q[3];
rz(1.833235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6824049) q[2];
sx q[2];
rz(-0.46899691) q[2];
sx q[2];
rz(-2.0501308) q[2];
rz(-1.5952716) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(-2.6923164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0794373) q[0];
sx q[0];
rz(-1.9154444) q[0];
sx q[0];
rz(0.011818258) q[0];
rz(2.2305409) q[1];
sx q[1];
rz(-1.3900737) q[1];
sx q[1];
rz(-1.3284838) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2171611) q[0];
sx q[0];
rz(-1.4086485) q[0];
sx q[0];
rz(-2.6171706) q[0];
rz(2.9507347) q[2];
sx q[2];
rz(-1.4578739) q[2];
sx q[2];
rz(2.0115122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3551533) q[1];
sx q[1];
rz(-0.82988534) q[1];
sx q[1];
rz(-2.607858) q[1];
rz(-1.4592917) q[3];
sx q[3];
rz(-1.498025) q[3];
sx q[3];
rz(-0.36117103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38982424) q[2];
sx q[2];
rz(-1.5247034) q[2];
sx q[2];
rz(-2.7218008) q[2];
rz(-1.3587562) q[3];
sx q[3];
rz(-0.32521453) q[3];
sx q[3];
rz(-1.9512008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457394) q[0];
sx q[0];
rz(-1.1379108) q[0];
sx q[0];
rz(2.6275291) q[0];
rz(-0.38814107) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(-1.2302037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82321754) q[0];
sx q[0];
rz(-1.8288695) q[0];
sx q[0];
rz(2.8283872) q[0];
rz(-pi) q[1];
rz(-2.7204334) q[2];
sx q[2];
rz(-0.76986662) q[2];
sx q[2];
rz(-1.0854967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88273222) q[1];
sx q[1];
rz(-0.89619918) q[1];
sx q[1];
rz(-1.098295) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7994653) q[3];
sx q[3];
rz(-0.91311753) q[3];
sx q[3];
rz(-2.2715037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5019044) q[2];
sx q[2];
rz(-1.014726) q[2];
sx q[2];
rz(0.91628966) q[2];
rz(-0.4942975) q[3];
sx q[3];
rz(-1.9435792) q[3];
sx q[3];
rz(-0.77427197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2275527) q[0];
sx q[0];
rz(-2.1844449) q[0];
sx q[0];
rz(0.48602948) q[0];
rz(-2.5388429) q[1];
sx q[1];
rz(-1.9662247) q[1];
sx q[1];
rz(-1.1154307) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19768342) q[0];
sx q[0];
rz(-2.3872774) q[0];
sx q[0];
rz(0.59342845) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79430842) q[2];
sx q[2];
rz(-1.8657078) q[2];
sx q[2];
rz(0.25943929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72157114) q[1];
sx q[1];
rz(-0.90769115) q[1];
sx q[1];
rz(2.1065191) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0904457) q[3];
sx q[3];
rz(-1.8060038) q[3];
sx q[3];
rz(3.1344828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9976161) q[2];
sx q[2];
rz(-2.82085) q[2];
sx q[2];
rz(-1.0616659) q[2];
rz(-1.6210506) q[3];
sx q[3];
rz(-1.1989667) q[3];
sx q[3];
rz(2.2475713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068037085) q[0];
sx q[0];
rz(-0.03980045) q[0];
sx q[0];
rz(-1.918248) q[0];
rz(2.6262737) q[1];
sx q[1];
rz(-0.66910187) q[1];
sx q[1];
rz(1.3656778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13785744) q[0];
sx q[0];
rz(-2.9107339) q[0];
sx q[0];
rz(2.0876711) q[0];
rz(2.0414822) q[2];
sx q[2];
rz(-2.9376279) q[2];
sx q[2];
rz(2.230913) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9357696) q[1];
sx q[1];
rz(-2.6403815) q[1];
sx q[1];
rz(-0.31822121) q[1];
x q[2];
rz(-1.6069769) q[3];
sx q[3];
rz(-1.2734969) q[3];
sx q[3];
rz(2.1373914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0627275) q[2];
sx q[2];
rz(-1.8026423) q[2];
sx q[2];
rz(-1.8143181) q[2];
rz(1.4158538) q[3];
sx q[3];
rz(-0.073181987) q[3];
sx q[3];
rz(-0.88254005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0075204) q[0];
sx q[0];
rz(-2.6044758) q[0];
sx q[0];
rz(-2.9448331) q[0];
rz(-0.21601954) q[1];
sx q[1];
rz(-2.1327503) q[1];
sx q[1];
rz(0.78741995) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741873) q[0];
sx q[0];
rz(-1.3816773) q[0];
sx q[0];
rz(-1.322079) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2741267) q[2];
sx q[2];
rz(-2.0624522) q[2];
sx q[2];
rz(1.9751538) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.710121) q[1];
sx q[1];
rz(-0.97968757) q[1];
sx q[1];
rz(1.6256871) q[1];
x q[2];
rz(2.9522252) q[3];
sx q[3];
rz(-1.8872617) q[3];
sx q[3];
rz(-2.6031983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54062033) q[2];
sx q[2];
rz(-1.5121907) q[2];
sx q[2];
rz(2.7579894) q[2];
rz(0.90841928) q[3];
sx q[3];
rz(-2.3288265) q[3];
sx q[3];
rz(-2.9878591) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28021321) q[0];
sx q[0];
rz(-2.6160243) q[0];
sx q[0];
rz(0.63661611) q[0];
rz(-2.5909297) q[1];
sx q[1];
rz(-1.6267585) q[1];
sx q[1];
rz(-2.6689463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7028684) q[0];
sx q[0];
rz(-1.4819549) q[0];
sx q[0];
rz(-0.6779365) q[0];
rz(-pi) q[1];
rz(2.6246895) q[2];
sx q[2];
rz(-3.0344525) q[2];
sx q[2];
rz(-0.60984367) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3797261) q[1];
sx q[1];
rz(-1.7032924) q[1];
sx q[1];
rz(-1.6220868) q[1];
rz(2.5100559) q[3];
sx q[3];
rz(-1.3270006) q[3];
sx q[3];
rz(2.7840028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4281281) q[2];
sx q[2];
rz(-1.2423923) q[2];
sx q[2];
rz(-3.0430651) q[2];
rz(-0.030700961) q[3];
sx q[3];
rz(-1.5701598) q[3];
sx q[3];
rz(-2.9876685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6105662) q[0];
sx q[0];
rz(-2.1838146) q[0];
sx q[0];
rz(0.73262334) q[0];
rz(0.0064119617) q[1];
sx q[1];
rz(-1.0259722) q[1];
sx q[1];
rz(-1.7195255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2910478) q[0];
sx q[0];
rz(-1.8191162) q[0];
sx q[0];
rz(0.19638176) q[0];
x q[1];
rz(1.2071916) q[2];
sx q[2];
rz(-2.0418613) q[2];
sx q[2];
rz(2.1944012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3141503) q[1];
sx q[1];
rz(-1.536751) q[1];
sx q[1];
rz(-2.9877547) q[1];
rz(-0.51958618) q[3];
sx q[3];
rz(-1.5056464) q[3];
sx q[3];
rz(1.5198073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0942568) q[2];
sx q[2];
rz(-1.3553268) q[2];
sx q[2];
rz(2.6195841) q[2];
rz(-0.020307288) q[3];
sx q[3];
rz(-2.4419407) q[3];
sx q[3];
rz(2.844753) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2688399) q[0];
sx q[0];
rz(-1.5248542) q[0];
sx q[0];
rz(2.010345) q[0];
rz(2.3616683) q[1];
sx q[1];
rz(-1.9386539) q[1];
sx q[1];
rz(-0.47526971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6492622) q[0];
sx q[0];
rz(-1.9340252) q[0];
sx q[0];
rz(0.80357768) q[0];
rz(-pi) q[1];
rz(-0.23481253) q[2];
sx q[2];
rz(-0.91660344) q[2];
sx q[2];
rz(1.7476817) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1058029) q[1];
sx q[1];
rz(-0.33638182) q[1];
sx q[1];
rz(0.85806429) q[1];
rz(-pi) q[2];
rz(0.17367878) q[3];
sx q[3];
rz(-2.4376737) q[3];
sx q[3];
rz(0.46797637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37107006) q[2];
sx q[2];
rz(-1.5658242) q[2];
sx q[2];
rz(-1.2062262) q[2];
rz(-1.7462339) q[3];
sx q[3];
rz(-2.5985056) q[3];
sx q[3];
rz(-3.0622845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4974925) q[0];
sx q[0];
rz(-2.946749) q[0];
sx q[0];
rz(-0.83644833) q[0];
rz(2.1438228) q[1];
sx q[1];
rz(-2.0496968) q[1];
sx q[1];
rz(-3.0539378) q[1];
rz(3.072425) q[2];
sx q[2];
rz(-1.3528878) q[2];
sx q[2];
rz(2.4895346) q[2];
rz(3.0042778) q[3];
sx q[3];
rz(-2.4581494) q[3];
sx q[3];
rz(0.70829151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
