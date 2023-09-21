OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4164299) q[0];
sx q[0];
rz(-0.13983146) q[0];
sx q[0];
rz(-2.5319985) q[0];
rz(-5.6145515) q[1];
sx q[1];
rz(0.86548391) q[1];
sx q[1];
rz(15.794985) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5370109) q[0];
sx q[0];
rz(-2.4049979) q[0];
sx q[0];
rz(-1.0010615) q[0];
x q[1];
rz(-0.76724648) q[2];
sx q[2];
rz(-1.478072) q[2];
sx q[2];
rz(-0.35284943) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.062292369) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(2.9855707) q[1];
x q[2];
rz(-1.5243516) q[3];
sx q[3];
rz(-1.159045) q[3];
sx q[3];
rz(-0.1842894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13880754) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(-2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2574629) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.9467547) q[0];
rz(3.0589814) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(3.1412178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21391348) q[0];
sx q[0];
rz(-1.1535026) q[0];
sx q[0];
rz(0.37950619) q[0];
x q[1];
rz(-1.8325009) q[2];
sx q[2];
rz(-1.3771025) q[2];
sx q[2];
rz(-1.9270093) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61987585) q[1];
sx q[1];
rz(-0.65345018) q[1];
sx q[1];
rz(-0.39342777) q[1];
rz(-pi) q[2];
rz(-0.46044465) q[3];
sx q[3];
rz(-1.5213955) q[3];
sx q[3];
rz(-1.3747665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80883819) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(2.3475032) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(-0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.2394543) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(0.40476558) q[0];
rz(1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(1.3084897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8264016) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(-2.1000923) q[0];
rz(-pi) q[1];
rz(1.8070418) q[2];
sx q[2];
rz(-1.0457195) q[2];
sx q[2];
rz(2.6785148) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32654999) q[1];
sx q[1];
rz(-2.3800441) q[1];
sx q[1];
rz(-1.9425068) q[1];
x q[2];
rz(-2.7029196) q[3];
sx q[3];
rz(-1.4891426) q[3];
sx q[3];
rz(-0.7338394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.219316) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(0.02040872) q[2];
rz(-1.071788) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(-0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999009) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(-0.91598696) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4823285) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(2.7042424) q[0];
rz(-pi) q[1];
rz(-2.3154503) q[2];
sx q[2];
rz(-2.38378) q[2];
sx q[2];
rz(-3.0405424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6009439) q[1];
sx q[1];
rz(-1.800866) q[1];
sx q[1];
rz(0.01407108) q[1];
rz(1.7980868) q[3];
sx q[3];
rz(-2.3152581) q[3];
sx q[3];
rz(-0.28971653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(2.7988953) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.3559568) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(3.0551531) q[0];
rz(1.3899639) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(-2.6729029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3062895) q[0];
sx q[0];
rz(-1.1790743) q[0];
sx q[0];
rz(2.3098371) q[0];
rz(-pi) q[1];
rz(-1.7519978) q[2];
sx q[2];
rz(-1.8743519) q[2];
sx q[2];
rz(0.62523491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9710755) q[1];
sx q[1];
rz(-1.7726388) q[1];
sx q[1];
rz(-1.7032196) q[1];
rz(3.0607871) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(-0.24085837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1420574) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-0.86432499) q[2];
rz(-0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.7140759) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(0.19700225) q[0];
rz(-1.7794094) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(-2.8053455) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7173548) q[0];
sx q[0];
rz(-0.80168085) q[0];
sx q[0];
rz(-0.63071155) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1543428) q[2];
sx q[2];
rz(-1.0580214) q[2];
sx q[2];
rz(-3.0836881) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0834004) q[1];
sx q[1];
rz(-2.5501049) q[1];
sx q[1];
rz(0.085043474) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.398596) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(-2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7476615) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(-2.8826662) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(2.0475725) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0703697) q[0];
sx q[0];
rz(-0.56814146) q[0];
sx q[0];
rz(-0.90026654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6254243) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(0.54578997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8158405) q[1];
sx q[1];
rz(-3.0430803) q[1];
sx q[1];
rz(-1.2835531) q[1];
x q[2];
rz(1.592698) q[3];
sx q[3];
rz(-1.4727482) q[3];
sx q[3];
rz(2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(0.32067498) q[2];
rz(-2.705412) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86673474) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(2.136769) q[0];
rz(2.5514305) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(0.80387962) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893387) q[0];
sx q[0];
rz(-0.58710557) q[0];
sx q[0];
rz(-3.0703074) q[0];
rz(-2.5875823) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(-2.6487034) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4623973) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(-1.4560844) q[1];
rz(-pi) q[2];
rz(1.6611093) q[3];
sx q[3];
rz(-0.90925018) q[3];
sx q[3];
rz(2.6722398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(-0.60339749) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(0.55019125) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(-0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-2.3908652) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.663561) q[0];
sx q[0];
rz(-0.41304092) q[0];
sx q[0];
rz(-0.77126276) q[0];
rz(1.6523916) q[2];
sx q[2];
rz(-1.2921234) q[2];
sx q[2];
rz(2.1747053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1322861) q[1];
sx q[1];
rz(-2.7704151) q[1];
sx q[1];
rz(2.1585141) q[1];
rz(-pi) q[2];
rz(0.084089355) q[3];
sx q[3];
rz(-0.30232271) q[3];
sx q[3];
rz(1.7526527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0728545) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(-1.1661952) q[2];
rz(1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-3.1197746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5584548) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(-1.3903842) q[0];
rz(-0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.4601382) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4936192) q[0];
sx q[0];
rz(-0.68035347) q[0];
sx q[0];
rz(-0.86662678) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0748126) q[2];
sx q[2];
rz(-1.5911615) q[2];
sx q[2];
rz(2.6838944) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32669386) q[1];
sx q[1];
rz(-1.9683451) q[1];
sx q[1];
rz(-1.7705998) q[1];
x q[2];
rz(2.4246033) q[3];
sx q[3];
rz(-0.90962142) q[3];
sx q[3];
rz(2.9644074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(-1.6798518) q[2];
rz(-1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5158952) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(1.760578) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(2.3537221) q[2];
sx q[2];
rz(-1.315016) q[2];
sx q[2];
rz(-1.6223326) q[2];
rz(-0.67301987) q[3];
sx q[3];
rz(-1.2613847) q[3];
sx q[3];
rz(-3.1223084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];