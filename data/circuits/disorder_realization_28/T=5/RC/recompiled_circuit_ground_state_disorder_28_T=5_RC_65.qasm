OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31792274) q[0];
sx q[0];
rz(1.6190785) q[0];
sx q[0];
rz(9.4879307) q[0];
rz(1.8237279) q[1];
sx q[1];
rz(-0.39931077) q[1];
sx q[1];
rz(-0.77594405) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6187817) q[0];
sx q[0];
rz(-1.543985) q[0];
sx q[0];
rz(0.045700886) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6444998) q[2];
sx q[2];
rz(-2.1462206) q[2];
sx q[2];
rz(-2.6445553) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9358238) q[1];
sx q[1];
rz(-0.73038126) q[1];
sx q[1];
rz(-1.8380866) q[1];
x q[2];
rz(-1.8100061) q[3];
sx q[3];
rz(-2.0612217) q[3];
sx q[3];
rz(1.3077206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45106384) q[2];
sx q[2];
rz(-1.3961926) q[2];
sx q[2];
rz(-0.53831354) q[2];
rz(-2.8484143) q[3];
sx q[3];
rz(-1.5436951) q[3];
sx q[3];
rz(1.8831016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9341768) q[0];
sx q[0];
rz(-2.2947312) q[0];
sx q[0];
rz(0.22879452) q[0];
rz(3.0711807) q[1];
sx q[1];
rz(-1.2733302) q[1];
sx q[1];
rz(-1.061903) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5492229) q[0];
sx q[0];
rz(-1.8350936) q[0];
sx q[0];
rz(2.6068674) q[0];
rz(-0.47669551) q[2];
sx q[2];
rz(-1.9856493) q[2];
sx q[2];
rz(1.7384993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3913566) q[1];
sx q[1];
rz(-1.8614381) q[1];
sx q[1];
rz(0.33927563) q[1];
x q[2];
rz(-2.7135647) q[3];
sx q[3];
rz(-1.7645538) q[3];
sx q[3];
rz(1.0953812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.215302) q[2];
sx q[2];
rz(-0.44507241) q[2];
sx q[2];
rz(-3.0253547) q[2];
rz(-0.91533533) q[3];
sx q[3];
rz(-1.6719619) q[3];
sx q[3];
rz(-0.62469971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6621458) q[0];
sx q[0];
rz(-1.5818469) q[0];
sx q[0];
rz(0.33552718) q[0];
rz(1.4115931) q[1];
sx q[1];
rz(-0.84452191) q[1];
sx q[1];
rz(1.1988877) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20366448) q[0];
sx q[0];
rz(-2.0988582) q[0];
sx q[0];
rz(0.33007921) q[0];
rz(2.6012318) q[2];
sx q[2];
rz(-2.4746404) q[2];
sx q[2];
rz(1.4948542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3366223) q[1];
sx q[1];
rz(-1.4545856) q[1];
sx q[1];
rz(2.3776965) q[1];
rz(-pi) q[2];
rz(1.1245804) q[3];
sx q[3];
rz(-0.5329105) q[3];
sx q[3];
rz(0.41934723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18232839) q[2];
sx q[2];
rz(-0.79572833) q[2];
sx q[2];
rz(1.1129334) q[2];
rz(0.5111323) q[3];
sx q[3];
rz(-1.7685578) q[3];
sx q[3];
rz(-2.6326877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48753259) q[0];
sx q[0];
rz(-3.0656116) q[0];
sx q[0];
rz(-2.7677166) q[0];
rz(-1.9591676) q[1];
sx q[1];
rz(-1.9805209) q[1];
sx q[1];
rz(-2.9808796) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99237011) q[0];
sx q[0];
rz(-0.1609896) q[0];
sx q[0];
rz(-2.5991254) q[0];
x q[1];
rz(-2.2437975) q[2];
sx q[2];
rz(-2.2429969) q[2];
sx q[2];
rz(0.77122818) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4121033) q[1];
sx q[1];
rz(-2.0827052) q[1];
sx q[1];
rz(-0.44513925) q[1];
rz(-pi) q[2];
rz(0.75467189) q[3];
sx q[3];
rz(-0.42413482) q[3];
sx q[3];
rz(-1.4671668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5339727) q[2];
sx q[2];
rz(-1.1090358) q[2];
sx q[2];
rz(0.71471659) q[2];
rz(-2.886582) q[3];
sx q[3];
rz(-1.8111572) q[3];
sx q[3];
rz(-0.79849201) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052499972) q[0];
sx q[0];
rz(-0.30464259) q[0];
sx q[0];
rz(2.3783045) q[0];
rz(0.25453645) q[1];
sx q[1];
rz(-1.3724047) q[1];
sx q[1];
rz(-1.6932142) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93212381) q[0];
sx q[0];
rz(-0.49977344) q[0];
sx q[0];
rz(0.80003663) q[0];
rz(0.82008597) q[2];
sx q[2];
rz(-1.4561355) q[2];
sx q[2];
rz(-1.3263758) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9336945) q[1];
sx q[1];
rz(-2.1842644) q[1];
sx q[1];
rz(1.3529569) q[1];
rz(1.9687327) q[3];
sx q[3];
rz(-2.2421561) q[3];
sx q[3];
rz(-1.9111247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.280507) q[2];
sx q[2];
rz(-2.0528767) q[2];
sx q[2];
rz(-2.1121934) q[2];
rz(-1.6086802) q[3];
sx q[3];
rz(-0.89919388) q[3];
sx q[3];
rz(-2.5077584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7480046) q[0];
sx q[0];
rz(-0.91359502) q[0];
sx q[0];
rz(-2.9314801) q[0];
rz(2.4402319) q[1];
sx q[1];
rz(-1.7751834) q[1];
sx q[1];
rz(2.0393541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022033545) q[0];
sx q[0];
rz(-2.4071472) q[0];
sx q[0];
rz(-0.84656113) q[0];
x q[1];
rz(2.4793825) q[2];
sx q[2];
rz(-1.5726316) q[2];
sx q[2];
rz(-2.3037132) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6893856) q[1];
sx q[1];
rz(-1.2579134) q[1];
sx q[1];
rz(-0.72231267) q[1];
x q[2];
rz(2.5999448) q[3];
sx q[3];
rz(-1.4337501) q[3];
sx q[3];
rz(0.26462072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29409757) q[2];
sx q[2];
rz(-1.5832486) q[2];
sx q[2];
rz(1.6451277) q[2];
rz(1.3930813) q[3];
sx q[3];
rz(-2.2025509) q[3];
sx q[3];
rz(2.4699672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74712718) q[0];
sx q[0];
rz(-1.324993) q[0];
sx q[0];
rz(-0.4766683) q[0];
rz(-1.3564302) q[1];
sx q[1];
rz(-1.2973659) q[1];
sx q[1];
rz(-2.2043601) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.40537) q[0];
sx q[0];
rz(-0.69044603) q[0];
sx q[0];
rz(1.5140945) q[0];
rz(0.71800443) q[2];
sx q[2];
rz(-1.6176301) q[2];
sx q[2];
rz(2.5281434) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8647178) q[1];
sx q[1];
rz(-1.9559033) q[1];
sx q[1];
rz(-0.44189806) q[1];
rz(-pi) q[2];
rz(-1.473533) q[3];
sx q[3];
rz(-1.124112) q[3];
sx q[3];
rz(2.6616055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.53200191) q[2];
sx q[2];
rz(-1.8014182) q[2];
sx q[2];
rz(1.0756005) q[2];
rz(-2.7514451) q[3];
sx q[3];
rz(-1.8307999) q[3];
sx q[3];
rz(-2.3384317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53668642) q[0];
sx q[0];
rz(-2.0669879) q[0];
sx q[0];
rz(2.5965776) q[0];
rz(0.99264985) q[1];
sx q[1];
rz(-2.1558709) q[1];
sx q[1];
rz(-1.6880021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3869868) q[0];
sx q[0];
rz(-2.0029066) q[0];
sx q[0];
rz(-2.291392) q[0];
rz(-pi) q[1];
rz(-3.0849669) q[2];
sx q[2];
rz(-1.7030962) q[2];
sx q[2];
rz(0.43435979) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0664898) q[1];
sx q[1];
rz(-0.76998341) q[1];
sx q[1];
rz(1.7335197) q[1];
rz(-0.66176118) q[3];
sx q[3];
rz(-0.076120928) q[3];
sx q[3];
rz(-2.1832328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7327026) q[2];
sx q[2];
rz(-1.3971034) q[2];
sx q[2];
rz(0.2962386) q[2];
rz(-0.57684165) q[3];
sx q[3];
rz(-2.7448765) q[3];
sx q[3];
rz(-1.860994) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7883564) q[0];
sx q[0];
rz(-2.3529973) q[0];
sx q[0];
rz(2.9175135) q[0];
rz(-0.60399404) q[1];
sx q[1];
rz(-0.9915587) q[1];
sx q[1];
rz(-1.4858861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1224269) q[0];
sx q[0];
rz(-2.6976524) q[0];
sx q[0];
rz(2.5371518) q[0];
rz(-pi) q[1];
rz(2.9913556) q[2];
sx q[2];
rz(-1.7662107) q[2];
sx q[2];
rz(3.0216211) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.09391245) q[1];
sx q[1];
rz(-0.96598066) q[1];
sx q[1];
rz(-1.259786) q[1];
rz(-pi) q[2];
rz(1.8698244) q[3];
sx q[3];
rz(-2.0439897) q[3];
sx q[3];
rz(-2.2132471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6622582) q[2];
sx q[2];
rz(-1.0666288) q[2];
sx q[2];
rz(2.7545641) q[2];
rz(-0.28338638) q[3];
sx q[3];
rz(-0.75205386) q[3];
sx q[3];
rz(0.40248218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1166444) q[0];
sx q[0];
rz(-2.8706757) q[0];
sx q[0];
rz(-1.1234294) q[0];
rz(-1.4943538) q[1];
sx q[1];
rz(-1.6554183) q[1];
sx q[1];
rz(-0.87596881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6982684) q[0];
sx q[0];
rz(-0.7900376) q[0];
sx q[0];
rz(2.7774071) q[0];
rz(-pi) q[1];
rz(0.89089762) q[2];
sx q[2];
rz(-0.65717893) q[2];
sx q[2];
rz(-2.3466722) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3800602) q[1];
sx q[1];
rz(-0.71543869) q[1];
sx q[1];
rz(-1.9081294) q[1];
rz(-0.42966162) q[3];
sx q[3];
rz(-1.4080099) q[3];
sx q[3];
rz(1.9229887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7312077) q[2];
sx q[2];
rz(-1.7816252) q[2];
sx q[2];
rz(0.13710055) q[2];
rz(-1.3827263) q[3];
sx q[3];
rz(-0.9226678) q[3];
sx q[3];
rz(1.2285129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41505861) q[0];
sx q[0];
rz(-2.6714323) q[0];
sx q[0];
rz(2.5191125) q[0];
rz(-2.7525735) q[1];
sx q[1];
rz(-2.7906489) q[1];
sx q[1];
rz(-0.53967478) q[1];
rz(2.6555003) q[2];
sx q[2];
rz(-1.9584502) q[2];
sx q[2];
rz(2.9050585) q[2];
rz(-2.9403654) q[3];
sx q[3];
rz(-1.9126127) q[3];
sx q[3];
rz(0.9780608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
